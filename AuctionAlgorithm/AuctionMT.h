/*
 * Auction.h
 *
 *  Created on: May 21, 2013
 *      Author: fb
 */

#ifndef AUCTIONMT_H_
#define AUCTIONMT_H_

#include "../internal.h"
#include "../ThreadBarrier.h"
#include "AuctionCommon.h"
#include <thread>

namespace LSAP
{

#define __AUCTION_DEBUG

/**
 *
 */
template<typename Scalar = double, typename MatrixType = Eigen::Matrix<double, -1, -1> >
class AuctionMT
{
public:

	typedef typename AuctionCommon<Scalar>::Scalars Scalars;
	typedef typename AuctionCommon<Scalar>::Locks Locks;
	typedef typename AuctionCommon<Scalar>::Indices Indices;
	typedef typename AuctionCommon<Scalar>::Edge Edge;
	typedef typename AuctionCommon<Scalar>::Edges Edges;
	typedef typename AuctionCommon<Scalar>::BidResults BidResults;

	AuctionMT() = delete;

	/**
	 * solve the assignment with the auction algorithm
	 * use real-types as coefficients, otherwise scaling will not work properly!
	 * @param a nxm weight matrix of type Scalar
	 * @return assignment matrix which has maximum value of objective function
	 */
	static Edges solve(const MatrixType & m, const size_t noOfThreads = 2) //, const Scalar e = __AUCTION_EPSILON_MULTIPLIER)
	{
		Edges edges;
		Scalars prices(m.cols(), 0.);
		Scalars profits(m.rows(), 1.);

		Scalar lambda = 0.;
		Scalar epsilon = (Scalar) __AUCTION_EPSILON_MULTIPLIER / m.cols();
		Locks lockedRows(m.rows(), false);
		Locks lockedCols(m.cols(), false);

		std::vector<std::thread> threadsForward, threadsReverse;

		// barriers for synchronization
		ThreadBarrier barrierParamsForward(noOfThreads);
		ThreadBarrier barrierResultsForward(noOfThreads);
		ThreadBarrier barrierParamsReverse(noOfThreads);
		ThreadBarrier barrierResultsReverse(noOfThreads);

		Indices iterationIntervalsForward(noOfThreads);

		Indices iterationIntervalsReverse (noOfThreads);

		// if matrix is dense eigen-type => pre-fill the iteration intervals for threads
		// otherwise do this within the forward/reverse iteration, because intervals depend on
		// column/row
		if ( std::is_same<MatrixType, Eigen::Matrix<Scalar, -1, -1>>::value )
			iterationIntervalsForward = AuctionCommon<Scalar>::splitToIntervals(noOfThreads, m.cols());

		// reverse intervals are the same for both matrix-types
		iterationIntervalsReverse = AuctionCommon<Scalar>::splitToIntervals(noOfThreads, m.rows());


		BidResults results(noOfThreads);

		size_t index = 0;	// representing row or column which is searched

		bool threadsActive = true;

		//bind threads to references and run them
		for (size_t t = 1; t < noOfThreads; t++)
		{
			std::function < void() > forwardThread = std::bind(threadForward, std::cref(m), t,
					std::cref(index), std::cref(iterationIntervalsForward[t-1]),
					std::cref(iterationIntervalsForward[t]),
					std::cref(prices), std::ref(results),
					std::ref(barrierParamsForward),
					std::ref(barrierResultsForward), std::ref(threadsActive));

			threadsForward.emplace_back(std::thread(forwardThread));

			std::function < void() > reverseThread = std::bind(threadReverse, std::cref(m), t,
					std::cref(index), std::cref(iterationIntervalsReverse[t-1]),
					std::cref(iterationIntervalsReverse[t]),
					std::cref(profits), std::ref(results),
					std::ref(barrierParamsReverse),
					std::ref(barrierResultsReverse), std::ref(threadsActive));

			threadsReverse.emplace_back(std::thread(reverseThread));
		}

		do
		{
			//		Step 1 (forward	auction cycle):
			//		Execute iterations of the forward auction algorithm until at least one
			//		more person becomes assigned. If there is an unassigned person left, go
			//		to step 2; else go to step 3.
			while (parallelForward(m, prices, profits, lockedRows, lockedCols,
					edges, lambda, epsilon, results, iterationIntervalsForward,
					noOfThreads, index, barrierParamsForward,
					barrierResultsForward))
				;

			if (!AuctionCommon<Scalar>::allPersonsAssigned(lockedRows))
			{
				//		Step 2 (reverse auction cycle):
				//		Execute several iterations of the reverse auction algorithm until at least
				//		one more object becomes assigned or until we have p_j <= lambda for all
				//		unassigned objects. If there is an unassigned person left, go to step 1
				//		else go to step 3
				while (!parallelReverse(m, prices, profits, lockedRows,
						lockedCols, edges, lambda, epsilon, results,
						iterationIntervalsReverse, noOfThreads, index,
						barrierParamsReverse, barrierResultsReverse)
						|| !AuctionCommon<Scalar>::unassignedObjectsLTlambda(
								lockedCols, prices, lambda))
					; // reverse auction
			}
			if (AuctionCommon<Scalar>::allPersonsAssigned(lockedRows))
			{
				//		Step 3 (reverse auction):
				//		Execute successive iterations of the reverse auction algorithm until the
				//		algorithm terminates with p_j <= lambda for all unassigned objects j
				while (true)
				{
					parallelReverse(m, prices, profits, lockedRows, lockedCols,
							edges, lambda, epsilon, results,
							iterationIntervalsReverse, noOfThreads, index,
							barrierParamsReverse, barrierResultsReverse);
					if (AuctionCommon<Scalar>::unassignedObjectsLTlambda(
							lockedCols, prices, lambda))
						break;
				}
				break;
			}
		} while (true);

		// kill threads hanging in barrier
		threadsActive = false;
		barrierParamsForward.kill();
		barrierParamsReverse.kill();
		barrierResultsForward.kill();
		barrierResultsReverse.kill();

		for (auto & t : threadsForward)
		{
			t.join();
		}
		for ( auto & t : threadsReverse)
		{
			t.join();
		}

		return edges;
	}

private:

	static void threadForward(const MatrixType & m, const size_t threadId, const size_t & i,
			const size_t & fromColIdx, const size_t & toColIdx, const Scalars & prices,
			BidResults & results, ThreadBarrier & barrierParams,
			ThreadBarrier & barrierResult, bool & threadActive)
	{
		while (threadActive)
		{
			if (barrierParams() ) return;
			AuctionCommon<Scalar>::forwardGS(m, threadId, i, fromColIdx, toColIdx, prices, results);
			if ( barrierResult() ) return;
		}
	}

	static void threadReverse(const MatrixType & m, const size_t threadId, const size_t & j,
			const size_t & fromRowIdx, const size_t & toRowIdx, const Scalars & profits,
			BidResults & results, ThreadBarrier & barrierParams,
			ThreadBarrier & barrierResult, bool & threadActive)
	{

		while (threadActive)
		{
			if (barrierParams()) return;
			AuctionCommon<Scalar>::reverseGS(m, threadId, j, fromRowIdx, toRowIdx,  profits, results);
			if (barrierResult()) return;
		}
	}

	static bool parallelReverse(const MatrixType & m, Scalars & prices,
			Scalars & profits, Locks & lockedRows, Locks & lockedCols,
			Edges & edges, Scalar & lambda, const Scalar epsilon,
			BidResults & results, Indices & iterationIntervals,
			const size_t noOfThreads, size_t & index,
			ThreadBarrier & barrierParams, ThreadBarrier & barrierResult)
	{

		bool assignmentFound = false;

		// wake up reverse threads
		barrierParams.wakeup();
		barrierResult.wakeup();

		// for all cols
		for (size_t j = 0; j < m.cols(); ++j)
		{
			//  column already locked ? p_j > lambda ?
			if (lockedCols[j])
				continue;

			if (!(prices[j] > lambda))
				continue;

			index = j; // global index for threads

			// sync and run local function for id 0
			barrierParams();
			AuctionCommon<Scalar>::reverseGS(m, 0, j, 0, iterationIntervals[0], profits, results);
			barrierResult();

			// values and indizes for the best assignment
			size_t i_j = 0;
			Scalar b_j = -__AUCTION_INF, g_j = -__AUCTION_INF;
			Scalar a_ij_j = 0.;
			bool foundBest = false; // found assignment?

			// find best and second best assignment in this row
			for (size_t t = 0; t < noOfThreads; t++)
			{
				// was an possible assignment found by the thread?
				if (results[t].assignmentFound)
				{
					bool found_b_j = false;

					// find best b_j and corresponding assignment
					if (results[t].bestBid > b_j)
					{
						g_j = b_j;
						a_ij_j = results[t].bestMatrixValue;
						i_j = results[t].bestIndex;
						b_j = results[t].bestBid;
						foundBest = found_b_j = true;
					}
					if (g_j < results[t].bestBid && !found_b_j)
						g_j = results[t].bestBid;
					if (g_j < results[t].secondBestBid)
						g_j = results[t].secondBestBid;
				}
			}

			if (foundBest) // best assignment found in row?
			{
				if (b_j >= lambda + epsilon)
				{
					const Scalar diff = g_j - epsilon; // G_j - E

					const Scalar max = lambda > diff ? lambda : diff; //  max { L, G_j - E}

					//	p_j = max { L, G_j - E}
					prices[j] = max;

					//	P_i_j = a_ij_j - max {L, G_j - E}
					profits[i_j] = a_ij_j - max;

					lockedRows[i_j] = true;
					lockedCols[j] = true;

					bool newEdgeFound = true;

					// if j_i was assigned to different i' to begin, remove (i', j_i) from S
					for (auto & e : edges)
						if (e.x == i_j) // change edge
						{
							lockedCols[e.y] = false; // unlock row i'
							newEdgeFound = false;
							e.y = j;
							e.v = a_ij_j;
							break;
						}
					if (newEdgeFound)
					{
						Edge newEdge;
						newEdge.x = i_j;
						newEdge.y = j;
						newEdge.v = a_ij_j;
						edges.push_back(newEdge);
					}

					assignmentFound = true;
				}
				else	// if B_j < L + E, case 2
				{
					//	p_j = B_j - E
					prices[j] = b_j - epsilon;

					/** standard lambda scaling **/
					unsigned int lowerThanLambda = 0;
					Scalar newLambda = lambda;

					// if the number of objectes k with p_k < lambda is bigger than (rows - cols)
					for (size_t k = 0; k < m.cols(); k++)
					{
						if (prices[k] < lambda)
						{
							lowerThanLambda++;
							if (prices[k] < newLambda)
								newLambda = prices[k];
						}
					}
					if (lowerThanLambda >= (m.cols() - m.rows()))
						lambda = newLambda;
				}
			}
		}

		barrierParams.sleep();
		barrierResult.sleep();

		return assignmentFound;

	}

	static bool parallelForward(const MatrixType & m, Scalars & prices,
			Scalars & profits, Locks & lockedRows, Locks & lockedCols,
			Edges & edges, const Scalar lambda, const Scalar epsilon,
			BidResults & results, Indices & iterationIntervals,
			const size_t noOfThreads, size_t & index,
			ThreadBarrier & barrierParams, ThreadBarrier & barrierResult)
	{
		bool assignmentFound = false;
		// for all rows

		barrierParams.wakeup();
		barrierResult.wakeup();

		for (size_t i = 0; i < m.rows(); ++i)
		{
			// if row is locked, there's already an assignment
			if (lockedRows[i])
				continue;

			index = i;

			// if sparse-matrix is used: split row/column dependent intervals for threads
			// offset is used for inner index offset
			size_t offset = 0;

			Indices it =
					AuctionCommon<Scalar>::splitToIntervalsByMatrixType(m, noOfThreads, i, true, offset );
			for ( size_t i = 0; i < it.size(); ++i ) iterationIntervals[i] = offset + it[i];


			// BARRIER
			barrierParams();
			AuctionCommon<Scalar>::forwardGS(m, 0, i, offset, iterationIntervals[0], prices, results);
			barrierResult();

			// values and indizes for the best assignment
			size_t j_i = 0;
			Scalar v_i = -__AUCTION_INF, w_i = -__AUCTION_INF;
			Scalar a_i_ji = 0.;
			bool foundBest = false; // found assignment?

			// find best and second best assignment in this row
			for (size_t t = 0; t < results.size(); t++)
			{
				// was an possible assignment found by the thread?
				if (results[t].assignmentFound)
				{
					bool found_v_i = false;
					// find best v_i and corresponding assignment
					if (results[t].bestBid > v_i)
					{
						w_i = v_i;
						a_i_ji = results[t].bestMatrixValue;
						j_i = results[t].bestIndex;
						v_i = results[t].bestBid;
						foundBest = found_v_i = true;
					}

					if (w_i < results[t].bestBid && !found_v_i)
						w_i = results[t].bestBid;
					if (w_i < results[t].secondBestBid)
						w_i = results[t].secondBestBid;
				}
			}

			if (foundBest) // best assignment found in row?
			{

				// P_i = w_i - epsilon
				profits[i] = w_i - epsilon;

				// diff = a(i,j_i) - w(i) + epsilon
				// prices(j_i) = max(lambda, a(i,j_i) - w(i) + epsilon)
				const Scalar diff = a_i_ji - w_i + epsilon;

				if (lambda <= diff)
				{
					prices[j_i] = diff;
					assignmentFound = true; // new assignment was found

					// assignment was made, so lock row and col
					lockedRows[i] = true;
					lockedCols[j_i] = true;

					bool newEdgeFound = true;

					// if j_i was assigned to different i' to begin, remove (i', j_i) from S
					for (auto & e : edges)
						if (e.y == j_i) // change edge to new assignment
						{
							newEdgeFound = false;
							lockedRows[e.x] = false; // unlock row
							e.v = a_i_ji;
							e.x = i;
							break;
						}

					if (newEdgeFound)	// add new edge to list
					{
						Edge newEdge;
						newEdge.x = i;
						newEdge.y = j_i;
						newEdge.v = a_i_ji;
						edges.push_back(newEdge);
//					std::cout << "adding edge (" << i << ", " << j_i << ")" << std::endl;
					}
				}
				else
					prices[j_i] = lambda;
			}
		}

		barrierParams.sleep();
		barrierResult.sleep();

		return assignmentFound;
	}

};

} /* namespace LSAP */
#endif /* AUCTIONMT_H_ */
