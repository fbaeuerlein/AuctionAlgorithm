/*
 * Auction.h
 *
 *  Created on: May 21, 2013
 *      Author: fb
 */

#ifndef AUCTION_H_
#define AUCTION_H_

#include "../internal.h"
#include "AuctionCommon.h"

//#define DEBUG_AUCTION

namespace LSAP
{

template<typename Scalar = double, typename MatrixType = Eigen::Matrix<double, -1, -1> >
class Auction
{
private:

	Auction()
	{
	}

	virtual ~Auction()
	{
	}
public:
	typedef typename AuctionCommon<Scalar>::Scalars Scalars;
	typedef typename AuctionCommon<Scalar>::Locks Locks;
	typedef typename AuctionCommon<Scalar>::Indizes Indizes;
	typedef typename AuctionCommon<Scalar>::Edge Edge;
	typedef typename AuctionCommon<Scalar>::Edges Edges;

	/**
	 * solve the assignment problem with the auction algorithm
	 * use real-types as coefficients, otherwise scaling will not work properly!
	 * @param a nxm weight matrix of type Scalar
	 * @return vector of Edges which represent the assignments
	 */
	static const Edges solve(const MatrixType & a)
	{
//		const Scalar e = __AUCTION_EPSILON_MULTIPLIER;
		const size_t rows = a.rows();
		const size_t cols = a.cols();

		Locks lockedRows(a.rows(), false);
		Locks lockedCols(a.cols(), false);
		Edges E;

		Scalar lambda = .0;
		Scalar epsilon = __AUCTION_EPSILON_MULTIPLIER / a.cols();

		// condition 3: initially set p_j >= lambda
		Scalars prices(cols, .0), profits(rows, 1.); // p-Vector  (1 to j) = p_j

		do
		{
			//		Step 1 (forward	auction cycle):
			//		Execute iterations of the forward auction algorithm until at least one
			//		more person becomes assigned. If there is an unassigned person left, go
			//		to step 2; else go to step 3.
			while (forward(a, E, prices, profits, lockedRows, lockedCols,
					lambda, epsilon))
				;

			if (!AuctionCommon<Scalar>::allPersonsAssigned(lockedRows))
			{

				//		Step 2 (reverse auction cycle):
				//		Execute several iterations of the reverse auction algorithm until at least
				//		one more object becomes assigned or until we have p_j <= lambda for all
				//		unassigned objects. If there is an unassigned person left, go to step 1
				//		else go to step 3
				while (!reverse(a, E, prices, profits, lockedRows, lockedCols,
						lambda, epsilon)
						|| !AuctionCommon<Scalar>::unassignedObjectsLTlambda(lockedCols, prices,
								lambda))
					; // reverse auction
			}

			if (AuctionCommon<Scalar>::allPersonsAssigned(lockedRows))
			{
				//		Step 3 (reverse auction):
				//		Execute successive iterations of the reverse auction algorithm until the
				//		algorithm terminates with p_j <= lambda for all unassigned objects j
				while (true)
				{
					reverse(a, E, prices, profits, lockedRows, lockedCols,
							lambda, epsilon);
					if (AuctionCommon<Scalar>::unassignedObjectsLTlambda(lockedCols, prices, lambda))
						break;
				}
				break;
			}

		} while (true);

#ifdef DEBUG_AUCTION
		for (auto i : E ) std::cout << i.x << ", " << i.y << " => " << prices[i.y] + profits[i.x]
		<< " a = " << a(i.x, i.y) << std::endl;
		PRINTVECTOR(prices)
		PRINTVECTOR(profits)
#endif
		return E;
	}
private:

	/**
	 * template specific implementation for finding the best entry in row
	 * used for dense matrix
	 * @param a	input matrix
	 * @param prices prices
	 * @param i row
	 * @param v_i best value
	 * @param w_i second best value
	 * @param a_i_ji value of entry
	 * @param j_i best entry
	 * @return true if assignment was found, otherwise false
	 */
	inline static const bool findForward(const Eigen::Matrix<Scalar, -1, -1> & a, const Scalars & prices, const size_t i,
			Scalar & v_i, Scalar & w_i, Scalar & a_i_ji, size_t & j_i)
	{
		const size_t cols = a.cols();

		bool assignmentFound = false;

		for (size_t j = 0; j < cols; j++) // for the j-th column
		{
			const Scalar aij = a(i,j);
			if ( aij == 0 ) continue;
			const Scalar diff = aij - prices[j];
			if (diff > v_i)
			{
				// if there already was an entry found, this is the second best
				if (assignmentFound)
					w_i = v_i;

				v_i = diff;
				j_i = j;
				a_i_ji = aij;
				assignmentFound = true;
			}
			if (diff > w_i && j_i != j)
				w_i = diff;
			// if no entry is bigger than v_i, check if there's still a bigger second best entry
		}
		return assignmentFound;
	}

	/**
	 * template specific implementation for finding the best entry in row
	 * using eigen sparse matrix with row major storage
	 * @param a	input matrix
	 * @param prices prices
	 * @param i row
	 * @param v_i best value
	 * @param w_i second best value
	 * @param a_i_ji value of entry
	 * @param j_i best entry
	 * @return true if assignment was found, otherwise false
	 */
	inline static const bool findForward(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & a, const Scalars & prices, const size_t i,
			Scalar & v_i, Scalar & w_i, Scalar & a_i_ji, size_t & j_i)
		{
//			assert(a.IsRowMajor);
			assert(a.isCompressed());

			bool assignmentFound = false;

			const auto outerStart = a.outerIndexPtr()[i];
			const auto outerEnd = a.outerIndexPtr()[i + 1];

			for (auto innerIdx = outerStart; innerIdx < outerEnd; innerIdx++) // for the j-th column
			{
				const auto j = a.innerIndexPtr()[innerIdx];
				const Scalar m_i_j = a.valuePtr()[innerIdx];
				const Scalar diff = m_i_j - prices[j];
				if (diff > v_i)
				{
					// if there already was an entry found, this is the second best
					if (assignmentFound)
						w_i = v_i;

					v_i = diff;
					j_i = j;
					a_i_ji = m_i_j;
					assignmentFound = true;
				}
				if (diff > w_i && j_i != j)
					w_i = diff;
				// if no entry is bigger than v_i, check if there's still a bigger second best entry
			}
			return assignmentFound;
		}

	/**
	 * forward cycle of auction algorithm
	 * @param a weight matrix (nxm)
	 * @param S assignment matrix (nxm)
	 * @param prices prices per object (m)
	 * @param profits profits per person (n)
	 * @param lambda bidding threshold lambda
	 * @param epsilon bidding increment
	 * @return true if assignment was made, false otherwise
	 */
	static bool forward(const MatrixType & a, Edges & E,
			Scalars & prices, Scalars & profits, Locks & lockedRows,
			Locks & lockedCols, Scalar & lambda, Scalar & epsilon)
	{
#ifdef DEBUG_AUCTION
		std::cout << "forwarding ..." << std::endl;
#endif
		const size_t rows = a.rows();
		bool assignmentFound = false;

		for (size_t i = 0; i < rows; i++) // for the i-th row/person
		{
			bool assignmentInThisIterationFound = false;

			// person already assigned?
			if (lockedRows[i])
				continue;

			// find an unassigned person i, its best object j_i
			// j_i = argmax {a_ij - p_j} for j in A(i) ( A(i) are the set of edges of the i-th row )
			// if a(i,j) = 0. it is not a valid edge
			size_t j_i = 0;

			//	v_i = max { a_ij - p_j} for j in A(i)				// maximum profit for person i
			//  v_i was already found = v_i
			//	w_i = max { a_ij - p_j} for j in A(i) and j != j_i  // second best profit
			//	if j_i is the only entry in A(i), w_i = - inf       // there's no second best profit
			Scalar w_i = -__AUCTION_INF, v_i = -__AUCTION_INF, a_i_ji = 0.;	// = max { a_ij - p_j}

			// find maximum profit i.e. j_i = arg max { a_ij - p_j} and second best
			// no possible assignment found?
			if (!findForward(a, prices, i, v_i, w_i, a_i_ji, j_i))
				continue;
#ifdef DEBUG_AUCTION
			std::cout << "i = " << i << " - j_i = " << j_i << " - v_i = " << v_i
			<< " - w_i = " << w_i << " - a_ij = " << a_i_ji << std::endl;
#endif
			assignmentInThisIterationFound = false;

//			std::cout << "assignment found .." << std::endl;
			const Scalar bid = a_i_ji - w_i + epsilon;

			//	P_i = w_i - E
			profits[i] = w_i - epsilon; // set new profit for person

			//	prices(j_i) = max(lambda, a(i,j_i) - w(i) + epsilon)
			// if lambda <= a_ij - w_i + E, add (i, j_i) to S
			if (lambda <= bid)
			{
				prices[j_i] = bid;
				// assignment was made, so lock row and col
				lockedRows[i] = true;
				lockedCols[j_i] = true;

				bool newEdge = true;

				// if j_i was assigned to different i' to begin, remove (i', j_i) from S
				for (auto & e : E)
					if (e.y == j_i) // change edge
					{
						lockedRows[e.x] = false; // unlock row i'
						newEdge = false;
						e.x = i;
						e.v = a_i_ji;
//						std::cout << " edges: " << E.size()
//								<< " changing edge ";
						break;
					}
				if (newEdge)
				{
					Edge e;
					e.x = i;
					e.y = j_i;
					e.v = a_i_ji;
					E.push_back(e);
//					std::cout << "edges: " << E.size() << " added edge ";
				}
				assignmentInThisIterationFound = true;
//				std::cout << "assignmentInThisIterationFound =  true"
//						<< std::endl;
			}
			else
			{
				prices[j_i] = lambda;
				assignmentInThisIterationFound = false;
//				std::cout << "assignmentInThisIterationFound = false"
//						<< std::endl;
			}
			if (assignmentInThisIterationFound)
				assignmentFound = true;
		}

//		std::cout << "returning " << ((assignmentFound) ? "true" : "false")
//				<< std::endl;
		return assignmentFound;

	}

	inline static const bool findReverse(const MatrixType & a, const Scalars & profits, const size_t j,
			Scalar & b_j, Scalar & g_j, Scalar & a_ij_j, size_t & i_j)
	{
		const size_t rows = a.rows();
		bool assignmentFound = false;

		for (size_t i = 0; i < rows; i++) // for the j-th column
		{
			const Scalar aij = a.coeff(i, j);
			if ( aij == 0 ) continue;
			const Scalar diff = aij - profits[i];
			if (diff > b_j)
			{
				// if there already was an entry found, this is the second best
				if (assignmentFound)
					g_j = b_j;

				b_j = diff;
				i_j = i;
				a_ij_j = aij;
				assignmentFound = true;
			}
			if (diff > g_j && i_j != i)
				g_j = diff;
		}
		return assignmentFound;
	}

	/**
	 * reverse cycle of auction algorithm
	 * @param a weight matrix (nxm)
	 * @param S assignment matrix (nxm)
	 * @param prices prices per object (m)
	 * @param profits profits per person (n)
	 * @param lambda bidding threshold lambda
	 * @param epsilon bidding increment
	 * @return true if assignment was made, false otherwise
	 */
	static bool reverse(const MatrixType & a, Edges & E,
			Scalars & prices, Scalars & profits, Locks & lockedRows,
			Locks & lockedCols, Scalar & lambda, const Scalar & epsilon)
	{
		const size_t rows = a.rows();
		const size_t cols = a.cols();

		bool assignmentFound = false;

		for (size_t j = 0; j < cols; j++) // for the j-th column (objects)
		{
			bool assignmentInThisIterationFound = false;

			// object already assigned,  p_j > lambda ?
			if (lockedCols[j])
				continue;

			if (!(prices[j] > lambda))
				continue;

			// Find an unassigned object j with p_j > lambda, its best person i_j
			// i_j = argmax {a_ij - profits[i]) für i aus B(j) (PI !!!)
			size_t i_j = 0;

			//g_j = max {a_ij - P_i} for i in B(j) and i != i_j
			// if j_i is the only entry in B(j), g_j = - inf ( g_j < b_j)
			//b_j = max {a_ij - P_i} for i in B(j)
			Scalar b_j = -__AUCTION_INF, g_j = -__AUCTION_INF, a_ij_j = 0.;

			// find maximum profit i.e. j_i = arg max { a_ij - p_j} and second best
			// no assignment found
			if (!findReverse(a, profits, j, b_j, g_j, a_ij_j, i_j))
				continue;

#ifdef DEBUG_AUCTION
			std::cout << "j = " << j << " i_j = " << i_j << " b_j = " << b_j << " g_j = " << g_j
			<< " a_ij_j = " << a_ij_j
			<< " p_j = " << prices[j] << " P_i = " << profits[i_j]<< std::endl;
#endif
			assignmentInThisIterationFound = false;

			//if b_j >= L + E, case 1:
			if (b_j >= (lambda + epsilon))
			{
#ifdef DEBUG_AUCTION
				std::cout << "  b_j >= lambda + epsilon" << std::endl;
#endif
				const Scalar diff = g_j - epsilon; // G_j - E

				const Scalar max = lambda > diff ? lambda : diff; //  max { L, G_j - E}

				//	p_j = max { L, G_j - E}
				prices[j] = max;

				//	P_i_j = a_i_jj - max {L, G_j - E}
				profits[i_j] = a_ij_j - max;

				lockedRows[i_j] = true;
				lockedCols[j] = true;

				bool newEdge = true;

				// if j_i was assigned to different i' to begin, remove (i', j_i) from S
				for (auto & e : E)
					if (e.x == i_j) // change edge
					{
						lockedCols[e.y] = false; // unlock row i'
						newEdge = false;
						e.y = j;
						e.v = a_ij_j;
#ifdef DEBUG_AUCTION
						std::cout << "  edges: " << E.size()
						<< "  changing edge ";
#endif
						break;

					}
				if (newEdge)
				{
					Edge e;
					e.x = i_j;
					e.y = j;
					e.v = a_ij_j;
					E.push_back(e);
#ifdef DEBUG_AUCTION
					std::cout << "  added edge " << E.size() << " ";
#endif
				}
				assignmentInThisIterationFound = true;
			}
			else	// if B_j < L + E, case 2
			{
				//	p_j = B_j - E
				prices[j] = b_j - epsilon;
#ifdef DEBUG_AUCTION
				std::cout << "  b_j < lambda + epsilon " << std::endl;
#endif
				/** standard lambda scaling **/
				size_t lowerThanLambda = 0;
				Scalar newLambda = lambda;

				// if the number of objectes k with p_k < lambda is bigger than (rows - cols)
				for (size_t k = 0; k < cols; k++)
				{
					if (prices[k] < lambda) // p_k < lambda
					{
						lowerThanLambda++;
						if (prices[k] < newLambda)
							newLambda = prices[k];
					}
				}
				// set new lambda
#ifdef DEBUG_AUCTION
				std::cerr << "  changing lambda from " << lambda << " to " << newLambda << std::endl;
#endif
				if (lowerThanLambda >= (cols - rows))
					lambda = newLambda;
				assignmentInThisIterationFound = false;
			}
			if (assignmentInThisIterationFound)
				assignmentFound = true;
		}
		return assignmentFound;
	}

private:

};

} /* namespace LSAP */
#endif /* AUCTION_H_ */
