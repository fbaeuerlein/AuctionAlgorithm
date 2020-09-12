/*
 * AuctionCommon.h
 *
 *  Created on: Jul 2, 2013
 *      Author: fb
 */

#ifndef AUCTIONCOMMON_H_
#define AUCTIONCOMMON_H_

#include "../internal.h"

namespace LSAP
{

template <typename Scalar = double>
class AuctionCommon
{
  public:
    struct Edge
    {
		Edge(size_t x_, size_t y_, Scalar v_ ) : x(x_), y(y_), v(v_) {}
        Edge() = default;
        size_t x;
        size_t y;
        Scalar v;
    };

    typedef typename std::vector<Scalar> Scalars;
    typedef typename std::vector<bool> Locks;
    typedef typename std::vector<size_t> Indices;
    typedef typename std::vector<Edge> Edges;
    typedef typename Eigen::Matrix<Scalar, -1, -1> WeightMatrix;

    class BidResult
    {

      public:
        BidResult() = default;

        BidResult(const size_t idx, const size_t bestIdx, const Scalar bestBid,
                  const Scalar secondBestBid, const Scalar bestMatrixValue)
            : index(idx)
            , bestIndex(bestIdx)
            , bestBid(bestBid)
            , secondBestBid(secondBestBid)
            , bestMatrixValue(bestMatrixValue)
            , assignmentFound(true)
        {
        }

        virtual ~BidResult() {}

        // set by thread
        size_t index{0};                     // index of examined row or column
        size_t bestIndex{0};                 // index of best entry
        Scalar bestBid{0}, secondBestBid{0}; // corresponding bids (w_i, v_i)
        Scalar bestMatrixValue{0};           // a_ij
        bool assignmentFound{false};         // assignment found ?
    };

    typedef typename std::vector<BidResult> BidResults;

    /**
     * gauss-seidel forward cycle for dense matrix
     * @param m
     * @param threadId
     * @param i
     * @param fromColIdx
     * @param toColIdx
     * @param prices
     * @param results
     */
    static void forwardGS(const WeightMatrix & m, const size_t threadId, const size_t & i,
                          const size_t & fromColIdx, const size_t & toColIdx,
                          const Scalars & prices, BidResults & results)
    {
        results[threadId].assignmentFound = false;
        size_t j_i = 0;
        Scalar v_i = -__AUCTION_INF, a_i_ji = 0., w_i = -__AUCTION_INF;

        // find maximum profit i.e. j_i = arg max { a_ij - p_j} and second best
        // iterate over columns
        for (size_t j = fromColIdx; j < toColIdx; j++) // for the j-th column
        {
            const Scalar m_i_j = m(i, j);
            if (m_i_j == 0)
                continue;

            const Scalar diff = m_i_j - prices[j];

            if (diff > v_i)
            {
                // if there already was an entry found, this is the second best
                if (results[threadId].assignmentFound)
                    w_i = v_i;

                a_i_ji = m_i_j;
                v_i = diff;
                j_i = j;
                results[threadId].assignmentFound = true;
            }
            if (diff > w_i && j_i != j)
                w_i = diff;
            // if no entry is bigger than v_i, check if there's still a bigger second
            // best entry
        }

        // best object for person i
        results[threadId].bestIndex = j_i;
        //	v_i = max { a_ij - p_j} for j in A(i)
        results[threadId].bestBid = v_i;
        //	w_i = max { a_ij - p_j} for j in A(i) and j != j_i
        results[threadId].secondBestBid = w_i;
        // a_i_ji = a(i, j_i), already found above
        results[threadId].bestMatrixValue = a_i_ji;
        return;
    }

    /**
     * gauss-seidel forward cycle for sparse matrix
     * @param m
     * @param threadId
     * @param i
     * @param fromColIdx
     * @param toColIdx
     * @param prices
     * @param results
     */
    static void forwardGS(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & m,
                          const size_t threadId, const size_t & i, const size_t & fromColIdx,
                          const size_t & toColIdx, const Scalars & prices, BidResults & results)
    {
        results[threadId].assignmentFound = false;

        size_t j_i = 0;
        Scalar v_i = -__AUCTION_INF, a_i_ji = 0., w_i = -__AUCTION_INF;

        //			std::cout << "threadForward: thread " << threadId << " iterating
        // from " << fromColIdx << " to " << toColIdx << std::endl;
        // find maximum profit i.e. j_i = arg max { a_ij - p_j} and second best
        // iterate over columns
        for (size_t innerIndex = fromColIdx; innerIndex < toColIdx;
             innerIndex++) // for the j-th column
        {

            const size_t j = m.innerIndexPtr()[innerIndex];
            const Scalar m_i_j = m.valuePtr()[innerIndex];
            const Scalar diff = m_i_j - prices[j];

            if (diff > v_i)
            {
                // if there already was an entry found, this is the second best
                if (results[threadId].assignmentFound)
                    w_i = v_i;

                a_i_ji = m_i_j;
                v_i = diff;
                j_i = j;
                results[threadId].assignmentFound = true;
            }
            if (diff > w_i && j_i != j)
                w_i = diff;
            // if no entry is bigger than v_i, check if there's still a bigger second
            // best entry
        }

        // best object for person i
        results[threadId].bestIndex = j_i;
        //	v_i = max { a_ij - p_j} for j in A(i)
        results[threadId].bestBid = v_i;
        //	w_i = max { a_ij - p_j} for j in A(i) and j != j_i
        results[threadId].secondBestBid = w_i;
        // a_i_ji = a(i, j_i), already found above
        results[threadId].bestMatrixValue = a_i_ji;
    }

    /**
     * gaus-seidel reverse cycle for sparse row major matrix
     * @param threadId
     * @param j
     * @param fromRowIdx
     * @param toRowIdx
     * @param m
     * @param profits
     * @param results
     */
    static void reverseGS(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & m,
                          const size_t threadId, const size_t & j, const size_t & fromRowIdx,
                          const size_t & toRowIdx, const Scalars & profits, BidResults & results)
    {

        // find maximum profit i.e. j_i = arg max { a_ij - p_j} and second best
        // iterate over columns
        size_t i_j = 0;
        Scalar b_j = -__AUCTION_INF, a_ij_j = 0.,
               g_j = -__AUCTION_INF; // = max { a_ij - P_i}

        // find maximum profit i.e. j_i = arg max { a_ij - p_j} and second best
        for (size_t i = fromRowIdx; i < toRowIdx; i++) // for the i-th row
        {
            const Scalar m_i_j = m.coeff(i, j);
            if (m_i_j == 0)
                continue;

            const Scalar diff = m_i_j - profits[i];

            if (diff > b_j)
            {
                // if there already was an entry found, this is the second best
                if (results[threadId].assignmentFound)
                    g_j = b_j;

                b_j = diff;
                i_j = i;
                a_ij_j = m_i_j;

                results[threadId].assignmentFound = true;
            }
            if (diff > g_j && i_j != i)
                g_j = diff;
        }

        // best person for object j
        results[threadId].bestIndex = i_j;
        //	b_i = max { a_ij - P_i} for i in A(j)
        results[threadId].bestBid = b_j;
        //	g_i = max { a_ij - P_i} for j in A(j) and i != i_j
        results[threadId].secondBestBid = g_j;
        // a_ij_j = a(i_j, j), already found above
        results[threadId].bestMatrixValue = a_ij_j;

        return;
    }

    /**
     * gaus-seidel reverse cycle for dense matrix
     * @param threadId
     * @param j
     * @param fromRowIdx
     * @param toRowIdx
     * @param m
     * @param profits
     * @param results
     */
    static void reverseGS(const WeightMatrix & m, const size_t threadId, const size_t & j,
                          const size_t & fromRowIdx, const size_t & toRowIdx,
                          const Scalars & profits, BidResults & results)
    {

        // find maximum profit i.e. j_i = arg max { a_ij - p_j} and second best
        // iterate over columns

        size_t i_j = 0;
        Scalar b_j = -__AUCTION_INF, a_ij_j = 0.,
               g_j = -__AUCTION_INF; // = max { a_ij - P_i}

        // find maximum profit i.e. j_i = arg max { a_ij - p_j} and second best
        for (size_t i = fromRowIdx; i < toRowIdx; i++) // for the i-th row
        {
            const Scalar m_i_j = m(i, j);
            if (m_i_j == 0)
                continue;

            const Scalar diff = m_i_j - profits[i];

            if (diff > b_j)
            {
                // if there already was an entry found, this is the second best
                if (results[threadId].assignmentFound)
                    g_j = b_j;

                b_j = diff;
                i_j = i;
                a_ij_j = m_i_j;

                results[threadId].assignmentFound = true;
            }
            if (diff > g_j && i_j != i)
                g_j = diff;
        }

        // best person for object j
        results[threadId].bestIndex = i_j;
        //	b_i = max { a_ij - P_i} for i in A(j)
        results[threadId].bestBid = b_j;
        //	g_i = max { a_ij - P_i} for j in A(j) and i != i_j
        results[threadId].secondBestBid = g_j;
        // a_ij_j = a(i_j, j), already found above
        results[threadId].bestMatrixValue = a_ij_j;

        return;
    }

    static Indices splitToIntervals(const size_t parts, const size_t width, const size_t offset = 0)
    {
        Indices intervals(parts, offset);
        const size_t base = width / parts;
        size_t reminder = width - (base * parts);
        size_t last = 0;
        for (size_t i = 0; i < parts; ++i)
        {
            intervals[i] = last + base;
            if (reminder > 0)
            {
                reminder--;
                intervals[i]++;
            }
            last = intervals[i];
        }
        return intervals;
    }

    static Indices
    splitToIntervalsByMatrixType(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & m,
                                 const size_t parts, const size_t index, bool major,
                                 size_t & offset)
    {
        const size_t entries = m.outerIndexPtr()[index + 1] - m.outerIndexPtr()[index];

        offset = m.outerIndexPtr()[index];

        return AuctionCommon<Scalar>::splitToIntervals(parts, entries, 0);
    }

    static Indices splitToIntervalsByMatrixType(const WeightMatrix & m, const size_t parts,
                                                const size_t index, bool major, size_t & offset)
    {
        const size_t entries = (major) ? m.rows() : m.cols();
        offset = 0;
        return AuctionCommon<Scalar>::splitToIntervals(parts, entries);
    }
    /**
     * returns true if p_j <= lambda for all unassigned objects.
     *
     * @param c locked columns
     * @param prices prices of objects
     * @param lambda bidding threshold
     * @return true if all prices of unassigned objects are below lambda,
     * otherwise false
     */
    static const bool unassignedObjectsLTlambda(const Locks & c, const Scalars & prices,
                                                const Scalar lambda)
    {
        for (size_t j = 0; j < c.size(); ++j)
            if (!c[j] && prices[j] > lambda)
                return false;
        return true;
    }

    /**
     * check if all persons are assigned
     * @return true if all persons are assigned, otherwise false
     */
    static const bool allPersonsAssigned(const Locks & r)
    {
        return std::all_of(r.begin(), r.end(), [](bool const & value) { return value; });
    }

  private:
    AuctionCommon();
    virtual ~AuctionCommon();
};

} /* namespace LSAP */
#endif /* AUCTIONCOMMON_H_ */
