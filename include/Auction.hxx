/*
 * Auction.h
 *
 *  Created on: May 21, 2013
 *      Author: fb
 */

#ifndef AUCTION_H_
#define AUCTION_H_

#include "../internal.h"
#include "AuctionCommon.hxx"
#include <cassert>

//#define DEBUG_AUCTION

namespace Auction
{

template <typename Scalar = double, typename MatrixType = Eigen::Matrix<Scalar, -1, -1>>
class Solver
{
  public:
    typedef typename AuctionCommon<Scalar>::Scalars Scalars;
    typedef typename AuctionCommon<Scalar>::Indices Indices;
    typedef typename AuctionCommon<Scalar>::Edge Edge;
    typedef typename AuctionCommon<Scalar>::Edges Edges;

    class Locks
    {
        public:
        Locks(size_t const & size) 
            : _locks(size, false)
        {}

        Locks() = delete;

        /**
         * @brief locks the i-th flag
         * 
         * @param index index of the flag
         */
        void lock(size_t const & index ) 
        { 
            assert(index < _locks.size());
            _locks[index] = true; 
            _locked++;
        }

        /**
         * @brief unlocks the i-th flag
         * 
         * @param index index of the flag
         */
        void unlock(size_t const & index ) 
        { 
            assert(index < _locks.size());
            _locks[index] = false; 
            _locked--;
        }

        /**
         * @brief indicates that the i-th flag is set locked
         * 
         * @param index index = i-th flag
         * @return true if flag is set locked
         * @return false otherwise
         */
        bool is_locked(size_t const & index) 
        {
            assert(index < _locks.size()); 
            return _locks[index]; 
        }

        /**
         * @brief indicates that all flags are locked
         * 
         * @return true if all flags are locked
         * @return false otherwise
         */
        bool all_locked() const { return _locked == _locks.size(); }

        /**
         * @brief returns the number of the lock flags
         * 
         * @return size_t number of the flags
         */
        size_t size() const { return _locks.size(); }
        private:
        std::vector<bool> _locks;   ///< storage of the lock flags
        size_t _locked{0};          ///< number of locked entries
    };

  private:
    Solver() = delete;

  protected:
    Solver(size_t const & rows, size_t const & cols)
        : _rows(rows)
        , _cols(cols)
        , _locked_rows(rows)
        , _locked_cols(cols)
        , _lambda(0.)
        , _epsilon(__AUCTION_EPSILON_MULTIPLIER / cols)
        , _prices(cols, .0)
        , _profits(rows, 1.) // condition 3: initially set p_j >= lambda
    {
    }

    size_t _rows, _cols;
    Locks _locked_rows, _locked_cols;
    Scalar _lambda, _epsilon;
    Scalars _prices, _profits;
    Edges _edges; // refactor to use vector index as row/col index



  public:
    /**
     * solve the assignment problem with the auction algorithm
     * use real-types as coefficients, otherwise scaling will not work properly!
     * @param a nxm weight matrix of type Scalar
     * @return vector of Edges which represent the assignments
     */
    static const Edges solve(const MatrixType & a)
    {
        Solver s(a.rows(), a.cols());

        do
        {
            //		Step 1 (forward	auction cycle):
            //		Execute iterations of the forward auction algorithm until at least one
            //		more person becomes assigned. If there is an unassigned person left, go
            //		to step 2; else go to step 3.
            while (s.forward(a))
                ;

            if (s.areAllPersonsAssigned())
            {
                //		Step 3 (reverse auction):
                //		Execute successive iterations of the reverse auction algorithm until
                // the 		algorithm terminates with p_j <= lambda for all unassigned objects j
                while (true)
                {
                    s.reverse(a);
                    if (s.unassignedObjectsLowerThanLambda())
                    {
                        break;
                    }
                }
                break;
            }
            else
            {

                //		Step 2 (reverse auction cycle):
                //		Execute several iterations of the reverse auction algorithm until at
                // least 		one more object becomes assigned or until we have p_j <=
                // lambda for all 		unassigned objects. If there is an unassigned person
                // left, go to step 1 else go to step 3
                while (!s.reverse(a) || !s.unassignedObjectsLowerThanLambda())
                    ; // reverse auction
            }

        } while (true);

        return s._edges;
    }

  private:
    /**
     * template specific implementation for finding the best entry in row
     * used for dense matrix
     * @param a	input matrix
     * @param i row
     * @param v_i best value
     * @param w_i second best value
     * @param a_i_ji value of entry
     * @param j_i best entry
     * @return true if assignment was found, otherwise false
     */
    bool findForward(const Eigen::Matrix<Scalar, -1, -1> & a, const size_t i, Scalar & v_i,
                     Scalar & w_i, Scalar & a_i_ji, size_t & j_i)
    {
        const size_t cols = a.cols();

        bool assignmentFound = false;

        for (size_t j = 0; j < cols; j++) // for the j-th column
        {
            const Scalar aij = a(i, j);
            if (aij == 0)
            {
                continue;
            }
            const Scalar diff = aij - _prices[j];
            if (diff > v_i)
            {
                // if there already was an entry found, this is the second best
                if (assignmentFound)
                {
                    w_i = v_i;
                }

                v_i = diff;
                j_i = j;
                a_i_ji = aij;
                assignmentFound = true;
            }
            if (diff > w_i && j_i != j)
            {
                w_i = diff;
            }
            // if no entry is bigger than v_i, check if there's still a bigger second best entry
        }
        return assignmentFound;
    }

    /**
     * template specific implementation for finding the best entry in row
     * using eigen sparse matrix with row major storage
     * @param a	input matrix
     * @param i row
     * @param v_i best value
     * @param w_i second best value
     * @param a_i_ji value of entry
     * @param j_i best entry
     * @return true if assignment was found, otherwise false
     */
    bool findForward(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & a, const size_t i,
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
            const Scalar diff = m_i_j - _prices[j];
            if (diff > v_i)
            {
                // if there already was an entry found, this is the second best
                if (assignmentFound)
                {
                    w_i = v_i;
                }

                v_i = diff;
                j_i = j;
                a_i_ji = m_i_j;
                assignmentFound = true;
            }
            if (diff > w_i && j_i != j)
            {
                w_i = diff;
            }
            // if no entry is bigger than v_i, check if there's still a bigger second best entry
        }
        return assignmentFound;
    }

    /**
     * @brief Forward cycle of auction algorithm
     * @param a weight matrix (nxm)
     * @return true if assignment was made, false otherwise
     */
    bool forward(const MatrixType & a)
    {

        const size_t rows = a.rows();
        bool assignmentFound = false;

        for (size_t i = 0; i < rows; i++) // for the i-th row/person
        {
            // person already assigned?
            if (_locked_rows.is_locked(i))
            {
                continue;
            }

            bool assignmentInThisIterationFound = false;

            // find an unassigned person i, its best object j_i
            // j_i = argmax {a_ij - p_j} for j in A(i) ( A(i) are the set of edges of the i-th row )
            // if a(i,j) = 0. it is not a valid edge
            size_t j_i = 0;

            //	v_i = max { a_ij - p_j} for j in A(i)				// maximum profit
            // for person
            // i
            //  v_i was already found = v_i
            //	w_i = max { a_ij - p_j} for j in A(i) and j != j_i  // second best profit
            //	if j_i is the only entry in A(i), w_i = - inf       // there's no second best profit
            Scalar w_i = -__AUCTION_INF, v_i = -__AUCTION_INF, a_i_ji = 0.; // = max { a_ij - p_j}

            // find maximum profit i.e. j_i = arg max { a_ij - p_j} and second best
            // no possible assignment found?
            if (!findForward(a, i, v_i, w_i, a_i_ji, j_i))
            {
                continue;
            }

            assignmentInThisIterationFound = false;

            const Scalar bid = a_i_ji - w_i + _epsilon;

            //	P_i = w_i - E
            _profits[i] = w_i - _epsilon; // set new profit for person

            //	prices(j_i) = max(lambda, a(i,j_i) - w(i) + epsilon)
            // if lambda <= a_ij - w_i + E, add (i, j_i) to S
            if (_lambda <= bid)
            {
                _prices[j_i] = bid;
                // assignment was made, so lock row and col
                _locked_rows.lock(i);
                _locked_cols.lock(j_i);

                bool newEdge = true;

                // if j_i was assigned to different i' to begin, remove (i', j_i) from S
                for (auto & e : _edges)
                {
                    if (e.y == j_i) // change edge
                    {
                        _locked_rows.unlock(e.x); // unlock row i'
                        newEdge = false;
                        e.x = i;
                        e.v = a_i_ji;
                        break;
                    }
                }

                if (newEdge)
                {
                    _edges.emplace_back(i, j_i, a_i_ji);
                }
                assignmentInThisIterationFound = true;
            }
            else
            {
                _prices[j_i] = _lambda;
                assignmentInThisIterationFound = false;
            }
            assignmentFound = assignmentInThisIterationFound;
        }

        return assignmentFound;
    }

    bool findReverse(const MatrixType & a, const size_t j, Scalar & b_j, Scalar & g_j,
                     Scalar & a_ij_j, size_t & i_j)
    {
        const size_t rows = a.rows();
        bool assignmentFound = false;

        for (size_t i = 0; i < rows; i++) // for the j-th column
        {
            const Scalar aij = a.coeff(i, j);
            if (aij == 0)
            {
                continue;
            }
            const Scalar diff = aij - _profits[i];
            if (diff > b_j)
            {
                // if there already was an entry found, this is the second best
                if (assignmentFound)
                {
                    g_j = b_j;
                }

                b_j = diff;
                i_j = i;
                a_ij_j = aij;
                assignmentFound = true;
            }
            if (diff > g_j && i_j != i)
            {
                g_j = diff;
            }
        }
        return assignmentFound;
    }

    /**
     * @brief Reverse cycle of auction algorithm
     * @param a weight matrix (nxm)
     * @return true if assignment was made, false otherwise
     */
    bool reverse(const MatrixType & a)
    {
        const size_t rows = a.rows();
        const size_t cols = a.cols();

        bool assignmentFound = false;

        for (size_t j = 0; j < cols; j++) // for the j-th column (objects)
        {
            bool assignmentInThisIterationFound = false;

            // object already assigned,  p_j > lambda ?
            if (_locked_cols.is_locked(j))
            {
                continue;
            }

            if (!(_prices[j] > _lambda))
            {
                continue;
            }

            // Find an unassigned object j with p_j > lambda, its best person i_j
            // i_j = argmax {a_ij - profits[i]) für i aus B(j) (PI !!!)
            size_t i_j = 0;

            // g_j = max {a_ij - P_i} for i in B(j) and i != i_j
            // if j_i is the only entry in B(j), g_j = - inf ( g_j < b_j)
            // b_j = max {a_ij - P_i} for i in B(j)
            Scalar b_j = -__AUCTION_INF, g_j = -__AUCTION_INF, a_ij_j = 0.;

            // find maximum profit i.e. j_i = arg max { a_ij - p_j} and second best
            // no assignment found
            if (!findReverse(a, j, b_j, g_j, a_ij_j, i_j))
            {
                continue;
            }

            assignmentInThisIterationFound = false;

            // if b_j >= L + E, case 1:
            if (b_j >= (_lambda + _epsilon))
            {
                const Scalar diff = g_j - _epsilon;         // G_j - E
                const Scalar max = std::max(_lambda, diff); //  max { L, G_j - E}

                //	p_j = max { L, G_j - E}
                _prices[j] = max;

                //	P_i_j = a_i_jj - max {L, G_j - E}
                _profits[i_j] = a_ij_j - max;

                _locked_rows.lock(i_j);
                _locked_cols.lock(j);

                bool newEdge = true;

                // if j_i was assigned to different i' to begin, remove (i', j_i) from S
                for (auto & e : _edges)
                {
                    if (e.x == i_j) // change edge
                    {
                        _locked_cols.unlock(e.y); // unlock col i'
                        newEdge = false;
                        e.y = j;
                        e.v = a_ij_j;

                        break;
                    }
                }
                if (newEdge)
                {
                    _edges.emplace_back(i_j, j, a_ij_j);
                }
                assignmentInThisIterationFound = true;
            }
            else // if B_j < L + E, case 2
            {
                //	p_j = B_j - E
                _prices[j] = b_j - _epsilon;

                /** standard lambda scaling **/
                size_t lowerThanLambda = 0;
                Scalar newLambda = _lambda;

                // if the number of objects k with p_k < lambda is bigger than (rows - cols)
                for (size_t k = 0; k < cols; k++)
                {
                    if (_prices[k] < _lambda) // p_k < lambda
                    {
                        lowerThanLambda++;
                        if (_prices[k] < newLambda)
                        {
                            newLambda = _prices[k];
                        }
                    }
                }
                if (lowerThanLambda >= (cols - rows))
                {
                    _lambda = newLambda;
                }
                assignmentInThisIterationFound = false;
            }
            assignmentFound = assignmentInThisIterationFound;
        }
        return assignmentFound;
    }

  protected: // enable testing of private methods via inheritance
    /**
     * check if all persons are assigned
     * @return true if all persons are assigned, otherwise false
     */
    bool areAllPersonsAssigned()
    {
        return _locked_rows.all_locked();
    }

    /**
     * returns true if p_j <= lambda for all unassigned objects.
     *
     * @return true if all prices of unassigned objects are below or equal lambda,
     * otherwise false
     */
    bool unassignedObjectsLowerThanLambda()
    {
        for (size_t j = 0; j < _locked_cols.size(); ++j)
        {
            if (!_locked_cols.is_locked(j) && _prices[j] > _lambda)
            {
                return false;
            }
        }
        return true;
    }
};

} // namespace Auction
#endif /* AUCTION_H_ */
