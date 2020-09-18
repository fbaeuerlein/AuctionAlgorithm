/*
 * Auction.h
 *
 *  Created on: May 21, 2013
 *      Author: fb
 * @see
 * https://dspace.mit.edu/bitstream/handle/1721.1/3302/P-2159-27790376.pdf?sequence=1&isAllowed=y
 */

#ifndef AUCTION_H_
#define AUCTION_H_

#define __AUCTION_EPSILON_MULTIPLIER 1e-5 // epsilon multiplier
#define __AUCTION_INF 1e6                 // infinity for setting second best match
#include "Edge.hxx"
#include "Matrix.hxx"
#include "Locks.hxx"
#include <tuple>

//#define DEBUG_AUCTION

namespace Auction
{

template <typename T>
using container_type = std::vector<T>;

template <typename Scalar>
using Edges = container_type<Edge<Scalar>>;

/**
 * @brief
 *
 * @warning cols must be lower or equal to rows!
 * @warning c_ij [0;1]
 *
 * @tparam DenseEigenMatrix<double>
 */
template <typename MatrixType = DenseEigenMatrix<double>>
class Solver
{
  private:
    typedef typename MatrixType::scalar_t Scalar;

  public:
    using Scalars = container_type<Scalar>;
    using Indices = container_type<size_t>;
    using Edges = Edges<Scalar>;
    using Edge = Edge<Scalar>;

    Solver() = delete;

    typedef detail::Locks Locks;
    Solver(MatrixType const & matrix)
        : _matrix(matrix)
        , _locked_rows(matrix.rows())
        , _locked_cols(matrix.cols())
        , _lambda(0.)
        , _epsilon(__AUCTION_EPSILON_MULTIPLIER / matrix.cols())
        , _prices(matrix.cols(), .0)
        , _profits(matrix.rows(), 1.) // condition 3: initially set p_j >= lambda
    {
        assert(matrix.cols() <= matrix.rows());
    }

  protected:
    MatrixType _matrix;
    Locks _locked_rows, _locked_cols;
    Scalar _lambda, _epsilon;
    Scalars _prices, _profits;
    Edges _edges; 

  public:
    Edges solve()
    {
        do
        {
            // Step 1 (forward	auction cycle):
            // Execute iterations of the forward auction algorithm until at least one
            // more person becomes assigned. If there is an unassigned person left, go
            // to step 2; else go to step 3.
            while (forward());

            if (areAllPersonsAssigned())
            {
                // Step 3 (reverse auction):
                // Execute successive iterations of the reverse auction algorithm until
                // the algorithm terminates with p_j <= lambda for all unassigned objects j
                while (true)
                {
                    reverse();
                    if (unassignedObjectsLowerThanLambda())
                    {
                        break;
                    }
                }
                break;
            }
            else
            {

                // Step 2 (reverse auction cycle):
                // Execute several iterations of the reverse auction algorithm until at
                // least one more object becomes assigned or until we have p_j <=
                // lambda for all unassigned objects. If there is an unassigned person
                // left, go to step 1 else go to step 3
                while (!reverse() || !unassignedObjectsLowerThanLambda());
            }

        } while (true);

        return _edges;
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
    bool findForward(size_t const & i, Scalar & v_i, Scalar & w_i, Scalar & a_i_ji, size_t & j_i)
    {
        bool assignmentFound = false;

        for (size_t j = 0; j < _matrix.cols(); j++) // for the j-th column
        {
            Scalar const aij = _matrix(i, j);
            if (aij == 0)
            {
                continue;
            }
            Scalar const diff = aij - _prices[j];
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

    void updateEdgeRowOrAddEdge(size_t const & j_i, size_t const & i)
    {
        auto it =
            std::find_if(_edges.begin(), _edges.end(), [&](Edge const & e) { return e.col == j_i; });
        if (it != _edges.end())
        {
            auto & edge = *it;
            _locked_rows.unlock(edge.row); // unlock previous row
            edge.row = i;
        }
        else
        {
            _edges.emplace_back(i, j_i);
        }
    }

    void updateEdgeColumnOrAddEdge(size_t const & i_j, size_t const & j)
    {
        auto it =
            std::find_if(_edges.begin(), _edges.end(), [&](Edge const & e) { return e.row == i_j; });

        // if j_i was assigned to different i' to begin, remove (i', j_i) from S
        if (it != _edges.end())
        {
            auto & edge = *it;
            _locked_cols.unlock(edge.col); // unlock col i'
            edge.col = j;
        }
        else
        {
            _edges.emplace_back(i_j, j);
        }
    }
    /**
     * @brief Forward cycle of auction algorithm
     * @param a weight matrix (nxm)
     * @return true if assignment was made, false otherwise
     */
    bool forward()
    {
        bool assignmentFound = false;

        for (size_t i = 0; i < _matrix.rows(); i++) // for the i-th row/person
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
            if (!findForward(i, v_i, w_i, a_i_ji, j_i))
            {
                continue;
            }

            assignmentInThisIterationFound = false;

            Scalar const bid = a_i_ji - w_i + _epsilon;

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

                updateEdgeRowOrAddEdge(j_i, i);

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

    bool findReverse(const size_t j, Scalar & b_j, Scalar & g_j, Scalar & a_ij_j, size_t & i_j)
    {
        bool assignmentFound = false;

        for (size_t i = 0; i < _matrix.rows(); i++) // for the j-th column
        {
            Scalar const aij = _matrix(i, j);
            if (aij == 0)
            {
                continue;
            }
            Scalar const diff = aij - _profits[i];
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

    void scaleLambda()
    {
        /** standard lambda scaling **/
        size_t lowerThanLambda = 0;
        Scalar newLambda = _lambda;

        // if the number of objects k with p_k < lambda is bigger than (rows - cols)
        for (size_t k = 0; k < _matrix.cols(); k++)
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
        if (lowerThanLambda >= (_matrix.rows() - _matrix.cols()))
        {
            _lambda = newLambda;
        }
    }

    /**
     * @brief Reverse cycle of auction algorithm
     * @param a weight matrix (nxm)
     * @return true if assignment was made, false otherwise
     */
    bool reverse()
    {
        bool assignmentFound = false;

        for (size_t j = 0; j < _matrix.cols(); j++) // for the j-th column (objects)
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
            if (!findReverse(j, b_j, g_j, a_ij_j, i_j))
            {
                continue;
            }

            assignmentInThisIterationFound = false;

            // if b_j >= L + E, case 1:
            if (b_j >= (_lambda + _epsilon))
            {
                Scalar const diff = g_j - _epsilon;         // G_j - E
                Scalar const max = std::max(_lambda, diff); //  max { L, G_j - E}

                //	p_j = max { L, G_j - E}
                _prices[j] = max;

                //	P_i_j = a_i_jj - max {L, G_j - E}
                _profits[i_j] = a_ij_j - max;

                _locked_rows.lock(i_j);
                _locked_cols.lock(j);

                updateEdgeColumnOrAddEdge(i_j, j);

                assignmentInThisIterationFound = true;
            }
            else // if B_j < L + E, case 2
            {
                //	p_j = B_j - E
                _prices[j] = b_j - _epsilon;

                assignmentInThisIterationFound = false;

                scaleLambda();
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
    bool areAllPersonsAssigned() { return _locked_rows.all_locked(); }

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

/**
 * solve the assignment problem with the auction algorithm
 * use real-types as coefficients, otherwise scaling will not work properly!
 * @param a nxm weight matrix of type Scalar
 * @return vector of Edges which represent the assignments
 */
template <typename Scalar = double, int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic,
          int Options = 0, int MaxRows = Eigen::Dynamic, int MaxCols = Eigen::Dynamic>
Edges<Scalar> solve(Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> const & matrix)
{
    Solver<DenseEigenMatrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> s(matrix);
    return s.solve();
}

/**
 * solve the assignment problem with the auction algorithm
 * use real-types as coefficients, otherwise scaling will not work properly!
 * @param a nxm weight matrix of type Scalar
 * @return vector of Edges which represent the assignments
 */
template <typename Scalar = double, int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic,
          int Options = 0, int MaxRows = Eigen::Dynamic, int MaxCols = Eigen::Dynamic>
Edges<Scalar>
solve(Eigen::Transpose<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> const & matrix)
{
    Solver<DenseEigenMatrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> s(matrix);
    return s.solve();
}
} // namespace Auction
#endif /* AUCTION_H_ */

    // TODO: move to matrix implementation!
    // /**
    //  * template specific implementation for finding the best entry in row
    //  * using eigen sparse matrix with row major storage
    //  * @param a	input matrix
    //  * @param i row
    //  * @param v_i best value
    //  * @param w_i second best value
    //  * @param a_i_ji value of entry
    //  * @param j_i best entry
    //  * @return true if assignment was found, otherwise false
    //  */
    // bool findForward(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & a, const size_t i,
    //                  Scalar & v_i, Scalar & w_i, Scalar & a_i_ji, size_t & j_i)
    // {
    //     //			assert(a.IsRowMajor);
    //     assert(a.isCompressed());

    //     bool assignmentFound = false;

    //     const auto outerStart = a.outerIndexPtr()[i];
    //     const auto outerEnd = a.outerIndexPtr()[i + 1];

    //     for (auto innerIdx = outerStart; innerIdx < outerEnd; innerIdx++) // for the j-th column
    //     {
    //         const auto j = a.innerIndexPtr()[innerIdx];
    //         const Scalar m_i_j = a.valuePtr()[innerIdx];
    //         const Scalar diff = m_i_j - _prices[j];
    //         if (diff > v_i)
    //         {
    //             // if there already was an entry found, this is the second best
    //             if (assignmentFound)
    //             {
    //                 w_i = v_i;
    //             }

    //             v_i = diff;
    //             j_i = j;
    //             a_i_ji = m_i_j;
    //             assignmentFound = true;
    //         }
    //         if (diff > w_i && j_i != j)
    //         {
    //             w_i = diff;
    //         }
    //         // if no entry is bigger than v_i, check if there's still a bigger second best entry
    //     }
    //     return assignmentFound;
    // }