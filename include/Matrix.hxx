/**
 * @file Matrix.hxx
 * @author Florian Baeuerlein (fbaeuerlein@gmail.com)
 * @brief Wrapper for different matrix types to simplify usage in algorithm
 * 
 * Uses static compile time polymorphism to avoid virtual function calls
 * 
 * @todo Use concepts in future to check interface
 * 
 * @date 2020-09-13
 */
#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace Auction
{
// use compile time polymorphism to have an abstract matrix
// type to make the concrete matrix type exchangable
// TODO: switch to concepts if there's time to

/**
 * @brief Wrapper type for dense Eigen matrices
 * 
 * @tparam Scalar 
 * @tparam Eigen::Dynamic 
 * @tparam Eigen::Dynamic 
 * @tparam Options
 * @tparam Eigen::Dynamic 
 * @tparam Eigen::Dynamic 
 */
template <typename Scalar = double, int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic, int Options = 0,
          int MaxRows = Eigen::Dynamic, int MaxCols = Eigen::Dynamic>
class DenseEigenMatrix
{
  public:
    typedef Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> matrix_t;
    typedef Scalar scalar_t;

    // use template to also support transpose return types

    template<typename MatrixType>
    DenseEigenMatrix(MatrixType && matrix) noexcept
        : _matrix(std::forward<MatrixType>(matrix))
    {
    }

    /**
     * @brief Returns the coefficient at row and column of the matrix
     * 
     * @param row row of the coefficient
     * @param column columnt of the coefficient
     * @return Scalar const& coefficient at (row, column)
     */
    Scalar const & operator()(size_t const & row, size_t const & column) const noexcept
    {
        return _matrix(row, column);
    }

    /**
     * @brief Returns the number of rows of the matrix
     * 
     * @return size_t number of rows
     */
    size_t rows() const noexcept { return _matrix.rows(); }

    /**
     * @brief returns the number of columns of the matrix
     * 
     * @return size_t number of cols
     */
    size_t cols() const noexcept { return _matrix.cols(); }

  private:
    matrix_t const _matrix;   ///< matrix type
};

} // namespace Auction