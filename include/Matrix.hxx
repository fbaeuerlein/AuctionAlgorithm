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
template <typename Scalar, int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic, int Options = 0,
          int MaxRows = Eigen::Dynamic, int MaxCols = Eigen::Dynamic>
class DenseEigenMatrix
{
  public:
    typedef Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> matrix_t;

    // use template to also support transpose return types
    template<typename T>
    DenseEigenMatrix(T const & matrix) noexcept
        : _matrix(matrix)
    {
    }

    template<typename T>
    DenseEigenMatrix(T && matrix) noexcept
        : _matrix(std::move(matrix))
    {
    }

    Scalar & operator()(size_t const & row, size_t const & column) noexcept
    {
        return _matrix(row, column);
    }

    Scalar const & operator()(size_t const & row, size_t const & column) const noexcept
    {
        return _matrix(row, column);
    }

    size_t rows() const noexcept { return _matrix.rows(); }
    size_t cols() const noexcept { return _matrix.cols(); }

  private:
    matrix_t _matrix;   ///< matrix type
};

} // namespace Auction