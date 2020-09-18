#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

namespace Auction
{
template <typename Scalar>
struct Edge
{
    Edge(size_t const & x_, size_t const & y_)
        : row(x_)
        , col(y_)
    {
    }
    Edge() = default;
    Edge(Edge const &) = default;
    Edge(Edge &&) = default;
    Edge & operator=(Edge const &) = default;
    Edge & operator=(Edge &&) = default;

    /**
     * @brief transpose the current edge (i.e. swap x and y)
     *
     */
    void transpose() { std::swap(row, col); }

    /**
     * @brief return an Edge that is the transposed one of this
     *
     * @return Edge transposed edge
     */
    Edge transposed() const { return Edge(col, row); }

    /**
     * @brief stream Edge to ostream
     *
     * @param os stream
     * @return std::ostream&
     */
    template<typename T>
    friend std::ostream & operator<<(std::ostream & os, Edge<T> const & e);

    size_t row{0};
    size_t col{0};
};

template<typename Scalar>
std::ostream & operator<<(std::ostream & os, Edge<Scalar> const & e)
{
    os << "(" << e.row << ", " << e.col << ")";
    return os;
}

} // namespace Auction