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
        : x(x_)
        , y(y_)
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
    void transpose() { std::swap(x, y); }

    /**
     * @brief return an Edge that is the transposed one of this
     *
     * @return Edge transposed edge
     */
    Edge transposed() const { return Edge(y, x); }

    /**
     * @brief stream Edge to ostream
     *
     * @param os stream
     * @return std::ostream&
     */
    template<typename T>
    friend std::ostream & operator<<(std::ostream & os, Edge<T> const & e);

    size_t x{0};
    size_t y{0};
};

template<typename Scalar>
std::ostream & operator<<(std::ostream & os, Edge<Scalar> const & e)
{
    os << "(" << e.x << ", " << e.y << ")";
    return os;
}

} // namespace Auction