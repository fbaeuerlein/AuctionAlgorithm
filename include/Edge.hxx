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
        // , v(v_)
    {
    }
    Edge() = default;
    Edge(Edge const &) = default;
    Edge(Edge &&) = default;
    Edge & operator=(Edge const &) = default;
    Edge & operator=(Edge &&) = default;

    size_t x{0};
    size_t y{0};
    // Scalar v{0};
};

} // namespace Auction