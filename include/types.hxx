#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

template <typename Scalar>
struct Edge
{
    Edge(size_t x_, size_t y_, Scalar v_)
        : x(x_)
        , y(y_)
        , v(v_)
    {
    }
    Edge() = default;
    Edge(Edge const &) = default;
    Edge(Edge &&) = default;
    Edge & operator=(Edge const &) = default;
    Edge & operator=(Edge &&) = default;

    size_t x{0};
    size_t y{0};
    Scalar v{0};
};

template <typename Scalar>
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

    virtual ~BidResult() = default;

    // set by thread
    size_t index{0};                     // index of examined row or column
    size_t bestIndex{0};                 // index of best entry
    Scalar bestBid{0}, secondBestBid{0}; // corresponding bids (w_i, v_i)
    Scalar bestMatrixValue{0};           // a_ij
    bool assignmentFound{false};         // assignment found ?
};