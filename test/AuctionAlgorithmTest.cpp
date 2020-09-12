#include <Auction.h>
#include <gtest/gtest.h>
#include <algorithm>

using namespace LSAP;

template <typename T>
T objectiveFunctionValue(typename Auction<T>::Edges const &edges) {
  T result = 0;
  for (auto const &e : edges) {
    result += e.v;
  }
  return result;
}

// template <typename T>
// T objectiveFunctionValue(typename Auction<T>::Edges const &edges) {
//   T result = 0;
//   for (auto const &e : edges) {
//     result += e.v;
//   }
//   return result;
// }

template <typename T> void printEdges(typename Auction<T>::Edges const &edges) {
  for (auto const &e : edges) {
    std::cout << "(" << e.x << ", " << e.y << ") = " << e.v << std::endl;
    ;
  }
}



// TODO: this example hangs infinitely if not transposed!
TEST(DISABLED_test_data, data_01) {

  // example from
  // https://documentation.sas.com/?docsetId=ornoaug&docsetTarget=ornoaug_optnet_examples04.htm&docsetVersion=14.3&locale=zh-CN
  Eigen::MatrixXd m = Eigen::Matrix<double, 5, 4>();

  // clang-format off
    m <<    1000, 36.7, 28.3, 36.1, 
         34.6,    1000,    1000, 26.2, 
         31.3,    1000, 27.1,    1000, 
         28.6,    1000, 29.1,    1000, 
         32.9,    1000, 26.6,    1000;
  // clang-format on

  // expected result - minimize weights
  // (0, 1) = 36.7, (1, 3) = 26.2, (2, 0) = 28.6, (3, 2) = 26.6

  auto const maxCoefficient = m.maxCoeff();

  auto const solution = Auction<double>::solve(m.transpose());
  auto x = objectiveFunctionValue<double>(solution);
  EXPECT_EQ(118.1, x);
  std::cout << (m) << std::endl;
  printEdges<double>(solution);
}

TEST(test_data, data_02) {

  // example from
  // https://developers.google.com/optimization/assignment/linear_assignment
  Eigen::MatrixXd m = Eigen::Matrix<double, 4, 4>();

  // clang-format off
    m <<   
        90, 76, 75, 70,
        35, 85, 55, 65,
        125, 95, 90, 105,
        45, 110, 95, 115;
  // clang-format on
  // expected result:
  // (0, 3) = 70, (1, 2) = 55, (2, 1) = 95, (3, 0) = 45
  // Total cost: 265
  m.normalize();

  for ( size_t i = 0; i < m.rows(); ++i )
    for ( size_t j = 0; j < m.cols(); ++j )
        m(i, j) = 1. - m(i, j);

  auto solution = Auction<double>::solve(m.normalized());
  auto x = objectiveFunctionValue<double>(solution);
  typedef typename Auction<double>::Edge Edge;
  std::sort(solution.begin(), solution.end(), []( Edge const & e1, Edge const & e2 ) { return e1.x < e2.x; });
  EXPECT_EQ(3, solution[0].y);
  EXPECT_EQ(2, solution[1].y);
  EXPECT_EQ(1, solution[2].y);
  EXPECT_EQ(0, solution[3].y);

}
