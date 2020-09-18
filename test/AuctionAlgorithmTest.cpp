#include <algorithm>
#include <auction>
#include <fstream>
#include <gtest/gtest.h>
#include <tuple>

using namespace Auction;

// template <typename T = DenseEigenMatrix<double>>
// auto objectiveFunctionValue(typename Solver<T>::Edges const &edges) -> typename T::scalar_t{
//   typename T::scalar_t result = 0;
//   for (auto const &e : edges) {
//     result += e.v;
//   }
//   return result;
// }

// template <typename T>
// T objectiveFunctionValue(typename Auction<T>::Edges const &edges) {
//   T result = 0;
//   for (auto const &e : edges) {
//     result += e.v;
//   }
//   return result;
// }
template <typename Scalar>
std::tuple<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Scalar>
readDataFromFile(std::string const & filename)
{
    std::ifstream file("test_data_01");
    unsigned int rows = 0, cols = 0;
    file >> rows;
    file >> cols;
    Scalar cost = 0;
    file >> cost;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> m(rows, cols);
    for (size_t row = 0; row < rows; ++row)
        for (size_t col = 0; col < cols; ++col)
            file >> m(row, col);

    return std::make_tuple(m, cost);
}

template <typename T>
void printEdges(typename Solver<T>::Edges const & edges)
{
    for (auto const & e : edges)
    {
        std::cout << e << std::endl;
        ;
    }
}

// TODO: this example hangs infinitely if not transposed!
TEST(DISABLED_test_data, data_01)
{

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

    // minimization problem to maximization problem
    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            m(i, j) = maxCoefficient - m(i, j);

    m.transposeInPlace();
    auto solution = Auction::solve<>(m);
    // auto x = objectiveFunctionValue<>(solution);
    // // EXPECT_EQ(118.1, x);
    // std::cout << (m) << std::endl;
    printEdges<DenseEigenMatrix<double>>(solution);
    typedef decltype(solution)::value_type Edge;

    std::sort(solution.begin(), solution.end(),
              [](Edge const & e1, Edge const & e2) { return e1.row < e2.row; });
    EXPECT_EQ(0, solution[1].row);
    EXPECT_EQ(1, solution[3].row);
    EXPECT_EQ(2, solution[0].row);
    EXPECT_EQ(3, solution[2].row);
}

TEST(test_data, data_02)
{

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

    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            m(i, j) = 1. - m(i, j);

    auto solution = Auction::solve<>(m.normalized());
    typedef typename Solver<>::Edge Edge;
    std::sort(solution.begin(), solution.end(),
              [](Edge const & e1, Edge const & e2) { return e1.row < e2.row; });
    EXPECT_EQ(3, solution[0].col);
    EXPECT_EQ(2, solution[1].col);
    EXPECT_EQ(1, solution[2].col);
    EXPECT_EQ(0, solution[3].col);
}

TEST(test_data, data_03)
{

    // example from
    // http://www.cse.psu.edu/~rtc12/CSE598C/comboptBlockICM.pdf
    Eigen::Matrix<double, 5, 5> m = Eigen::Matrix<double, 5, 5>();

    // clang-format off
    m <<   
      0.95,   0.76,   0.62,   0.41,   0.06,
      0.23,   0.46,   0.79,   0.94,   0.35,
      0.61,   0.02,   0.92,   0.92,   0.81,   
      0.49,   0.82,   0.74,   0.41,   0.01,   
      0.89,   0.44,   0.18,   0.89,   0.14;
    // clang-format on

    m.normalize();

    // for ( size_t i = 0; i < m.rows(); ++i )
    //   for ( size_t j = 0; j < m.cols(); ++j )
    //       m(i, j) = 1. - m(i, j);

    auto solution = Auction::solve<>(m);
    // auto x = objectiveFunctionValue<>(solution);
    typedef typename Solver<>::Edge Edge;
    std::sort(solution.begin(), solution.end(),
              [](Edge const & e1, Edge const & e2) { return e1.row < e2.row; });
    EXPECT_EQ(0, solution[0].col);
    EXPECT_EQ(2, solution[1].col);
    EXPECT_EQ(4, solution[2].col);
    EXPECT_EQ(1, solution[3].col);
    EXPECT_EQ(3, solution[4].col);
}

// https://towardsdatascience.com/operations-research-in-r-assignment-problem-4a1f92a09ab
TEST(test_data, data_04)
{

    // example from
    // https://developers.google.com/optimization/assignment/linear_assignment
    Eigen::MatrixXd m = Eigen::Matrix<double, 3, 3>();

    // clang-format off
    m << 0.15, 0.10, 0.09, 
         0.09, 0.15, 0.10, 
         0.10, 0.12, 0.08;
    // clang-format on
    // expected result:
    // (0, 1) = .10, (1, 0) = .09, (2,2) = 0.08
    // Total cost: 265
    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            m(i, j) = 1. - m(i, j);

    auto solution = Auction::solve<>(m);
    typedef typename Solver<>::Edge Edge;
    std::sort(solution.begin(), solution.end(),
              [](Edge const & e1, Edge const & e2) { return e1.row < e2.row; });
    EXPECT_EQ(1, solution[0].col);
    EXPECT_EQ(0, solution[1].col);
    EXPECT_EQ(2, solution[2].col);
}
typedef ::testing::TestWithParam<std::string> ParametrizedMaximizationTest;

// Testcases from https://github.com/bmc/munkres/blob/master/test/test_munkres.py
TEST_P(ParametrizedMaximizationTest, result_has_expected_cost)
{
    auto const n = GetParam(); // this is why we need a fixture
    Eigen::MatrixXd m;
    double cost;
    std::tie(m, cost) = readDataFromFile<double>("test_data_01");

    auto _m = m.normalized();
    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            _m(i, j) = 1. - _m(i, j);

    auto solution = Auction::solve<>(_m);

    double sum = 0.;
    for (auto e : solution)
        sum += m(e.row, e.col);

    ASSERT_EQ(cost, sum);
}

// prefix (i.e. instantiation name), fixture name, parameters
INSTANTIATE_TEST_SUITE_P(something, ParametrizedMaximizationTest,
                        ::testing::Values("test_data_01", "test_data_02", "test_data_03"));
