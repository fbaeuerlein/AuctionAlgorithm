
#include "AuctionAlgorithm/Auction.h"
#include "AuctionAlgorithm/AuctionMT.h"

#include <random>
using namespace LSAP;

typedef double Scalar;

typedef Eigen::Matrix<Scalar, -1, -1> WeightMatrix;

const Scalar objFuncValue(const typename Auction<Scalar>::Edges & edges )
{
	Scalar val = 0.;
	for ( auto & e : edges)
		val += e.v;

	return val;
}

int main(int argc, char **argv)
{
	const size_t rows = 5000, cols = 5000;

	// assert that rows <= cols and coefficients are between 0 and 1!
	Eigen::MatrixXd m = Eigen::MatrixXd::Random(rows, cols);

	// shift coefficients to get positive values
	for ( size_t i = 0; i < rows; ++i)
		for ( size_t j = 0; j < cols; ++j )
			m(i, j) += 1.;

	m /= m.maxCoeff(); // normalize to 0..1

	// create sparse matrix from dense
	Eigen::SparseMatrix<double, Eigen::RowMajor> s(rows, cols);
	s = m.sparseView(); // fill sparse from dense
	s.makeCompressed(); // be sure to compress!

	// result type
	Auction<double>::Edges solution;

	// single threaded computation and some time measurement
	MEASURE_DURATION_SINGLE((solution = Auction<double>::solve(m)));
	std::cout << "objective function value: " << objFuncValue(solution) << std::endl;

//	for ( auto & e: solution )
//		std::cout << "(" << e.x << ", " << e.y << ") ";
//	std::cout << std::endl;

	MEASURE_DURATION_SINGLE((solution = Auction<double, Eigen::SparseMatrix<double, Eigen::RowMajor> >::solve(s)));
	std::cout << "objective function value: " << objFuncValue(solution) << std::endl;

//
//	// multi threaded computation (2 threads) and time measurement
	MEASURE_DURATION_SINGLE((solution = AuctionMT<double>::solve(m, 2)));
	std::cout << "objective function value: " << objFuncValue(solution) << std::endl;

	MEASURE_DURATION_SINGLE((solution = AuctionMT<double, Eigen::SparseMatrix<double, Eigen::RowMajor> >::solve(s, 2)));
	std::cout << "objective function value: " << objFuncValue(solution) << std::endl;


}
