
#include "AuctionAlgorithm/Auction.h"
#include "AuctionAlgorithm/AuctionMT.h"
#include "Matrix/SparseMatrix.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
using namespace LSAP;

typedef double Scalar;

typedef Eigen::Matrix<Scalar, -1, -1> WeightMatrix;

int main(int argc, char **argv)
{
	const size_t rows = 2000, cols = 5000;

	// assert that rows <= cols and coefficients are between 0 and 1!
	Eigen::MatrixXd m = Eigen::MatrixXd::Random(rows, cols);

	// shift coefficients to get positive values
	for ( size_t i = 0; i < rows; ++i)
		for ( size_t j = 0; j < cols; ++j )
			m(i, j) += 1.;

	m /= m.maxCoeff(); // normalize to 0..1

	// create sparse matrix from dense
	// this could take some time!
	// currently there's an own sparse class, might be changed to eigen sparse type
	// this matrix type stores the matrix in rowMajor AND colMajor format, so it uses
	// a lot of space! (this was due to testing purposes)
	SparseMatrix<double> s(m);

	// result type
	Edges solution;

	// single threaded computation and some time measurement
	MEASURE_DURATION_SINGLE((solution = Auction<double>::solve(m)));
	MEASURE_DURATION_SINGLE((solution = Auction<double, SparseMatrix<double> >::solve(s)));

	// multi threaded computation (2 threads) and time measurement
	MEASURE_DURATION_SINGLE((solution = AuctionMT<double>::solve(m, 2)));
	MEASURE_DURATION_SINGLE((solution = AuctionMT<double, SparseMatrix<double> >::solve(s, 2)));

}
