
#include <auction>

#include <random>
using namespace Auction;

typedef double Scalar;

typedef Eigen::Matrix<Scalar, -1, -1> WeightMatrix;

const Scalar objFuncValue(const typename Solver<Scalar>::Edges & edges )
{
	Scalar val = 0.;
	for ( auto & e : edges)
		val += e.v;

	return val;
}

int main(int argc, char **argv)
{
	const size_t rows = 1000, cols = 1000;

	// assert that rows <= cols and coefficients are between 0 and 1!
	Eigen::MatrixXd m = Eigen::MatrixXd::Random(rows, cols);

	// shift coefficients to get positive values
	for ( size_t i = 0; i < rows; ++i)
		for ( size_t j = 0; j < cols; ++j )
			m(i, j) = std::abs(m(i, j));

	m /= m.maxCoeff(); // normalize to 0..1

	// std::cout << m << std::endl;

	// create sparse matrix from dense
	Eigen::SparseMatrix<double, Eigen::RowMajor> s(rows, cols);
	s = m.sparseView(); // fill sparse from dense
	s.makeCompressed(); // be sure to compress!

	// single threaded computation and some time measurement
	// TODO check if implicit conversion is not possible!
	auto solution = Solver<double>::solve(m);
	std::cout << "objective function value: " << objFuncValue(solution) << std::endl;

	// for ( auto & e: solution )
	// 	std::cout << "(" << e.x << ", " << e.y << ") ";
	// std::cout << std::endl;


}
