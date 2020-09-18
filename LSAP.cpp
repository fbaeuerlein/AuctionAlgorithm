
#include <auction>

#include <random>

template<typename Scalar>
const Scalar objFuncValue(const typename Auction::Edges<Scalar> & edges, Eigen::MatrixXd const & m )
{
	Scalar val = 0.;
	for ( auto & e : edges)
		val += m(e.row, e.col);

	return val;
}

int main(int argc, char **argv)
{
	const size_t rows = 100, cols = 100;

	// assert that rows <= cols and coefficients are between 0 and 1!
	Eigen::MatrixXd m = Eigen::MatrixXd::Random(rows, cols);

	// shift coefficients to get positive values
	for ( size_t i = 0; i < rows; ++i)
		for ( size_t j = 0; j < cols; ++j )
			m(i, j) = std::abs(m(i, j));

	m.normalize();
	// std::cout << m << std::endl;

	// create sparse matrix from dense
	Eigen::SparseMatrix<double, Eigen::RowMajor> s(rows, cols);
	s = m.sparseView(); // fill sparse from dense
	s.makeCompressed(); // be sure to compress!

	// single threaded computation and some time measurement
	// TODO check if implicit conversion is not possible!
	auto solution = Auction::solve<>(m);
	std::cout << "objective function value: " << objFuncValue(solution, m) << std::endl;

	for ( auto & e: solution )
		std::cout << e << " ";
	std::cout << std::endl;


}
