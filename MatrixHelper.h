/*
 * MatrixHelper.h
 *
 *  Created on: Jun 30, 2013
 *      Author: fb
 */

#ifndef MATRIXHELPER_H_
#define MATRIXHELPER_H_

#include "internal.h"
namespace LSAP
{

template<typename Scalar = double>
class MatrixHelper
{
private:
	MatrixHelper() {}
	virtual ~MatrixHelper() {}
public:

	typedef typename Eigen::Matrix<Scalar, -1, -1> WeightMatrix;

	/**
	 * calculate objective function value (sum over x_ij * c_ij)
	 * @param w weight matrix
	 * @param s assignment matrix
	 * @return value of objective function
	 */
	static const Scalar objFuncValue(const WeightMatrix & w, const AssignmentMatrix & s) {
		Scalar ret = 0.;
		for (unsigned int i = 0; i < w.rows(); ++i)
			for (unsigned int j = 0; j < w.cols(); ++j)
				ret += w(i, j) * s(i, j);
		return ret;
	}

	static const bool validAssignment(const WeightMatrix & w, const AssignmentMatrix & s)
	{
		for (unsigned int i = 0; i < w.rows(); ++i)
			for (unsigned int j = 0; j < w.cols(); ++j)
				if ( w(i, j) == 0. &&  s(i, j) == 1 ) return false;
		return true;
	}

	/**
	 * calculate objective function value (sum over x_ij * c_ij)
	 * @param w weight matrix
	 * @param e array of edges (assumed to be of capacity rows)
	 * @return value of objective function
	 */
	static const Scalar objFuncValue(const WeightMatrix & w, Edge * e) {
		Scalar ret = 0.;
		for (unsigned int i = 0; i < w.rows(); ++i)
				ret += w(e[i].x, e[i].y);
		return ret;
	}

	/**
	 * calculate objective function value (sum over x_ij * c_ij)
	 * @param w weight matrix
	 * @param e vector of edges
	 * @return value of objective function
	 */
	static const Scalar objFuncValue(const WeightMatrix & w, const Edges & e) {
		Scalar ret = 0.;
		for (unsigned int i = 0; i < e.size(); ++i)
			ret += w(e[i].x, e[i].y);

		return ret;
	}

	/**
	 * calculate objective function value (sum over x_ij * c_ij)
	 * @param e vector of edges with contained values
	 * @return value of objective function
	 */
	static const Scalar objFuncValue(const Edges & edges) {
		Scalar ret = 0.;
		for (auto & e : edges) ret += e.v;

		return ret;
	}
};

} /* namespace LSAP */
#endif /* MATRIXHELPER_H_ */
