/*
 * CRSMatrix.h
 *
 *  Created on: May 8, 2013
 *      Author: fb
 */

#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

#include "../internal.h"

#define __SPARSE_MATRIX_EPSILON 0
#define LOG_SPARSEMATRIX __PROJECT_LOG << "[SparseMatrix] "
//#define DEBUG

namespace LSAP
{

template<typename Scalar = double>
class SparseMatrix
{
public:

	enum CompressionType
	{
		rowMajor = true, colMajor = false
	};

	typedef Eigen::Matrix<Scalar, -1, -1> weightMatrix;

	typedef typename std::vector<Scalar> Scalars;
	typedef std::vector<size_t> Indizes;

	SparseMatrix() :
			_rmOuterSize(0), _rmInnerSize(0), _cmOuterSize(0), _cmInnerSize(0), _NumberOfValues(
					0), _rows(0), _cols(0), _zero(0.)
	{
		_CompressionType = rowMajor;
	}

	SparseMatrix(const weightMatrix & m) :
			_rmOuterSize(0), _rmInnerSize(0), _cmOuterSize(0), _cmInnerSize(0), _NumberOfValues(
					0), _rows(0), _cols(0), _zero(0.)
	{
		_CompressionType = rowMajor; // compatibility

		constructFromEigenMatrix(m);

#ifdef DEBUG
		LOG_SPARSEMATRIX << "input matrix: " << std::endl << m << std::endl;

		LOG_SPARSEMATRIX << "constructing from Eigen matrix ..." << std::endl;

		LOG_SPARSEMATRIX << "using " << (_CompressionType ? "row" : "col") << " major representation! " << std::endl;

		LOG_SPARSEMATRIX << "values = ";
		for ( size_t i = 0; i < _NumberOfValues; ++i) __PROJECT_LOG << _rmValues[i] << " ";
		__PROJECT_LOG << "\n";

		LOG_SPARSEMATRIX << "InnerIndizes = ";
		for ( size_t i = 0; i < _NumberOfValues; ++i) __PROJECT_LOG << _rmInnerIndizes[i] << " ";
		__PROJECT_LOG << "\n";

		LOG_SPARSEMATRIX << "OuterIndizes = ";
		for ( size_t i = 0; i < _rmOuterSize; ++i) __PROJECT_LOG << _rmOuterIndizes[i] << " ";
		__PROJECT_LOG << "\n";

#endif
	}

	virtual ~SparseMatrix()
	{

	}

	/**
	 * eigen like operator for accessing the matrix-values
	 * @param r row
	 * @param c column
	 * @return value of (r,c)
	 */
	const Scalar & operator()(const size_t r, const size_t c) const
	{
		assert(r >= 0 && r < _rows && c >= 0 && c < _cols);

		// determine outer and inner index values
		const size_t outerIdx = (_CompressionType) ? r : c;
		const size_t innerIdx = (_CompressionType) ? c : r;

		// get end-index of outer indizes
		for (size_t i = _rmOuterIndizes[outerIdx];
				i < _rmOuterIndizes[outerIdx + 1]; ++i)
		{
			if (_rmInnerIndizes[i] > innerIdx)
				return _zero;
			if (_rmInnerIndizes[i] == innerIdx)
				return _rmValues[i];
		}

		return _zero;
	}

private:

	void constructFromEigenMatrix(const weightMatrix & m)
	{
		_rows = m.rows();
		_cols = m.cols();

		// count number of non-zero values
		for (size_t r = 0; r < _rows; ++r)
			for (size_t c = 0; c < _cols; ++c)
				if (abs(m(r, c) > __SPARSE_MATRIX_EPSILON))
					_NumberOfValues++;

		// determine array sizes dependend on compression type
		_rmOuterSize = _rows + 1;
		_rmInnerSize = _cols;

		// determine array sizes dependend on compression type
		_cmOuterSize = _cols + 1;
		_cmInnerSize = _rows;

		// allocate  values and indizes
		_rmOuterIndizes = Indizes(_rmOuterSize);
		_rmInnerIndizes = Indizes(_NumberOfValues);
		_rmValues = Scalars(_NumberOfValues);

		_cmOuterIndizes = Indizes(_cmOuterSize);
		_cmInnerIndizes = Indizes(_NumberOfValues);
		_cmValues = Scalars(_NumberOfValues);

		// set allocated index-array to 0
		for (auto & i : _rmOuterIndizes)
			i = 0;
		for (auto & i : _cmOuterIndizes)
			i = 0;

		// last position holds the size of values/innerindizes
		_rmOuterIndizes[_rmOuterSize - 1] = _NumberOfValues;
		_cmOuterIndizes[_cmOuterSize - 1] = _NumberOfValues;

		// read matrix data for row-major representation
		size_t valuesSoFar = 0;

		for (size_t o = 0; o < _rmOuterSize - 1; ++o)
		{
			size_t valuesPerOuter = 0;
			for (size_t i = 0; i < _rmInnerSize; ++i)
			{
				const Scalar mVal = m(o, i);
				if (mVal != 0)
				{
					_rmValues[valuesSoFar] = mVal;
					_rmInnerIndizes[valuesSoFar] = i;
					valuesPerOuter++;
					valuesSoFar++;
				}
			}
			_rmOuterIndizes[o + 1] = _rmOuterIndizes[o] + valuesPerOuter;
		}

		// read matrix data for col-major representation
		valuesSoFar = 0;

		for (size_t o = 0; o < _cmOuterSize - 1; ++o)
		{
			size_t valuesPerOuter = 0;
			for (size_t i = 0; i < _cmInnerSize; ++i)
			{
				const Scalar mVal = m(i, o);
				if (mVal != 0)
				{
					_cmValues[valuesSoFar] = mVal;
					_cmInnerIndizes[valuesSoFar] = i;
					valuesPerOuter++;
					valuesSoFar++;
				}
			}
			_cmOuterIndizes[o + 1] = _cmOuterIndizes[o] + valuesPerOuter;
		}
	}

public:

	inline const Scalar maxCoeff() const
	{
		Scalar c = 0.;
		for (size_t i = 0; i < _NumberOfValues; ++i)
			c = (_rmValues[i] > c) ? _rmValues[i] : c;
		return c;
	}

	inline const Scalar minNonZeroCoeff() const
	{
		assert(_rmValues.size() > 0);
		Scalar c = _rmValues[0];
		for (size_t i = 0; i < _NumberOfValues; ++i)
			c = (_rmValues[i] < c) ? _rmValues[i] : c;
		return c;
	}

	inline Scalars & valuesRM()
	{
		return _rmValues;
	}
	inline Indizes & outerIndizesRM()
	{
		return _rmOuterIndizes;
	}
	inline Indizes & innerIndizesRM()
	{
		return _rmInnerIndizes;
	}

	inline const size_t & outerIndexRM(const size_t o) const
	{
		return _rmOuterIndizes[o];
	}
	inline const size_t & innerIndexRM(const size_t i) const
	{
		return _rmInnerIndizes[i];
	}
	inline const Scalar & valueRM(const size_t i) const
	{
		return _rmValues[i];
	}

	inline const size_t & outerIndexCM(const size_t o) const
	{
		return _cmOuterIndizes[o];
	}
	inline const size_t & innerIndexCM(const size_t i) const
	{
		return _cmInnerIndizes[i];
	}
	inline const Scalar & valueCM(const size_t i) const
	{
		return _cmValues[i];
	}

	inline Scalars & valuesCM()
	{
		return _cmValues;
	}
	inline Indizes & outerIndizesCM()
	{
		return _cmOuterIndizes;
	}
	inline Indizes & innerIndizesCM()
	{
		return _cmInnerIndizes;
	}

	inline const size_t & rows() const
	{
		return _rows;
	}
	inline const size_t & cols() const
	{
		return _cols;
	}

	inline const size_t size() const
	{
		return _NumberOfValues;
	}

	inline const bool compressionType() const
	{
		return _CompressionType;
	}

	// compatibility to old version
	inline Scalars & values()
	{
		return _rmValues;
	}
	inline Indizes & outerIndizes()
	{
		return _rmOuterIndizes;
	}
	inline Indizes & innerIndizes()
	{
		return _rmInnerIndizes;
	}

private:
	// row-major representation
	Scalars _rmValues;
	Indizes _rmOuterIndizes; // outer index (holds rowIdx if CRS, colIdx if CCS)
	Indizes _rmInnerIndizes; // inner index
	size_t _rmOuterSize;
	size_t _rmInnerSize;

	// col-major representation
	Scalars _cmValues;
	Indizes _cmOuterIndizes; // outer index (holds rowIdx if CRS, colIdx if CCS)
	Indizes _cmInnerIndizes; // inner index
	size_t _cmOuterSize;
	size_t _cmInnerSize;

	// common
	size_t _NumberOfValues;

	size_t _rows, _cols;

	CompressionType _CompressionType;

	Scalar _zero;	// for returning reference to values

};

template<typename Scalar = double>
class SparseMatrixRM
{
public:

	typedef Eigen::Matrix<Scalar, -1, -1> weightMatrix;

	typedef typename std::vector<Scalar> Scalars;
	typedef std::vector<size_t> Indizes;

	SparseMatrixRM() :
			_rmOuterSize(0), _rmInnerSize(0), _NumberOfValues(0), _rows(0), _cols(0), _zero(0.)
	{}

	SparseMatrixRM(const weightMatrix & m) :
			_rmOuterSize(0), _rmInnerSize(0), _NumberOfValues(
					0), _rows(0), _cols(0), _zero(0.)
	{
		constructFromEigenMatrix(m);

#ifdef DEBUG
		LOG_SPARSEMATRIX << "input matrix: " << std::endl << m << std::endl;

		LOG_SPARSEMATRIX << "constructing from Eigen matrix ..." << std::endl;

		LOG_SPARSEMATRIX << "using " << (_CompressionType ? "row" : "col") << " major representation! " << std::endl;

		LOG_SPARSEMATRIX << "values = ";
		for ( size_t i = 0; i < _NumberOfValues; ++i) __PROJECT_LOG << _rmValues[i] << " ";
		__PROJECT_LOG << "\n";

		LOG_SPARSEMATRIX << "InnerIndizes = ";
		for ( size_t i = 0; i < _NumberOfValues; ++i) __PROJECT_LOG << _rmInnerIndizes[i] << " ";
		__PROJECT_LOG << "\n";

		LOG_SPARSEMATRIX << "OuterIndizes = ";
		for ( size_t i = 0; i < _rmOuterSize; ++i) __PROJECT_LOG << _rmOuterIndizes[i] << " ";
		__PROJECT_LOG << "\n";

#endif
	}

	virtual ~SparseMatrixRM()
	{

	}

	/**
	 * eigen like operator for accessing the matrix-values
	 * @param r row
	 * @param c column
	 * @return value of (r,c)
	 */
	const Scalar & operator()(const size_t r, const size_t c) const
	{
		assert(r >= 0 && r < _rows && c >= 0 && c < _cols);

		// get end-index of outer indizes
		for (size_t i = _rmOuterIndizes[r];
				i < _rmOuterIndizes[r + 1]; ++i)
		{
			if (_rmInnerIndizes[i] < c) continue;
			if (_rmInnerIndizes[i] == c)
				return _rmValues[i];
			if (_rmInnerIndizes[i] > c)
				return _zero;
		}

		return _zero;
	}

private:

	void constructFromEigenMatrix(const weightMatrix & m)
	{
		_rows = m.rows();
		_cols = m.cols();

		// count number of non-zero values
		for (size_t r = 0; r < _rows; ++r)
			for (size_t c = 0; c < _cols; ++c)
				if (abs(m(r, c) > __SPARSE_MATRIX_EPSILON))
					_NumberOfValues++;

		// determine array sizes dependend on compression type
		_rmOuterSize = _rows + 1;
		_rmInnerSize = _cols;

		// allocate  values and indizes
		_rmOuterIndizes = Indizes(_rmOuterSize);
		_rmInnerIndizes = Indizes(_NumberOfValues);
		_rmValues = Scalars(_NumberOfValues);

		// set allocated index-array to 0
		for (auto & i : _rmOuterIndizes)
			i = 0;

		// last position holds the size of values/innerindizes
		_rmOuterIndizes[_rmOuterSize - 1] = _NumberOfValues;

		// read matrix data for row-major representation
		size_t valuesSoFar = 0;

		for (size_t o = 0; o < _rmOuterSize - 1; ++o)
		{
			size_t valuesPerOuter = 0;
			for (size_t i = 0; i < _rmInnerSize; ++i)
			{
				const Scalar mVal = m(o, i);
				if (mVal != 0)
				{
					_rmValues[valuesSoFar] = mVal;
					_rmInnerIndizes[valuesSoFar] = i;
					valuesPerOuter++;
					valuesSoFar++;
				}
			}
			_rmOuterIndizes[o + 1] = _rmOuterIndizes[o] + valuesPerOuter;
		}

	}

public:

	inline const Scalar maxCoeff() const
	{
		Scalar c = 0.;
		for (size_t i = 0; i < _NumberOfValues; ++i)
			c = (_rmValues[i] > c) ? _rmValues[i] : c;
		return c;
	}

	inline const Scalar minNonZeroCoeff() const
	{
		assert(_rmValues.size() > 0);
		Scalar c = _rmValues[0];
		for (size_t i = 0; i < _NumberOfValues; ++i)
			c = (_rmValues[i] < c) ? _rmValues[i] : c;
		return c;
	}

	inline Scalars & values()
	{
		return _rmValues;
	}
	inline Indizes & outerIndizes()
	{
		return _rmOuterIndizes;
	}
	inline Indizes & innerIndizes()
	{
		return _rmInnerIndizes;
	}

	inline const size_t & outerIndex(const size_t o) const
	{
		return _rmOuterIndizes[o];
	}
	inline const size_t & innerIndex(const size_t i) const
	{
		return _rmInnerIndizes[i];
	}
	inline const Scalar & value(const size_t i) const
	{
		return _rmValues[i];
	}

	inline const size_t & rows() const
	{
		return _rows;
	}
	inline const size_t & cols() const
	{
		return _cols;
	}

	inline const size_t size() const
	{
		return _NumberOfValues;
	}

private:
	// row-major representation
	Scalars _rmValues;
	Indizes _rmOuterIndizes; // outer index (holds rowIdx if CRS, colIdx if CCS)
	Indizes _rmInnerIndizes; // inner index
	size_t _rmOuterSize;
	size_t _rmInnerSize;

	// common
	size_t _NumberOfValues;

	size_t _rows, _cols;

	Scalar _zero;	// for returning reference to values

};

template<typename Scalar = double>
class SparseMatrixCM
{
public:

	typedef Eigen::Matrix<Scalar, -1, -1> weightMatrix;

	typedef typename std::vector<Scalar> Scalars;
	typedef std::vector<size_t> Indizes;

	SparseMatrixCM() :
			_cmOuterSize(0), _cmInnerSize(0), _NumberOfValues(0), _rows(0), _cols(0), _zero(0.)
	{
	}

	SparseMatrixCM(const weightMatrix & m) :
			_cmOuterSize(0), _cmInnerSize(0), _NumberOfValues(0), _rows(0), _cols(0), _zero(0.)
	{
		constructFromEigenMatrix(m);

#ifdef DEBUG
		LOG_SPARSEMATRIX << "input matrix: " << std::endl << m << std::endl;

		LOG_SPARSEMATRIX << "constructing from Eigen matrix ..." << std::endl;

		LOG_SPARSEMATRIX << "using " << (_CompressionType ? "row" : "col") << " major representation! " << std::endl;

		LOG_SPARSEMATRIX << "values = ";
		for ( size_t i = 0; i < _NumberOfValues; ++i) __PROJECT_LOG << _rmValues[i] << " ";
		__PROJECT_LOG << "\n";

		LOG_SPARSEMATRIX << "InnerIndizes = ";
		for ( size_t i = 0; i < _NumberOfValues; ++i) __PROJECT_LOG << _rmInnerIndizes[i] << " ";
		__PROJECT_LOG << "\n";

		LOG_SPARSEMATRIX << "OuterIndizes = ";
		for ( size_t i = 0; i < _rmOuterSize; ++i) __PROJECT_LOG << _rmOuterIndizes[i] << " ";
		__PROJECT_LOG << "\n";

#endif
	}

	virtual ~SparseMatrixCM()
	{

	}

	/**
	 * eigen like operator for accessing the matrix-values
	 * @param r row
	 * @param c column
	 * @return value of (r,c)
	 */
	const Scalar & operator()(const size_t r, const size_t c) const
	{
		assert(r >= 0 && r < _rows && c >= 0 && c < _cols);

		// get end-index of outer indizes
		for (size_t i = _cmOuterIndizes[c];
				i < _cmOuterIndizes[c + 1]; ++i)
		{
			if (_cmInnerIndizes[i] < r)
				continue;
			if (_cmInnerIndizes[i] == r)
				return _cmValues[i];
			if (_cmInnerIndizes[i] > r)
				return _zero;
		}

		return _zero;
	}

private:

	void constructFromEigenMatrix(const weightMatrix & m)
	{
		_rows = m.rows();
		_cols = m.cols();

		// count number of non-zero values
		for (size_t r = 0; r < _rows; ++r)
			for (size_t c = 0; c < _cols; ++c)
				if (abs(m(r, c) > __SPARSE_MATRIX_EPSILON))
					_NumberOfValues++;

		// determine array sizes dependend on compression type
		_cmOuterSize = _cols + 1;
		_cmInnerSize = _rows;

		// allocate  values and indizes
		_cmOuterIndizes = Indizes(_cmOuterSize);
		_cmInnerIndizes = Indizes(_NumberOfValues);
		_cmValues = Scalars(_NumberOfValues);

		// set allocated index-array to 0
		for (auto & i : _cmOuterIndizes)
			i = 0;

		// last position holds the size of values/innerindizes
		_cmOuterIndizes[_cmOuterSize - 1] = _NumberOfValues;

		// read matrix data for row-major representation
		size_t valuesSoFar = 0;

		for (size_t o = 0; o < _cmOuterSize - 1; ++o)
		{
			size_t valuesPerOuter = 0;
			for (size_t i = 0; i < _cmInnerSize; ++i)
			{
				const Scalar mVal = m(i, o);
				if (mVal != 0)
				{
					_cmValues[valuesSoFar] = mVal;
					_cmInnerIndizes[valuesSoFar] = i;
					valuesPerOuter++;
					valuesSoFar++;
				}
			}
			_cmOuterIndizes[o + 1] = _cmOuterIndizes[o] + valuesPerOuter;
		}
	}

public:

	inline const Scalar maxCoeff() const
	{
		Scalar c = 0.;
		for (size_t i = 0; i < _NumberOfValues; ++i)
			c = (_cmValues[i] > c) ? _cmValues[i] : c;
		return c;
	}

	inline const Scalar minNonZeroCoeff() const
	{
		assert(_cmValues.size() > 0);
		Scalar c = _cmValues[0];
		for (size_t i = 0; i < _NumberOfValues; ++i)
			c = (_cmValues[i] < c) ? _cmValues[i] : c;
		return c;
	}

	inline const size_t & outerIndex(const size_t o) const
	{
		return _cmOuterIndizes[o];
	}
	inline const size_t & innerIndex(const size_t i) const
	{
		return _cmInnerIndizes[i];
	}
	inline const Scalar & value(const size_t i) const
	{
		return _cmValues[i];
	}

	inline Scalars & values()
	{
		return _cmValues;
	}
	inline Indizes & outerIndizes()
	{
		return _cmOuterIndizes;
	}
	inline Indizes & innerIndizes()
	{
		return _cmInnerIndizes;
	}

	inline const size_t & rows() const
	{
		return _rows;
	}
	inline const size_t & cols() const
	{
		return _cols;
	}

	inline const size_t size() const
	{
		return _NumberOfValues;
	}

private:
	// col-major representation
	Scalars _cmValues;
	Indizes _cmOuterIndizes; // outer index (holds rowIdx if CRS, colIdx if CCS)
	Indizes _cmInnerIndizes; // inner index
	size_t _cmOuterSize;
	size_t _cmInnerSize;

	// common
	size_t _NumberOfValues;

	size_t _rows, _cols;

	Scalar _zero;	// for returning reference to values

};

}
#endif /* SPARSEMATRIX_H_ */
