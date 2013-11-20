/*
 * internal.h
 *
 *  Created on: Apr 25, 2013
 *      Author: fb
 */

#ifndef INTERNAL_H_
#define INTERNAL_H_

#define __PROJECT_LOG std::cerr
#include <vector>
#include <eigen3/Eigen/Core>
#include <boost/chrono/chrono_io.hpp>
#include <eigen3/Eigen/Sparse>
#include "Matrix/SparseMatrix.h"
#define NO_CL
#define __AUCTION_EPSILON_MULTIPLIER 1e-5 	// epsilon multiplier
#define __AUCTION_INF 1e6 					// infinity for setting second best match

typedef Eigen::Matrix<size_t, -1, -1> AssignmentMatrix;

// which type of coefficient in weight matrix
typedef double Scalar;

// represents edges/assignments
typedef std::pair<size_t, size_t> Assignment;
typedef std::vector<Assignment> Assignments;

struct Edge
{
	size_t x;
	size_t y;
	Scalar v;
};

typedef std::vector<Edge> Edges;

// locking of rows/cols
typedef std::vector<bool> LockVector;

enum APSolvingMode { MAXIMIZATION, MINIMIZATION };


#define MEASURE_DURATION_SINGLE(__MD_COMMAND__) \
{\
double __md_accumulator = 0;\
std::cout << "command  <" << #__MD_COMMAND__ << "> took "; \
boost::chrono::high_resolution_clock::time_point __md_start = boost::chrono::high_resolution_clock::now();\
__MD_COMMAND__ ;\
boost::chrono::nanoseconds __md_sec = boost::chrono::high_resolution_clock::now() - __md_start;\
__md_accumulator += __md_sec.count();\
std::cout << ( __md_accumulator / 1000. )  << " us " << std::endl;} \

#define MEASURE_DURATION_SINGLE_WITHOUT_PRINTING_COMMAND(__MD_COMMAND__) \
{\
double __md_accumulator = 0;\
boost::chrono::high_resolution_clock::time_point __md_start = boost::chrono::high_resolution_clock::now();\
__MD_COMMAND__ ;\
boost::chrono::nanoseconds __md_sec = boost::chrono::high_resolution_clock::now() - __md_start;\
__md_accumulator += __md_sec.count();\
std::cout << ( __md_accumulator / 1000. )  << " us " << std::endl;} \


#define MEASURE_DURATION_SINGLE_STORED(__MD_COMMAND__, __MD_STORE__) \
{\
double __md_accumulator = 0;\
boost::chrono::high_resolution_clock::time_point __md_start = boost::chrono::high_resolution_clock::now();\
__MD_COMMAND__ ;\
boost::chrono::nanoseconds __md_sec = boost::chrono::high_resolution_clock::now() - __md_start;\
__md_accumulator += __md_sec.count();\
__MD_STORE__ = (double) ( __md_accumulator / 1000. );} \

#endif /* INTERNAL_H_ */
