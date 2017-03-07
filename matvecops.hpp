#ifndef MATVECOPS_HPP
#define MATVECOPS_HPP

#include <vector>

/* Function that gives a dot product of two vectors  */
double dot_pro(std::vector<double> x,
	       std::vector<double> y);

/* Function that computes the norm of a vector  */
double vec_norm(std::vector<double> x);

/* Function that add two vectors  */
std::vector<double> vec_add(std::vector<double> x,
			    std::vector<double> y);

/* Function that mutliplies a vector by a scalar  */
std::vector<double> vec_sca(std::vector<double> x,
			    double c);

/* Function that subtracts one vector from another  */
std::vector<double> vec_sub(std::vector<double> x,
			    std::vector<double> y);

/* Function that multiplies sparse matrix of the CSR form to a vector  */
std::vector<double> CSR_mat_vec(std::vector<double> val,
				std::vector<int>    row_ptr,
				std::vector<int>    col_idx,
				std::vector<double> x);

#endif /* MATVECOPS_HPP */
