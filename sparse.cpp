#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "COO2CSR.hpp"
#include "matvecops.hpp"
#include "sparse.hpp"

/* Modifying the matrix dimensions */
void SparseMatrix::Resize(int nrows, int ncols)
{
  this->nrows = nrows;
  this->ncols = ncols;
}

/* Method to add entru to matrix in COO format */
void SparseMatrix::AddEntry(int i, int j, double val)
{
  i_idx.push_back(i);
  j_idx.push_back(j);
  a.push_back(val);
}


/* Method to convert COO matrix to CSR format using provided function */
void SparseMatrix::ConvertToCSR()
{
  COO2CSR(a, i_idx, j_idx);

}

/* Method to perform sparse matrix vector multiplication using CSR formatter matrix */
std::vector<double> SparseMatrix::MulVec(std::vector<double> &vec)
{
  std::vector<double> vec1 = CSR_mat_vec(a, i_idx, j_idx, vec);
  return vec1;
}

