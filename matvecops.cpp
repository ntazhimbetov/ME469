#include <vector>
#include <cmath>

#include "matvecops.hpp"

/* Multiplying elementwise and taking the sum*/
double dot_pro(std::vector<double> x,
               std::vector<double> y)
{
  int len = (int) x.size();
  double z = 0;
  for (int i = 0; i < len; i++)
    {
      z += x[i]*y[i];
    }
  return z;
}

/* Taking the square root of a dot product of one vector  */
double vec_norm(std::vector<double> x)
{
  double l2_norm = sqrt(dot_pro(x, x));
  return l2_norm;
}

/* Elementwise addition  */
std::vector<double> vec_add(std::vector<double> x,
                            std::vector<double> y)
{
  int len = (int) x.size();
  std::vector<double> z;
  for (int i = 0; i < len; i++)
    {
      z.push_back(x[i] + y[i]);
    }
  return z;
}

/* Elementwise scaling  */
std::vector<double> vec_sca(std::vector<double> x,
			    double c)
{
  int len = (int) x.size();
  std::vector<double> y;
  for (int i = 0; i < len; i++)
    {
      y.push_back(c*x[i]);
    }
  return y;
}

/* Elementwise subtraction  */
std::vector<double> vec_sub(std::vector<double> x,
                            std::vector<double> y)
{
  int len = (int) x.size();
  std::vector<double> z;
  for (int i = 0; i < len; i++)
    {
      z.push_back(x[i] - y[i]);
    }
  return z;
}

/* Sparse matrix-vector multiplication  */
std::vector<double> CSR_mat_vec(std::vector<double> val,
                                std::vector<int>    row_ptr,
                                std::vector<int>    col_idx,
                                std::vector<double> x)
{
  int len = (int) row_ptr.size();
  std::vector<double> y;
  int col_num = 0;
  for (int i = 0; i < len-1; i++)
    {
      int new_col_num = row_ptr[i+1] - row_ptr[i];
      double row_sum = 0;
      for (int j = col_num; j < col_num + new_col_num; j++)
	{
	  row_sum += val[j]*x[col_idx[j]];
	}
      col_num = new_col_num + col_num;
      y.push_back(row_sum);
    }

  return y;
}
