#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CGSolver.hpp"
#include "matvecops.hpp"

int CGSolver(std::vector<double> &val,
             std::vector<int>    &i_idx,
             std::vector<int>    &j_idx,
	           std::vector<double> &b,
	           std::vector<double> &x,
	           double              tol)
{
  /* r_0 = b - Ax  */
  std::vector<double> r0 = vec_sub(b, CSR_mat_vec(val, i_idx, j_idx, x));
  double l2norm_r0 = vec_norm(r0);
  /* p_0 = r_0  */
  std::vector<double> p0 = r0;
  int niter = 0;

  int nitermax = (int) x.size();
  while (niter < nitermax)
    {
      niter += 1;
      /* Finding the alpha  */
      std::vector<double> Ap = CSR_mat_vec(val, i_idx, j_idx, p0);
      double l2norm_r_n_2 = dot_pro(r0, r0);
      double alpha = l2norm_r_n_2/dot_pro(p0, Ap);

      /* Updating x */
      x = vec_add(x, vec_sca(p0, alpha));

      /* Updating r */
      std::vector<double > rn_1 = vec_sub(r0, vec_sca(Ap, alpha));
      double l2norm_r = vec_norm(rn_1);

      /* checking for threshold */
      if (l2norm_r/l2norm_r0 < tol)
        {
          break;
        }

      double beta = pow(l2norm_r, 2)/l2norm_r_n_2;
      p0 = vec_add(rn_1, vec_sca(p0, beta));
      r0 = rn_1;
    }  
    
  return niter;
}


//--design_0
//--very good! i like your decomposition into writing/solving fns
//--END
