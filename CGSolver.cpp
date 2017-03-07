#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CGSolver.hpp"
#include "matvecops.hpp"
#include "sparse.hpp"


/* Writing into solution function */
void Writing_Solution(std::vector<double> x,
		      std::vector<double> b,
		      std::string         soln_prefix,
		      int                 niter,
		      int                 colid)
{
  /* Setting up the output file name */
  std::stringstream output0;
  output0 << soln_prefix;
  output0 << std::setfill('0') << std::setw(3) << niter;
  output0 << ".txt";
  std::ofstream g(output0.str());
  if (g.is_open())
    {
      /* Writing the isothermic hot boundary condition */
      for (int i = 0; i < colid; i++)
        {
          g << std::scientific;
          g << std::setprecision(4);
          g << b[i] << std::endl;
        }
      g << b[0] << std::endl; // periodic point at the top right
      for (int i = 0; i < (int) x.size(); i++)
        {
          g << x[i] << std::endl;
	  /* Writing the right side periodic boundary */
	  if (i % colid == colid - 1)
	    {
	      g << x[i + 1 - colid] << std::endl;
	    }
        }
      /* Writing the bottom cold condition */
      for (int i = (int) b.size() - colid; i < (int) b.size(); i++)
        {
          g << b[i] << std::endl;
        }
      // perdiodic point at the bottom right
      g << b[(int) b.size() - colid + 1] << std::endl;
      g.close();
    }
}

int CGSolver(SparseMatrix         A,
	     std::vector<double> &b,
	     std::vector<double> &x,
	     double              tol,
	     int                 colid,
	     std::string         soln_prefix)
{
  /* r_0 = b - Ax  */
  std::vector<double> r0 = vec_sub(b, A.MulVec(x));
  double l2norm_r0 = vec_norm(r0);
  /* p_0 = r_0  */
  std::vector<double> p0 = r0;
  int niter = 0;

  /* Printing out the initial solution */
  Writing_Solution(x, b, soln_prefix, niter, colid);

  int nitermax = (int) x.size();

  while (niter < nitermax)
    {
      niter += 1;
      /* Finding the alpha  */
      std::vector<double> Ap = A.MulVec(p0);
      double l2norm_r_n_2 = dot_pro(r0, r0);
      double alpha = l2norm_r_n_2/dot_pro(p0, Ap);

      /* Updating x */
      x = vec_add(x, vec_sca(p0, alpha));


      /* Printing the 10x solution */
      if (niter % 10 == 0)
	{
	  Writing_Solution(x, b, soln_prefix, niter, colid);
	}

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

  /* Printing the final solution */
  Writing_Solution(x, b, soln_prefix, niter, colid);
  
  return niter;
}


//--design_0
//--very good! i like your decomposition into writing/solving fns
//--END
