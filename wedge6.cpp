/*****************************************************************************
*
* This file is part of GwMove hydrogeological software developed by
* Dr. M. A. Sbai
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*
*****************************************************************************/

// Project headers
#include <wedge6.h>

// C++ headers
#include <cmath>
#include <iostream> 
using namespace std;


namespace happ
{
   Wedge6::Wedge6(double K[3], double X[3][6], unsigned int conn[6])
   {
      unsigned int i, j;

      for (i = 0; i < 3; i++)
         tensor[i] = K[i];

      // sets nodal coordinates matrix
      for (j = 0; j < 3; j++) {
         for (i = 0; i < 6; i++) {
            x[j][i] = X[j][i];
         }
      }

      // sets element connectivity
      for (i = 0; i < 6; i++) {
         connectivity[i] = conn[i];
      }

      // initialize the local FE matrices
      for (j = 0; j < 6; j++) {
         mass[j] = 0.;
         for (i = 0; i < 6; i++)
            matrix[j][i] = 0.;
      }

      // other initializations
      Wedge6::sqrt3 = sqrt(3.);
      Wedge6::w = sqrt3 / 4.;
      Wedge6::SetGaussCoord();
   }

   Wedge6::~Wedge6()
   {
   }

   void Wedge6::SetGaussCoord()
   {
      // xi coordinates
      gp[0][0] = 1.0; gp[0][1] = 1.0; gp[0][2] = -0.5;
      gp[0][3] = -0.5; gp[0][4] = -0.5; gp[0][5] = -0.5;

      // eta coordinates 
      gp[1][0] = 0.0;      gp[1][1] = 0.0;      gp[1][2] = sqrt3 / 2.;
      gp[1][3] = -sqrt3 / 2.; gp[1][4] = sqrt3 / 2.; gp[1][5] = -sqrt3 / 2.;

      // zeta coordinates 
      gp[2][0] = 1.0; gp[2][1] = -1.0; gp[2][2] = 1.0;
      gp[2][3] = 1.0; gp[2][4] = -1.0; gp[2][5] = -1.0;
   }

   void Wedge6::SetGaussNumber(unsigned int gauss)
   {
      ngauss = gauss;
   }

   /**
    * Calculate local shape functions at gauss integration
    * point number 'ngauss'.
    *
    * These shape functions are given according to Zienkiewicz (1977):
    *
    * b1(xi,eta,zeta) = 1/6 * (1 + 2*xi)*(1 - zeta)
    * b2(xi,eta,zeta) = 1/6 * (1 - xi + sqrt(3)*eta)*(1 - zeta)
    * b3(xi,eta,zeta) = 1/6 * (1 - xi - sqrt(3)*eta)*(1 - zeta)
    * b4(xi,eta,zeta) = 1/6 * (1 + 2*xi)*(1 + zeta)
    * b5(xi,eta,zeta) = 1/6 * (1 - xi + sqrt(3)*eta)*(1 + zeta)
    * b6(xi,eta,zeta) = 1/6 * (1 - xi - sqrt(3)*eta)*(1 + zeta)
    */
   void Wedge6::ComputeBasis()
   {
      unsigned int n;
      const double sqrt3 = sqrt(3.);

      double a[6];
      for (n = 0; n < 3; n++) {
         a[n] = (1 - gp[2][ngauss]) / 6.;
      }
      for (n = 3; n < 6; n++) {
         a[n] = (1 + gp[2][ngauss]) / 6.;
      }

      basis[0] = (1 + 2 * gp[0][ngauss])*a[0];
      basis[1] = (1 - gp[0][ngauss] + sqrt3*gp[1][ngauss])*a[1];
      basis[2] = (1 - gp[0][ngauss] - sqrt3*gp[1][ngauss])*a[2];
      basis[3] = (1 + 2 * gp[0][ngauss])*a[3];
      basis[4] = (1 - gp[0][ngauss] + sqrt3*gp[1][ngauss])*a[4];
      basis[5] = (1 - gp[0][ngauss] - sqrt3*gp[1][ngauss])*a[5];
   }

   void Wedge6::ComputeLocalGradients()
   {
      unsigned int n, i, j;

      // calculates derivatives of basis functions versus local coordinates
      double a[6];
      for (n = 0; n < 3; n++) {
         a[n] = (1 - gp[2][ngauss]) / 6.;
      }
      for (n = 3; n < 6; n++) {
         a[n] = (1 + gp[2][ngauss]) / 6.;
      }

      // d(b_i)/d(xi)
      by[0][0] = 2 * a[0];
      by[0][1] = -a[1];
      by[0][2] = -a[2];
      by[0][3] = 2 * a[3];
      by[0][4] = -a[4];
      by[0][5] = -a[5];

      // d(b_i)/d(eta)
      by[1][0] = 0;
      by[1][1] = +sqrt3*a[1];
      by[1][2] = -sqrt3*a[2];
      by[1][3] = 0;
      by[1][4] = +sqrt3*a[4];
      by[1][5] = -sqrt3*a[5];

      // d(b_i)/d(zeta)
      by[2][0] = -(1 + 2 * gp[0][ngauss]) / 6.;
      by[2][1] = -(1 - gp[0][ngauss] + sqrt3*gp[1][ngauss]) / 6.;
      by[2][2] = -(1 - gp[0][ngauss] - sqrt3*gp[1][ngauss]) / 6.;
      by[2][3] = -by[2][0];
      by[2][4] = -by[2][1];
      by[2][5] = -by[2][2];

      for (j = 0; j < 3; j++) {
         for (i = 0; i < 6; i++) {
            by[j][i] *= w;
         }
      }
   }

   void Wedge6::ComputeGradients()
   {
      unsigned int n;

      // compute local basis gradients
      Wedge6::ComputeLocalGradients();

      // construct jacobian matrix
      Wedge6::ComputeJacobian();

      // calculate determinant and inverse of the jacobian matrix
      Wedge6::ComputeInverseJacobian();

      // calculate derivatives of basis functions versus global coordinates
      for (n = 0; n < 6; n++) {
         for (unsigned int i = 0; i < 3; i++) {
            bx[i][n] = 0.;
            for (unsigned int j = 0; j < 3; j++) {
               bx[i][n] += by[j][n] * jacinv[j][i];
            }
         }
      }

   }

   void Wedge6::ComputeJacobian()
   {
      unsigned int i, j, n;
      for (i = 0; i < 3; i++) {
         for (j = 0; j < 3; j++) {
            jac[j][i] = 0.;
            for (n = 0; n < 6; n++) {
               jac[j][i] += by[i][n] * x[j][n];
            }
         }
      }
   }

   void Wedge6::ComputeInverseJacobian()
   {
      // first compute determinant if not done yet
      Wedge6::ComputeJacobianDet();

      // calculate inverse matrix coefficients
      jacinv[0][0] = (jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2]) / det;
      jacinv[1][0] = (jac[2][0] * jac[1][2] - jac[1][0] * jac[2][2]) / det;
      jacinv[2][0] = (jac[1][0] * jac[2][1] - jac[2][0] * jac[1][1]) / det;
      jacinv[0][1] = (jac[2][1] * jac[0][2] - jac[0][1] * jac[2][2]) / det;
      jacinv[1][1] = (jac[0][0] * jac[2][2] - jac[2][0] * jac[0][2]) / det;
      jacinv[2][1] = (jac[2][0] * jac[0][1] - jac[0][0] * jac[2][1]) / det;
      jacinv[0][2] = (jac[0][1] * jac[1][2] - jac[1][1] * jac[0][2]) / det;
      jacinv[1][2] = (jac[1][0] * jac[0][2] - jac[0][0] * jac[1][2]) / det;
      jacinv[2][2] = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) / det;
   }

   void Wedge6::ComputeJacobianDet()
   {
      det = jac[0][0] * (jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2]) -
         jac[1][0] * (jac[0][1] * jac[2][2] - jac[2][1] * jac[0][2]) +
         jac[2][0] * (jac[0][1] * jac[1][2] - jac[1][1] * jac[0][2]);
   }

   void Wedge6::ComputeStiffness()
   {
      unsigned int l, k;

      // sum contributions over gauss points
      for (ngauss = 0; ngauss < 6; ngauss++) {

         // calculate basis functions gradients
         Wedge6::ComputeGradients();
         // for unsaturated flow: ComputeGradients() must compute or return the
         // current average pressure head at the cell center ==> but could be
         // directly implemented in a separate function !

         // add contributions for all node-pairs
         for (l = 0; l < 6; l++) {
            for (k = 0; k < 6; k++) {
               // we compute only lower triangular matrix
               if (connectivity[k] <= connectivity[l])
                  matrix[k][l] += det*(tensor[0] * bx[0][l] * bx[0][k] +
                  tensor[1] * bx[1][l] * bx[1][k] +
                  tensor[2] * bx[2][l] * bx[2][k]);
               // in case of unsaturated flow: all stiffness matrix terms
               // must be multiplied by the relative permeability function
               // of the soil type that fills this cell ==> depends on the
               // type of model used (analytical or tabulated). Analytical
               // models are for example the well established Van-Genuchten
               // Mualem-VG and Corey curves. Tabulated data could be either
               // interpolated linearly or using 3rd order splines to get
               // smooth 2nd order differentiables functions. 
            }
         }

      }
   }

   void Wedge6::ComputeMass()
   {
      Wedge6::ComputeBasis();
      Wedge6::ComputeLocalGradients();
      Wedge6::ComputeJacobian();
      Wedge6::ComputeJacobianDet();

      // sum contributions over gauss points
      for (ngauss = 0; ngauss < 6; ngauss++) {
         for (unsigned int l = 0; l < 6; l++)
            mass[l] += basis[l] * det;
      }
   }

}
