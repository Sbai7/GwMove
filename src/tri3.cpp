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
#include <tri3.h>

// C++ headers
#include <cmath>
#include <iostream> 
using namespace std;


namespace happ
{
   Tri3::Tri3(double K[2], double X[2][3], unsigned int conn[3])
   {
      unsigned int i, j;

      for (i = 0; i < 2; i++)
         tensor[i] = K[i];

      // sets local FE coordinates of nodes
      for (j = 0; j < 2; j++) {
         for (i = 0; i < 3; i++) {
            x[j][i] = X[j][i];
         }
      }

      // sets element connectivity
      for (i = 0; i < 3; i++) {
         connectivity[i] = conn[i];
      }

      // initialize the local FE matrices
      for (j = 0; j < 3; j++) {
         mass[j] = 0.;
         for (i = 0; i < 3; i++)
            matrix[j][i] = 0.;
      }

      // other initializations 
      Tri3::sqrt3 = sqrt(3.);
      Tri3::w[0] = 9.0*sqrt3 / 16.;
      for (i = 1; i < 4; i++) {
         Tri3::w[i] = sqrt3 / 16.;
      }
      Tri3::SetGaussCoord();
      Tri3::SetLocalCoord();
   }

   Tri3::~Tri3()
   {
   }

   void Tri3::SetGaussCoord()
   {
      // xi coordinates
      gp[0][0] = 0.0; gp[0][1] = 1.0;
      gp[0][2] = -0.5; gp[0][3] = -0.5;

      // eta coordinates 
      gp[1][0] = 0.0;      gp[1][1] = 0.0;
      gp[1][2] = sqrt3 / 2.; gp[1][3] = -sqrt3 / 2.;
   }

   void Tri3::SetLocalCoord()
   {
      unsigned int i, j;
      for (j = 0; j < 2; j++) {
         for (i = 0; i < 3; i++) {
            lc[j][i] = gp[j][i + 1];
         }
      }
   }

   void Tri3::SetGaussNumber(unsigned int gauss)
   {
      ngauss = gauss;
   }

   /**
     Calculate local shape functions at gauss integration
     point number 'ngauss'.

     These shape functions are given according to Zienkiewicz (1977):
     \f[
     b1(\xi,\eta) = \frac{1}{3} (1+2\xi)
     \f]
     \f[
     b2(\xi,\eta) = \frac{1}{3} (1-\xi+\sqrt{3}\eta)
     \f]
     \f[
     b3(\xi,\eta) = \frac{1}{3} (1-\xi-\sqrt{3}\eta)
     \f]
     */
   void Tri3::ComputeBasis()
   {
      basis[0] = (1 + 2 * gp[0][ngauss]) / 3.;
      basis[1] = (1 - gp[0][ngauss] + sqrt3*gp[1][ngauss]) / 3.;
      basis[2] = (1 - gp[0][ngauss] - sqrt3*gp[1][ngauss]) / 3.;
   }

   void Tri3::ComputeLocalGradients()
   {
      // d(b_i)/d(xi)
      by[0][0] = +2. / 3.;
      by[1][1] = -1. / 3.;
      by[1][2] = -1. / 3.;

      // d(b_i)/d(eta)
      by[1][0] = 0.;
      by[1][1] = +sqrt3 / 3.;
      by[1][2] = -sqrt3 / 3.;

      // multiply by the weight at current gauss point 
      unsigned int i, j;
      for (j = 0; j < 2; j++) {
         for (i = 0; i < 3; i++) {
            by[j][i] *= w[ngauss];
         }
      }
   }

   void Tri3::ComputeGradients()
   {
      unsigned int n;

      // compute local basis gradients
      Tri3::ComputeLocalGradients();

      // construct jacobian matrix
      Tri3::ComputeJacobian();

      // calculate determinant and inverse of the jacobian matrix
      Tri3::ComputeInverseJacobian();

      // calculate derivatives of basis functions versus global coordinates
      for (n = 0; n < 3; n++) {
         for (unsigned int i = 0; i < 2; i++) {
            bx[i][n] = 0.;
            for (unsigned int j = 0; j < 2; j++) {
               bx[i][n] += by[j][n] * jacinv[j][i];
            }
         }
      }

   }

   void Tri3::ComputeJacobian()
   {
      unsigned int i, j, n;
      for (i = 0; i < 2; i++) {
         for (j = 0; j < 2; j++) {
            jac[j][i] = 0.;
            for (n = 0; n < 3; n++) {
               jac[j][i] += by[i][n] * x[j][n];
            }
         }
      }
   }

   void Tri3::ComputeInverseJacobian()
   {
      // first compute determinant if not done yet
      Tri3::ComputeJacobianDet();

      // calculate inverse matrix coefficients
      jacinv[0][0] = +jac[1][1] / det;
      jacinv[1][0] = -jac[1][0] / det;
      jacinv[0][1] = -jac[0][1] / det;
      jacinv[1][1] = +jac[0][0] / det;
   }

   void Tri3::ComputeJacobianDet()
   {
      det = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1];
   }

   void Tri3::ComputeStiffness()
   {
      unsigned int l, k;

      // sum contributions over gauss points
      for (ngauss = 0; ngauss < 4; ngauss++) {

         // calculate basis functions gradients
         Tri3::ComputeGradients();
         // for unsaturated flow: ComputeGradients() must compute or return the
         // current average pressure head at the cell center ==> but could be
         // directly implemented in a separate function !

         // add contributions for all node-pairs
         for (l = 0; l < 3; l++) {
            for (k = 0; k < 3; k++) {
               // we compute only lower triangular matrix
               if (connectivity[k] <= connectivity[l])
                  matrix[k][l] += det*(tensor[0] * bx[0][l] * bx[0][k] +
                  tensor[1] * bx[1][l] * bx[1][k]);
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

   void Tri3::ComputeMass()
   {
      Tri3::ComputeBasis();
      Tri3::ComputeLocalGradients();
      Tri3::ComputeJacobian();
      Tri3::ComputeJacobianDet();

      // sum contributions over gauss points
      for (ngauss = 0; ngauss < 4; ngauss++) {
         for (unsigned int l = 0; l < 3; l++)
            mass[l] += basis[l] * det;
      }
   }

}
