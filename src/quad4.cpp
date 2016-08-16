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
#include <quad4.h>

// C++ headers
#include <cmath>
#include <iostream> 
using namespace std;


namespace happ
{
   Quad4::Quad4(double K[2], double X[2][4], unsigned int conn[4])
   {
      unsigned int i, j;

      for (i = 0; i < 2; i++)
         tensor[i] = K[i];

      // sets local FE coordinates of nodes
      for (j = 0; j < 2; j++) {
         for (i = 0; i < 4; i++) {
            x[j][i] = X[j][i];
         }
      }

      // sets element connectivity
      for (i = 0; i < 4; i++) {
         connectivity[i] = conn[i];
      }

      // initialize the local FE matrices
      for (j = 0; j < 4; j++) {
         mass[j] = 0.;
         for (i = 0; i < 4; i++)
            matrix[j][i] = 0.;
      }

      // other initializations 
      Quad4::invSqrt3 = 1. / sqrt(3.);
      Quad4::w = 1.0;
      Quad4::SetGaussCoord();
      Quad4::SetLocalCoord();
   }

   Quad4::~Quad4()
   {
   }

   void Quad4::SetGaussCoord()
   {
      // xi coordinates
      gp[0][0] = -invSqrt3; gp[0][1] = +invSqrt3;
      gp[0][2] = +invSqrt3; gp[0][3] = -invSqrt3;

      // eta coordinates 
      gp[1][0] = -invSqrt3; gp[1][1] = -invSqrt3;
      gp[1][2] = +invSqrt3; gp[1][3] = +invSqrt3;
   }

   void Quad4::SetLocalCoord()
   {
      // xi coordinates
      lc[0][0] = -1.0; lc[0][1] = +1.0;
      lc[0][2] = +1.0; lc[0][3] = -1.0;

      // eta coordinates 
      lc[1][0] = -1.0; lc[1][1] = -1.0;
      lc[1][2] = +1.0; lc[1][3] = +1.0;
   }

   void Quad4::SetGaussNumber(unsigned int gauss)
   {
      ngauss = gauss;
   }

   /**
     Calculate local shape functions at gauss integration
     point number 'ngauss'.

     These shape functions are given according to Zienkiewicz (1977):
     \f[
     b1(\xi,\eta) = \frac{1}{4} (1-\xi)(1-\eta)
     \f]
     \f[
     b2(\xi,\eta) = \frac{1}{4} (1+\xi)(1-\eta)
     \f]
     \f[
     b3(\xi,\eta) = \frac{1}{4} (1+\xi)(1+\eta)
     \f]
     \f[
     b4(\xi,\eta) = \frac{1}{4} (1-\xi)(1+\eta)
     \f]
     */
   void Quad4::ComputeBasis()
   {
      unsigned int n;

      for (n = 0; n < 4; n++) {
         basis[n] = (1 + gp[0][ngauss] * lc[0][n])*
            (1 + gp[1][ngauss] * lc[1][n]) / 4.;
      }
   }

   void Quad4::ComputeLocalGradients()
   {
      unsigned int n;


      // calculates derivatives of basis functions versus local coordinates
      for (n = 0; n < 4; n++) {
         by[0][n] = w*lc[0][n] * (1 + gp[1][ngauss] * lc[1][n]) / 4.;
         by[1][n] = w*lc[1][n] * (1 + gp[0][ngauss] * lc[0][n]) / 4.;
      }
   }

   void Quad4::ComputeGradients()
   {
      unsigned int n;

      // compute local basis gradients
      Quad4::ComputeLocalGradients();

      // construct jacobian matrix
      Quad4::ComputeJacobian();

      // calculate determinant and inverse of the jacobian matrix
      Quad4::ComputeInverseJacobian();

      // calculate derivatives of basis functions versus global coordinates
      for (n = 0; n < 4; n++) {
         for (unsigned int i = 0; i < 2; i++) {
            bx[i][n] = 0.;
            for (unsigned int j = 0; j < 2; j++) {
               bx[i][n] += by[j][n] * jacinv[j][i];
            }
         }
      }

   }

   void Quad4::ComputeJacobian()
   {
      unsigned int i, j, n;
      for (i = 0; i < 2; i++) {
         for (j = 0; j < 2; j++) {
            jac[j][i] = 0.;
            for (n = 0; n < 4; n++) {
               jac[j][i] += by[i][n] * x[j][n];
            }
         }
      }
   }

   void Quad4::ComputeInverseJacobian()
   {
      // first compute determinant if not done yet
      Quad4::ComputeJacobianDet();

      // calculate inverse matrix coefficients
      jacinv[0][0] = +jac[1][1] / det;
      jacinv[1][0] = -jac[1][0] / det;
      jacinv[0][1] = -jac[0][1] / det;
      jacinv[1][1] = +jac[0][0] / det;
   }

   void Quad4::ComputeJacobianDet()
   {
      det = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1];
   }

   void Quad4::ComputeStiffness()
   {
      unsigned int l, k;

      // sum contributions over gauss points
      for (ngauss = 0; ngauss < 4; ngauss++) {

         // calculate basis functions gradients
         Quad4::ComputeGradients();
         // for unsaturated flow: ComputeGradients() must compute or return the
         // current average pressure head at the cell center ==> but could be
         // directly implemented in a separate function !

         // add contributions for all node-pairs
         for (l = 0; l < 4; l++) {
            for (k = 0; k < 4; k++) {
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

   void Quad4::ComputeMass()
   {
      Quad4::ComputeBasis();
      Quad4::ComputeLocalGradients();
      Quad4::ComputeJacobian();
      Quad4::ComputeJacobianDet();

      // sum contributions over gauss points
      for (ngauss = 0; ngauss < 4; ngauss++) {
         for (unsigned int l = 0; l < 4; l++)
            mass[l] += basis[l] * det;
      }
   }
}
