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
#include <hex8.h>

// C++ headers
#include <cmath>
using namespace std;


namespace happ
{
   Hex8::Hex8(const double K[3], const double X[3][8], const unsigned int conn[8])
   {
      unsigned int i, j;

      for (i = 0; i < 3; i++)
         tensor[i] = K[i];

      // sets local FE coordinates of nodes
      for (j = 0; j < 3; j++) {
         for (i = 0; i < 8; i++) {
            x[j][i] = X[j][i];
         }
      }

      // sets element connectivity
      for (i = 0; i < 8; i++) {
         connectivity[i] = conn[i];
      }

      // initialize the local FE matrices
      for (j = 0; j < 8; j++) {
         mass[j] = 0.;
         for (i = 0; i < 8; i++)
            matrix[j][i] = 0.;
      }

      // other initializations 
      Hex8::invSqrt3 = 1. / sqrt(3.);
      Hex8::w = 1.0;
      Hex8::SetLocalCoord();
      Hex8::SetGaussCoord();
   }

   Hex8::~Hex8()
   {
   }

   void Hex8::SetGaussCoord()
   {
      unsigned int i, j;
      for (j = 0; j < 3; j++) {
         for (i = 0; i < 8; i++) {
            gp[j][i] = invSqrt3*lc[j][i];
         }
      }
   }

   void Hex8::SetLocalCoord()
   {
      //lc[3][8] = {
      //   { -1, +1, +1, -1, -1, +1, +1, -1 },
      //   { -1, -1, +1, +1, -1, -1, +1, +1 },
      //   { -1, -1, -1, -1, +1, +1, +1, +1 }
      //};

      // xi
      lc[0][0] = -1; lc[0][1] = +1; lc[0][2] = +1; lc[0][3] = -1;
      lc[0][4] = -1; lc[0][5] = +1; lc[0][6] = +1; lc[0][7] = -1;

      // eta
      lc[1][0] = -1; lc[1][1] = -1; lc[1][2] = +1; lc[1][3] = +1;
      lc[1][4] = -1; lc[1][5] = -1; lc[1][6] = +1; lc[1][7] = +1;

      // zeta 
      lc[2][0] = -1; lc[2][1] = -1; lc[2][2] = -1; lc[2][3] = -1;
      lc[2][4] = +1; lc[2][5] = +1; lc[2][6] = +1; lc[2][7] = +1;
   }

   void Hex8::SetGaussNumber(const unsigned int gauss)
   {
      ngauss = gauss;
   }

   void Hex8::ComputeBasis()
   {
      unsigned int n;
      for (n = 0; n < 8; n++) {
         basis[n] = (1 + lc[0][ngauss] * lc[0][n])*
                    (1 + lc[1][ngauss] * lc[1][n])*
                    (1 + lc[2][ngauss] * lc[2][n]) / 8.;
      }
   }

   void Hex8::ComputeLocalGradients()
   {
      // calculates derivatives of basis functions versus local coordinates
      unsigned int n;
      for (n = 0; n < 8; n++) {
         by[0][n] = w*lc[0][n] * (1 + gp[1][ngauss] * lc[1][n])*(1 + gp[2][ngauss] * lc[2][n]) / 8.;
         by[1][n] = w*lc[1][n] * (1 + gp[0][ngauss] * lc[0][n])*(1 + gp[2][ngauss] * lc[2][n]) / 8.;
         by[2][n] = w*lc[2][n] * (1 + gp[0][ngauss] * lc[0][n])*(1 + gp[1][ngauss] * lc[1][n]) / 8.;
      }
   }

   void Hex8::ComputeGradients()
   {
      unsigned int n;

      // compute local basis gradients
      Hex8::ComputeLocalGradients();

      // construct jacobian matrix
      Hex8::ComputeJacobian();

      // calculate determinant and inverse of the jacobian matrix
      Hex8::ComputeInverseJacobian();

      // calculate derivatives of basis functions versus global coordinates
      for (n = 0; n < 8; n++) {
         for (unsigned int i = 0; i < 3; i++) {
            bx[i][n] = 0.;
            for (unsigned int j = 0; j < 3; j++) {
               bx[i][n] += by[j][n] * jacinv[j][i];
            }
         }
      }

   }

   void Hex8::ComputeJacobian()
   {
      unsigned int i, j, n;
      for (i = 0; i < 3; i++) {
         for (j = 0; j < 3; j++) {
            jac[j][i] = 0.;
            for (n = 0; n < 8; n++) {
               jac[j][i] += by[i][n] * x[j][n];
            }
         }
      }
   }

   void Hex8::ComputeInverseJacobian()
   {
      // first compute determinant if not done yet
      Hex8::ComputeJacobianDet();

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

   void Hex8::ComputeJacobianDet()
   {
      det = jac[0][0] * (jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2]) -
         jac[1][0] * (jac[0][1] * jac[2][2] - jac[2][1] * jac[0][2]) +
         jac[2][0] * (jac[0][1] * jac[1][2] - jac[1][1] * jac[0][2]);
   }

   void Hex8::ComputeStiffness(const eFlowProcess process)
   {
      unsigned int l, k;

      // sum contributions over gauss points
      for (ngauss = 0; ngauss < 8; ngauss++) {

         // calculate basis functions gradients
         Hex8::ComputeGradients();
         // for unsaturated flow: ComputeGradients() must compute or return the
         // current average pressure head at the cell center ==> but could be
         // directly implemented in a separate function !

         // add contributions for all node-pairs
         for (l = 0; l < 8; l++) {
            for (k = 0; k < 8; k++) {
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

   void Hex8::ComputeStiffRowsSum()
   {
      unsigned int row, col;
      double sum;

      // work on a temporary array whose lower entries are filled 
      // since there're zero in 'matrix' (due to symmetry) 
      double buffer[8][8]; 
      for (row = 0; row < 8; row++) {
         for (col = 0; col < 8; col++) {
            if (col >= row)
               buffer[row][col] = matrix[row][col];
            else 
               buffer[row][col] = matrix[col][row];
         }
      }

      for (row = 0; row < 8; row++) {
         sum = 0.0;
         for (col = 0; col < 8; col++) sum += buffer[row][col];
         rowSum[row] = sum;
      }
   }

   bool Hex8::IsPositiveDefinite()
   {
      unsigned int row; 

      /* 
       * positive definitess is verified when:
       * - The sum of row entries is positive 
       */
      Hex8::ComputeStiffRowsSum(); 

      bool positive = true;
      for (row = 0; row < 8; row++) {
         if (rowSum[row] < -std::numeric_limits<double>::epsilon()*10.0) {
            positive = false;
            break;
         }
      }

      return positive;
   }

   bool Hex8::IsMmatrix() 
   {
      unsigned int row, col;

      /*
       * According to Axelsson and Barker (1984) an M-matrix satisfy three
       * conditions:
       * - All diagonals are strictly positive
       * - All off-diagonals are null or negative 
       * - Is positive definite (i.e. invertible) 
       */
      
      // verify first condition
      bool diagPositive = true; 
      for (row = 0; row < 8; row++) {
         if (matrix[row][row] <= -std::numeric_limits<double>::epsilon()*10.0) {
            diagPositive = false;
            break;
         }
      }
      if (diagPositive == false) return false;

      // verify second condition 
      bool offDiags = true; 
      for (row = 0; row < 8; row++) {
         for (col = 0; col < 8; col++) {
            if (col == row) continue;
            if (matrix[row][col] > -std::numeric_limits<double>::epsilon()*10.0) {
               offDiags = false;
               break;
            }
         }
      }
      if (offDiags == false) return false;

      // verify third condition 
      bool posdef = Hex8::IsPositiveDefinite();
      if (posdef == false) return false;

      // All tests passed
      return true;
   }

   void Hex8::ComputeMass(const eFlowProcess process)
   {
      // sum contributions over gauss points
      for (ngauss = 0; ngauss < 8; ngauss++) {
         Hex8::ComputeBasis();
         Hex8::ComputeLocalGradients();
         Hex8::ComputeJacobian();
         Hex8::ComputeJacobianDet();

         for (unsigned int l = 0; l < 8; l++)
            mass[l] += basis[l] * det;
      }
   }

   bool Hex8::Interpolate(const Tpoint3D& pt, const double var[8], double& value)
   {
      // chech whether the point 'pt' is inside the element otherwise return false !
      // Need to implement the method Hex8::IsInside(Tpoint3D&)

      unsigned int i;

      value = 0.0;
      for (ngauss = 0; ngauss < 8; ngauss++) {
         Hex8::ComputeBasis();

         for (i = 0; i < 8; i++) 
            value += var[i] * basis[i];
      }

      return true;
   }

   bool Hex8::Interpolate(const std::vector<Tpoint3D> pt,
      const double var[8],
      std::vector<double> values)
   {
      unsigned int s1 = pt.size(); 
      values.resize(s1);

      unsigned int i, j;

      for (j = 0; j < s1; j++) values[j] = 0.0;

      for (ngauss = 0; ngauss < 8; ngauss++) {
         Hex8::ComputeBasis();

         for (j = 0; j < s1; j++)
            for (i = 0; i < 8; i++) values[j] += var[i] * basis[i];
      }

      return true;
   }

   bool Hex8::IsInside(const Tpoint3D& pt) const
   {
      // see implementation in 'ptracker.cpp' and copy it here
      return true; 
   }

   Tpoint3D Hex8::CellCentroidVelocity(const double head[8])
   {
      double velocity[3];
      unsigned int i, j;
      for (i = 0; i < 3; i++) velocity[i] = 0.0;

      for (ngauss = 0; ngauss < 8; ngauss++) {

         // Calculate global derivatives of basis functions (if not yet)
         Hex8::ComputeGradients();

         // calculate directional gradients of groundwater heads
         for (j = 0; j < 3; j++) {
            for (i = 0; i < 8; i++) {
               velocity[j] += -head[i] * bx[j][i];
            }
         }

         // Eventually call unsat(ie,p,wc,uc) when dealing 
         // with unsaturated flow processes ... 

      }

      // now calculate cell centered velocity vector
      // by simple application of Darcy's law locally ... 
      for (i = 0; i < 3; i++) velocity[i] *= tensor[i];

      // copy to Tpoint3D type 
      Tpoint3D v(velocity[0], velocity[1], velocity[2]);
      return v;
   }

   void Hex8::Print(ostream& os) const 
   {
      unsigned int row, col;

      // 'by' array 
      os << "Deriv vs. local coordinates array 'by' =" << endl;
      os << setiosflags(ios::scientific);
      os << setprecision(3);
      for (row = 0; row < 3; row++) {
         for (col = 0; col < 8; col++)
            os << setw(15) << by[row][col];
         os << endl;
      }
      os << endl;

      // Jacobian matrix
      os << "Jacobian matrix 'jac' =" << endl;
      for (row = 0; row < 3; row++) {
         for (col = 0; col < 3; col++)
            os << setw(15) << jac[row][col];
         os << endl;
      }
      os << endl;

      // Jacobian determinant 
      os << "Jacobian matrix determinant 'det' = "
         << det << endl << endl;

      // Inverse of the Jacobian matrix 
      os << "Inverse Jacobian matrix 'jacinv' =" << endl;
      for (row = 0; row < 3; row++) {
         for (col = 0; col < 3; col++)
            os << setw(15) << jacinv[row][col];
         os << endl;
      }
      os << endl;

      // 'bx' derivatives array
      os << "Deriv vs. global coordinates array 'bx' =" << endl;
      for (row = 0; row < 3; row++) {
         for (col = 0; col < 8; col++)
            os << setw(15) << bx[row][col];
         os << endl;
      }
      os << endl;
   }

   void Hex8::PrintMatrix(ostream& os) const 
   {
      unsigned int row, col;

      os << "Local finite element stiffness 'matrix' =" << endl;
      os << setiosflags(ios::scientific);
      os << setprecision(3);
      for (row = 0; row < 8; row++) {
         for (col = 0; col < 8; col++)
            os << setw(15) << matrix[row][col];
         os << endl;
      }
      os << endl;
   }

   void Hex8::PrintRowSum(ostream& os) const
   {
      unsigned int col;

      os << "Rows sum of the local stiffness 'matrix' =" << endl;
      os << setiosflags(ios::scientific);
      os << setprecision(3);
      for (col = 0; col < 8; col++)
         os << setw(15) << rowSum[col];
      os << endl;
   }

   void Hex8::PrintMassMatrix(ostream& os) const
   {
      unsigned int col;

      os << "Local finite element 'mass' matrix =" << endl;
      os << setiosflags(ios::scientific);
      os << setprecision(3);
      for (col = 0; col < 8; col++)
         os << setw(15) << mass[col];
      os << endl;
   }

}
