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

#ifndef HAPP_TRI3
#define HAPP_TRI3


namespace happ
{
   /**
    * Three nodes triangular finite element class.
    * It does numerical integration and setup of shape functions,
    * gradients, and Jacobian calculations for a linear triangle
    * with 3 corner nodes.
    */
   class Tri3
   {
   private:
      /// number of gauss point
      unsigned int ngauss;

      double sqrt3;

      /// weight of four-point Gaussian quadrature rule 
      double w[4];

      /// coordinates of gauss points
      double gp[2][4];

      /// sets gauss coordinates
      void SetGaussCoord();

      /// local coordinates of corner nodes in (xi,eta) space
      double lc[2][3];

      /// sets local coordinates of quad nodes
      void SetLocalCoord();

   public:
      /// local matrix dimension
      unsigned int dimension;

      /// local finite element matrix
      double matrix[3][3];

      /// mass-lumped storativity matrix
      double mass[3];

      /// components of hydraulic conductivity
      double tensor[2];

      /// nodal coordinates matrix
      double x[2][3];

      /// nodal connectivity of this element
      unsigned int connectivity[3];

      /// basis functions
      double basis[3];

      /// derivatives of basis functions versus local coordinates
      double by[2][3];

      /// derivatives of basis functions versus global coordinates
      double bx[2][3];

      /// Jacobian matrix
      double jac[2][2];

      /// inverse of the jacobian matrix
      double jacinv[2][2];

      /// determinant of Jacobian matrix
      double det;


      /// constructor
      Tri3(double K[2], double X[2][3], unsigned int conn[3]);

      /// destructor
      ~Tri3();

      /// Sets gauss point number
      void SetGaussNumber(unsigned int gauss);

      /// Compute basis functions
      void ComputeBasis();

      /// Compute local basis functions gradients
      void ComputeLocalGradients();

      /// Compute basis functions gradients
      void ComputeGradients();

      /// Compute Jacobian matrix
      void ComputeJacobian();

      /// Compute jacobian determinant
      void ComputeJacobianDet();

      /// Compute inverse Jacobian matrix
      void ComputeInverseJacobian();

      /// Compute local stiffness FE matrix
      void ComputeStiffness();

      /// Compute local mass FE matrix
      void ComputeMass();

   };
}

#endif // HAPP_TRI3
