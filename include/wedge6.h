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

#ifndef HAPP_WEDGE6
#define HAPP_WEDGE6


namespace happ
{
   /**
    * Linear Wedge or prism finite element class.
    * It does numerical integration and setup of shape functions,
    * gradients, and Jacobian calculations for one wedge element.
    */
   class Wedge6
   {
   private:
      /// number of gauss point
      unsigned int ngauss;

      double sqrt3;

      // weight of six-points Gaussian quadrature rule 
      double w;

      // coordinates of gauss points
      double gp[3][6];

      /// sets gauss coordinates 
      void SetGaussCoord();

   public:
      /// local matrix dimension
      unsigned int dimension;

      /// local finite element matrix
      double matrix[6][6];

      /// mass-lumped storativity matrix
      double mass[6];

      /// components of hydraulic conductivity
      double tensor[3];

      /// nodal coordinates matrix
      double x[3][6];

      /// nodal connectivity of this element
      unsigned int connectivity[6];

      /// basis functions
      double basis[6];

      /// derivatives of basis functions versus local coordinates
      double by[3][6];

      /// derivatives of basis functions versus global coordinates
      double bx[3][6];

      /// Jacobian matrix
      double jac[3][3];

      /// inverse of the jacobian matrix
      double jacinv[3][3];

      /// determinant of Jacobian matrix
      double det;


      /// constructor
      Wedge6(double K[3], double X[3][6], unsigned int conn[6]);

      /// destructor
      ~Wedge6();

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

#endif // HAPP_WEDGE6
