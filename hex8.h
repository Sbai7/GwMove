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

#ifndef HAPP_HEX8
#define HAPP_HEX8

#include <point3d.h> 
#include <soil.h>

#include <vector>
#include <iostream> 
using namespace std;


namespace happ
{
   /**
    * Linear Hexahedral finite element class.
    * It does numerical integration and setup of shape functions,
    * gradients, and Jacobian calculations for one hexahedral.
    */
   class Hex8
   {
   private:
      /// number of gauss point
      unsigned int ngauss;

      double invSqrt3;

      /// weight of the eight-point Gaussian quadrature rule
      double w;

      /// coordinates of gauss points
      double gp[3][8];

      /// sets coordinates of Gauss points 
      void SetGaussCoord();

      /// array of local coordinates of nodes in reference cube
      int lc[3][8];

      /// sets local coordinates values
      void SetLocalCoord();

      /// local flow process enum type 
      enum eFlowProcess {
         confined, 
         unconfined, 
         unconfined_lam, 
         unsaturated
      };

   public:
      /// local matrix dimension
      unsigned int dimension;

      /// local finite element matrix
      double matrix[8][8];

      /// mass-lumped storativity matrix
      double mass[8];

      /// components of hydraulic conductivity
      double tensor[3];

      /// nodal coordinates matrix
      double x[3][8];

      /// nodal connectivity of this element
      unsigned int connectivity[8];

      /// basis functions
      double basis[8];

      /// derivatives of basis functions versus local coordinates
      double by[3][8];

      /// derivatives of basis functions versus global coordinates
      double bx[3][8];

      /// Jacobian matrix
      double jac[3][3];

      /// inverse of the jacobian matrix
      double jacinv[3][3];

      /// determinant of Jacobian matrix
      double det;

      /// rowise sum of the stiffness matrix 
      double rowSum[8];


      /// constructor
      Hex8(const double K[3], const double X[3][8], const unsigned int conn[8]);

      /// destructor
      virtual ~Hex8();

      /// Return local coordinates 
      inline void GetLocalCoord(int Lc[3][8]) const {
         unsigned int i, j;
         for (j = 0; j < 3; j++) {
            for (i = 0; i < 8; i++) {
               Lc[j][i] = lc[j][i];
            }
         }
      }

      /// Sets gauss point number
      virtual void SetGaussNumber(const unsigned int gauss);

      /// Compute basis functions contribution from nth gauss point.
      virtual void ComputeBasis();

      /// Compute local basis functions gradients
      virtual void ComputeLocalGradients();

      /// Compute basis functions gradients
      virtual void ComputeGradients();

      /// Compute Jacobian matrix
      virtual void ComputeJacobian();

      /// Compute jacobian determinant
      virtual void ComputeJacobianDet();

      /// Compute inverse Jacobian matrix
      virtual void ComputeInverseJacobian();

      /// Compute local stiffness FE matrix
      virtual void ComputeStiffness(const eFlowProcess process = eFlowProcess::confined);

      /// Compute rows sum of the local stiffness FE matrix 
      virtual void ComputeStiffRowsSum();

      /// Is the local stiffness matrix positive definite? 
      virtual bool IsPositiveDefinite(); 

      /// Is the local stiffness matrix an M-matrix? 
      virtual bool IsMmatrix(); 

      /// Compute local mass FE matrix
      virtual void ComputeMass(const eFlowProcess process = eFlowProcess::confined);

      /// Interpolate variable at a position inside the element 
      virtual bool Interpolate(const Tpoint3D& pt, const double var[8], double& value);
      virtual bool Interpolate(const std::vector<Tpoint3D> pt,
                               const double var[8],
                               std::vector<double> values);
      /* Finite element interpolation is a non-trivial task since it requires 
       * calculation of the inverse mapping for the hexahedral cell.
       * This could be performed analytically or numerically.
       * - For an analytic derivation see the following paper:
       *   * Yuan et al., 1994. The inverse mapping and distorsion measures 
       *     for 8-node hexahedral isoparametric elements. Comput. Mech., 
       *     14, 189-199.
       *   The authors use the thory of geodesics in differential geometry  
       *   to do such thing. 
       *   The principal limitation of the analytic method is the non 
       *   invertibility for highly distorted elements. However, it is the 
       *   fastest.
       * - For a numerical inversion the trick consists on iterative solution 
       *   of a nonlinear system from an initial trial position. For a 
       *   description along with formulas derivations see the papers: 
       *   * Silva et al., 2007. Exact and efficient interpolation using 
       *     finite elements shape functions.
       *   * Yogev et al., 2010. A novel energy-based approach for merging 
       *     finite elements. IJNME, 1010.
       */
      
      /// Detect if a given point is inside this cell
      virtual bool IsInside(const Tpoint3D& pt) const;
      /* This is a non trivial test.
       * See the algorithms developed in the following paper:
       *    Sila et al., 2007. Exact and efficient interpolation using 
       *    finite elements shape functions.
       * Three search techniques algorithms were developed and compared:
       * ST1: sequential scan of the mesh using the cross-product test.
       *      the latter consists in a series of cross and dot products.
       * ST2: sequential scan of the mesh using a combination of cross-
       *      product and bounding-box tests. 
       * ST3: uses cross-product and bounding-box tests to search a 
       *      small list of elements obtained from a virtual mesh. 
       */

      /// Calculate cell centered velocity. 
      // In all rigor this method should be inerited from an upper
      // base class desribing general 'Element' data & methods. 
      virtual Tpoint3D CellCentroidVelocity(const double head[8]);

      /// prints calculated FE arrays of the hexahedral element 
      virtual void Print(ostream& os) const;

      /// Prints locally assembled stiffness FE matrix
      virtual void PrintMatrix(ostream& os) const; 

      /// Prints rows sums of the locally assembled stiffness FE matrix
      virtual void PrintRowSum(ostream& os) const;

      /// Prints locally assembled mass matrix 
      virtual void PrintMassMatrix(ostream& os) const;

   };
}

#endif // HAPP_HEX8 