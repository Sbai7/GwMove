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

#ifndef HAPP_POINT
#define HAPP_POINT

#define _USE_MATH_DEFINES // mathematical constants

// Project headers
#include <platform.h>

// C++ headers
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;


const double pointPrecision = 1.E-8; 


namespace happ
{
   /**
    * X, Y, Z space coordinates data structure.
    */

   struct Tpoint3D
   {
      /// X-coordinate
      double x;

      /// Y-coordinate
      double y;

      /// Z-ccordinate
      double z;

      /// Default constructor
      inline Tpoint3D(double a = 0., double b = 0., double c = 0.) {
         x = a;
         y = b;
         z = c;
      };

      /// overloaded [] operator
      inline double operator[](unsigned int i) {
         double a = 0.;
         switch (i) {
         case 0:
            a = x;
            break;
         case 1:
            a = y;
            break;
         case 2:
            a = z;
            break;
         default:
            break;
         }
         return a;
      };

      /// Truncates coordinates if they smaller than machine's precision
      void	Truncate() {
         PlatformConstants platfm;
         if (fabs(x) < platfm.epsilon)
            x = 0.;
         if (fabs(y) < platfm.epsilon)
            y = 0.;
         if (fabs(z) < platfm.epsilon)
            z = 0.;
      };

      /// overloaded = operator for points
      Tpoint3D& operator =(const Tpoint3D& a) {
         this->x = a.x;
         this->y = a.y;
         this->z = a.z;

         return *this;
      };

   };

   /// overloaded + operator for points
   Tpoint3D operator+(Tpoint3D& a, Tpoint3D& b);

   /// overloaded - operator for points
   Tpoint3D operator-(Tpoint3D& a, Tpoint3D& b);

   /// overloaded / operator for points
   Tpoint3D operator/(Tpoint3D& a, double b);

   /// overloaded * operator for points
   Tpoint3D operator*(Tpoint3D& a, double b);

   /// overloaded == operator for points
   bool operator ==(Tpoint3D& a, Tpoint3D& b);

   /// overloaded < operator for points
   bool operator <(Tpoint3D& a, Tpoint3D& b);

   /// overloaded <= operator for points
   bool operator <=(Tpoint3D& a, Tpoint3D& b);

   /// Calculates distance between two points
   double Distance(Tpoint3D& a, Tpoint3D& b);

   /// Calculates square of distance between two points
   double SquareDistance(Tpoint3D& a, Tpoint3D& b);

   /// Calculates maximal coordinates from two points
   Tpoint3D pmax(Tpoint3D& a, Tpoint3D& b);

   /// Calculates minimal coordinates from two points
   Tpoint3D pmin(Tpoint3D& a, Tpoint3D& b);

   /// Generates an initial point for randomized ray tracing
   Tpoint3D RandomRay(double ray);

   /// Tests if a given 2D point is inside a polygon 
   bool IsInsidePolygon(Tpoint3D& pt, Tpoint3D polygon[], const unsigned int np);

   /// overloaded output operator for points
   ostream& operator <<(ostream& os, Tpoint3D& pt);

   /// overloaded input operator for points
   istream& operator>>(istream& is, Tpoint3D& pt);

   /**
    * Point3D is a base class for geometric type classes.
    * All members are public for fast computational access.
    */
   class Point3D
   {
   public:

      /// x, y, z coordinates of the points
      Tpoint3D	pt;

      /// Indicates whether the point is calculated from intersections or not
      bool isCalculated;

      /// Default constructor
      Point3D(Tpoint3D point = Tpoint3D(), bool cal = false);

      /// Copy constructor
      Point3D(const Point3D& point);

      /// virtual destructor
      virtual ~Point3D();

      // overloaded output operator for 3D points 
      friend ostream& operator <<(ostream& os, Point3D* point);

      /// overloaded input operator for 3D points
      friend istream& operator >>(istream& is, Point3D* point);
   };
}

#endif // HAPP_POINT
