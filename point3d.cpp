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
#include <point3d.h>
#include <algorithm>

// C++ headers
#include <ctime>
#include <cstdlib>


namespace happ
{
   Tpoint3D operator+(Tpoint3D& a, Tpoint3D& b)
   {
      Tpoint3D pt;

      pt.x = a.x + b.x;
      pt.y = a.y + b.y;
      pt.z = a.z + b.z;

      return pt;
   }

   Tpoint3D operator-(Tpoint3D& a, Tpoint3D& b)
   {
      Tpoint3D pt;

      pt.x = a.x - b.x;
      pt.y = a.y - b.y;
      pt.z = a.z - b.z;

      return pt;
   }

   Tpoint3D operator/(Tpoint3D& a, double b)
   {
      Tpoint3D pt;
      if (b) {
         pt.x = a.x / b;
         pt.y = a.y / b;
         pt.z = a.z / b;
      }
      return pt;
   }

   Tpoint3D operator*(Tpoint3D& a, double b)
   {
      Tpoint3D pt;
      pt.x = a.x*b;
      pt.y = a.y*b;
      pt.z = a.z*b;

      return pt;
   }

   bool operator ==(Tpoint3D& a, Tpoint3D& b)
   {
      PlatformConstants platfm;

      if (Distance(a, b) < platfm.epsilon)
         return true;
      return false;
   }

   bool operator < (Tpoint3D& a, Tpoint3D& b)
   {
      // first tests for equality
      if (a == b) return false;

      if (a.x < b.x && a.y < b.y && a.z < b.z)
         return true;

      return false;
   }

   bool operator <=(Tpoint3D& a, Tpoint3D& b)
   {
      if (a == b) return true;

      if (a.x <= b.x && a.y <= b.y && a.z <= b.z)
         return true;

      return false;
   }

   double Distance(Tpoint3D& a, Tpoint3D& b)
   {
      double dist;

      dist = sqrt(SquareDistance(a, b));
      return dist;
   }

   double SquareDistance(Tpoint3D& a, Tpoint3D& b)
   {
      double d1, d2, d3, dist;

      d1 = a.x - b.x;
      d2 = a.y - b.y;
      d3 = a.z - b.z;
      dist = d1*d1 + d2*d2 + d3*d3;

      return dist;
   }

   Tpoint3D pmax(Tpoint3D& a, Tpoint3D& b)
   {
      Tpoint3D pt;

      pt.x = max(a.x, b.x);
      pt.y = max(a.y, b.y);
      pt.z = max(a.z, b.z);

      return pt;
   }

   Tpoint3D pmin(Tpoint3D& a, Tpoint3D& b)
   {
      Tpoint3D pt;

      pt.x = min(a.x, b.x);
      pt.y = min(a.y, b.y);
      pt.z = min(a.z, b.z);

      return pt;
   }

   Tpoint3D RandomRay(double r)
   {
      srand((unsigned)time(NULL));

      /*
       * Generate a random point on a sphere of radius 1.
       * the sphere is sliced at z, and a random point at angle t
       * generated on the circle of intersection.
       */
      double z = 2.0 *(double)rand() / RAND_MAX - 1.0,
         t = 2.0 * M_PI *(double)rand() / RAND_MAX,
         w = sqrt(1. - z*z),
         x = r * w * cos(t),
         y = r * w * sin(t);

      z *= r;

      return Tpoint3D(x, y, z);
   }

   bool IsInsidePolygon(Tpoint3D& pt, Tpoint3D polygon[], const unsigned int np)
   {
      unsigned int nIntersect = 0; // number of intersections
      unsigned int i;
      Tpoint3D pt1, pt2;
      double x;

      pt1 = polygon[0] - pt;
      for (i = 1; i < np; i++) {
         pt2 = polygon[0] - pt;
         if (fabs(pt1.x) < pointPrecision && fabs(pt1.y) < pointPrecision) return true;
         if ((pt1.y > 0.0) != (pt2.y > 0.0)) {
            x = (pt1.x*pt2.y - pt2.x*pt1.y) / (pt2.y - pt1.y);
            if (x > 0.0) nIntersect++;
         }
         pt1 = pt2;
      }

      if ((nIntersect % 2) == 1) return true;
      else return false;
   }

   ostream& operator <<(ostream& os, Tpoint3D& pt)
   {
      os << setiosflags(ios::fixed) << setprecision(8) << setiosflags(ios::showpoint);
      os << setw(20) << pt.x << setw(20) << pt.y << setw(20) << pt.z;
      return os;
   }

   istream& operator>>(istream& is, Tpoint3D& pt)
   {
      is >> pt.x >> pt.y >> pt.z;
      return is;
   }

   Point3D::Point3D(Tpoint3D point, bool cal)
   {
      pt = point;
      isCalculated = cal;
   }

   Point3D::Point3D(const Point3D& point)
   {
      pt = point.pt;
      isCalculated = point.isCalculated;
   }

   Point3D::~Point3D()
   {
   }

   ostream& operator <<(ostream& os, Point3D* point)
   {
      os << point->pt << " " << point->isCalculated;
      return os;
   }

   istream& operator >>(istream& is, Point3D* point)
   {
      bool cal;
      is >> point->pt >> cal;
      point->isCalculated = cal;

      return is;
   }
}