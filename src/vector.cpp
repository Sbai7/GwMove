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

/// Implementation of data type vector

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>

#include <vector.h>

namespace happ
{

   Vector::Vector(const Vector &v)
   {
      int s = v.Size();
      if (s > 0)
      {
         allocsize = size = s;
         data = new double[s];
         for (int i = 0; i < s; i++)
         {
            data[i] = v(i);
         }
      }
      else
      {
         allocsize = size = 0;
         data = NULL;
      }
   }

   void Vector::Load(std::istream **in, int np, int *dim)
   {
      int i, j, s;

      s = 0;
      for (i = 0; i < np; i++)
      {
         s += dim[i];
      }

      SetSize(s);

      int p = 0;
      for (i = 0; i < np; i++)
         for (j = 0; j < dim[i]; j++)
         {
            *in[i] >> data[p++];
         }
   }

   void Vector::Load(std::istream &in, int Size)
   {
      SetSize(Size);

      for (int i = 0; i < size; i++)
      {
         in >> data[i];
      }
   }

   double &Vector::Elem(int i)
   {
      return operator()(i);
   }

   const double &Vector::Elem(int i) const
   {
      return operator()(i);
   }

   double Vector::operator*(const double *v) const
   {
      int s = size;
      const double *d = data;
      double prod = 0.0;
#ifdef HAPP_USE_OPENMP
#pragma omp parallel for reduction(+:prod)
#endif
      for (int i = 0; i < s; i++)
      {
         prod += d[i] * v[i];
      }
      return prod;
   }

   double Vector::operator*(const Vector &v) const
   {
#ifdef HAPP_DEBUG
      if (v.size != size)
      {
         happ_error("Vector::operator*(const Vector &) const");
      }
#endif

      return operator*(v.data);
   }

   Vector &Vector::operator=(const double *v)
   {
      for (int i = 0; i < size; i++)
      {
         data[i] = v[i];
      }
      return *this;
   }

   Vector &Vector::operator=(const Vector &v)
   {
      SetSize(v.Size());
      for (int i = 0; i < size; i++)
      {
         data[i] = v.data[i];
      }
      return *this;
   }

   Vector &Vector::operator=(double value)
   {
      int i, s = size;
      double *p = data, v = value;
      for (i = 0; i < s; i++)
      {
         *(p++) = v;
      }
      return *this;
   }

   Vector &Vector::operator*=(double c)
   {
      for (int i = 0; i < size; i++)
      {
         data[i] *= c;
      }
      return *this;
   }

   Vector &Vector::operator/=(double c)
   {
      double m = 1.0 / c;
      for (int i = 0; i < size; i++)
      {
         data[i] *= m;
      }
      return *this;
   }

   Vector &Vector::operator-=(double c)
   {
      for (int i = 0; i < size; i++)
      {
         data[i] -= c;
      }
      return *this;
   }

   Vector &Vector::operator-=(const Vector &v)
   {
#ifdef HAPP_DEBUG
      if (size != v.size)
      {
         happ_error("Vector::operator-=(const Vector &)");
      }
#endif
      for (int i = 0; i < size; i++)
      {
         data[i] -= v(i);
      }
      return *this;
   }

   Vector &Vector::operator+=(const Vector &v)
   {
#ifdef HAPP_DEBUG
      if (size != v.size)
      {
         happ_error("Vector::operator+=(const Vector &)");
      }
#endif
      for (int i = 0; i < size; i++)
      {
         data[i] += v(i);
      }
      return *this;
   }

   Vector &Vector::Add(const double a, const Vector &Va)
   {
#ifdef HAPP_DEBUG
      if (size != Va.size)
      {
         happ_error("Vector::Add(const double, const Vector &)");
      }
#endif
      if (a != 0.0)
      {
         for (int i = 0; i < size; i++)
         {
            data[i] += a * Va(i);
         }
      }
      return *this;
   }

   Vector &Vector::Set(const double a, const Vector &Va)
   {
#ifdef HAPP_DEBUG
      if (size != Va.size)
      {
         happ_error("Vector::Set(const double, const Vector &)");
      }
#endif
      for (int i = 0; i < size; i++)
      {
         data[i] = a * Va(i);
      }
      return *this;
   }

   void Vector::SetVector(const Vector &v, int offset)
   {
      int vs = v.Size();
      double *vp = v.data, *p = data + offset;

#ifdef HAPP_DEBUG
      if (offset + vs > size)
      {
         happ_error("Vector::SetVector(const Vector &, int)");
      }
#endif

      for (int i = 0; i < vs; i++)
      {
         p[i] = vp[i];
      }
   }

   void Vector::Neg()
   {
      for (int i = 0; i < size; i++)
      {
         data[i] = -data[i];
      }
   }

   void add(const Vector &v1, const Vector &v2, Vector &v)
   {
#ifdef HAPP_DEBUG
      if (v.size != v1.size || v.size != v2.size)
      {
         happ_error("add(Vector &v1, Vector &v2, Vector &v)");
      }
#endif

#ifdef HAPP_USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < v.size; i++)
      {
         v.data[i] = v1.data[i] + v2.data[i];
      }
   }

   void add(const Vector &v1, double alpha, const Vector &v2, Vector &v)
   {
#ifdef HAPP_DEBUG
      if (v.size != v1.size || v.size != v2.size)
      {
         happ_error("add(Vector &v1, double alpha, Vector &v2, Vector &v)");
      }
#endif
      if (alpha == 0.0)
      {
         v = v1;
      }
      else if (alpha == 1.0)
      {
         add(v1, v2, v);
      }
      else
      {
         const double *v1p = v1.data, *v2p = v2.data;
         double *vp = v.data;
         int s = v.size;
#ifdef HAPP_USE_OPENMP
#pragma omp parallel for
#endif
         for (int i = 0; i < s; i++)
         {
            vp[i] = v1p[i] + alpha*v2p[i];
         }
      }
   }

   void add(const double a, const Vector &x, const Vector &y, Vector &z)
   {
#ifdef HAPP_DEBUG
      if (x.size != y.size || x.size != z.size)
         happ_error("add(const double a, const Vector &x, const Vector &y,"
         " Vector &z)");
#endif
      if (a == 0.0)
      {
         z = 0.0;
      }
      else if (a == 1.0)
      {
         add(x, y, z);
      }
      else
      {
         const double *xp = x.data;
         const double *yp = y.data;
         double       *zp = z.data;
         int            s = x.size;

#ifdef HAPP_USE_OPENMP
#pragma omp parallel for
#endif
         for (int i = 0; i < s; i++)
         {
            zp[i] = a * (xp[i] + yp[i]);
         }
      }
   }

   void add(const double a, const Vector &x,
      const double b, const Vector &y, Vector &z)
   {
#ifdef HAPP_DEBUG
      if (x.size != y.size || x.size != z.size)
         happ_error("add(const double a, const Vector &x,\n"
         "    const double b, const Vector &y, Vector &z)");
#endif
      if (a == 0.0)
      {
         z.Set(b, y);
      }
      else if (b == 0.0)
      {
         z.Set(a, x);
      }
      else if (a == 1.0)
      {
         add(x, b, y, z);
      }
      else if (b == 1.0)
      {
         add(y, a, x, z);
      }
      else if (a == b)
      {
         add(a, x, y, z);
      }
      else
      {
         const double *xp = x.data;
         const double *yp = y.data;
         double       *zp = z.data;
         int            s = x.size;

#ifdef HAPP_USE_OPENMP
#pragma omp parallel for
#endif
         for (int i = 0; i < s; i++)
         {
            zp[i] = a * xp[i] + b * yp[i];
         }
      }
   }

   void subtract(const Vector &x, const Vector &y, Vector &z)
   {
#ifdef HAPP_DEBUG
      if (x.size != y.size || x.size != z.size)
      {
         happ_error("subtract(const Vector &, const Vector &, Vector &)");
      }
#endif
      const double *xp = x.data;
      const double *yp = y.data;
      double       *zp = z.data;
      int            s = x.size;

#ifdef HAPP_USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < s; i++)
      {
         zp[i] = xp[i] - yp[i];
      }
   }

   void subtract(const double a, const Vector &x, const Vector &y, Vector &z)
   {
#ifdef HAPP_DEBUG
      if (x.size != y.size || x.size != z.size)
         happ_error("subtract(const double a, const Vector &x,"
         " const Vector &y, Vector &z)");
#endif

      if (a == 0.)
      {
         z = 0.;
      }
      else if (a == 1.)
      {
         subtract(x, y, z);
      }
      else
      {
         const double *xp = x.data;
         const double *yp = y.data;
         double       *zp = z.data;
         int            s = x.size;

#ifdef HAPP_USE_OPENMP
#pragma omp parallel for
#endif
         for (int i = 0; i < s; i++)
         {
            zp[i] = a * (xp[i] - yp[i]);
         }
      }
   }

   void Vector::median(const Vector &lo, const Vector &hi)
   {
      double *v = data;

      for (int i = 0; i < size; i++)
      {
         if (v[i] < lo[i])
         {
            v[i] = lo[i];
         }
         else if (v[i] > hi[i])
         {
            v[i] = hi[i];
         }
      }
   }

   void Vector::GetSubVector(const Array<int> &dofs, Vector &elemvect) const
   {
      int i, j, n = dofs.Size();

      elemvect.SetSize(n);

      for (i = 0; i < n; i++)
      {
         if ((j = dofs[i]) >= 0)
         {
            elemvect(i) = data[j];
         }
         else
         {
            elemvect(i) = -data[-1 - j];
         }
      }
   }

   void Vector::GetSubVector(const Array<int> &dofs, double *elem_data) const
   {
      int i, j, n = dofs.Size();

      for (i = 0; i < n; i++)
      {
         if ((j = dofs[i]) >= 0)
         {
            elem_data[i] = data[j];
         }
         else
         {
            elem_data[i] = -data[-1 - j];
         }
      }
   }

   void Vector::SetSubVector(const Array<int> &dofs, const Vector &elemvect)
   {
      int i, j, n = dofs.Size();

      for (i = 0; i < n; i++)
      {
         if ((j = dofs[i]) >= 0)
         {
            data[j] = elemvect(i);
         }
         else
         {
            data[-1 - j] = -elemvect(i);
         }
      }
   }

   void Vector::SetSubVector(const Array<int> &dofs, double *elem_data)
   {
      int i, j, n = dofs.Size();

      for (i = 0; i < n; i++)
      {
         if ((j = dofs[i]) >= 0)
         {
            data[j] = elem_data[i];
         }
         else
         {
            data[-1 - j] = -elem_data[i];
         }
      }
   }

   void Vector::AddElementVector(const Array<int> &dofs, const Vector &elemvect)
   {
      int i, j, n = dofs.Size();

      for (i = 0; i < n; i++)
         if ((j = dofs[i]) >= 0)
         {
            data[j] += elemvect(i);
         }
         else
         {
            data[-1 - j] -= elemvect(i);
         }
   }

   void Vector::AddElementVector(const Array<int> &dofs, double *elem_data)
   {
      int i, j, n = dofs.Size();

      for (i = 0; i < n; i++)
      {
         if ((j = dofs[i]) >= 0)
         {
            data[j] += elem_data[i];
         }
         else
         {
            data[-1 - j] -= elem_data[i];
         }
      }
   }

   void Vector::AddElementVector(const Array<int> &dofs, const double a,
      const Vector &elemvect)
   {
      int i, j, n = dofs.Size();

      for (i = 0; i < n; i++)
         if ((j = dofs[i]) >= 0)
         {
            data[j] += a * elemvect(i);
         }
         else
         {
            data[-1 - j] -= a * elemvect(i);
         }
   }

   void Vector::SetSubVectorComplement(const Array<int> &dofs, const double val)
   {
      Vector dofs_vals;
      GetSubVector(dofs, dofs_vals);
      operator=(val);
      SetSubVector(dofs, dofs_vals);
   }

   void Vector::Print(std::ostream &out, int width) const
   {
      if (!size) { return; }

      for (int i = 0; 1;)
      {
         out << data[i];
         i++;
         if (i == size)
         {
            break;
         }
         if (i % width == 0)
         {
            out << '\n';
         }
         else
         {
            out << ' ';
         }
      }
      out << '\n';
   }

   void Vector::Print_HYPRE(std::ostream &out) const
   {
      int i;
      std::ios::fmtflags old_fmt = out.flags();
      out.setf(std::ios::scientific);
      std::streamsize old_prec = out.precision(14);

      out << size << '\n';  // number of rows

      for (i = 0; i < size; i++)
      {
         out << data[i] << '\n';
      }

      out.precision(old_prec);
      out.flags(old_fmt);
   }

   void Vector::Randomize(int seed)
   {
      // static unsigned int seed = time(0);
      const double max = (double)(RAND_MAX)+1.;

      if (seed == 0)
      {
         seed = (int)time(0);
      }

      // srand(seed++);
      srand((unsigned)seed);

      for (int i = 0; i < size; i++)
      {
         data[i] = fabs(rand() / max);
      }
   }

   double Vector::Norml2() const
   {
      return sqrt((*this)*(*this));
   }

   double Vector::Normlinf() const
   {
      double max = fabs(data[0]);

      for (int i = 1; i < size; i++)
         if (fabs(data[i]) > max)
         {
            max = fabs(data[i]);
         }

      return max;
   }

   double Vector::Norml1() const
   {
      double sum = 0.0;

      for (int i = 0; i < size; i++)
      {
         sum += fabs(data[i]);
      }

      return sum;
   }

   double Vector::Normlp(double p) const
   {
      HAPP_ASSERT(p > 0.0, "Vector::Normlp");
      if (p == 1.0)
      {
         return Norml1();
      }
      if (p == 2.0)
      {
         return Norml2();
      }
      if (p < std::numeric_limits<double>::infinity())
      {
         double sum = 0.0;
         for (int i = 0; i < size; i++)
         {
            sum += pow(fabs(data[i]), p);
         }
         return pow(sum, 1.0 / p);
      }
      else
      {
         return Normlinf();
      }
   }

   double Vector::Max() const
   {
      double max = data[0];

      for (int i = 1; i < size; i++)
         if (data[i] > max)
         {
            max = data[i];
         }

      return max;
   }

   double Vector::Min() const
   {
      double min = data[0];

      for (int i = 1; i < size; i++)
         if (data[i] < min)
         {
            min = data[i];
         }

      return min;
   }

   double Vector::Sum() const
   {
      double sum = 0.0;

      for (int i = 0; i < size; i++)
         sum += data[i];

      return sum;
   }

   double Vector::DistanceTo(const double *p) const
   {
      return Distance(data, p, size);
   }

   double Vector::Mean() const
   {
      return Sum()/size; 
   }

   double Vector::Variance() const 
   {
      double mu = Mean();
      double sum = 0.0; 

      for (int i = 0; i < size; i++)
         sum += (data[i] - mu)*(data[i] - mu);

      return sum/size;
   }

   double Vector::StandardDev() const
   {
      return sqrt(Variance());
   }

   double Vector::StandardDev2() const
   {
      double mu = Mean();
      double sum = 0.0;
      unsigned int i;

      for (i = 0; i < size; i++)
         sum += (data[i] - mu)*(data[i] - mu);

      return sqrt(sum / (size-1));
   }

   double Vector::Skewness() const
   {
      double mu = Mean();
      double sum = 0.0;

      for (int i = 0; i < size; i++)
         sum += std::pow((data[i] - mu),3);

      sum /= (size * std::pow(StandardDev(),3) );

      return sum;
   }

   double Vector::Skewness2() const
   {
      double mu = Mean();
      double sum = 0.0;

      for (int i = 0; i < size; i++)
         sum += std::pow((data[i] - mu), 3);

      sum /= ((size-1) * std::pow(StandardDev2(),3) );

      return sum;
   }

   double Vector::Kurtosis() const
   {
      double mu = Mean();
      double sum = 0.0;

      for (int i = 0; i < size; i++)
         sum += std::pow((data[i] - mu), 4);

      sum /= (size * std::pow(StandardDev(), 4));

      return sum;
   }

   double Vector::Kurtosis2() const
   {
      double mu = Mean();
      double sum = 0.0;

      for (int i = 0; i < size; i++)
         sum += std::pow((data[i] - mu), 4);

      sum /= ((size-1) * std::pow(StandardDev2(), 4));

      return sum;
   }

   void Vector::Percentiles(double values[11]) const
   {
      unsigned int i,j; 

      unsigned int* count = new unsigned int[size];
      double* quant = new double[size]; 

      for (i = 0; i < size; i++) count[i] = 0; 
      for (i = 0; i < size; i++) quant[i] = 0.0;

      // compute quantiles of each entry  
      for (i = 0; i < size; i++)
      {
         for (j = 0; j < size; j++)
         {
            if (i == j) continue; 
            if (data[i] > data[j]) count[i]++; 
         }
         quant[i]  = (double)count[i] / (size - 1);
         quant[i] *= 100.0;
      }

      // search for each percentile 
      double percent[11] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0};
      for (i = 0; i < 11; i++)
      {
         double eps = 100.0;
         double old_eps = eps;
         for (j = 0; j < size; j++)
         {
            eps = std::abs(quant[j] - percent[i]);
            if (eps < old_eps) {
               values[i] = data[j];
               old_eps = eps;
            }
         }
      }

      delete[] count; 
      delete[] quant; 
   }

   double Vector::Covariance(const Vector& y) const
   {
      if (size != y.size) return 0.0;

      double mu_x = this->Mean();
      double mu_y = y.Mean();

      double sum = 0.0;

      for (int i = 0; i < size; i++)
      {
         sum += (this->data[i] - mu_x)*(y.data[i] - mu_y);
      }

      return sum / size;
   }

   double Vector::Covariance2(const Vector& y) const
   {
      return size * Covariance(y) / (size - 1);
   }

   double Vector::Pearson(const Vector& y) const
   {
      // sum of x (this vector)
      double sx = this->Sum();

      // sum of y vector 
      double sy = y.Sum(); 

      // sum of x*y vector 
      double sxy = 0.0; 
      unsigned int i; 
      for (i = 0; i < size; i++)
      {
         sxy += data[i] * y.data[i];
      }

      // sum of x^2
      double sx2 = 0.0;
      for (i = 0; i < size; i++)
      {
         sx2 += data[i] * data[i];
      }

      // sum of y^2
      double sy2 = 0.0;
      for (i = 0; i < size; i++)
      {
         sy2 += y.data[i] * y.data[i];
      }

      // return Pearson correlation rank coefficient 
      double sp = (size * sxy) - (sx * sy);
      sp /= std::sqrt( (size * sx2) - (sx * sx) );
      sp /= std::sqrt( (size * sy2) - (sy * sy) );
      return sp;
   }
}
