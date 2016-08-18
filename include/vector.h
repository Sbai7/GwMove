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

#ifndef HAPP_VECTOR
#define HAPP_VECTOR

// Data type vector

#include <array.h>

#include <cmath>
#include <iostream>
#if defined(_MSC_VER) && (_MSC_VER < 1800)
#include <float.h>
#define isfinite _finite
#endif

namespace happ
{

   /** Count the number of entries in an array of doubles for which isfinite
   is false, i.e. the entry is a NaN or +/-Inf. */
   inline int CheckFinite(const double *v, const int n);

   /// Vector data type.
   class Vector
   {
   protected:

      int size, allocsize;
      double * data;

   public:

      /// Default constructor for Vector. Sets size = 0 and data = NULL
      Vector() { allocsize = size = 0; data = 0; }

      /// Copy constructor
      Vector(const Vector &);

      /// Creates vector of size s.
      explicit Vector(int s);

      /// Creates a vector referencing an array of doubles, owned by someone else.
      Vector(double *_data, int _size)
      {
         data = _data; size = _size; allocsize = -size;
      }

      /// Reads a vector from multiple files
      void Load(std::istream ** in, int np, int * dim);

      /// Load a vector from an input stream.
      void Load(std::istream &in, int Size);

      /// Load a vector from an input stream.
      void Load(std::istream &in) { int s; in >> s; Load(in, s); }

      /// Resizes the vector if the new size is different
      void SetSize(int s);

      void SetData(double *d) { data = d; }

      void SetDataAndSize(double *d, int s)
      {
         data = d; size = s; allocsize = -s;
      }

      void NewDataAndSize(double *d, int s)
      {
         if (allocsize > 0) { delete[] data; }
         SetDataAndSize(d, s);
      }

      void MakeDataOwner() { allocsize = abs(allocsize); }

      /// Destroy a vector
      void Destroy();

      /// Returns the size of the vector.
      inline int Size() const { return size; }

      // double *GetData() { return data; }

      inline double *GetData() const { return data; }

      inline operator double *() { return data; }

      inline operator const double *() const { return data; }

      inline bool OwnsData() const { return (allocsize > 0); }

      /// Changes the ownership of the data; after the call the Vector is empty
      inline void StealData(double **p)
      {
         *p = data; data = 0; size = allocsize = 0;
      }

      /// Changes the ownership of the data; after the call the Vector is empty
      inline double *StealData() { double *p; StealData(&p); return p; }

      /// Sets value in vector. Index i = 0 .. size-1
      double & Elem(int i);

      /// Sets value in vector. Index i = 0 .. size-1
      const double & Elem(int i) const;

      /// Sets value in vector. Index i = 0 .. size-1
      inline double & operator() (int i);

      /// Sets value in vector. Index i = 0 .. size-1
      inline const double & operator() (int i) const;

      double operator*(const double *) const;

      /// Return the inner-product.
      double operator*(const Vector &v) const;

      Vector & operator=(const double *v);

      /// Redefine '=' for vector = vector.
      Vector & operator=(const Vector &v);

      /// Redefine '=' for vector = constant.
      Vector & operator=(double value);

      Vector & operator*=(double c);

      Vector & operator/=(double c);

      Vector & operator-=(double c);

      Vector & operator-=(const Vector &v);

      Vector & operator+=(const Vector &v);

      /// (*this) += a * Va
      Vector & Add(const double a, const Vector &Va);

      /// (*this) = a * x
      Vector & Set(const double a, const Vector &x);

      void SetVector(const Vector &v, int offset);

      /// (*this) = -(*this)
      void Neg();

      /// Swap the contents of two Vectors
      inline void Swap(Vector &other);

      /// Do v = v1 + v2.
      friend void add(const Vector &v1, const Vector &v2, Vector &v);

      /// Do v = v1 + alpha * v2.
      friend void add(const Vector &v1, double alpha, const Vector &v2, Vector &v);

      /// z = a * (x + y)
      friend void add(const double a, const Vector &x, const Vector &y, Vector &z);

      /// z = a * x + b * y
      friend void add(const double a, const Vector &x,
         const double b, const Vector &y, Vector &z);

      /// Do v = v1 - v2.
      friend void subtract(const Vector &v1, const Vector &v2, Vector &v);

      /// z = a * (x - y)
      friend void subtract(const double a, const Vector &x,
         const Vector &y, Vector &z);

      /// v = median(v,lo,hi) entrywise.  Implementation assumes lo <= hi.
      void median(const Vector &lo, const Vector &hi);

      void GetSubVector(const Array<int> &dofs, Vector &elemvect) const;
      void GetSubVector(const Array<int> &dofs, double *elem_data) const;

      void SetSubVector(const Array<int> &dofs, const Vector &elemvect);
      void SetSubVector(const Array<int> &dofs, double *elem_data);

      /// Add (element) subvector to the vector.
      void AddElementVector(const Array<int> & dofs, const Vector & elemvect);
      void AddElementVector(const Array<int> & dofs, double *elem_data);
      void AddElementVector(const Array<int> & dofs, const double a,
         const Vector & elemvect);

      /// Set all vector entries NOT in the 'dofs' array to the given 'val'.
      void SetSubVectorComplement(const Array<int> &dofs, const double val);

      /// Prints vector to stream out.
      void Print(std::ostream & out = std::cout, int width = 8) const;

      /// Prints vector to stream out in HYPRE_Vector format.
      void Print_HYPRE(std::ostream &out) const;

      /// Set random values in the vector.
      void Randomize(int seed = 0);
      /// Returns the l2 norm of the vector.
      double Norml2() const;
      /// Returns the l_infinity norm of the vector.
      double Normlinf() const;
      /// Returns the l_1 norm of the vector.
      double Norml1() const;
      /// Returns the l_p norm of the vector.
      double Normlp(double p) const;
      /// Returns the maximal element of the vector.
      double Max() const;
      /// Returns the minimal element of the vector.
      double Min() const;
      /// Return the sum of the vector entries
      double Sum() const;
      /// Compute the Euclidean distance to another vector.
      double DistanceTo(const double *p) const;

      /// Returns the mean of the vector. 
      double Mean() const; 
      /// Returns the variance of the vector. 
      double Variance() const;
      /// Returns the standard deviation of the vector.
      double StandardDev() const; 
      /// Returns Bessel's standard deviation of the Vector.
      double StandardDev2() const; 
      /// Returns skewness (third moment) of vector distribution.
      double Skewness() const; 
      /// Returns sample skewness of vector distribution.
      double Skewness2() const; 
      /// Returns kurtosis (fourth moment) of vector distribution. 
      double Kurtosis() const; 
      /// Returns sample kurtosis of vector distribution. 
      double Kurtosis2() const; 
      /// Returns percentiles of the vector. 
      void Percentiles(double values[11]) const; 
      /// Returns the covariance with another vector. 
      double Covariance(const Vector& y) const; 
      /// Returns sample covarinace with another vector. 
      double Covariance2(const Vector& y) const;
      /// Calculate Pearson's Correlation meaasure between two vectors. 
      double Pearson(const Vector& y) const;

      /** Count the number of entries in the Vector for which isfinite
      is false, i.e. the entry is a NaN or +/-Inf. */
      int CheckFinite() const { return happ::CheckFinite(data, size); }

      /// Destroys vector.
      virtual ~Vector();
   };

   // Inline methods

   inline int CheckFinite(const double *v, const int n)
   {
      // isfinite didn't appear in a standard until C99, and later C++11
      // It wasn't standard in C89 or C++98.  PGI as of 14.7 still defines
      // it as a macro, which sort of screws up everybody else.
      int bad = 0;
      for (int i = 0; i < n; i++)
      {
#ifdef isfinite
         if (!isfinite(v[i]))
#else
         if (!std::isfinite(v[i]))
#endif
         {
            bad++;
         }
      }
      return bad;
   }

   inline Vector::Vector(int s)
   {
      if (s > 0)
      {
         allocsize = size = s;
         data = new double[s];
      }
      else
      {
         allocsize = size = 0;
         data = NULL;
      }
   }

   inline void Vector::SetSize(int s)
   {
      if (s == size)
      {
         return;
      }
      if (s <= abs(allocsize))
      {
         size = s;
         return;
      }
      if (allocsize > 0)
      {
         delete[] data;
      }
      allocsize = size = s;
      data = new double[s];
   }

   inline void Vector::Destroy()
   {
      if (allocsize > 0)
      {
         delete[] data;
      }
      allocsize = size = 0;
      data = NULL;
   }

   inline double & Vector::operator() (int i)
   {
      HAPP_ASSERT(data && i >= 0 && i < size,
         "index [" << i << "] is out of range [0," << size << ")");

      return data[i];
   }

   inline const double & Vector::operator() (int i) const
   {
      HAPP_ASSERT(data && i >= 0 && i < size,
         "index [" << i << "] is out of range [0," << size << ")");

      return data[i];
   }

   inline void Vector::Swap(Vector &other)
   {
      happ::Swap(size, other.size);
      happ::Swap(allocsize, other.allocsize);
      happ::Swap(data, other.data);
   }

   /// Specialization of the template function Swap<> for class Vector
   template<> inline void Swap<Vector>(Vector &a, Vector &b)
   {
      a.Swap(b);
   }

   inline Vector::~Vector()
   {
      if (allocsize > 0)
      {
         delete[] data;
      }
   }

   inline double Distance(const double *x, const double *y, const int n)
   {
      using namespace std;
      double d = 0.0;

      for (int i = 0; i < n; i++)
      {
         d += (x[i] - y[i])*(x[i] - y[i]);
      }

      return sqrt(d);
   }

}

#endif // #ifndef HAPP_VECTOR
