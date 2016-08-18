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

// Abstract array data type

#include <array.h>

namespace happ
{

   BaseArray::BaseArray(int asize, int ainc, int elementsize)
   {
      if (asize > 0)
      {
         data = new char[asize * elementsize];
         size = allocsize = asize;
      }
      else
      {
         data = 0;
         size = allocsize = 0;
      }
      inc = ainc;
   }

   BaseArray::~BaseArray()
   {
      if (allocsize > 0)
      {
         delete[](char*)data;
      }
   }

   void BaseArray::GrowSize(int minsize, int elementsize)
   {
      void *p;
      int nsize = (inc > 0) ? abs(allocsize) + inc : 2 * abs(allocsize);
      if (nsize < minsize) { nsize = minsize; }

      p = new char[nsize * elementsize];
      if (size > 0)
      {
         memcpy(p, data, size * elementsize);
      }
      if (allocsize > 0)
      {
         delete[](char*)data;
      }
      data = p;
      allocsize = nsize;
   }

   template <class T>
   void Array<T>::Print(std::ostream &out, int width)
   {
      for (int i = 0; i < size; i++)
      {
         out << ((T*)data)[i];
         if (!((i + 1) % width) || i + 1 == size)
         {
            out << '\n';
         }
         else
         {
            out << " ";
         }
      }
   }

   template <class T>
   void Array<T>::Save(std::ostream &out)
   {
      out << size << '\n';
      for (int i = 0; i < size; i++)
      {
         out << operator[](i) << '\n';
      }
   }

   template <class T>
   T Array<T>::Max() const
   {
      HAPP_ASSERT(size > 0, "Array is empty with size " << size);

      T max = operator[](0);
      for (int i = 1; i < size; i++)
         if (max < operator[](i))
         {
            max = operator[](i);
         }

      return max;
   }

   template <class T>
   T Array<T>::Min() const
   {
      HAPP_ASSERT(size > 0, "Array is empty with size " << size);

      T min = operator[](0);
      for (int i = 1; i < size; i++)
         if (operator[](i) < min)
         {
            min = operator[](i);
         }

      return min;
   }

   // Partial Sum
   template <class T>
   void Array<T>::PartialSum()
   {
      T sum = static_cast<T>(0);
      for (int i = 0; i < size; i++)
      {
         sum += operator[](i);
         operator[](i) = sum;
      }
   }

   // Sum
   template <class T>
   T Array<T>::Sum()
   {
      T sum = static_cast<T>(0);
      for (int i = 0; i < size; i++)
      {
         sum += operator[](i);
      }

      return sum;
   }

   template <class T>
   int Array<T>::IsSorted()
   {
      T val_prev = operator[](0), val;
      for (int i = 1; i < size; i++)
      {
         val = operator[](i);
         if (val < val_prev)
         {
            return 0;
         }
         val_prev = val;
      }

      return 1;
   }

   template class Array<int>;
   template class Array<double>;

}
