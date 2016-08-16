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

#ifndef HAPP_PLATFORM
#define HAPP_PLATFORM


namespace happ
{
   /**
    * Class for calculation of platform dependent constants.
    */
   class PlatformConstants
   {
   protected:
      /**
       * Calculates the smallest possible floating point number.
       * Adapated from the book: "Scientific C++: building numerical
       * libraries the object-oriented way", by G. Buzzi-Ferraris (1993).
       */
      double GetPlatformEpsilon() {
         double epsilon = 1.;
         while ((1. + epsilon) != 1.) epsilon /= 2.;
         epsilon *= 2.;
         return epsilon;
      };

   public:
      /// Machine's smallest floating point value
      double epsilon;

      /// Default constructor
      PlatformConstants() {
         epsilon = GetPlatformEpsilon();
      };

      /// Destructor
      virtual ~PlatformConstants() {};

   };
}

#endif // HAPP_PLATFORM
