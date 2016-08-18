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

#ifndef HydroApp_ERROR
#define HydroApp_ERROR

//comm #include "../config/config.hpp"
#include <iomanip>
#include <sstream>

namespace happ
{

   void happ_error(const char *msg = NULL);

   void happ_warning(const char *msg = NULL);

}

#ifndef _HAPP_FUNC_NAME
#ifndef _WIN32
// This is nice because it shows the class and method name
#define _HAPP_FUNC_NAME __PRETTY_FUNCTION__
// This one is C99 standard.
//#define _HAPP_FUNC_NAME __func__
#else
// for Visual Studio C++
#define _HAPP_FUNC_NAME __FUNCSIG__
#endif
#endif

// Common error message and abort macro
#define _HAPP_MESSAGE(msg, warn)                                        \
   {                                                                    \
      std::ostringstream happMsgStream;                                 \
      happMsgStream << std::setprecision(16);                           \
      happMsgStream << std::setiosflags(std::ios_base::scientific);     \
      happMsgStream << msg << '\n';                                     \
      happMsgStream << " ... at line " << __LINE__;                     \
      happMsgStream << " in " << _HAPP_FUNC_NAME << " of file ";        \
      happMsgStream << __FILE__ << ".";                                 \
      happMsgStream << std::ends;                                       \
      if (!(warn))                                                      \
         happ::happ_error(happMsgStream.str().c_str());                 \
            else                                                              \
         happ::happ_warning(happMsgStream.str().c_str());               \
   }

// Outputs lots of useful information and aborts.
// For all of these functions, "msg" is pushed to an ostream, so you can
// write useful (if complicated) error messages instead of writing
// out to the screen first, then calling abort.  For example:
// HAPP_ABORT( "Unknown geometry type: " << type );
#define HAPP_ABORT(msg) _HAPP_MESSAGE("HAPP abort: " << msg, 0)

// Does a check, and then outputs lots of useful information if the test fails
#define HAPP_VERIFY(x, msg)                             \
   if (!(x))                                            \
      {                                                    \
      _HAPP_MESSAGE("Verification failed: ("            \
                    << #x << ") is false: " << msg, 0); \
      }

// Use this if the only place your variable is used is in ASSERTs
// For example, this code snippet:
//   int err = MPI_Reduce(ldata, maxdata, 5, MPI_INT, MPI_MAX, 0, MyComm);
//   HAPP_CONTRACT_VAR(err);
//   HAPP_ASSERT( err == 0, "MPI_Reduce gave an error with length "
//                       << ldata );
#define HAPP_CONTRACT_VAR(x) if (0 && &x == &x){}

// Now set up some optional checks, but only if the right flags are on
#ifdef HAPP_DEBUG

#define HAPP_ASSERT(x, msg)                             \
   if (!(x))                                            \
      {                                                    \
      _HAPP_MESSAGE("Assertion failed: ("               \
                    << #x << ") is false: " << msg, 0); \
      }

#else

// Get rid of all this code, since we're not checking.
#define HAPP_ASSERT(x, msg)

#endif

// Generate a warning message - always generated, regardless of MFEM_DEBUG.
#define HAPP_WARNING(msg) _HAPP_MESSAGE("HAPP Warning: " << msg, 1)

#endif // #ifndef HydroApp_ERROR
