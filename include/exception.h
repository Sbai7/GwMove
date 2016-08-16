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

#ifndef _EXCEPTION_H
#define _EXCEPTION_H

// C style headers
#include <stdio.h>


/**
 * Base class for all error handling or exception classes.
 */
class Error 
{
private: 
   /// submitted error message.
   const char *msg; 

   /// indicate the filename where the error has occured.
   char *filename; 

   /// indicate the line number where the error has occured.
   int linenum;

public: 
   /**
    * Default constructor.
    */
   Error()
   {
      msg = NULL; 
      filename = NULL; 
      linenum = 0; 
   };

   /**
    * Another constructor.
    * \param *m pointer to the null etrminated error message string.
    * \param *fn pointer to the filename.
    * \param ln a line number variable.
    */
   Error(const char *m,char *fn = NULL,int ln = 0):msg(m),filename(fn),linenum(ln)
   {
      // Nothing to do here.
   };

   /**
    * Class destructor.
    */
   ~Error()
   {
      linenum = 0;
   };

   /** 
    * Write the error to the debugging console 
    * and return the error message.
    * \return pointer to the error message.
    */
   const char *message()
   {
      if (msg != NULL) {
#ifdef _WINDOWS
         TRACE0(msg);
#endif // _WINDOWS
         return msg;
      }
      else {
#ifdef _WINDOWS
         TRACE0("Excepts error thrown\n");
#endif // _WINDOWS
         return "Excepts error thrown"; 
      }
   };

   /** 
    * Return full report with filename and line number.
    * \return pointer to the report null terminated string.
    */
   char *report()
   {
      // A buffer for the error message. It is dangerous to have a
      // fixed length buffer for storing the error message!
      static char str[1024];
        
      if (linenum == 0) {
         if (filename == NULL) {
#ifdef _WINDOWS
            TRACE0("Excepts error thrown\n");
#endif // _WINDOWS
            return "Excepts error thrown"; 
         }
         else
            sprintf(str,"%s in %s\n",msg,filename); 
      }
      else {  
         if (filename == NULL) {
#ifdef _WINDOWS
            TRACE0("Excepts error thrown\n");
#endif //  _WINDOWS
            return "Excepts error thrown"; 
         }
         else 
            sprintf(str,"%s in %s at line %d\n", msg, filename, linenum);
      }

#ifdef _WINDOWS
      TRACE0(str);
#endif // _WINDOWS

      return str;
   };
};


/**
 * Input/Output error handling class.
 */
class IO_Error : public Error
{
public:
   IO_Error();
   IO_Error(char *m,char *fn = NULL,int ln = 0) : Error(m,fn,ln){  };
 };


/**
 * Message Passing error handling class.
 */
class MP_Error : public Error
{ 
public:
   MP_Error();
   MP_Error(char *m,char *fn = NULL,int ln = 0) : Error(m,fn,ln){  };
};


/**
 * Memory error handling class.
 */
class Mem_Error : public Error
{
public:
   Mem_Error();
   Mem_Error(char *m,char *fn = NULL,int ln = 0) : Error(m,fn,ln){  };
};


/**
 * Array bounds error handling class.
 */
class Bounds_Error : public Error 
{
public:
   Bounds_Error();
   Bounds_Error(char *m,char *fn = NULL,int ln = 0) : Error(m,fn,ln){  };
};

#endif // #ifndef _EXCEPTION_H

