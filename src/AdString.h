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

/*
 * Adil's CString like class.
 * For an easiear and portable CString objects
 *
 * Created 14 April 2004
 * Last update of 14 April 2004
 *
 * Author: A. Sbai, Ph.D. (a.sbai@brgm.fr)
 */

#ifndef HAPP_STRING
#define HAPP_STRING

// C headers
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <exception.h>

// C++ headers
#include <strstream>
#include <iostream>
using namespace std;

namespace happ
{
#ifndef NO_BOOL
#ifndef bool
#define bool int
#define true 1
#define false 0
#endif
#endif

#define PLUS_CHARGE  '+'
#define MINUS_CHARGE '-'


   /**
    * A high level, and very useful, String class.
    * This class mimic the nice MFC CString class, but this one
    * is portable, and is enhanced for processing of chemical strings.
    *
    * \author Adil Sbai, Ph.D. (a.sbai@brgm.fr)
    * \version 1.0
    */
   class String
   {
      /// name of species or item considered
      char* str;

      /// private copy of memcopy:
      char* copy(const char* s) {
         return (char*)((s == NULL) ? NULL : strcpy(new char[strlen(s) + 1], s));
      }

      /// private copy of strcat:
      char* cat(const char *s1, const char* s2) {
         if (s1 == NULL)
            return copy(s2);
         else if (s2 == NULL)
            return copy(s1);
         else {
            char *s = new char[strlen(s1) + strlen(s2) + 1];
            strcpy(s, s1);
            return strcat(s, s2);
         }
      }

      /// another version of our private copy of strcat:
      char* cat(const char *s1, const char c) {
         if (s1 == NULL) {
            if (c == '\0')
               return NULL;
            char *s = new char[2];
            s[0] = c;
            s[1] = '\0';
            return s;
         }
         else if (c == '\0')
            return copy(s1);
         else {
            int l = (int)strlen(s1);
            char *s = new char[l + 2];
            strcpy(s, s1);
            s[l] = c;
            s[l + 1] = '\0';
            return s;
         }
      }

      // another version of our private copy of strcat:
      char* cat(const char c, const char *s1) {
         if (s1 == NULL) {
            if (c == '\0')
               return NULL;
            char *s = new char[2];
            s[0] = c;
            s[1] = '\0';
            return s;
         }
         else if (c == '\0')
            return copy(s1);
         else {
            int l = (int)strlen(s1);
            char *s = new char[l + 2];
            s[0] = c;
            strcat(s, s1);
            s[l + 1] = '\0';
            return s;
         }
      }

   public:
      /**
       * Default constructor
       */
      String(const char* s = NULL) {
         str = copy(s);
      }

      /**
       * Another constructor.
       */
      String(const String& s) {
         str = copy(s.GetBuffer());
      }

      /**
       * Returns the number of characters in the String object.
       */
      int GetLength() const {
         if (IsEmpty())
            return 0;
         return (int)strlen(str);
      }

      /**
       * Converts a specified range of characters
       * in this String to lowercase characters.
       */
      void MakeLower(const int from, int to) {
         if (!IsEmpty()) {
            int l = GetLength();
            if (to > l)
               to = l;
            for (int i = from; i < to; i++)
               str[i] = tolower(str[i]);
         }
      }

      /**
       * Converts all the characters in this
       * string to lowercase characters.
       */
      void MakeLower() {
         if (!IsEmpty()) {
            int l = GetLength();
            for (int i = 0; i < l; i++)
               str[i] = tolower(str[i]);
         }
      }

      /**
       * Converts a specified range of characters
       * in this string to uppercase characters.
       */
      void MakeUpper(const int from, int to) {
         if (!IsEmpty()) {
            int l = GetLength();
            if (to > l)
               to = l;
            for (int i = from; i < to; i++)
               str[i] = toupper(str[i]);
         }
      }

      /**
       * Converts all the characters in this
       * string to uppercase characters.
       */
      void MakeUpper() {
         if (!IsEmpty()) {
            int l = GetLength();
            for (int i = 0; i < l; i++)
               str[i] = toupper(str[i]);
         }
      }

      /**
       * Forces a string to have 0 length
       */
      void Empty() {
         if (str != NULL) {
            delete[] str;
            str = NULL;
         }
      }

      // destructor:
      ~String() {
         Empty();
      }

      /**
       * Returns a pointer to the characters in the string.
       */
      char* GetBuffer() const {
         return str;
      }

      /**
       * Overloads the + operator, such that we can use
       * string = string + c.
       */
      friend String operator+(String& L, const char c) {
         String tmp;
         tmp.str = L.cat(L.GetBuffer(), c);
         return tmp;
      }

      /**
       * Overloads the + operator, such that we can use
       * string = c + string.
       */
      friend String operator+(const char c, String& L) {
         String tmp;
         tmp.str = L.cat(c, L.GetBuffer());
         return tmp;
      }

      /**
       * Overloads the + operator, such that we can use
       * String = String + "some text".
       */
      friend String operator+(String& L, const char* s) {
         String tmp;
         tmp.str = L.cat(L.GetBuffer(), s);
         return tmp;
      }

      /**
       * Overloads the + operator, such that we can use
       * String = "some text" + String.
       */
      friend String operator+(const char* s, String& L) {
         String tmp;
         tmp.str = L.cat(s, L.GetBuffer());
         return tmp;
      }

      /**
       * Overloads the + operator, such that we can use
       * String = String + String.
       */
      friend String operator+(String& s, String& L) {
         String tmp;
         tmp.str = L.cat(s.GetBuffer(), L.GetBuffer());
         return tmp;
      }

      /**
       * Overloaded =operator, case of a String at the right.
       * This function allows something like String = char, building a
       * String of the character.
       */
      String& operator=(const char c) {
         Empty();
         str = new char[1];
         str[0] = c;
         return *this;
      }

      /**
       * Overloaded =operator, case of a const char* at the right.
       * This function allows an asignment like String = char *.
       */
      String& operator=(const char* s) {
         Empty();
         str = copy(s);
         return *this;
      }

      /**
       * Overloaded =operator, case of a String at the right.
       * This function allows something like String1 = String2.
       */
      String & operator=(const String & L) {
         if (this == &L)
            return *this;
         Empty();
         str = copy(L.str);
         return *this;
      }

      /**
       * Overloaded =operator, case of a float value at the right.
       * This function allows something like String = 3.
       */
      String& operator=(const int N) {
         strstream s;
         s << N << '\0';
         Empty();
         str = copy(s.str());
         return *this;
      }

      /**
       * Overloaded =operator, case of a float value at the right.
       * This function allows something like String = 0.32.
       */
      String& operator=(const float N) {
         strstream s;
         s << N << '\0';
         Empty();
         str = copy(s.str());
         return *this;
      }

      /**
       * Overloaded =operator, case of a float value at the right.
       * This function allows something like String = 0.32.
       */
      String& operator=(const double N) {
         strstream s;
         s << N << '\0';
         Empty();
         str = copy(s.str());
         return *this;
      }

      /**
       * Overloaded += operator. Adds a character to the current String.
       */
      String& operator+=(const char c) {
         char* tmp = cat(str, c);
         Empty();
         str = tmp;
         return *this;
      }

      /**
       * Overloaded += operator. Adds a String to the current String.
       */
      String& operator+=(const char* s) {
         char* tmp = cat(str, s);
         Empty();
         str = tmp;
         return *this;
      }

      /**
       * Overloaded += operator. Adds an integer to the current String
       */
      String& operator+=(const int i) {
         strstream s;
         s << str << i << '\0';
         Empty();
         str = copy(s.str());
         return *this;
      }

      /**
       * Overloaded += operator. Adds an integer to the current String
       */
      String& operator+=(const float f) {
         strstream s;
         s << str << f << '\0';
         Empty();
         str = copy(s.str());
         return *this;
      }

      /**
       * Overloaded += operator. Adds a String to the current String.
       */
      String& operator+=(const String& s) {
         char* tmp = cat(str, s.GetBuffer());
         Empty();
         str = tmp;
         return *this;
      }

      /**
       * Overloaded == operator. Compares a String with the current String.
       */
      bool operator==(String& s) const {
         return Compare(s);
      }

      /**
       * Overloaded == operator. Compares the current string with a
       * pointer to char.
       */
      bool operator==(const char* s) const {
         return Compare(s);
      }

      /**
       * Overloaded != operator. Inequality of two strings.
       */
      bool operator!=(String& s) const {
         return !Compare(s);
      }

      /**
       * Overloaded != operator. Compares the current string with a
       * pointer to char for inequality.
       */
      bool operator!=(const char* s) const {
         return !Compare(s);
      }

      /**
       * Compares the current String to another String object.
       * Compares, character per character, the contents of a
       * String with the current String.
       */
      bool Compare(String& L) const {
         return (str == NULL || L.str == NULL || strcmp(L.str, str)) ? false : true;
      }

      /**
      * Compares the current String with a pointer to char.
      * Compares, character per character, the contents of a
      * String (char*) with the current String.
      */
      bool Compare(const char* s) const {
         return strcmp(s, str) ? false : true;
      }

      /**
       * Sets a character at a given position.
       */
      void SetAt(int index, char c) const {
         if (str != NULL) {
            str[index] = c;
         }
      }

      /**
       * Returns the character at a given position.
       */
      char GetAt(const int index) const {
         if (IsEmpty())
            return '\0';
         return str[index];
      }

      /**
       * Finds the first match of character c.in the String
       * k = 0 if starting from the beginning, -1 if
       * starting from the end, -2 if starting from the end-1 etc.
       * it returns the position of character c, or -1 if not found.
       */
      int Find(const char c, const int k = 0) const {
         if (!IsEmpty()) {
            int i = (k < 0) ? GetLength() + k + 1 : k;
            int end = (k < 0) ? 0 : GetLength();
            int inc = (k < 0) ? -1 : 1;
            int a = (i > end) ? end : i;
            int b = (i > end) ? i : end;
            while (i >= a && i <= b) {
               if (str[i] == c)
                  return i;
               i += inc;
            }
         }
         return -1;
      }

      /**
       * Inserts a character at the given index.
       */
      void InsertAt(const int index, const char c) {
         if (!IsEmpty()) {
            int l = GetLength() + 1;
            char* newStr = new char[l + 1];
            int i, j;
            for (i = 0, j = 0; i < l; i++) {
               if (i == index)
                  newStr[j++] = c;
               newStr[j++] = str[i];
            }
            delete[] str;
            str = newStr;
         }
      }

      /**
       * Returns a subString of the String, starting at a specific index.
       */
      String& SubString(const int i) const {
         int l = GetLength();
         if (i < l)
            return *new String(str + i);
         else
            return *new String;
      }

      /**
       * Returns a substring between starting and ending indices.
       */
      String& SubString(const int from, const int to) const {
         int l = GetLength();
         if (from < 0 || to > l || to < from)
            return *new String;
         else {
            char* s = new char[to - from + 1];
            int i, k;
            for (k = 0, i = from; i < to; i++, k++) s[k] = str[i];
            s[k] = '\0';
            return *new String(s);
         }
      }

      /**
       * Appends a string to the current string.
       * returns new string which is the sum of the current string and the argument.
       */
      void Append(String& s) {
         char* newStr = cat(str, s.GetBuffer());
         if (!IsEmpty())
            delete[] str;
         str = newStr;
      }

      /**
       * Appends a string to the current string.
       * returns new string which is the sum of the current string and the argument.
       */
      void Append(const char* s) {
         char* newStr = cat(str, s);
         if (!IsEmpty())
            delete[] str;
         str = newStr;
      }

      /**
       * Returns true if the current string starts with the specified character.
       */
      bool IsFirst(char c) const {
         return (str != NULL && str[0] == c) ? true : false;
      }

      /**
       * Returns true if the current string starts with the specified string.
       */
      bool IsFirst(String& s) const {
         return (str != NULL && !strncmp(s.GetBuffer(), str, s.GetLength())) ? true : false;
      }

      /**
       * Returns true if the current string starts with the specified string.
       */
      bool IsFirst(const char* s) const {
         return (str != NULL && !strncmp(s, str, strlen(s))) ? true : false;
      }

      /**
       * Returns true if the current string ends with the specified string.
       */
      bool IsLast(const char* s) const {
         if (str != NULL) {
            int len = (int)strlen(s);
            if (len <= GetLength() && !strcmp(str + GetLength() - len, s))
               return true;
         }
         return false;
      }

      /**
       * Returns true if the current string ends with the specified string.
       */
      bool IsLast(String& s) const {
         return IsLast(s.GetBuffer());
      }

      /**
       * Returns true if the current string starts with a digit.
       */
      bool IsFirstDigit() const {
         const char digits[] = "0123456789";
         return (str != NULL && !strcspn(str, digits)) ? true : false;
      }

      /**
       *
       */
      bool IsDigit() const {
         const char digits[] = "0123456789";
         for (int i = 0; i < GetLength(); i++) {
            char *s = &str[i];
            if (strcspn(s, digits)) return false;
         }
         return true;
      }

      /**
       * Returns the index of the last valid stochiometric digit with
       * which the actual string starts. The stochiometric fraction,
       * either integer or float is given in f.
       */
      int GetFirstDigits(float &f) const {
         const char validChar[] = "0123456789.";
         int i;

         if (str != NULL) {
            int len = (int)strspn(str, validChar);
            char* s = new char[len];
            for (i = 0; i < len; i++)
               s[i] = str[i];
            // append null-terminating string char
            s[i] = '\0';
            f = (float)atof(s);
            //delete [] s;
            return len;
         }
         else
            return -1;
      }

      /**
       * Tests whether the string contains anything else than NULL.
       * returns 1 if the current string is NULL, 0 otherwise.
       */
      bool IsEmpty() const {
         return (str == NULL) ? true : false;
      }

      /**
       * Overloaded output stream.
       */
      friend ostream& operator<<(ostream& out, const String& L) {
         if (!L.IsEmpty())
            out << L.str;
         else
            out << "(null)";
         return out;
      }

      /**
       * Trims the string by deleting leading and trailing spaces.
       * returns a string with spaces at the start and end.
       */
      void Trim() {
         if (!IsEmpty()) {
            int i, k, start, end, l = GetLength();
            start = 0;
            while (isspace(str[start])) start++;
            end = l - 1;
            while (isspace(str[end])) end--;
            end++;
            if (end > start) {
               char* newStr = new char[end - start + 1];
               for (i = start, k = 0; i < end; i++, k++) newStr[k] = str[i];
               newStr[k] = '\0';
               delete[] str;
               str = newStr;
            }
            else Empty();
         }
      }

      /**
       * returns the charge of a chemical formula string.
       */
      float GetCharge()
      {
         float charge;
         char chargeSign;
         bool bad = false;
         /*
         if (IsEmpty()) {
         IO_Error err = IO_Error("Warning: Chemical formula string is empty!");
         cout << err.message() << endl;
         bad = true;
         return false;
         }
         */
         // cleanup
         Trim();

         // find the first charge sign in this string
         int firstCharge = Find(PLUS_CHARGE);
         if (firstCharge < 0) { // not found
            firstCharge = Find(MINUS_CHARGE);
            if (firstCharge < 0) {
               charge = 0.0;
               return charge;
            }
            else
               chargeSign = MINUS_CHARGE;
         }
         else
            chargeSign = PLUS_CHARGE;

         if (chargeSign == MINUS_CHARGE || chargeSign == PLUS_CHARGE) {
            bool natural = false;  // e.g. Ca+2

            // search for the last occurance of the charge
            int lastCharge = Find(chargeSign, -1);
            if (lastCharge < 0) {
               bad = true;
               return 0; // something is going wrong
            }

            // specie is given in compact notation, e.g. Ca+2
            else if (lastCharge == firstCharge) {
               String val = SubString(firstCharge);
               if (val.GetLength() == 1)
                  charge = (chargeSign == PLUS_CHARGE) ? 1.0f : -1.0f;
               else if (!val.IsEmpty()) {
                  char *str;
                  if (val.IsDigit())
                     charge = (float)strtod(val.GetBuffer(), &str);
                  else
                     return 0; // bad !
               }
            }

            // specie is given in natural notation, and
            // the detected position is the last character
            else if (lastCharge == GetLength() - 1){
               natural = true;
               for (int i = firstCharge; i < lastCharge + 1; i++){
                  char symbl = GetAt(i);
                  if (symbl != chargeSign) {
                     char* msg;
                     sprintf(msg, "chemical formula %s is invalid", str);
                     IO_Error err = IO_Error(msg);
                     cout << err.message() << endl;
                     bad = true;
                     return 0; // bad name
                  }
               }
               charge = (chargeSign == PLUS_CHARGE) ?
                  (float)(lastCharge - firstCharge + 1) :
                  (float)(-lastCharge + firstCharge - 1);
            }

            else {
               bad = true;
               return 0; // bad
            }
         }
         return charge;
      }

      /**
       * return another string after excluding the charge fraction
       * of this current string.
       */
      String& ExcludeCharge()
      {
         if (IsEmpty())
            return *this;

         else {
            Trim();

            // find the first charge sign in this string
            int firstCharge = Find(PLUS_CHARGE);
            if (firstCharge < 0) { // not found
               firstCharge = Find(MINUS_CHARGE);
               if (firstCharge < 0)
                  return *this;
            }

            // now extracts what we need
            return (SubString(0, firstCharge));
         }
      }

      /**
       * makes more compact the charge radical of a chemical symbol.
       */
      void CompactCharge()
      {
         int charge = (int)GetCharge();
         if (charge == 0) return;

         *this = ExcludeCharge();
         if (charge > 0)
            *this += '+';
         else if (charge < 0)
            *this += '-';
         if (abs(charge) > 1) {
            char buffer[2];
            sprintf(buffer, "%d", abs(charge));
            *this += String(buffer);
            Trim();
         }
      }

      /**
       *
       */
      bool RemoveTrailingBrackets()
      {
         int firstBracket = Find('(', 0);
         if (firstBracket == -1) return false;
         if (!IsLast(")")) return false;
         *this = SubString(0, firstBracket);
         return true;
      }

      /**
       * Overloads the > operator, such that we can compare two strings.
       */
      friend int operator>(String & a, String & b) {
         if (a.IsEmpty())
            return 0;
         if (b.IsEmpty())
            return 1;
         if (a.GetAt(0) > b.GetAt(0))
            return 1;
         if (a.GetAt(0) < b.GetAt(0))
            return 0;
         return (a.SubString(1) > b.SubString(1));
      }

      /**
       * Overloads the < operator, such that we can compare two strings.
       */
      friend int operator<(String & a, String & b) {
         return (b > a);
      }

      /**
       * Overloads the >= operator, such that we can compare two strings.
       */
      friend int operator>=(String & a, String & b) {
         return (!(a < b));
      }

      /**
       * Overloads the <= operator, such that we can compare two strings.
       */
      friend int operator<=(String & a, String & b) {
         return (!(a > b));
      }

      /**
       * operator[] returns character at index i
       */
      char & operator[] (const int i) {
#ifdef _DEBUG
         if (i < 0 || i >= GetLength() ) {
            cerr << "Warning: trying to access beyond String's buffer (1): " << str
               << "\n\t" << i << " >= " << GetLength() << endl;
            exit(0);
         }
#endif
         return str[i];
      }

      /**
       * operator[] returns character at index i
       */
      const char & operator[] (const int i) const {
#ifdef _DEBUG
         if (i < 0 || i >= GetLength() ) {
            cerr << "Warning: trying to access beyond String's buffer (1): " << str
               << "\n\t" << i << " >= " << GetLength() << endl;
            exit(0);
         }
#endif
         return str[i];
      }

   };

}

#endif // HAPP_STRING
