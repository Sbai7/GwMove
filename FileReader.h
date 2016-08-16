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

#ifndef HAPP_FILE_READER
#define HAPP_FILE_READER

// Project headers
#include <AdString.h>

// C style headers
#include <stdlib.h>

// C++ headers 
#include <fstream>
#include <iostream>
#include <iomanip>
#include <list>
using namespace std;


namespace happ
{
#define COMMENT_CHAR '#'
#define TOKENSEP " \t,"

#define STRLEN 10000
#define MAXTOKENS 500
#define SLEN 4000
#define MTOKENS 200

   /**
    * Defines a base class for input streams objects processing.
    * This class has all funcationnalities needed for high level
    * processing of input streams. It tokenizes strings automatically,
    * and much more.
    *
    * \author Adil Sbai, Ph.D. (a.sbai@brgm.fr)
    * \version 1.0
    */
   class FileReader
   {
   private:
      /// file input stream definition
      ifstream inf;

      /// general input stream (can be std, or file), must be pointer
      istream *in;

      /// this string contains the line of input, will have length 'mslen'
      char* str;

      /// list of tokens of the line, will have length 'mtokens'
      char** tokenlist;

      /// list of token separators
      char* tokensep;

      /// number of tokens
      int Ntok;

      /// shift number, used to skip one or more tokens
      int shift;

      /// actual linenumber
      long int linenr;

      /// maximum length of the string
      int mslen;

      /// maximum number of tokens allowed
      int mtokens;

      /// character used to indicate comments in the file
      char comment;

      /// ALL the lines read during the object's life
      long int linesread;

      /// is TRUE if something went wrong
      int failbit;

      /// is false if we're working from file, true if it's the standard input
      int fromstd;

      /// name of the file
      String filename;

      /**
       * Private tokenize function for this class.
       * \param s pointer of chars to tokenize.
       * \param b if true we tokenize till the first comment,
       * otherwise, we continue.
       */
      int Tokenize(char* s, const bool b = true)
      {
         int i = 0;
         if (s == NULL)
            return 0;
         char* thistoken = strtok(s, tokensep);
         while (i < mtokens && thistoken != NULL) {
            if (b && thistoken[0] == COMMENT_CHAR)
               break;
            tokenlist[i] = thistoken;
            thistoken = strtok(NULL, tokensep);
            i++;
         }
         tokenlist[i] = NULL;
         return i;
      }

      /**
       * Reads one line, puts the content in a string and tokenizes.
       * \return 1 if successful, 0 if EOF.
       */
      int ReadLine()
      {
         int i;
         do {
            linenr++;
            if (in->eof() || !in->getline(str, mslen))
               return 0;
            while (str[(i = (int)strlen(str)) - 1] == '\\') {
               linenr++;
               if (in->eof() || !in->getline(str + i - 1, mslen - i + 1))
                  return 0;
            }
            i = 0;
            while (str[i] == ' ') i++;
         } while (str[i] == comment || str[i] == '\0');
         return 1;
      }

      /**
       * Reads a line, even if it is commented out.
       * \return 1 if successful, 0 if EOF.
       */
      int ReadAllLine()
      {
         int i;
         do {
            linenr++;
            if (in->eof() || !in->getline(str, mslen))
               return 0;
            while (str[(i = (int)strlen(str)) - 1] == '\\') {
               linenr++;
               if (in->eof() || !in->getline(str + i - 1, mslen - i + 1))
                  return 0;
            }
            i = 0;
            while (str[i] == ' ') i++;
         } while (str[i] == '\0');
         return 1;
      }

      /**
       * reads a line, even if it is commented out and even if it is blank
       * and even if it ends with a backslash. Literally.
       * \return 1 if successful, 0 if EOF.
       */
      int ReadAllLine2()
      {
         linenr++;
         if (in->eof() || !in->getline(str, mslen))
            return 0;
         return 1;
      }

      /**
       * Copies a string into a newly memory-allocated string
       */
      char* MemCopy(const char* s)
      {
         return (char*)((s == NULL) ? NULL : strcpy(new char[strlen(s) + 1], s));
      }

      void Init(const int l, const int mtok, const char* ts, const char c)
      {
         comment = c;
         tokenlist = NULL;
         str = NULL;
         if ((mslen = l) != 0)
            str = new char[mslen];
         else
            failbit = 1;
         if ((mtokens = mtok) != 0)
            tokenlist = new char*[mtokens];
         else
            failbit = 1;
         tokensep = MemCopy(ts);
      }

   public:

      /**
       * Opens an inputstream.
       * The inputstream is given by a name in the argument.
       * \param f name in terms of a char-pointer.
       * \return true if OK, false otherwise.
       */
      bool Open(const char* f)
      {
         filename = f;
         fromstd = failbit = shift = 0;
         linesread = linenr = (long)0;
         if (!filename.IsEmpty()) {
            //inf.open(f,ios::in|ios::nocreate);
            inf.open(f, ios::in);
            in = &inf;
         }
         if (!failbit) failbit = in->fail();
         return (failbit) ? false : true;
      }

      /**
       * Opens an inputstream.
       * The inputstream is initialized before hand.
       * \param is valid (i.e., opened) inputstream.
       * \return true if OK, false otherwise.
       */
      bool Open(istream& is)
      {
         fromstd = 1;
         failbit = shift = 0;
         linesread = linenr = (long)0;
         in = &is;
         return true;
      }

      /**
       * Opens an inputstream.
       * The inputstream is specified by a name in a String instance.
       * \param s name of the input stream.
       * \return true if OK, false otherwise.
       */
      bool Open(String& s)
      {
         return Open(s.GetBuffer());
      }

      /**
       * Constructor, builds the FReader.
       * \param f pointer to the name of the input stream.
       * \param l maximum length of one input line (defaults to SLEN).
       * \param mtok maximum number of tokens (defaults to MTOKENS).
       * \param ts pointer to token-separation characters (defaults to TOKENSEP).
       * \param c character used to specify comments lines (skipped).
       */
      FileReader(const char* f, const int l = SLEN, const int mtok = MTOKENS,
         const char* ts = TOKENSEP, const char c = COMMENT_CHAR)
      {
         Init(l, mtok, ts, c);
         Open(f);
      }

      /**
       * Constructor, builds the FReader.
       * \param s instance of a String object, name of the input stream.
       * \param l maximum length of one input line (defaults to SLEN).
       * \param mtok maximum number of tokens (defaults to MTOKENS).
       * \param ts pointer to token-separation characters (defaults to TOKENSEP).
       * \param c character used to specify comments lines (skipped).
       */
      FileReader(String& s, const int l = SLEN, const int mtok = MTOKENS,
         const char* ts = TOKENSEP, const char c = COMMENT_CHAR)
      {
         Init(l, mtok, ts, c);
         Open(s);
      }

      /**
       * Constructor, builds the FReader.
       * \param is input stream, must be opened (and validated) before hand.
       * \param l maximum length of one input line (defaults to SLEN).
       * \param mtok maximum number of tokens (defaults to MTOKENS).
       * \param ts pointer to token-separation characters (defaults to TOKENSEP).
       * \param c character used to specify comments lines (skipped).
       */
      FileReader(istream& is, const int l = SLEN, const int mtok = MTOKENS,
         const char* ts = TOKENSEP, const char c = COMMENT_CHAR)
      {
         Init(l, mtok, ts, c);
         Open(is);
      }

      /**
       * Constructor, builds the FReader.
       * With no arguments, a default file reader with the standard input stream
       * (keyboard) is instantiated.
       */
      FileReader()
      {
         Init(SLEN, MTOKENS, TOKENSEP, COMMENT_CHAR);
         Open(cin);
      }

      /**
       * Class destructor.
       */
      virtual ~FileReader()
      {
         if (str != NULL)
            delete[] str;
         if (tokenlist != NULL)
            delete[] tokenlist;
         if (tokensep != NULL)
            delete[] tokensep;
      }

      /**
       * Return token number i from the current input line.
       */
      char* GetToken(const int i)
      {
         return tokenlist[i + shift];
      }

      /**
       * Returns true if the token with index i equals the string in argument.
       * \param index index of the token to compare with.
       * \param s string to compare with.
       * \return true if the token with index i equals the string in argument.
       */
      bool TokenCompare(const int index, const char* s)
      {
         return strcmp(s, GetToken(index)) ? false : true;
      }

      /**
       * Returns true if the token with index i equals the string in argument.
       * \param index index of the token to compare with.
       * \param s string to compare with.
       * \return true if the token with index i equals the string in argument.
       */
      bool TokenCompare(const int index, String& s)
      {
         return s.Compare(GetToken(index));
      }

      /**
       * Returns true if the token with index i starts with the string in argument.
       * \param index index of the token to compare with.
       * \param s string to compare with.
       * \return true if the token with index i starts with the string in argument.
       */
      bool TokenStartsWith(const int index, String& s)
      {
         String tmp(GetToken(index));
         return tmp.IsFirst(s);
      }

      /**
       * Returns true if the token with index i starts with the string in argument.
       * \param index index of the token to compare with.
       * \param s string to compare with.
       * \return true if the token with index i starts with the string in argument.
       */
      bool TokenStartsWith(const int index, const char* s)
      {
         String tmp(GetToken(index));
         return tmp.IsFirst(s);
      }

      /**
       * Returns true if the token with index i ends with the string in argument.
       * \param index index of the token to compare with.
       * \param s string to compare with.
       * \return true if the token with index i ends with the string in argument.
       */
      bool TokenEndsWith(const int index, String& s)
      {
         String tmp(GetToken(index));
         return tmp.IsLast(s);
      }

      /**
       * Returns true if the token with index i ends with the string in argument.
       * \param index index of the token to compare with.
       * \param s string to compare with.
       * \return true if the token with index i ends with the string in argument.
       */
      bool TokenEndsWith(const int index, const char* s)
      {
         String tmp(GetToken(index));
         return tmp.IsLast(s);
      }

      /**
       * \return the number of tokens of the current line.
       */
      int GetNtokens()
      {
         return Ntok - shift;
      }

      /**
       * \return the entire current input line.
       */
      char* GetString()
      {
         return str;
      }

      /**
       * \return the current filename.
       */
      String & GetFilename()
      {
         return filename;
      }

      /**
       * gets a new line (updates str, tokenlist, Ntok, status...).
       * \return 0 if EOF, else the number of tokens read.
       */
      int GetLine()
      {
         shift = Ntok = 0;
         if (ReadLine())
            Ntok = Tokenize(MemCopy(str));
         return Ntok;
      }

      /**
       * as getline function, but gets all lines.
       * Hence lines starting with the comment-character are
       * not disregarded, but blank lines are
       * \return 0 if EOF, the number of tokens of this line otherwise.
       */
      int GetAllLines()
      {
         shift = Ntok = 0;
         if (ReadAllLine())
            Ntok = Tokenize(MemCopy(str), false);
         return Ntok;
      }

      /**
       * as getline function, but returns the line exactly as is.
       * \return 0 if EOF, 1 otherwise.
       */
      int GetLiteral()
      {
         shift = Ntok = 0;
         if (ReadAllLine2() && str)
            Ntok = Tokenize(MemCopy(str), false);
         return 1;
      }

      /**
       * \return the current line number.
       */
      int ThisLine()
      {
         return linenr;
      }

      /**
       * \return the total number of lines read since this object was created.
       */
      int All()
      {
         return linesread + linenr;
      }

      /**
       * \return true if something went wrong, false otherwise.
       */
      int Fail()
      {
         return failbit;
      }

      /**
       * rewind the file, clear the flags. This is obtained by clearing
       * the stream, and setting the stream flags to the beginning of
       * of the file (seekg(ios::beg).
       */
      void Rewind()
      {
         //in->seekg(0L,ios::beg); // put pointer at begin of file
         in->clear();              // erase eof flag    
         in->seekg(ios::beg);      // put pointer at begin of file
         linesread += long(linenr);
         Ntok = shift = 0;
         linenr = (long)0;
      }

      /**
       * shift all items to the left by the specified number.
       * The result is as if the actual line started one or several
       * tokens later.
       * \param i number of tokens you want to shift.
       */
      void LeftShift(const int i)
      {
         shift += i;
      }

      /**
       * \return true if we are reading from the standard-input, false otherwise
       */
      int Interactive()
      {
         return fromstd;
      }

      /**
       * returns true if the current line starts with the specified character.
       * \param c character which is searched for.
       * \return true if the first character of the current line starts with c.
       */
      bool StartsWith(const char c)
      {
         return (str[0] == c) ? true : false;
      }

      /**
       * Sets the character used to indicate a comment-line.
       * \param c character which indicates a comment.
       */
      void SetComment(char c)
      {
         comment = c;
      }

      /**
       * Sets the characters used to separate tokens.
       * \param s string which indicates a comment.
       */
      void SetTokenSep(char* s)
      {
         tokensep = MemCopy(s);
      }

      /**
       * This function sets the string of the filereader.
       * This function allows to use FReader without a file, but
       * merely as a string-parser.
       * \param s string to be parsed.
       */
      void PushString(String s)
      {
         PushString(s.GetBuffer());
      }

      /**
       * This function sets the string of the filereader.
       * This function allows to use FReader without a file, but
       * merely as a string-parser.
       * \param s string to be parsed.
       */
      void PushString(const char* s)
      {
         str = MemCopy(s);
         Ntok = Tokenize(MemCopy(str));
         shift = 0;
      }

      /**
       * This function rewind the opnened file, and reads
       * the file until the first keyword 'key' is found.
       * The string 'key' must be the first string in the line.
       * \param key the keyword string to which we'll skip.
       * \return true if successful, false if this section doesn't exist.
       */
      bool GotoSection(String key)
      {
         // rewind the opened stream
         Rewind();

         bool ok = false; // ok if the keyword is found
         int notEOF;

         do {
            if ((notEOF = GetLine()) != 0) {
               // does the current line starts with the given key ?
               if (GetNtokens() > 0)
                  ok = TokenCompare(0, key);
            }
         } while (ok == false && notEOF);

         return ok;
      }

      /**
       * This function returns in a list of strings, all lines
       * of text between the actual key and one in a list of
       * keywords. Comment lines are excluded.
       * Useful to get quickly a file section to be post
       * processed otherwise.
       * \param key the keyword string to which we'll skip.
       * \param keywords list of other possible keywords.
       * \param text list of lines read through.
       * \return true if successful, false if this section doesn't exist.
       */
      bool GetSectionStrings(String& key, std::list<String>& keywords, std::list<String>& text)
      {
         // Init fill in the list of strings
         text.clear();

         // the key must exist in the list of keywords
         bool found = false;
         std::list <String>::iterator iter;
         for (iter = keywords.begin(); iter != keywords.end(); iter++){
            if (*iter == key) {
               found = true;
               break;
            }
         }
         if (!found)
            return false;

         if (GotoSection(key)) { // check if this data block exists ...
            bool HitKey = false;
            int NotEOF;
            do {
               if ((NotEOF = GetLine()) > 0) {
                  for (iter = keywords.begin(); iter != keywords.end(); iter++){
                     if (TokenCompare(0, *iter))
                        HitKey = true;
                  }
                  if (!HitKey)
                     text.push_back(String(GetString()));
               }
            } while (HitKey == false && NotEOF);
            return true;
         }
         return false;
      }
   };

}

#endif // HAPP_FILE_READER
