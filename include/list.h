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
 * Adil's list template class.
 * For a fast and light management of doubly linked lists.
 *
 * Created 13 April 2004
 * Last update of 13 April 2004
 *
 * Author: A. Sbai, Ph.D. (a.sbai@brgm.fr)
 */

#ifndef HAPP_LIST
#define HAPP_LIST

#include <exception.h>
//using namespace std;


namespace happ
{
   /// A forward definition for friend status
   template <class C> class ListIterator;

   /**
    * Private class for internal use by the list class.
    * Any item has two pointers to the previous and
    * next elements in the same list.
    *
    * \author Adil Sbai, Ph.D. (a.sbai@brgm.fr)
    * \version 1.0
    */
   template <class A> class listelem
   {
   public:
      listelem *next, *prev;
      A value;

      listelem()
      {
         next = prev = NULL;
      }

      // Copy constructor
      listelem(listelem<A>& copy)
      {
         next = copy.next;
         prev = copy.prev;
         value = copy.value;
      }
   };


   /**
    * Implements a light weight, but optimized, doubly linked list.
    */
   template <class B> class list{
   private:
      /// pointer to top element in the list.
      listelem<B> *head;

      /// length of the list.
      int length;

      /*
       * Take advantage of the fact that the
       * list may often be accessed sequentially.
       */
      /// index of last accessed element.
      int last_index;

      /// pointer to last accessed item in the list.
      listelem<B> *last_ptr;

   public:
      friend class ListIterator<B>;

      /**
       * Perform all default initialization.
       */
      void Init()
      {
         while (!IsEmpty())
            Pop();

         length = 0;
         last_index = 0;
         last_ptr = NULL;
      }

      /// Default constructor.
      list()
      {
         head = NULL;
         length = 0;
         last_index = 0;
         last_ptr = NULL;
      }

      /// Default destructor.
      ~list()
      {
         Init();
      }

      /// Copy constructor.
      list(list<B>& copy)
      {
         Init();
         for (int i = 0; i < copy.length; i++)
            Push(copy[i]);
      }

      /**
       * Overloading of addition operator.
       * Insert another list of the same type at the end of this list.
       */
      list<B>& operator+=(list<B>& ls)
      {
         for (int i = 0; i < ls.length; i++)
            Push(ls[i]);

         return (*this);
      }

      /**
       * Overloading of equal = operator.
       */
      list<B>& operator=(list<B>& copy)
      {
         // Remove any existing elements
         while (!IsEmpty())
            Pop();

         for (int i = 0; i < copy.length; i++)
            Push(copy[i]);

         return (*this);
      }

      /**
       * Checks if a given element is inside the list.
       * \return the index of the searched element, -1 otherwise.
       */
      int IsInside(B& elem)
      {
         for (int i = 0; i < length; i++){
            if ((*this)[i] == elem){
               return i;
            }
         }
         return -1;
      }

      /**
       * Move the list elements from the copy to this list.
       * The memory is physically transfered from the copy list to
       * the tail of the current list.
       */
      void Move(list<B>& copy)
      {
         // Remove any existing elements
         //while (!IsEmpty())
         // Pop();

         listelem<B> *cur = head;

         if (cur == NULL) {
            // We are starting with a NULL list
            head = copy.head;
         }
         else {
            // Find the end of the list
            while (cur->next)
               cur = cur->next;

            cur->next = copy.head;

            if (cur->next)
               cur->next->prev = cur;
         }

         length += copy.length;

         last_ptr = head;
         last_index = 0;

         copy.head = NULL;
         copy.last_ptr = NULL;
         copy.length = 0;
         copy.last_index = 0;
      }

      /**
       * Make i the first item in the list.
       */
      void BringToFront(int i)
      {
         if ((i < 0) || (i >= length))
            throw Error("Error! List index out of bounds!.",
            __FILE__, __LINE__);

         listelem<B> *cur = last_ptr;
         int temp_index = i;

         if (i >= last_index) {
            i -= last_index;
            // Go to the i^{th} element of the list.
            while (i != 0) {
               cur = cur->next;
               i--;
            }
         }
         else /* (i < last_index) */ {
            i = last_index - i;
            // Go to the i^{th} element of the list.
            while (i != 0) {
               cur = cur->prev;
               i--;
            }
         }

         // cur now points to element i. We need to send it to the 
         // head of the list.
         if (cur == head) {
            // Do nothing
            last_ptr = head;
            last_index = 0;

            // All done
            return;
         }

         // Pull cur out of the list
         cur->prev->next = cur->next;

         // If cur->next != NULL
         if (cur->next) {
            cur->next->prev = cur->prev;
         }

         // Make head the second element and 
         // make cur the first element
         head->prev = cur;
         cur->next = head;
         cur->prev = NULL;

         head = cur;

         last_ptr = head;
         last_index = 0;
      }

      /**
       * Make i the last item in the list.
       */
      void SendToBack(int i)
      {
         if ((i < 0) || (i >= length))
            throw Error("Error! List index out of bounds!.",
            __FILE__, __LINE__);

         listelem<B> *cur = last_ptr;
         int temp_index = i;

         if (i >= last_index) {
            i -= last_index;
            // Go to the i^{th} element of the list.
            while (i != 0) {
               cur = cur->next;
               i--;
            }
         }
         else /* (i < last_index) */ {
            i = last_index - i;
            // Go to the i^{th} element of the list.
            while (i != 0) {
               cur = cur->prev;
               i--;
            }
         }

         // cur now points to element i. We need to send it to the 
         // back of the list.
         if (cur->next == NULL) {
            // Do nothing
            last_ptr = cur;
            last_index = temp_index;

            // All done
            return;
         }

         // Pull cur out of the list
         if (cur->prev) {
            cur->prev->next = cur->next;
         }

         // Find the terminal pointer
         listelem<B> *tmp_ptr = cur;
         while (tmp_ptr->next)
            tmp_ptr = tmp_ptr->next;

         // Finish pulling cur out of the list
         cur->next->prev = cur->prev;

         // Put cur at the tail of the list
         cur->prev = tmp_ptr;
         tmp_ptr->next = cur;
         cur->next = NULL;

         last_ptr = cur;
         last_index = length - 1;
      }

      /**
       * Make "thing" the i^{th} element of the list.
       */
      void Insert(B &thing, int i)
      {
         listelem<B> *cur = last_ptr;

         if ((i < 0) || (i > length))
            throw Error("Error! Attempt to add an "
            "element beyond the bounds of the list.",
            __FILE__, __LINE__);

         // Create the first element of the list?
         if (length == 0) {
            head = new listelem<B>;
            if (head == NULL)
               throw Error("Error! Can't allocate memory for first "
               "element of list", __FILE__, __LINE__);

            head->value = thing;
            last_ptr = head;
            last_index = 0;
         }
         // Add to the end of the list
         else if (i == length) {
            /* i must be greater than last_index */
            i -= last_index;
            while (i != 1) {
               cur = cur->next;
               i--;
            }

            cur->next = new listelem<B>;
            cur->next->value = thing;
            cur->next->prev = cur;
            last_ptr = cur->next;
            last_index = length;
         }
         // Add to the head of the list
         else if (i == 0) {
            last_ptr = new listelem<B>;
            last_ptr->next = head;
            last_ptr->prev = NULL;
            last_ptr->value = thing;
            head->prev = last_ptr;
            head = last_ptr;
            last_index = 0;
         }
         else {
            listelem<B> *temp;
            int temp_index = i;

            if (i >= last_index) {
               i -= last_index;
               while (i != 0) {
                  cur = cur->next;
                  i--;
               }
            }
            else /*(i < last_index) */ {
               i = last_index - i;
               while (i != 0) {
                  cur = cur->prev;
                  i--;
               }
            }

            temp = new listelem<B>;
            temp->value = thing;

            temp->next = cur;
            temp->prev = cur->prev;
            cur->prev = temp;
            temp->prev->next = temp;

            last_ptr = temp;
            last_index = temp_index;
         }

         length++;
      }

      /**
       * Make "thing" the last elememt of the list.
       */
      void Push(B &thing)
      {
         Insert(thing, length);
      }

      /**
       * Remove the i^{th} element of the list.
       */
      void Remove(int i)
      {
         listelem<B> *cur = last_ptr;
         int temp_index = i;

         if ((i < 0) || (i >= length))
            throw Error("Error! Attempt to remove an non-existant"
            "element the list.", __FILE__, __LINE__);

         if (i >= last_index) {
            i -= last_index;
            while (i > 0) {
               cur = cur->next;
               i--;
            }
         }
         else /* (i < last_index) */ {
            i = last_index - i;
            while (i > 0) {
               cur = cur->prev;
               i--;
            }
         }

         // Delete the i^{th} element from ...
         // ... the middle if the list
         if ((cur->prev != NULL) && (cur->next != NULL)) {
            cur->prev->next = cur->next;
            cur->next->prev = cur->prev;
            last_ptr = cur->prev;
            last_index = temp_index - 1;
         }
         // ... the end of the list
         else if (cur->prev != NULL) {
            cur->prev->next = NULL;
            last_ptr = cur->prev;
            last_index = temp_index - 1;
         }
         // ... the begining of the list.
         else if (cur->next != NULL) {
            head = cur->next;
            last_ptr = cur->next;
            last_index = 0;
            cur->next->prev = NULL;
         }
         else /* Last element */
            last_index = 0;

         length--;
         delete cur;
      }

      /**
       * Remove the last element of the list.
       * DOES NOT return the value of the popped element.
       */
      void Pop()
      {
         Remove(length - 1);
      }

      /**
       * returns the number of items in the list.
       */
      int GetLength() const
      {
         return length;
      };

      int IsEmpty() const
      {
         return !length;
      };

      /**
       * overloaded access by index operator [].
       */
      B& operator [] (int i)
      {
         listelem<B> *cur = last_ptr;
         int temp_index = i;

         if ((i < 0) || (i >= length))
            throw Error("Error! List index out of bounds!.",
            __FILE__, __LINE__);

         if (i >= last_index){
            i -= last_index;
            // Go to the i^{th} element of the list.
            while (i != 0){
               cur = cur->next;
               i--;
            }
         }
         else /* (i < last_index) */{
            i = last_index - i;
            // Go to the i^{th} element of the list.
            while (i != 0){
               cur = cur->prev;
               i--;
            }
         }

         last_ptr = cur;
         last_index = temp_index;

         return cur->value;
      }

      /**
       * overloaded output stream operator << .

       friend ostream& operator<<(ostream& out, list<B> &L)
       {
       for (int i=0; i<L.GetLength(); i++) {
       out << L[i] << "  ";
       }
       out << endl;
       return out;
       }
       */
#ifdef MFC
      /**
       * List serialization method.
       */
      void Serialize(CArchive& ar)
      {
         if(ar.IsStoring()){
            ar << length;

            for(int i = 0; i < length; i++)
               (*this)[i].Serialize(ar);
         }
         else{
            int len;
            ar >> len;

            // A little sanity check for the list length ...
            if (len < 0)
               AfxThrowArchiveException(CArchiveException::badIndex,
               ar.m_pDocument->GetTitle());

            for (int i = 0; i < len; i++){
               B tmp;
               tmp.Serialize(ar);
               Push(tmp);
            }
         }
      }
#endif // WIN32
   };


   /// Provide sequential list access
   template <class C> class ListIterator{
   private:
      int index;
      listelem<C> *cur;
   public:
      ListIterator(list<C> &lpList){
         cur = lpList.head;
         index = 0;
      };

      C* operator()(int &m_index){
         if (cur == NULL){
            return NULL;
         }

         // Save and increment the index
         m_index = index;

         index++;

         listelem<C> *ret_value = cur;

         cur = cur->next;

         return &(ret_value->value);
      }

      C* operator()(){
         if (cur == NULL){
            return NULL;
         }

         // Increment the index
         index++;

         listelem<C> *ret_value = cur;

         cur = cur->next;

         return &(ret_value->value);
      }
   };
}

#endif // HAPP_LIST
