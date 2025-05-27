#ifndef _LIST_H_
#define _LIST_H_

#include "../types.h"

#include <stdlib.h>

List new_List(delfunc pDeleteFunction, copyfunc pCopyFunction, strfunc pToStringFunction);
void del_List(List *self);
List copy_List(List self);
String str_List(List self);

size_t List_GetSize(List self);

void *List_Get(List self, size_t pIndex);

void List_Append(List self, void *pValue);
void List_Remove(List self, size_t pIndex);

#endif // !_LIST_H_