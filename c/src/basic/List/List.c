#include "List.h"

#include <stdio.h>

#include "../String/String.h"

struct cNode;
typedef struct cNode *Node;

Node new_Node(void *pValue, delfunc pDelete, copyfunc pCopy, strfunc pToString);
void del_Node(Node *self);
Node copy_Node(Node self);
String str_Node(Node self);

struct cNode
{
    void *_value;
    delfunc _delete;
    copyfunc _copy;
    strfunc _toString;
};

Node new_Node(void *pValue, delfunc pDelete, copyfunc pCopy, strfunc pToString)
{
    Node result = (Node)malloc(sizeof(struct cNode));

    result->_value = pValue;
    result->_delete = pDelete;
    result->_copy = pCopy;
    result->_toString = pToString;

    return result;
}

void del_Node(Node *self)
{
    if (*self == NULL)
        return;

    (*self)->_delete(&((*self)->_value));

    free(*self);

    *self = NULL;

    return;
}

Node copy_Node(Node self)
{
    return new_Node(self->_copy(self->_value), self->_delete, self->_copy, self->_toString);
}

String str_Node(Node self)
{
    return str("Node::ToString() not implemented yet.");
}

struct cList
{
    size_t _size;
    Node *_values;
    delfunc _delete;
    copyfunc _copy;
    strfunc _toString;
};

List new_List(delfunc pDeleteFunction, copyfunc pCopyFunction, strfunc pToStringFunction)
{
    List result = (List)malloc(sizeof(struct cList));

    result->_size = 0;
    result->_values = (Node *)malloc(sizeof(Node) * result->_size);
    result->_delete = pDeleteFunction;
    result->_copy = pCopyFunction;
    result->_toString = pToStringFunction;

    return result;
}

void del_List(List *self)
{
    if (*self == NULL)
        return;

    for (size_t i = 0; i < (*self)->_size; i++)
        del_Node(&((*self)->_values[i]));
    free((*self)->_values);
    free(*self);

    *self = NULL;

    return;
}

List copy_List(List self)
{
    List result = new_List(self->_delete, self->_copy, self->_toString);

    result->_size = self->_size;
    result->_values = realloc(result->_values, sizeof(Node) * result->_size);

    for (size_t i = 0; i < result->_size; i++)
        result->_values[i] = copy_Node(self->_values[i]);

    return result;
}

String str_List(List self)
{
    return str("List::ToString() not implemented yet.");
}

size_t List_GetSize(List self)
{
    return self->_size;
}

void *List_Get(List self, size_t pIndex)
{
    if (!(0 <= pIndex && pIndex < self->_size))
    {
        printf("index out of bounds");
        exit(-1);
    }

    return self->_values[pIndex]->_value;
}

void List_Append(List self, void *pValue)
{
    self->_values = realloc(self->_values, sizeof(Node) * (self->_size + 1));
    self->_values[self->_size] = new_Node(pValue, self->_delete, self->_copy, self->_toString);
    self->_size++;
}

void List_Remove(List self, size_t pIndex)
{
    if (!(0 <= pIndex && pIndex < self->_size))
    {
        printf("index out of bounds");
        exit(-1);
    }

    self->_size--;
    del_Node(&(self->_values[pIndex]));
    for (size_t i = pIndex; i < self->_size; i++)
        self->_values[i] = self->_values[i + 1];
    self->_values = realloc(self->_values, sizeof(Node) * self->_size);
}