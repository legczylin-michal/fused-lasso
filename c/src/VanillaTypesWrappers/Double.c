#include "Double.h"

#include <stdlib.h>

Double new_Double(double value)
{
    Double result = (Double)malloc(sizeof(struct cDouble));

    result->value = value;

    return result;
}

void del_Double(Double *self)
{
    if (*self == NULL)
        return;

    free(*self);

    *self = NULL;

    return;
}

Double copy_Double(Double self)
{
    return new_Double(self->value);
}