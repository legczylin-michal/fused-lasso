#include "SizeT.h"

#include "../basic/String/String.h"

SizeT new_SizeT(size_t value)
{
    SizeT result = (SizeT)malloc(sizeof(struct cSizeT));

    result->value = value;

    return result;
}

void del_SizeT(SizeT *self)
{
    if (*self == NULL)
        return;

    free(*self);

    *self = NULL;

    return;
}

SizeT copy_SizeT(SizeT self)
{
    return new_SizeT(self->value);
}

String str_SizeT(SizeT self)
{
    return str(self->value);
}