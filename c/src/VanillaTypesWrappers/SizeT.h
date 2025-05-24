#ifndef _SIZE_T_H_
#define _SIZE_T_H_

#include <stdlib.h>

struct cSizeT
{
    size_t value;
};
typedef struct cSizeT *SizeT;

SizeT new_SizeT(size_t value);
void del_SizeT(SizeT *self);
SizeT copy_SizeT(SizeT self);

#endif // !_SIZE_T_H_