#ifndef _FLSA_2D_H_
#define _FLSA_2D_H_

#include <stdlib.h>

#include "..\matrix\Matrix.h"

struct cFLSA_2D;
typedef struct cFLSA_2D *FLSA_2D;

FLSA_2D new_FLSA_2D(double pLambda1, double pLambda2, size_t pMaxIter, double pTol);
void del_FLSA_2D(FLSA_2D *self);

void FLSA_2D_Fit(FLSA_2D self, Matrix pY);

Matrix FLSA_2D_GetCoef(FLSA_2D self);

#endif //!_FLSA_2D_H_