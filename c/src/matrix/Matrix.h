#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdlib.h>

struct cMatrix;
typedef struct cMatrix *Matrix;

Matrix new_Matrix(size_t pNumRows, size_t pNumCols);
void del_Matrix(Matrix *self);

size_t Matrix_GetNumRows(Matrix self);
size_t Matrix_GetNumCols(Matrix self);

double Matrix_Get(Matrix self, size_t pRow, size_t pCol);
void Matrix_Set(Matrix self, size_t pRow, size_t pCol, double pNewValue);

Matrix Matrix_Add(Matrix A, Matrix B);
Matrix Matrix_Multiply(Matrix A, Matrix B);

Matrix Matrix_Repeat(Matrix A, size_t n, size_t axis);

void Matrix_Print(Matrix self);

Matrix Matrix_ReadCsv(const char *filepath, const char *separator);
void Matrix_SaveCsv(Matrix self, const char *filepath, const char *separator);

Matrix RandomBinomial(size_t n, double p, size_t rows, size_t cols);
Matrix RandomUniform(double low, double high, size_t rows, size_t cols);
Matrix RandomNormal(double mean, double std, size_t rows, size_t cols);

#endif // !_MATRIX_H_