#include "Matrix.h"

#include <stdio.h>
#include <math.h>

struct cMatrix
{
    size_t _numRows;
    size_t _numCols;
    size_t _length;
    double *_values;
};

Matrix new_Matrix(size_t pNumRows, size_t pNumCols)
{
    Matrix result = (Matrix)malloc(sizeof(struct cMatrix));

    result->_numRows = pNumRows;
    result->_numCols = pNumCols;
    result->_length = result->_numRows * result->_numCols;
    result->_values = (double *)malloc(sizeof(double) * result->_length);
    for (size_t i = 0; i < result->_length; i++)
        result->_values[i] = 0.0f;

    return result;
}

void del_Matrix(Matrix *self)
{
    if (*self == NULL)
        return;

    free((*self)->_values);
    free(*self);

    *self = NULL;

    return;
}

size_t Matrix_GetNumRows(Matrix self)
{
    return self->_numRows;
}

size_t Matrix_GetNumCols(Matrix self)
{
    return self->_numCols;
}

size_t Matrix_Coord2Index(Matrix self, size_t pRow, size_t pCol)
{
    return pRow * self->_numCols + pCol;
}

double Matrix_Get(Matrix self, size_t pRow, size_t pCol)
{
    return self->_values[Matrix_Coord2Index(self, pRow, pCol)];
}

void Matrix_Set(Matrix self, size_t pRow, size_t pCol, double pNewValue)
{
    self->_values[Matrix_Coord2Index(self, pRow, pCol)] = pNewValue;
}

Matrix Matrix_Add(Matrix A, Matrix B)
{
    if (!(A->_numRows == B->_numRows && A->_numCols == B->_numCols))
    {
        printf("unmatching dimensions");
        exit(-1);
    }

    Matrix result = new_Matrix(A->_numRows, A->_numCols);

    for (size_t i = 0; i < A->_length; i++)
        result->_values[i] = A->_values[i] + B->_values[i];

    return result;
}

Matrix Matrix_Multiply(Matrix A, Matrix B)
{
    if (!(A->_numRows == B->_numRows && A->_numCols == B->_numCols))
    {
        printf("unmatching dimensions");
        exit(-1);
    }

    Matrix result = new_Matrix(A->_numRows, A->_numCols);

    for (size_t i = 0; i < A->_length; i++)
        result->_values[i] = A->_values[i] * B->_values[i];

    return result;
}

Matrix Matrix_Repeat(Matrix A, size_t n, size_t axis)
{
    Matrix result = new_Matrix(A->_numRows * (axis == 0 ? n : 1), A->_numCols * (axis == 1 ? n : 1));

    for (size_t i = 0; i < result->_numRows; i++)
        for (size_t j = 0; j < result->_numCols; j++)
            result->_values[Matrix_Coord2Index(result, i, j)] = A->_values[Matrix_Coord2Index(A, i / (axis == 0 ? n : 1), j / (axis == 1 ? n : 1))];

    return result;
}

void Matrix_Print(Matrix self)
{
    for (size_t row = 0; row < self->_numRows; row++)
    {
        for (size_t col = 0; col < self->_numCols; col++)
            printf("%.7lf ", self->_values[Matrix_Coord2Index(self, row, col)]);

        printf("\n");
    }
}

char *doubleAsString(double value)
{
    int bs = snprintf(NULL, 0, "%lf", value);
    char *result = (char *)malloc(sizeof(char) * (bs + 1));
    snprintf(result, bs + 1, "%lf", value);

    return result;
}

void Matrix_SaveCsv(Matrix self, const char *filepath, const char *separator)
{
    FILE *fptr = fopen(filepath, "w");

    if (fptr == NULL)
    {
        perror("Error opening file");
        exit(-1);
    }

    for (size_t row = 0; row < self->_numRows; row++)
    {
        for (size_t col = 0; col < self->_numCols; col++)
        {
            char *str = doubleAsString(self->_values[Matrix_Coord2Index(self, row, col)]);
            fputs(str, fptr);
            free(str);
            if (col == self->_numCols - 1)
                break;
            fputs(separator, fptr);
        }
        if (row == self->_numRows - 1)
            break;
        fputc('\n', fptr);
    }

    fclose(fptr);
}

size_t binomial(size_t n, double p)
{
    size_t successes = 0;

    for (size_t i = 0; i < n; i++)
    {
        double random_value = (double)rand() / RAND_MAX;

        if (random_value < p)
            successes++;
    }

    return successes;
}

double uniform(double low, double high)
{
    return low + ((double)rand() / RAND_MAX) * (high - low);
}

double standard_normal()
{
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;

    // apply Box-Muller transform
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);

    return z0;
}

double normal(double mean, double std)
{
    return mean + standard_normal() * std;
}

Matrix RandomBinomial(size_t n, double p, size_t rows, size_t cols)
{
    Matrix result = new_Matrix(rows, cols);

    for (size_t i = 0; i < result->_length; i++)
        result->_values[i] = binomial(n, p);

    return result;
}

Matrix RandomUniform(double low, double high, size_t rows, size_t cols)
{
    Matrix result = new_Matrix(rows, cols);

    for (size_t i = 0; i < result->_length; i++)
        result->_values[i] = uniform(low, high);

    return result;
}

Matrix RandomNormal(double mean, double std, size_t rows, size_t cols)
{
    Matrix result = new_Matrix(rows, cols);

    for (size_t i = 0; i < result->_length; i++)
        result->_values[i] = normal(mean, std);

    return result;
}