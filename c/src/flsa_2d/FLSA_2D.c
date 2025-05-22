#include "FLSA_2D.h"

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

int sign(double a)
{
    if (a < 0)
        return -1;

    return a > 0 ? 1 : 0;
}

double maximum(double a, double b)
{
    return a > b ? a : b;
}

double soft_threshold(double a, double b)
{
    return (double)sign(a) * maximum(abs(a) - b, 0);
}

struct cFLSA_2D
{
    double _lambda1;
    double _lambda2;
    size_t _maxIter;
    double _tol;

    size_t _n1;
    size_t _n2;
    size_t _numGroups;
    size_t *_N;
    size_t **_G;
    size_t *_N_neighbours;
    size_t **_neighbours;
    size_t **_w;
    double *_y_bar;
    double *_gamma;
};

FLSA_2D new_FLSA_2D(double pLambda1, double pLambda2, size_t pMaxIter, double pTol)
{
    FLSA_2D result = (FLSA_2D)malloc(sizeof(struct cFLSA_2D));

    result->_lambda1 = pLambda1;
    result->_lambda2 = pLambda2;
    result->_maxIter = pMaxIter;
    result->_tol = pTol;

    result->_n1 = 0;
    result->_n2 = 0;
    result->_numGroups = 0;
    result->_N = (size_t *)malloc(sizeof(size_t) * result->_numGroups);
    result->_G = (size_t **)malloc(sizeof(size_t *) * result->_numGroups);
    result->_N_neighbours = (size_t *)malloc(sizeof(size_t) * result->_numGroups);
    result->_neighbours = (size_t **)malloc(sizeof(size_t *) * result->_numGroups);
    result->_w = (size_t **)malloc(sizeof(size_t *) * result->_numGroups);
    result->_y_bar = (double *)malloc(sizeof(double) * result->_numGroups);
    result->_gamma = (double *)malloc(sizeof(double) * result->_numGroups);

    return result;
}

void del_FLSA_2D(FLSA_2D *self)
{
    if (*self == NULL)
        return;

    free((*self)->_gamma);

    free((*self)->_y_bar);

    for (size_t i = 0; i < (*self)->_numGroups; i++)
        free((*self)->_w[i]);
    free((*self)->_w);

    for (size_t i = 0; i < (*self)->_numGroups; i++)
        free((*self)->_neighbours[i]);
    free((*self)->_neighbours);

    free((*self)->_N_neighbours);

    for (size_t i = 0; i < (*self)->_numGroups; i++)
        free((*self)->_G[i]);
    free((*self)->_G);

    free((*self)->_N);

    free(*self);

    *self = NULL;

    return;
}

// since it is hard to assert equality of floating-point numbers
// we use this function instead
bool FLSA_2D_RoughlyEqual(FLSA_2D self, double a, double b)
{
    return abs(a - b) < self->_tol;
}

// function to check whether coordinates (i, j) are reasonable
// to access value of matrix of n1 by n2 dimensions
bool FLSA_2D_InBounds(FLSA_2D self, int i, int j)
{
    return 0 <= i && i < self->_n1 && 0 <= j && j < self->_n2;
}

// convert 2D position into respective flattened 1D position
size_t FLSA_2D_Coord2Index(FLSA_2D self, size_t i, size_t j)
{
    return i * self->_n2 + j;
}

void FLSA_2D_Fit(FLSA_2D self, Matrix pY)
{
    // set increment for lambda2
    double delta = 1e-3;
    // read shape of incoming data
    self->_n1 = Matrix_GetNumRows(pY);
    self->_n2 = Matrix_GetNumCols(pY);
    // number of initial groups
    self->_numGroups = self->_n1 * self->_n2;
    // list that stores size of every group
    self->_N = realloc(self->_N, sizeof(size_t) * self->_numGroups);
    // list of groups (lists) containing indexes of points that belong to said group (ints)
    // every pixel is its own group at the beginning
    self->_G = realloc(self->_G, sizeof(size_t *) * self->_numGroups);
    for (size_t i = 0; i < self->_numGroups; i++)
    {
        self->_N[i] = 1;
        self->_G[i] = (size_t *)malloc(sizeof(size_t) * self->_N[i]);
        self->_G[i][0] = i;
    }
    // list that stores number of neighbours of every group
    self->_N_neighbours = realloc(self->_N_neighbours, sizeof(size_t) * self->_numGroups);
    // list that stores neighbours for every group
    self->_neighbours = realloc(self->_neighbours, sizeof(size_t *) * self->_numGroups);
    // list that stores shared border lengths of given group with its respective neighbours
    self->_w = realloc(self->_w, sizeof(size_t *) * self->_numGroups);
    for (size_t i = 0; i < self->_n1; i++)
    {
        for (size_t j = 0; j < self->_n2; j++)
        {
            // initiate list of neighbours
            self->_neighbours[FLSA_2D_Coord2Index(self, i, j)] = (size_t *)malloc(sizeof(size_t) * 4);
            self->_w[FLSA_2D_Coord2Index(self, i, j)] = (size_t *)malloc(sizeof(size_t) * 4);
            size_t size = 0;
            // if there is neighbour above
            if (FLSA_2D_InBounds(self, (int)i - 1, j))
            {
                // add it with shared border length 1
                self->_neighbours[FLSA_2D_Coord2Index(self, i, j)][size] = FLSA_2D_Coord2Index(self, i - 1, j);
                self->_w[FLSA_2D_Coord2Index(self, i, j)][size] = 1;
                size++;
            }

            if (FLSA_2D_InBounds(self, i + 1, j))
            {
                self->_neighbours[FLSA_2D_Coord2Index(self, i, j)][size] = FLSA_2D_Coord2Index(self, i + 1, j);
                self->_w[FLSA_2D_Coord2Index(self, i, j)][size] = 1;
                size++;
            }

            if (FLSA_2D_InBounds(self, i, (int)j - 1))
            {
                self->_neighbours[FLSA_2D_Coord2Index(self, i, j)][size] = FLSA_2D_Coord2Index(self, i, j - 1);
                self->_w[FLSA_2D_Coord2Index(self, i, j)][size] = 1;
                size++;
            }

            if (FLSA_2D_InBounds(self, i, j + 1))
            {
                self->_neighbours[FLSA_2D_Coord2Index(self, i, j)][size] = FLSA_2D_Coord2Index(self, i, j + 1);
                self->_w[FLSA_2D_Coord2Index(self, i, j)][size] = 1;
                size++;
            }
            self->_N_neighbours[FLSA_2D_Coord2Index(self, i, j)] = size;
            self->_neighbours[FLSA_2D_Coord2Index(self, i, j)] = realloc(self->_neighbours[FLSA_2D_Coord2Index(self, i, j)], sizeof(size_t) * size);
            self->_w[FLSA_2D_Coord2Index(self, i, j)] = realloc(self->_w[FLSA_2D_Coord2Index(self, i, j)], sizeof(size_t) * size);
        }
    }
    // averga value of data at every group
    self->_y_bar = realloc(self->_y_bar, sizeof(double) * self->_numGroups);
    for (size_t i = 0; i < self->_n1; i++)
        for (size_t j = 0; j < self->_n2; j++)
            self->_y_bar[FLSA_2D_Coord2Index(self, i, j)] = Matrix_Get(pY, i, j);
    // initial values of gamma for lambda2=0
    self->_gamma = realloc(self->_gamma, sizeof(double) * self->_numGroups);
    for (size_t i = 0; i < self->_numGroups; i++)
        self->_gamma[i] = soft_threshold(self->_y_bar[i], self->_lambda1);

    // checking if initiation works properly
    // for (size_t i = 0; i < self->_numGroups; i++)
    // {
    //     printf("%d-th group has %d elements:\n", i, self->_N[i]);
    //     for (size_t j = 0; j < self->_N[i]; j++)
    //         printf("\t%d", self->_G[i][j]);
    //     printf("\n%d-th group has %d neighbours:\n", i, self->_N_neighbours[i]);
    //     for (size_t j = 0; j < self->_N_neighbours[i]; j++)
    //         printf("\t%d with shared border length %d", self->_neighbours[i][j], self->_w[i][j]);
    //     printf("\n%d-th group's y_bar=%.7lf and gamma=%.7lf", i, self->_y_bar[i], self->_gamma[i]);
    //     printf("\n");
    // }

    // smooth cycle — increment lambda2 from 0 up to self._lambda2 by ~delta
    double step = (self->_lambda2 - delta) / ceil((self->_lambda2 - delta) / delta);
    for (float lambda2 = delta; lambda2 <= self->_lambda2; lambda2 += step)
    {
        // convergence key
        bool converged = false;
        // retry cycles until converged
        for (size_t _ = 0; _ < self->_maxIter; _++)
        {
        }
        // if algorithm did not converge
        if (!converged)
        {
            // show log
            // printf("algorithm did not converge");
        }
        // iterate over every group
        size_t k = 0;
        while (k < self->_numGroups)
        {
            // ...
            k++;
        }
    }
}

Matrix FLSA_2D_GetCoef(FLSA_2D self)
{
    return new_Matrix(self->_n1, self->_n2);
}