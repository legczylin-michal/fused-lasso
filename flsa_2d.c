#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void flsa_2d(double lambda_1, double lambda_2, size_t max_iter, double tol, double **y, size_t n1, size_t n2);

int main(int argc, char const *argv[])
{
    int n1 = 128, n2 = 128;
    double **y = (double **)malloc(sizeof(double *) * n1 * n2);
    for (int i = 0; i < n1; i++)
    {
        y[i] = (double *)malloc(sizeof(double) * n2);
        for (int j = 0; j < n2; j++)
            y[i][j] = 1;
    }

    flsa_2d(1.0, 1.0, 1000, 0.0001, y, n1, n2);

    for (int i = 0; i < n1; i++)
        free(y[i]);
    free(y);

    printf("\nSUCCESSFULLY COMPILED\n");

    return 0;
}

bool in_bounds(int i, int j, int n1, int n2)
{
    return 0 <= i && i < n1 && 0 <= j && j < n2;
}

int coord2index(int i, int j, int n1, int n2)
{
    return i * n1 + j;
}

int sign(double a)
{
    if (a < 0)
        return -1;
    else if (a > 0)
        return 1;

    return 0;
}

double maximum(double a, double b)
{
    return a > b ? a : b;
}

double soft_threshold(double a, double b)
{
    return sign(a) * maximum(abs(a) - b, 0);
}

void flsa_2d(double lambda_1, double lambda_2, size_t max_iter, double tol, double **y, size_t n1, size_t n2)
{
    double delta = 1e-9;

    // initiate set of groups
    // every pixel is its own group
    // store indexes of pixels — pixel p=(i, j) has index i * n1 + j
    int **G = (int **)malloc(sizeof(int *) * n1 * n2);
    for (int i = 0; i < n1 * n2; i++)
    {
        G[i] = (int *)malloc(sizeof(int) * 1);
        G[i][0] = i;
    }
    // set of neighbours
    int *Nneighbours = (int *)malloc(sizeof(int) * n1 * n2);
    int **neighbours = (int **)malloc(sizeof(int *) * n1 * n2);
    for (int i = 0; i < n1; i++)
    {
        for (int j = 0; j < n2; j++)
        {
            int size = 0;
            neighbours[coord2index(i, j, n1, n2)] = (int *)malloc(sizeof(int) * size);

            if (in_bounds(i - 1, j, n1, n2))
            {
                size++;
                neighbours[coord2index(i, j, n1, n2)] = realloc(neighbours[coord2index(i, j, n1, n2)], sizeof(int) * size);
                neighbours[coord2index(i, j, n1, n2)][size - 1] = coord2index(i - 1, j, n1, n2);
            }

            if (in_bounds(i + 1, j, n1, n2))
            {
                size++;
                neighbours[coord2index(i, j, n1, n2)] = realloc(neighbours[coord2index(i, j, n1, n2)], sizeof(int) * size);
                neighbours[coord2index(i, j, n1, n2)][size - 1] = coord2index(i + 1, j, n1, n2);
            }

            if (in_bounds(i, j - 1, n1, n2))
            {
                size++;
                neighbours[coord2index(i, j, n1, n2)] = realloc(neighbours[coord2index(i, j, n1, n2)], sizeof(int) * size);
                neighbours[coord2index(i, j, n1, n2)][size - 1] = coord2index(i, j - 1, n1, n2);
            }

            if (in_bounds(i, j + 1, n1, n2))
            {
                size++;
                neighbours[coord2index(i, j, n1, n2)] = realloc(neighbours[coord2index(i, j, n1, n2)], sizeof(int) * size);
                neighbours[coord2index(i, j, n1, n2)][size - 1] = coord2index(i, j + 1, n1, n2);
            }

            Nneighbours[coord2index(i, j, n1, n2)] = size;
        }
    }

    int *N = (int *)malloc(sizeof(int) * n1 * n2);
    for (int i = 0; i < n1 * n2; i++)
        N[i] = 1;

    double *y_bar = (double *)malloc(sizeof(double) * n1 * n2);
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++)
            y_bar[coord2index(i, j, n1, n2)] = y[i][j];

    double *gamma = (double *)malloc(sizeof(double) * n1 * n2);

    for (int i = 0; i < n1 * n2; i++)
        gamma[i] = soft_threshold(y_bar[i], lambda_1);

    // for (int i = 0; i < n1 * n2; i++)
    // {
    //     printf("%d-th group has %d elements:\n", i, N[i]);

    //     for (int j = 0; j < N[i]; j++)
    //     {
    //         printf("\t%d", G[i][j]);
    //     }

    //     printf("\n");
    // }

    // for (int i = 0; i < n1 * n2; i++)
    // {
    //     printf("%d-th group has %d neighbours:\n", i, Nneighbours[i]);

    //     for (int j = 0; j < Nneighbours[i]; j++)
    //     {
    //         printf("\t%d", neighbours[i][j]);
    //     }

    //     printf("\n");
    // }

    free(gamma);
    free(y_bar);
    free(N);
    for (int i = 0; i < n1 * n2; i++)
        free(neighbours[i]);
    free(neighbours);
    free(Nneighbours);
    for (int i = 0; i < n1 * n2; i++)
        free(G[i]);
    free(G);
}