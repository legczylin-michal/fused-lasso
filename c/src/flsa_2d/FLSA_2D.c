#include "FLSA_2D.h"

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "..\List\List.h"
#include "..\VanillaTypesWrappers\SizeT.h"
#include "..\VanillaTypesWrappers\Double.h"

double absolute(double a)
{
    return a < 0 ? -a : a;
}

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
    return (double)sign(a) * maximum(absolute(a) - b, 0);
}

void swap(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

void sort(List values)
{
    for (size_t i = 0; i < List_GetSize(values); i++)
    {
        size_t indexOfMin = i;
        for (size_t j = i; j < List_GetSize(values); j++)
            if (((SizeT)List_Get(values, j))->value < ((SizeT)List_Get(values, indexOfMin))->value)
                indexOfMin = j;

        size_t tmp = ((SizeT)List_Get(values, i))->value;
        ((SizeT)List_Get(values, i))->value = ((SizeT)List_Get(values, indexOfMin))->value;
        ((SizeT)List_Get(values, indexOfMin))->value = tmp;
    }
}

void sortTwoArraysAccordingToFirstOne(List values1, List values2)
{
    if (!(List_GetSize(values1) == List_GetSize(values2)))
    {
        printf("unmatched sizes");
        exit(-1);
    }

    for (size_t i = 0; i < List_GetSize(values1); i++)
    {
        size_t indexOfMin = i;
        for (size_t j = i; j < List_GetSize(values1); j++)
            if (((Double)List_Get(values1, j))->value < ((Double)List_Get(values1, indexOfMin))->value)
                indexOfMin = j;

        double tmp = ((Double)List_Get(values1, i))->value;
        ((Double)List_Get(values1, i))->value = ((Double)List_Get(values1, indexOfMin))->value;
        ((Double)List_Get(values1, indexOfMin))->value = tmp;

        tmp = ((Double)List_Get(values2, i))->value;
        ((Double)List_Get(values2, i))->value = ((Double)List_Get(values2, indexOfMin))->value;
        ((Double)List_Get(values2, indexOfMin))->value = tmp;
    }
}

bool solve(List breakpoints, List penalties, double bias, double *buffer)
{
    if (!(List_GetSize(breakpoints) == List_GetSize(penalties)))
    {
        printf("unmatched sizes");
        exit(-1);
    }

    double v = bias;
    for (size_t i = 0; i < List_GetSize(penalties); i++)
        v += ((Double)List_Get(penalties, i))->value;

    if (v <= ((Double)List_Get(breakpoints, 0))->value)
    {
        *buffer = v;
        return true;
    }

    v -= 2 * ((Double)List_Get(penalties, 0))->value;

    for (size_t i = 1; i < List_GetSize(breakpoints); i++)
    {
        if (((Double)List_Get(breakpoints, i - 1))->value < v && v <= ((Double)List_Get(breakpoints, i))->value)
        {
            *buffer = v;
            return true;
        }

        v -= 2 * ((Double)List_Get(penalties, i))->value;
    }

    if (((Double)List_Get(breakpoints, List_GetSize(breakpoints) - 1))->value < v)
    {
        *buffer = v;
        return true;
    }

    *buffer = 0.0f;
    return false;
}

size_t IndexAfterRemoval(size_t index, List removed)
{
    size_t i = 0;
    for (size_t j = 0; j < List_GetSize(removed); j++)
        if (index > ((SizeT)List_Get(removed, j))->value)
            i++;

    return index - i;
}

struct cFLSA_2D
{
    double _lambda1;
    double _lambda2;
    size_t _maxIter;
    double _tol;

    size_t _n1;
    size_t _n2;
    List _N;
    List _G;
    List _neighbours;
    List _w;
    List _y_bar;
    List _gamma;
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
    result->_N = new_List((delfunc)del_SizeT, (copyfunc)copy_SizeT);
    result->_G = new_List((delfunc)del_List, (copyfunc)copy_List);
    result->_neighbours = new_List((delfunc)del_List, (copyfunc)copy_List);
    result->_w = new_List((delfunc)del_List, (copyfunc)copy_List);
    result->_y_bar = new_List((delfunc)del_Double, (copyfunc)copy_Double);
    result->_gamma = new_List((delfunc)del_Double, (copyfunc)copy_Double);

    return result;
}

void del_FLSA_2D(FLSA_2D *self)
{
    if (*self == NULL)
        return;

    del_List(&((*self)->_N));
    del_List(&((*self)->_G));
    del_List(&((*self)->_neighbours));
    del_List(&((*self)->_w));
    del_List(&((*self)->_y_bar));
    del_List(&((*self)->_gamma));

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

double FLSA_2D_Loss(FLSA_2D self, List gamma)
{
    double firstTerm = 0.0f;
    double secondTerm = 0.0f;
    double thirdTerm = 0.0f;

    for (size_t k = 0; k < List_GetSize(self->_G); k++)
    {
        firstTerm += ((SizeT)List_Get(self->_N, k))->value * pow(((Double)List_Get(self->_y_bar, k))->value - ((Double)List_Get(gamma, k))->value, 2);
        secondTerm += ((SizeT)List_Get(self->_N, k))->value * abs(((Double)List_Get(gamma, k))->value);
        for (size_t i = 0; i < List_GetSize((List)List_Get(self->_neighbours, k)); i++)
        {
            size_t k_prime = ((SizeT)List_Get((List)List_Get(self->_neighbours, k), i))->value;
            thirdTerm += ((SizeT)List_Get((List)List_Get(self->_w, k), i))->value * abs(((Double)List_Get(gamma, k))->value - ((Double)List_Get(gamma, k_prime))->value);
        }
    }

    return 0.5 * firstTerm + self->_lambda1 * secondTerm + self->_lambda2 / 2 * thirdTerm;
}

void FLSA_2D_SetParametersBasedOnY(FLSA_2D self, Matrix pY)
{
    // read shape of incoming data
    self->_n1 = Matrix_GetNumRows(pY);
    self->_n2 = Matrix_GetNumCols(pY);
    for (size_t k = 0; k < self->_n1 * self->_n2; k++)
    {
        List kThGroup = new_List((delfunc)del_SizeT, (copyfunc)copy_SizeT);
        List_Append(kThGroup, new_SizeT(k));
        List_Append(self->_G, kThGroup);
        List_Append(self->_N, new_SizeT(1));
    }
    for (size_t i = 0; i < self->_n1; i++)
    {
        for (size_t j = 0; j < self->_n2; j++)
        {
            List_Append(self->_neighbours, new_List((delfunc)del_SizeT, (copyfunc)copy_SizeT));
            List_Append(self->_w, new_List((delfunc)del_SizeT, (copyfunc)copy_SizeT));
            // if there is neighbour above
            if (FLSA_2D_InBounds(self, (int)i - 1, j))
            {
                List_Append((List)List_Get(self->_neighbours, FLSA_2D_Coord2Index(self, i, j)), new_SizeT(FLSA_2D_Coord2Index(self, i - 1, j)));
                List_Append((List)List_Get(self->_w, FLSA_2D_Coord2Index(self, i, j)), new_SizeT(1));
            }
            if (FLSA_2D_InBounds(self, i + 1, j))
            {
                List_Append((List)List_Get(self->_neighbours, FLSA_2D_Coord2Index(self, i, j)), new_SizeT(FLSA_2D_Coord2Index(self, i + 1, j)));
                List_Append((List)List_Get(self->_w, FLSA_2D_Coord2Index(self, i, j)), new_SizeT(1));
            }
            if (FLSA_2D_InBounds(self, i, (int)j - 1))
            {
                List_Append((List)List_Get(self->_neighbours, FLSA_2D_Coord2Index(self, i, j)), new_SizeT(FLSA_2D_Coord2Index(self, i, j - 1)));
                List_Append((List)List_Get(self->_w, FLSA_2D_Coord2Index(self, i, j)), new_SizeT(1));
            }
            if (FLSA_2D_InBounds(self, i, j + 1))
            {
                List_Append((List)List_Get(self->_neighbours, FLSA_2D_Coord2Index(self, i, j)), new_SizeT(FLSA_2D_Coord2Index(self, i, j + 1)));
                List_Append((List)List_Get(self->_w, FLSA_2D_Coord2Index(self, i, j)), new_SizeT(1));
            }
        }
    }
    // average value of data at every group
    for (size_t i = 0; i < self->_n1; i++)
        for (size_t j = 0; j < self->_n2; j++)
            List_Append(self->_y_bar, new_Double(Matrix_Get(pY, i, j)));
    // initial values of gamma for lambda2=0
    for (size_t k = 0; k < List_GetSize(self->_G); k++)
        List_Append(self->_gamma, new_Double(soft_threshold(((Double)List_Get(self->_y_bar, k))->value, self->_lambda1)));
}

void FLSA_2D_Print(FLSA_2D self)
{
    printf("lambda1=%.5lf, lambda2=%.5lf, max_iter=%d, tol=%.5lf\n", self->_lambda1, self->_lambda2, self->_maxIter, self->_tol);
    printf("shape=(%d, %d)\n", self->_n1, self->_n2);
    printf("there are %d groups\n", List_GetSize(self->_G));
    for (size_t k = 0; k < List_GetSize(self->_G); k++)
    {
        printf("%d-th group => [", k);
        for (size_t i = 0; i < ((SizeT)List_Get(self->_N, k))->value; i++)
            printf("\t%d, ", ((SizeT)List_Get((List)List_Get(self->_G, k), i))->value);
        printf("], its neighbours => [");
        for (size_t i = 0; i < List_GetSize((List)List_Get(self->_neighbours, k)); i++)
        {
            printf("\t%d => %d, ", ((SizeT)List_Get((List)List_Get(self->_neighbours, k), i))->value, ((SizeT)List_Get((List)List_Get(self->_w, k), i))->value);
        }
        printf("], y_bar=%lf, gamma=%lf\n", ((Double)List_Get(self->_y_bar, k))->value, ((Double)List_Get(self->_gamma, k))->value);
    }
    printf("Example loss value: %lf\n", FLSA_2D_Loss(self, self->_gamma));
}

void FLSA_2D_Fit(FLSA_2D self, Matrix pY)
{

    // set increment for lambda2
    double delta = 1e-1;
    // initiate parameters based on y data matrix
    FLSA_2D_SetParametersBasedOnY(self, pY);
    // assert correctness of initialisation
    // FLSA_2D_Print(self);
    // smooth cycle — increment lambda2 from 0 up to self._lambda2 by ~delta
    double step = (self->_lambda2 - delta) / ceil((self->_lambda2 - delta) / delta);
    for (double lambda2 = delta; lambda2 <= self->_lambda2; lambda2 += step)
    {
        // convergence key
        bool converged = false;
        // retry cycles until converged
        for (size_t _ = 0; _ < self->_maxIter; _++)
        {
            // copy current value of gamma so to compare if there was a successful update
            List currentGamma = copy_List(self->_gamma);
            // iterate over every coordinate
            for (size_t k = 0; k < List_GetSize(self->_G); k++)
            {
                // descent part — try to optimise only in gamma[k] direction
                // derivative is piece-wise linear at breakpoints
                List breakpoints = new_List((delfunc)del_Double, (copyfunc)copy_Double);
                List_Append(breakpoints, new_Double(0.0f));
                for (size_t i = 0; i < List_GetSize((List)List_Get(self->_neighbours, k)); i++)
                    List_Append(breakpoints, new_Double(((Double)List_Get(self->_gamma, ((SizeT)List_Get((List)List_Get(self->_neighbours, k), i))->value))->value));
                // with scaling parameters stored in penalties
                List penalties = new_List((delfunc)del_Double, (copyfunc)copy_Double);
                List_Append(penalties, new_Double(self->_lambda1));
                for (size_t i = 0; i < List_GetSize((List)List_Get(self->_neighbours, k)); i++)
                    List_Append(penalties, new_Double(lambda2 * ((SizeT)List_Get((List)List_Get(self->_w, k), i))->value / List_GetSize((List)List_Get(self->_G, k))));
                // we need to sort by breakpoints
                sortTwoArraysAccordingToFirstOne(breakpoints, penalties);
                // and since it is piece-wise linear, solution is unique if one exists
                double solution;
                bool solutionFound = solve(breakpoints, penalties, ((Double)List_Get(self->_y_bar, k))->value, &solution);
                // copy current value of gamma to compare if update with found solution for gamma[k] is better than before
                List tmpGamma = copy_List(self->_gamma);
                // if there is no solution (0 falls in one of the breaks)
                if (!solutionFound)
                {
                    double minValue = INFINITY;
                    // we check what value gives loss function with gamma[k] set to breakpoint
                    // and select such value that gives the smallest value of loss
                    for (size_t i = 0; i < List_GetSize(breakpoints); i++)
                    {
                        ((Double)List_Get(tmpGamma, k))->value = ((Double)List_Get(breakpoints, i))->value;
                        // evaluate
                        double value = FLSA_2D_Loss(self, tmpGamma);
                        if (value < minValue)
                        {
                            minValue = value;
                            solution = ((Double)List_Get(breakpoints, i))->value;
                        }
                    }
                }
                // so we substitute to check if such an update is good
                ((Double)List_Get(tmpGamma, k))->value = solution;
                // if it is better
                bool betterSolFound = false;
                if (FLSA_2D_Loss(self, tmpGamma) < FLSA_2D_Loss(self, self->_gamma))
                {
                    // then use this value
                    del_List(&self->_gamma);
                    self->_gamma = copy_List(tmpGamma);
                    betterSolFound = true;
                }
                // deallocate memory
                del_List(&tmpGamma);
                del_List(&penalties);
                del_List(&breakpoints);
                // and continue to the next coordinate
                if (betterSolFound)
                    continue;
                // if there was no successful one-parameter-at-a-time update we consider fusion part
                // iterate over every neighbouring group of k-th group and consider fusion
                for (size_t i = 0; i < List_GetSize((List)List_Get(self->_neighbours, k)); i++)
                {
                    size_t k_prime = ((SizeT)List_Get((List)List_Get(self->_neighbours, k), i))->value;

                    // simulate size of potential group
                    size_t N_m = ((SizeT)List_Get(self->_N, k))->value + ((SizeT)List_Get(self->_N, k_prime))->value;
                    // simulate average data value at potential group
                    double y_bar_m = (((SizeT)List_Get(self->_N, k))->value * ((Double)List_Get(self->_y_bar, k))->value + ((SizeT)List_Get(self->_N, k_prime))->value * ((Double)List_Get(self->_y_bar, k_prime))->value) / N_m;
                    //
                    List neighbours_m = new_List((delfunc)del_SizeT, (copyfunc)copy_SizeT);
                    List w_m = new_List((delfunc)del_SizeT, (copyfunc)copy_SizeT);
                    // take neighbours of k-th group (exclude k_prime-th group)
                    for (size_t j = 0; j < List_GetSize((List)List_Get(self->_neighbours, k)); j++)
                    {
                        size_t curNeighbour = ((SizeT)List_Get((List)List_Get(self->_neighbours, k), j))->value;

                        if (curNeighbour != k_prime)
                        {
                            List_Append(neighbours_m, new_SizeT(curNeighbour));
                            List_Append(w_m, new_SizeT(((SizeT)List_Get((List)List_Get(self->_w, k), j))->value));
                        }
                    }
                    // take neighbours of k_prime-th group (exclude k-th group)
                    for (size_t j = 0; j < List_GetSize((List)List_Get(self->_neighbours, k_prime)); j++)
                    {
                        size_t kPrimeCurNeighbour = ((SizeT)List_Get((List)List_Get(self->_neighbours, k_prime), j))->value;

                        if (kPrimeCurNeighbour != k)
                        {
                            bool kPrimeCurNeighbourIsNeighbourOfK = false;
                            size_t kPrimeCurNeighbourIndex;
                            for (size_t l = 0; l < List_GetSize(neighbours_m); l++)
                            {
                                // if group with such index was already added
                                if (((SizeT)List_Get(neighbours_m, l))->value == kPrimeCurNeighbour)
                                {
                                    kPrimeCurNeighbourIsNeighbourOfK = true;
                                    kPrimeCurNeighbourIndex = l;
                                    break;
                                }
                            }

                            if (kPrimeCurNeighbourIsNeighbourOfK)
                                // then just update shared border length
                                ((SizeT)List_Get(w_m, kPrimeCurNeighbourIndex))->value += ((SizeT)List_Get((List)List_Get(self->_w, k_prime), j))->value;
                            // otherwise
                            else
                            {
                                // add as new neighbour
                                List_Append(neighbours_m, new_SizeT(kPrimeCurNeighbour));
                                List_Append(w_m, new_SizeT(((SizeT)List_Get((List)List_Get(self->_w, k_prime), j))->value));
                            }
                        }
                    }

                    // descent part — try to optimise only in gamma[k] direction
                    // derivative is piece-wise linear at breakpoints
                    List breakpoints = new_List((delfunc)del_Double, (copyfunc)copy_Double);
                    List_Append(breakpoints, new_Double(0.0f));
                    for (size_t j = 0; j < List_GetSize(neighbours_m); j++)
                        List_Append(breakpoints, new_Double(((Double)List_Get(self->_gamma, ((SizeT)List_Get(neighbours_m, j))->value))->value));
                    // with scaling parameters stored in penalties
                    List penalties = new_List((delfunc)del_Double, (copyfunc)copy_Double);
                    List_Append(penalties, new_Double(self->_lambda1));
                    for (size_t j = 0; j < List_GetSize(w_m); j++)
                        List_Append(penalties, new_Double(lambda2 * ((SizeT)List_Get(w_m, j))->value / N_m));
                    // we need to sort by breakpoints
                    sortTwoArraysAccordingToFirstOne(breakpoints, penalties);
                    // and since it is piece-wise linear, solution is unique if one exists
                    double solution;
                    bool solutionFound = solve(breakpoints, penalties, y_bar_m, &solution);
                    // copy current value of gamma to compare if update with found solution for gamma[k] is better than before
                    List tmpGamma = copy_List(self->_gamma);
                    // if there is no solution (0 falls in one of the breaks)
                    if (!solutionFound)
                    {
                        double minValue = INFINITY;
                        // we check what value gives loss function with gamma[k] set to breakpoint
                        // and select such value that gives the smallest value of loss
                        for (size_t j = 0; j < List_GetSize(breakpoints); j++)
                        {
                            ((Double)List_Get(tmpGamma, k))->value = ((Double)List_Get(breakpoints, j))->value;
                            ((Double)List_Get(tmpGamma, k_prime))->value = ((Double)List_Get(breakpoints, j))->value;
                            // evaluate
                            double value = FLSA_2D_Loss(self, tmpGamma);
                            if (value < minValue)
                            {
                                minValue = value;
                                solution = ((Double)List_Get(breakpoints, j))->value;
                            }
                        }
                    }
                    // so we substitute to check if such an update is good
                    ((Double)List_Get(tmpGamma, k))->value = solution;
                    ((Double)List_Get(tmpGamma, k_prime))->value = solution;
                    // if it is better
                    if (FLSA_2D_Loss(self, tmpGamma) < FLSA_2D_Loss(self, self->_gamma))
                    {
                        // then use this value
                        del_List(&self->_gamma);
                        self->_gamma = copy_List(tmpGamma);
                    }
                    // deallocate memory
                    del_List(&tmpGamma);
                    del_List(&penalties);
                    del_List(&breakpoints);
                    del_List(&w_m);
                    del_List(&neighbours_m);
                }
            }
            // when previous value of gamma is not substantailly diffferent from newly acquired one
            double norm = 0.0f;
            for (size_t k = 0; k < List_GetSize(self->_G); k++)
                norm += pow(((Double)List_Get(currentGamma, k))->value - ((Double)List_Get(self->_gamma, k))->value, 2);
            if (sqrt(norm) < self->_tol)
                // then we found solution for given lambda2
                converged = true;
            // deallocate memory
            del_List(&currentGamma);
            // since solution was found, we can move to fusion
            if (converged)
                break;
        }
        // if algorithm did not converge
        if (!converged)
            // show log
            printf("algorithm did not converge\n");
        // iterate over every group
        size_t k = 0;
        while (k < List_GetSize(self->_G))
        {
            // store indexes of neighbours which are to be fused
            List indexesOfGroupsToFuse = new_List((delfunc)del_SizeT, (copyfunc)copy_SizeT);
            // iterate over every neighbour
            for (size_t i = 0; i < List_GetSize((List)List_Get(self->_neighbours, k)); i++)
            {
                size_t k_prime = ((SizeT)List_Get((List)List_Get(self->_neighbours, k), i))->value;

                // if calculated values of gamma for group k and respective neighbour are close enough and non-zero
                if (FLSA_2D_RoughlyEqual(self, ((Double)List_Get(self->_gamma, k))->value, ((Double)List_Get(self->_gamma, k_prime))->value) && !FLSA_2D_RoughlyEqual(self, ((Double)List_Get(self->_gamma, k))->value, 0))
                {
                    // then we want to fuse these two together
                    List_Append(indexesOfGroupsToFuse, new_SizeT(k_prime));
                    // we fuse neighbour into k-th group!
                    // update sizes of groups
                    ((SizeT)List_Get(self->_N, k))->value += ((SizeT)List_Get(self->_N, k_prime))->value;
                    // update groups — union of points
                    for (size_t j = 0; j < List_GetSize((List)List_Get(self->_G, k_prime)); j++)
                        List_Append((List)List_Get(self->_G, k), new_SizeT(((SizeT)List_Get((List)List_Get(self->_G, k_prime), j))->value));
                    // update average values of data of groups
                    ((Double)List_Get(self->_y_bar, k))->value = (((Double)List_Get(self->_y_bar, k))->value * (((SizeT)List_Get(self->_N, k))->value - ((SizeT)List_Get(self->_N, k_prime))->value) + ((Double)List_Get(self->_y_bar, k_prime))->value * ((SizeT)List_Get(self->_N, k_prime))->value) / ((SizeT)List_Get(self->_N, k))->value;
                    // take neighbours of k_prime-th group (exclude k-th group)
                    for (size_t j = 0; j < List_GetSize((List)List_Get(self->_neighbours, k_prime)); j++)
                    {
                        size_t kPrimeCurNeighbour = ((SizeT)List_Get((List)List_Get(self->_neighbours, k_prime), j))->value;

                        if (kPrimeCurNeighbour != k)
                        {
                            bool kPrimeCurNeighbourIsNeighbourOfK = false;
                            size_t kPrimeCurNeighbourIndex;
                            for (size_t l = 0; l < List_GetSize((List)List_Get(self->_neighbours, k)); l++)
                            {
                                // if group with such index was already added
                                if (((SizeT)List_Get((List)List_Get(self->_neighbours, k), l))->value == kPrimeCurNeighbour)
                                {
                                    kPrimeCurNeighbourIsNeighbourOfK = true;
                                    kPrimeCurNeighbourIndex = l;
                                    break;
                                }
                            }

                            if (kPrimeCurNeighbourIsNeighbourOfK)
                            {
                                ((SizeT)List_Get((List)List_Get(self->_w, k), kPrimeCurNeighbourIndex))->value += ((SizeT)List_Get((List)List_Get(self->_w, k_prime), j))->value;
                            }
                            else
                            {
                                List_Append((List)List_Get(self->_neighbours, k), new_SizeT(kPrimeCurNeighbour));
                                List_Append((List)List_Get(self->_w, k), new_SizeT(((SizeT)List_Get((List)List_Get(self->_w, k_prime), j))->value));
                            }
                        }
                    }
                    size_t kPrimeGroupIndexAsNeighbourOfK;
                    for (size_t j = 0; j < List_GetSize((List)List_Get(self->_neighbours, k)); j++)
                        if (((SizeT)List_Get((List)List_Get(self->_neighbours, k), j))->value == k_prime)
                        {
                            kPrimeGroupIndexAsNeighbourOfK = j;
                            break;
                        }
                    List_Remove((List)List_Get(self->_neighbours, k), kPrimeGroupIndexAsNeighbourOfK);
                    List_Remove((List)List_Get(self->_w, k), kPrimeGroupIndexAsNeighbourOfK);
                    // update such values for neighbours too
                    for (size_t j = 0; j < List_GetSize((List)List_Get(self->_neighbours, k)); j++)
                    {
                        size_t neighbourIndex = ((SizeT)List_Get((List)List_Get(self->_neighbours, k), j))->value;
                        size_t neighbourW = ((SizeT)List_Get((List)List_Get(self->_w, k), j))->value;

                        // iterate over neighbours of neighbour in order to remove k_prime
                        for (size_t l = 0; l < List_GetSize((List)List_Get(self->_neighbours, neighbourIndex)); l++)
                        {
                            if (k == ((SizeT)List_Get((List)List_Get(self->_neighbours, neighbourIndex), l))->value)
                            {
                                ((SizeT)List_Get((List)List_Get(self->_w, neighbourIndex), l))->value = neighbourW;
                            }
                            else if (k_prime == ((SizeT)List_Get((List)List_Get(self->_neighbours, neighbourIndex), l))->value)
                            {
                                List_Remove((List)List_Get(self->_neighbours, neighbourIndex), l);
                                List_Remove((List)List_Get(self->_w, neighbourIndex), l);
                            }
                        }
                    }
                }
            }
            sort(indexesOfGroupsToFuse);
            List map = new_List((delfunc)del_SizeT, (copyfunc)copy_SizeT);
            for (size_t j = 0; j < List_GetSize(self->_G); j++)
                List_Append(map, new_SizeT(IndexAfterRemoval(j, indexesOfGroupsToFuse)));

            for (size_t j = 0; j < List_GetSize(indexesOfGroupsToFuse); j++)
            {
                List_Remove(self->_N, ((SizeT)List_Get(indexesOfGroupsToFuse, List_GetSize(indexesOfGroupsToFuse) - j - 1))->value);
                List_Remove(self->_G, ((SizeT)List_Get(indexesOfGroupsToFuse, List_GetSize(indexesOfGroupsToFuse) - j - 1))->value);
                List_Remove(self->_y_bar, ((SizeT)List_Get(indexesOfGroupsToFuse, List_GetSize(indexesOfGroupsToFuse) - j - 1))->value);
                List_Remove(self->_neighbours, ((SizeT)List_Get(indexesOfGroupsToFuse, List_GetSize(indexesOfGroupsToFuse) - j - 1))->value);
                List_Remove(self->_w, ((SizeT)List_Get(indexesOfGroupsToFuse, List_GetSize(indexesOfGroupsToFuse) - j - 1))->value);
                List_Remove(self->_gamma, ((SizeT)List_Get(indexesOfGroupsToFuse, List_GetSize(indexesOfGroupsToFuse) - j - 1))->value);
            }

            for (size_t j = 0; j < List_GetSize(self->_neighbours); j++)
            {
                for (size_t l = 0; l < List_GetSize((List)List_Get(self->_neighbours, j)); l++)
                {
                    ((SizeT)List_Get((List)List_Get(self->_neighbours, j), l))->value = ((SizeT)List_Get(map, ((SizeT)List_Get((List)List_Get(self->_neighbours, j), l))->value))->value;
                }
            }

            del_List(&map);
            del_List(&indexesOfGroupsToFuse);

            k++;
        }
    }
}

Matrix FLSA_2D_GetCoef(FLSA_2D self)
{
    Matrix result = new_Matrix(self->_n1, self->_n2);

    for (size_t k = 0; k < List_GetSize(self->_G); k++)
    {
        for (size_t i = 0; i < List_GetSize((List)List_Get(self->_G, k)); i++)
        {
            size_t p = ((SizeT)List_Get((List)List_Get(self->_G, k), i))->value;

            size_t col = p % self->_n2;
            size_t row = (p - col) / self->_n2;

            Matrix_Set(result, row, col, ((Double)List_Get(self->_gamma, k))->value);
        }
    }

    return result;
}