#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix\Matrix.h"
#include "flsa_2d\FLSA_2D.h"
#include "list\List.h"
#include "VanillaTypesWrappers\SizeT.h"
#include "VanillaTypesWrappers\Double.h"

int main(int argc, char const *argv[])
{
    // set seed
    srand(666);

    // set original dimension of picture as well as scale ratio
    size_t o_dim = 8, scale = 4;

    // create matrix of zeros with ones placed in random spots. This decides where fields will be
    Matrix random_binomial = RandomBinomial(1, 0.3, o_dim, o_dim);
    // create matrix with random values for every field
    Matrix random_uniform = RandomUniform(0, 1, o_dim, o_dim);
    // apply, so now every field has its random value, but 0 everywhere else
    Matrix y_original_small = Matrix_Multiply(random_binomial, random_uniform);

    // repeat by now created matrix along 0-th axis (repeat rows)
    Matrix y_original_small_repeated_axis_0 = Matrix_Repeat(y_original_small, scale, 0);
    // repeat by now created matrix along 1-st axis (repeat columns)
    Matrix y_original = Matrix_Repeat(y_original_small_repeated_axis_0, scale, 1);

    // create matrix of noise
    Matrix random_normal = RandomNormal(0, 0.1, Matrix_GetNumRows(y_original), Matrix_GetNumCols(y_original));
    // apply noise to original data matrix
    Matrix y = Matrix_Add(y_original, random_normal);

    // create FLSA 2D model and specify parameters
    FLSA_2D model = new_FLSA_2D(1.0f, 1.0f, 1000, 0.0001f);
    // fit model
    FLSA_2D_Fit(model, y);

    // save original matrix as csv for future visualisation
    Matrix_SaveCsv(y_original, "data/y_original.csv", ",");
    // save noised matrix as csv for future visualisation
    Matrix_SaveCsv(y, "data/y.csv", ",");
    // save estimated parameters by model as csv fo future visualisation
    Matrix_SaveCsv(FLSA_2D_GetCoef(model), "data/gamma.csv", ",");

    // deallocate memory
    del_FLSA_2D(&model);
    del_Matrix(&y);
    del_Matrix(&random_normal);
    del_Matrix(&y_original);
    del_Matrix(&y_original_small_repeated_axis_0);
    del_Matrix(&y_original_small);
    del_Matrix(&random_uniform);
    del_Matrix(&random_binomial);

    // assess that program completed successfully
    printf("\nSUCCESSFULLY COMPILED\n");

    // return success compilation code
    return 0;
}