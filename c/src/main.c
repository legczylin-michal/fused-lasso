#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "matrix\Matrix.h"
#include "flsa_2d\FLSA_2D.h"

#include "basic/String/String.h"

char *prettifyDurationInSeconds(time_t pDurationInSeconds)
{
    time_t seconds = pDurationInSeconds % 60;
    time_t durationInMinutes = (pDurationInSeconds - seconds) / 60;
    time_t minutes = durationInMinutes % 60;
    time_t hours = (durationInMinutes - minutes) / 60;

    int bs = snprintf(NULL, 0, "%dh %dm %ds", hours, minutes, seconds);
    char *result = (char *)malloc(sizeof(char) * (bs + 1));
    snprintf(result, bs + 1, "%dh %dm %ds", hours, minutes, seconds);

    return result;
}

int main(int argc, char const *argv[])
{
    // set seed
    srand(time(NULL));

    bool loadY = false;
    String filepath;
    double lambda1 = 1.0f, lambda2 = 1.0f, tolerance = 0.0001f;
    size_t originalSize = 8, scaleFactor = 4, maximumIterations = 1000;
    if (argc != 1)
    {
        for (size_t i = 1; i < argc; i++)

        {
            if (!strcmp(argv[i], "--load-y"))
            {
                loadY = true;
                filepath = str(argv[i + 1]);
                i++;
            }
            else if (!strcmp(argv[i], "--size"))
            {
                originalSize = atoi(argv[i + 1]);
                i++;
            }
            else if (!strcmp(argv[i], "--scale"))
            {
                scaleFactor = atoi(argv[i + 1]);
                i++;
            }
            else if (!strcmp(argv[i], "--lambda1"))
            {
                lambda1 = atof(argv[i + 1]);
                i++;
            }
            else if (!strcmp(argv[i], "--lambda2"))
            {
                lambda2 = atof(argv[i + 1]);
                i++;
            }
            else if (!strcmp(argv[i], "--max-iter"))
            {
                maximumIterations = atoi(argv[i + 1]);
                i++;
            }
            else if (!strcmp(argv[i], "--tol"))
            {
                tolerance = atof(argv[i + 1]);
                i++;
            }
            else if (!strcmp(argv[i], "--help"))
            {
                printf("Usage: main.exe [--option value]\n");
                printf("Options:\n");
                printf("%-10s %-10s load specified data matrix. Otherwise randomly generate one (default).\n", "--load-y", "value");
                printf("%-10s %-10s random generation only. Specify original data matrix dimensions. (default 8)\n", "--size", "value");
                printf("%-10s %-10s random generation only. Speocify scaling factor of data matrix (default 4).\n", "--scale", "value");
                printf("%-10s %-10s penalty parameter for coefficients (default 1).\n", "--lambda1", "value");
                printf("%-10s %-10s penalty parameter for differences in coefficients (default 1).\n", "--lambda2", "value");
                printf("%-10s %-10s maximum number of iterations to look for convergence (default 1000).\n", "--max-iter", "value");
                printf("%-10s %-10s numbers different by less then this value are considered equal (default 0.0001).\n", "--tol", "value");
                printf("%-10s %-10s calls this menu.\n", "--help", "");

                return 0;
            }
        }
    }

    Matrix y;

    if (loadY)
    {
        y = Matrix_ReadCsv(String_GetCharacters(filepath), ",");
    }
    else
    {
        // create matrix of zeros with ones placed in random spots. This decides where fields will be
        Matrix random_binomial = RandomBinomial(1, 0.3, originalSize, originalSize);
        // create matrix with random values for every field
        Matrix random_uniform = RandomUniform(0, 1, originalSize, originalSize);
        // apply, so now every field has its random value, but 0 everywhere else
        Matrix y_original_small = Matrix_Multiply(random_binomial, random_uniform);

        // repeat by now created matrix along 0-th axis (repeat rows)
        Matrix y_original_small_repeated_axis_0 = Matrix_Repeat(y_original_small, scaleFactor, 0);
        // repeat by now created matrix along 1-st axis (repeat columns)
        Matrix y_original = Matrix_Repeat(y_original_small_repeated_axis_0, scaleFactor, 1);

        // create matrix of noise
        Matrix random_normal = RandomNormal(0, 0.1, Matrix_GetNumRows(y_original), Matrix_GetNumCols(y_original));
        // apply noise to original data matrix
        y = Matrix_Add(y_original, random_normal);

        // save original matrix as csv for future visualisation
        Matrix_SaveCsv(y_original, "data/y_original.csv", ",");
        // save noised matrix as csv for future visualisation
        Matrix_SaveCsv(y, "data/y.csv", ",");

        del_Matrix(&random_normal);
        del_Matrix(&y_original);
        del_Matrix(&y_original_small_repeated_axis_0);
        del_Matrix(&y_original_small);
        del_Matrix(&random_uniform);
        del_Matrix(&random_binomial);
    }

    // create FLSA 2D model and specify parameters
    FLSA_2D model = new_FLSA_2D(lambda1, lambda2, maximumIterations, tolerance);
    // fit model
    time_t start = time(NULL);
    FLSA_2D_Fit(model, y);
    time_t end = time(NULL);
    char *str = prettifyDurationInSeconds(end - start);
    printf("evaluation took %s\n", str);
    free(str);

    // save estimated parameters by model as csv fo future visualisation
    Matrix_SaveCsv(FLSA_2D_GetCoef(model), "data/gamma_c.csv", ",");

    // deallocate memory
    del_FLSA_2D(&model);
    del_Matrix(&y);

    // assess that program completed successfully
    printf("\nSUCCESSFULLY COMPILED\n");

    // return success compilation code
    return 0;
}