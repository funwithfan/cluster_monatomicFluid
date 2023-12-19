#include <stdio.h>
#include <stdbool.h>

void histogram(int *array, int size, int *count, int offset);
void initializeIntArray(int size, int *array, int value);
void initializeDoubleArray(int size, double *array, double value);
void initializeBoolArray(int size, bool *array, bool value);
void initializeIntMatrix(int rows, int cols, int **array, int value);
void initializeDouble2DMatrix(int rows, int cols, double **array, double value);
int findMaxValueIndex(int *array, int size);
double meanDoubleArray(double *array, int length);
double maxDouble2DMatrix(int rows, int cols, double **array);
double minDouble2DMatrix(int rows, int cols, double **array);
void sumRowsInt2DMatrix(int rows, int cols, int **array, int *sumRows);
void sumRowsIntSymmetricMatrix(int *array, int size, int *sumRows);
double sumAllDouble2DMatrix(int rows, int cols, double **array);
void doubleArraysDotProduct(int size, double *array1, double *array2, double *product);
void intArraysMutiplyConstant(int size, int *array, double constant, double *product);
void doubleArraysDotDivision(int size, double *array1, double *array2, double *product);
void saveDouble2DMatrix(char *filename, int rows, int cols, double **array);
void saveInt2DMatrix(char *filename, int rows, int cols, int **array);
void saveInt1DArray(char *filename, int length, int *array);

int SymmetricIndex(int size, int i, int j);
void setDoubleSymmetricMatrixValue(double *array, int size, int i, int j, double value);
double getDoubleSymmetricMatrixValue(double *array, int size, int i, int j);
void setIntSymmetricMatrixValue(int *array, int size, int i, int j, int value);
int getIntSymmetricMatrixValue(int *array, int size, int i, int j);