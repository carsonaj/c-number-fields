#include <inttypes.h>

#ifndef NUMBER_FIELD_H
#define NUMBER_FIELD_H

#define MAX_DEG 100
#define MAX_SIZE 50

typedef struct Polynomial Polynomial;
typedef struct Matrix Matrix;
typedef struct PolyMatrix PolyMatrix;
typedef struct NFNumber NFNumber;

struct Polynomial {
    int deg;
    double coefs[MAX_DEG+1];
};

struct Matrix {
    int rows;
    int cols;
    double data[MAX_SIZE];
};

struct PolyMatrix {
    int rows;
    int cols;
    Polynomial data[MAX_SIZE];
};

struct NFNumber {
    Polynomial min_poly;
    Polynomial number;
};

// array
int arr_min(double a, double b);

void arr_dubswap(double *arr, int i, int j);
void arr_intswap(int *arr, int i, int j);
double arr_sum(double *arr, int length);

// merge-sort an array
void arr_merge_sort(double *arr, int l, int r);

//quick-sort an array
void arr_quick_sort(double *arr, int l, int r);

// binary-search an array
int arr_binary_search(double *arr, int l, int r, int x);
//-------------------------------------------------
//-------------------------------------------------

// counting

int64_t cnt_factorial(int n, int k);
int64_t cnt_permutation(int n, int k);
int64_t cnt_combination(int n, int k);


// polynomial




// structure
Polynomial ply_create(int deg);
int ply_get_deg(Polynomial poly);
double ply_get_coef(Polynomial poly, int i);
void ply_set_coef(Polynomial *poly, int i, double val);
Polynomial ply_copy(Polynomial p);
void ply_print(Polynomial poly);

// mathematics

// algebra
int ply_equal(Polynomial p, Polynomial q);
Polynomial ply_zero();
int ply_is_zero(Polynomial p);
Polynomial ply_sum(Polynomial poly1, Polynomial poly2);
Polynomial ply_product(Polynomial poly1, Polynomial poly2);
Polynomial ply_scale(double s, Polynomial p);
Polynomial ply_neg(Polynomial p);
PolyMatrix ply_division(Polynomial f, Polynomial g);
PolyMatrix ply_gcd(Polynomial p, Polynomial q);
Polynomial pymod_reduce(Polynomial p, Polynomial m);
Polynomial pymod_inv(Polynomial p, Polynomial m);
Polynomial pymod_sum(Polynomial p, Polynomial q, Polynomial m);
Polynomial pymod_product(Polynomial p, Polynomial q, Polynomial m);

// analysis
double ply_evaluate(Polynomial poly, double x);
Polynomial ply_differentiate(Polynomial poly, int n);

// families of polynomials
Polynomial ply_monomial(int n);
Polynomial ply_legendre(int n);
//----------------------------------------------------
//-----------------------------------------------------

//matrix

// Matrix
// structure
Matrix mat_create(int rows, int cols);
Matrix mat_zero(int rows, int cols);
double mat_get_element(Matrix mat, int row, int col);
void mat_set_element(Matrix *mat, int row, int col, double element);
Matrix mat_get_rows(Matrix mat, int rows, int *rows_arr);
Matrix mat_get_cols(Matrix mat, int cols, int *cols_arr);
Matrix mat_join(Matrix A, Matrix B, int axis);
void mat_print(Matrix mat);
Matrix mat_copy(Matrix mat);

// mathematics
void mat_row_op1(Matrix *mat, int i, int j);
void mat_row_op2(Matrix *mat, int i, double k);
void mat_row_op3(Matrix *mat, int i, int j, double k);
int mat_equal(Matrix A, Matrix B);
Matrix mat_product(Matrix A, Matrix B);
Matrix mat_had_product(Matrix A, Matrix B);
Matrix mat_scale(double c, Matrix mat);
Matrix mat_sum(Matrix A, Matrix B);
Matrix mat_transpose(Matrix mat);

// algorithms
void mat_ref(Matrix *mat);
void mat_rref(Matrix *mat);
Matrix mat_solve_system(Matrix A, Matrix b);


//------------------------------------------------------------

// PolyMatrix
//structure
PolyMatrix pymat_create(int rows, int cols);
PolyMatrix pymat_zero(int rows, int cols);
Polynomial pymat_get_element(PolyMatrix mat, int row, int col);
void pymat_set_element(PolyMatrix *mat, int row, int col, Polynomial element);
PolyMatrix pymat_get_rows(PolyMatrix mat, int rows, int *rows_arr);
PolyMatrix pymat_get_cols(PolyMatrix mat, int cols, int *cols_arr);
PolyMatrix pymat_join(PolyMatrix A, PolyMatrix B, int axis);
//void pymat_print(PolyMatrix mat);
PolyMatrix pymat_copy(PolyMatrix mat);

// mathematics
int pymat_equal(PolyMatrix A, PolyMatrix B);
PolyMatrix pymat_product(PolyMatrix A, PolyMatrix B);
PolyMatrix pymat_had_product(PolyMatrix A, PolyMatrix B);
PolyMatrix pymat_scale(double c, PolyMatrix mat);
PolyMatrix pymat_poly_scale(Polynomial p, PolyMatrix mat);
PolyMatrix pymat_sum(PolyMatrix A, PolyMatrix B);
PolyMatrix pymat_transpose(PolyMatrix mat);










//-----------------------------------------------------
//-----------------------------------------------------

// number field

// structure
NFNumber nf_create(Polynomial min_poly, Polynomial number);
void nf_print(NFNumber x);

// algebra
NFNumber nf_sum(NFNumber x, NFNumber y);
NFNumber nf_neg(NFNumber x);
NFNumber nf_product(NFNumber x, NFNumber y);
NFNumber nf_inv(NFNumber x);



#endif
