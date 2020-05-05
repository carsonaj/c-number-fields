#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include "number_field.h"

#define TRUE 1
#define FALSE 0

// arrays

// helpful functions
static int min(double a, double b) {
    if (a < b)
        return a;
    else
        return b;
}
//----------------------------------------------------------------------------

// array algorithms

void arr_dubswap(double *arr, int i, int j) {
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

void arr_intswap(int *arr, int i, int j) {
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

double arr_sum(double *arr, int length) {
    double sum = 0;
    int i;
    for (i=0; i<length; i=i+1)
        sum = sum + arr[i];
    return sum;
}

// merge-sort an array

// merge two sorted subarrays arr[l,l+1,...m] arr[m+1,m+2,...,r]
static void merge(double *arr, int l, int m, int r) {
    int i, j, k;
    int n1 = m-l+1;
    int n2 = r-m;

    //create temporary left and right arrays
    double la[n1], ra[n2];

    // copy data to temp arrays
    for (i=0; i<n1; i=i+1)
        la[i] = arr[l+i];

    for (j=0; j<n2; j=j+1)
        ra[j] = arr[m+1+j];

    // merge the temp arrays
    i = 0;
    j = 0;
    k = l;

    while (i<n1 && j<n2) {
        if (la[i] <= ra[j]){
            arr[k] = la[i];
            i = i+1;
        }
        else {
            arr[k] = ra[j];
            j = j+1;
        }
        k = k+1;
    }

    while (i<n1) {
        arr[k] = la[i];
        i = i+1;
        k = k+1;
    }

    while (j<n2) {
        arr[k] = ra[j];
        j = j+1;
        k = k+1;
    }
}

void arr_merge_sort(double *arr, int l, int r) {
    if (l < r) {
        int m = l + (r-l)/2;
        arr_merge_sort(arr, l, m);
        arr_merge_sort(arr, m+1, r);

        merge(arr, l, m, r);
    }
}

//quick-sort an array
typedef struct Tuple_ Tuple;

struct Tuple_ {
    int values[2];
};

// parttion arr with pivot as last element
static Tuple partition(double *arr, int l, int r) {
    int i = l, j = l, k = 0;
    while (j < r-k) {
        if (arr[j] == arr[r-k]) {
            if (arr[j] == arr[r-k-1]) {
                k = k+1;
                continue;
            }
            else {
                arr_dubswap(arr, j, r-k-1);
                k = k+1;
            }
        }
        if (arr[j] < arr[r-k]) {
            if (i != j)
                arr_dubswap(arr, i, j);
            i = i+1;
        }
        j = j+1;
    }

    int n, m = min(j-i-1, k);
    for (n=0; n<=m; n=n+1)
        arr_dubswap(arr, i+n, r-n);

    Tuple x;
    if (k < r-l) {
        x.values[0] = i-1;
        x.values[1] = i+k+1;
    }
    else {
        x.values[0] = l;
        x.values[1] = r;
    }

    return x;
}

void arr_quick_sort(double *arr, int l, int r) {
    if (l < r) {
        Tuple x = partition(arr, l, r);
        arr_quick_sort(arr, l, x.values[0]);
        arr_quick_sort(arr, x.values[1], r);
    }
}

// binary-search an array

//search

int arr_binary_search(double *arr, int l, int r, int x) {
    if (r > l) {
        int m = l + (r-l)/2;

        if (arr[m] == x)
            return m;
        else if (arr[m]<x)
            return arr_binary_search(arr, m+1, r, x);
        else
            return arr_binary_search(arr, l, m, x);
    }
    else {
        if (arr[l] == x)
            return l;
        else
            return -1;
    }
}
//------------------------------------------------------------------
//-----------------------------------------------------------------















//counting

int64_t cnt_factorial(int n, int k) {
    assert(k<=n);
    assert(0<=k);
    if (n==0) {
        return 1;
    }
    else {
        assert(k>=1);
        int64_t prod = 1;
        int i;
        for (i=k; i<=n; i++) {
            prod = prod*i;
        }
        return prod;
    }
}

int64_t cnt_permutation(int n, int k) {
    assert(0<=k);
    if (k==0) {
        return 1;
    }
    if (n<k) {
        return 0;
    }
    else {
        return cnt_factorial(n,n-k+1);
    }
}

int64_t cnt_combination(int n, int k) {
    assert(0<=k);
    if (k==0) {
        return 1;
    }
    if (n<k) {
        return 0;
    }
    else {
        return cnt_permutation(n,k)/cnt_factorial(k,1);
    }
}
//-----------------------------------------------------
//-----------------------------------------------------













//polynomial

#define TRUE 1
#define FALSE 0

// helpful functions
static Polynomial max_deg(Polynomial p1, Polynomial p2) {
    int deg1 = p1.deg;
    int deg2 = p2.deg;
    if (deg1>deg2)
        return p1;
    else if (deg2>deg1)
        return p2;
    else {
        if (ply_is_zero(p1))
            return p2;
        else
            return p1;
    }
}

static Polynomial min_deg(Polynomial p1, Polynomial p2) {
    int deg1 = p1.deg;
    int deg2 = p2.deg;
    if (deg1>deg2)
        return p2;
    else if (deg2>deg1)
        return p1;
    else {
        if (ply_is_zero(p1))
            return p1;
        else
            return p2;
    }
}

// data structure
Polynomial ply_create(int deg) {
    assert(deg >= 0 && deg <= MAX_DEG);
    Polynomial poly;
    poly.deg = deg;

    return poly;
}

int ply_get_deg(Polynomial poly) {
    int deg = poly.deg;
    return deg;
}

double ply_get_coef(Polynomial poly, int i) {
    return poly.coefs[i];
}

void ply_set_coef(Polynomial *poly, int i, double val) {
    poly->coefs[i] = val;
}

Polynomial ply_copy(Polynomial p) {
    Polynomial copy = ply_create(p.deg);
    int i;
    for (i=0; i<=copy.deg; i++) {
        double coef = ply_get_coef(p, i);
        ply_set_coef(&copy, i, coef);
    }

    return copy;
}

void ply_print(Polynomial poly) {
    int deg = poly.deg;
    double coef0 = ply_get_coef(poly, 0);
    double coefn = ply_get_coef(poly, deg); //PROBLEM
    if (deg!=0) {
        assert(coefn!=0);
        int k;
        int i;
        for (i=0; i<=deg; i++) {
            double coef = ply_get_coef(poly, i);
            if (coef != 0) {
                k = i;
                break;
            }
        }
        if (k<deg) {
            double coef = ply_get_coef(poly, k);
            printf("%.2fx^%d ", coef, k);
            for (i=k+1; i<=deg-1; i++) {
                double coef = ply_get_coef(poly, i);
                if (coef>0) {
                    printf("+ %.2fx^%d ", coef, i);
                }
                else if (coef<0)
                    printf("- %.2fx^%d ", -coef, i);
                }
            if (coefn>0)
                printf("+ %.2fx^%d\n", coefn, deg);
            else
                printf("- %.2fx^%d\n", -coefn, deg);
        }

        else if (k==deg) {
            printf("%.3fx^%d\n", coefn, deg);
        }
    }

    else if (deg==0) {
        printf("%.2f\n", coef0);
    }
}
//----------------------------------------------------------------------------

//mathematics

// algebra

int ply_equal(Polynomial p, Polynomial q) {
    int m = p.deg;
    int n = q.deg;

    if (m != n)
        return FALSE;
    else {
        int i;
        for (i=0; i<=m; i++){
            double co_p = ply_get_coef(p, i);
            double co_q = ply_get_coef(q, i);
            if (co_p != co_q)
                return FALSE;
        }

        return TRUE;
    }
}

Polynomial ply_zero() {
    Polynomial z = ply_create(0);
    z.coefs[0] = 0.0;
    return z;
    }

int ply_is_zero(Polynomial p) {
    int deg = p.deg;
    if (ply_get_coef(p, deg)==0)
        return TRUE;
    else
        return FALSE;

}

Polynomial ply_sum(Polynomial poly1, Polynomial poly2) {

    Polynomial p1 = max_deg(poly1, poly2);
    Polynomial p2 = min_deg(poly1, poly2);
    int deg1 = p1.deg;
    int deg2 = p2.deg;

    double sum_coefs[deg1+1];
    int i;
    for (i=0; i<=deg2; i++) {
        sum_coefs[i] = ply_get_coef(poly1, i) + ply_get_coef(poly2, i);
    }
    for (i=deg2+1; i<=deg1; i++) {
        sum_coefs[i] = ply_get_coef(p1, i);
    }
    int deg = deg1;
    while (deg>=0) {
        if (sum_coefs[deg]!=0)
            break;
        if (deg==0)
            break;
        deg = deg-1;
    }

    Polynomial sum_poly = ply_create(deg);
    for (i=0; i<=deg; i++) {
        double coef = sum_coefs[i];
        ply_set_coef(&sum_poly, i, coef);
    }

    return sum_poly;

}

Polynomial ply_product(Polynomial poly1, Polynomial poly2) {
    if (ply_is_zero(poly1)==TRUE)
        return poly1;
    else if (ply_is_zero(poly2)==TRUE)
        return poly2;
    else {
        Polynomial p1 = max_deg(poly1, poly2);
        Polynomial p2 = min_deg(poly1, poly2);
        int deg1 = p1.deg;
        int deg2 = p2.deg;
        int deg = deg1+deg2;

        int k;
        Polynomial prod_poly = ply_create(deg);
        for (k=0; k<=deg; k++) {
            double sum_k = 0;
                int l;
                for (l=0; l<=k; l++) {
                    if ((l<=deg1)&&(k-l<=deg2))
                        sum_k = sum_k+(ply_get_coef(p1, l)*ply_get_coef(p2, k-l));
                }
            ply_set_coef(&prod_poly, k, sum_k);
        }

        return prod_poly;
    }
}

Polynomial ply_scale(double s, Polynomial p) {
    if (s==0) {
        Polynomial z = ply_zero();
        return z;
    }
    else {
        int n = p.deg;
        Polynomial sp = ply_create(n);

        int i;
        for(i=0; i<=n; i++) {
            double coef = ply_get_coef(p, i);
            ply_set_coef(&sp, i, s*coef);
        }

        return sp;
    }
}

Polynomial ply_neg(Polynomial p) {
    Polynomial neg = ply_scale(-1.0, p);
    return neg;
}









//---------------------------------------------------------------
static PolyMatrix sub_division(Polynomial f, Polynomial g) {
    Polynomial q;
    Polynomial r;

    int deg_q = f.deg - g.deg;
    q = ply_monomial(deg_q);

    double coef = ply_get_coef(f, f.deg)/ply_get_coef(g, g.deg);
    ply_set_coef(&q, deg_q, coef);

    Polynomial prod = ply_product(g, q);
    Polynomial scaled = ply_scale(-1.0, prod);
    r = ply_sum(f, scaled);

    PolyMatrix pair = pymat_create(1, 2);
    pymat_set_element(&pair, 0, 0, q);
    pymat_set_element(&pair, 0, 1, r);

    return pair;
}

// division algorithm: divides f by g and
// returns struct containing quotient, remainder

// remember to free pair and the polys it contains after use
PolyMatrix ply_division(Polynomial f, Polynomial g) {
    assert(ply_is_zero(g)==FALSE);
    Polynomial q0 = ply_zero();
    Polynomial r0 = ply_copy(f);

    Polynomial q = q0;
    Polynomial r = r0;
    PolyMatrix pair;

    if (f.deg < g.deg || ply_is_zero(f)==TRUE) {
        pair = pymat_create(1,2);
        pymat_set_element(&pair, 0, 0, q);
        pymat_set_element(&pair, 0, 1, r);

        return pair;
    }
    else {
        while (r.deg >= g.deg && ply_is_zero(r)==FALSE) {
            pair = sub_division(r, g);
            Polynomial temp_q = pymat_get_element(pair, 0, 0);
            Polynomial temp_r = pymat_get_element(pair, 0, 1);

            q = ply_sum(q, temp_q);
            r = ply_copy(temp_r);

        }

        pair = pymat_create(1,2);
        pymat_set_element(&pair, 0, 0, q);
        pymat_set_element(&pair, 0, 1, r);

        return pair;
    }
}

// return gcd and Bezout coefficients for f,g
PolyMatrix ply_gcd(Polynomial f, Polynomial g) {
    assert((ply_is_zero(f)==TRUE && ply_is_zero(f)==TRUE)==FALSE);
    PolyMatrix result = pymat_create(1, 3);

    if (ply_is_zero(f)==TRUE) {
        pymat_set_element(&result, 0, 0, g);
        pymat_set_element(&result, 0, 1, ply_zero());
        pymat_set_element(&result, 0, 2, ply_monomial(0));

        return result;
    }

    else if (ply_is_zero(g)==TRUE) {
        pymat_set_element(&result, 0, 0, f);
        pymat_set_element(&result, 0, 1, ply_monomial(0));
        pymat_set_element(&result, 0, 2, ply_zero());

        return result;
    }

    else {
        int flip = FALSE;
        Polynomial p1;
        Polynomial p2;
        if (f.deg > g.deg) {
            p1 = f;
            p2 = g;
        }
        else {
            p1 = g;
            p2 = f;
            flip = TRUE;
        }

        PolyMatrix a = pymat_create(1, 2);
        PolyMatrix b = pymat_create(1, 2);
        PolyMatrix c = pymat_create(1, 2);

        pymat_set_element(&a, 0, 0, ply_monomial(0));
        pymat_set_element(&a, 0, 1, ply_zero());
        pymat_set_element(&b, 0, 0, ply_zero());
        pymat_set_element(&b, 0, 1, ply_monomial(0));
        pymat_set_element(&c, 0, 0, ply_zero());
        pymat_set_element(&c, 0, 1, ply_monomial(0));

        while(1) {
            PolyMatrix qr = ply_division(p1, p2);
            Polynomial q = pymat_get_element(qr, 0, 0);
            Polynomial r = pymat_get_element(qr, 0, 1);

            if (ply_is_zero(r))
                break;

            PolyMatrix scaled = pymat_poly_scale(q, b);
            PolyMatrix neg_qb = pymat_scale(-1.0, scaled);
            c = pymat_sum(a, neg_qb);

            p1 = ply_copy(p2);
            p2 = ply_copy(r);

            a = b;
            b = c;

        }
        double coef = ply_get_coef(p2, p2.deg);
        if (coef != 1.0) {
            Polynomial t_p2 = p2;
            PolyMatrix t_c = c;

            p2 = ply_scale(1.0/coef, t_p2);
            c = pymat_scale(1.0/coef, t_c);

        }

        Polynomial c00 = pymat_get_element(c, 0, 0);
        Polynomial c01 = pymat_get_element(c, 0, 1);

        pymat_set_element(&result, 0, 0 ,p2);
        if (flip == TRUE) {
            pymat_set_element(&result, 0, 1, c01);
            pymat_set_element(&result, 0, 2, c00);
        }
        else {
            pymat_set_element(&result, 0, 1, c00);
            pymat_set_element(&result, 0, 2, c01);
        }


        return result;
    }
}

//polynomial arithmetic modulo m(x)
//--------------------------------------------------
Polynomial pymod_reduce(Polynomial p, Polynomial m) {
    assert(ply_is_zero(m)==FALSE);
    PolyMatrix pair = ply_division(p, m);
    Polynomial pmodm = pymat_get_element(pair, 0, 1);
    return pmodm;
}

Polynomial pymod_inv(Polynomial p, Polynomial m) {
    PolyMatrix res = ply_gcd(p,m);
    Polynomial pinv_modm = pymat_get_element(res, 0, 1);
    if (pinv_modm.deg >= m.deg) {
        pinv_modm = pymod_reduce(p, m);
    }

    return pinv_modm;
}

Polynomial pymod_sum(Polynomial p, Polynomial q, Polynomial m) {
    Polynomial sum = ply_sum(p, q);
    sum = pymod_reduce(sum, m);

    return sum;
}

Polynomial pymod_product(Polynomial p, Polynomial q, Polynomial m) {
    Polynomial prod = ply_product(p, q);
    prod = pymod_reduce(prod, m);

    return prod;
}






// analysis

// Horner's method
double ply_evaluate(Polynomial poly, double x) {
    int deg = poly.deg;
    double coefn = ply_get_coef(poly, deg);
    assert(coefn!=0);

    double val = coefn;
    if (deg==0) {
        return val;
    }
    else if (deg>0) {
        while (deg>0) {
            val = val*x + ply_get_coef(poly, deg-1);
            deg = deg-1;
        }
    }

    return val;
}

Polynomial ply_differentiate(Polynomial poly, int n) {
    int deg = poly.deg;
    if (deg<n) {
        Polynomial deriv = ply_create(0);
        ply_set_coef(&deriv, 0, 0.0);
        return deriv;
    }
    else {
        int dn_deg = deg-n;
        double dn_coefs[dn_deg+1];
        int i;
        for (i=0; i<=dn_deg; i++) {
            dn_coefs[i] = cnt_factorial(n+i, i+1)*ply_get_coef(poly, n+i);
        }

        while (dn_deg>=0) {
            if (dn_coefs[dn_deg]!=0)
                break;
            if (dn_deg==0)
                break;
            dn_deg = dn_deg-1;
        }

        Polynomial deriv = ply_create(dn_deg);
        for (i=0; i<=dn_deg; i++) {
            double coef = dn_coefs[i];
            ply_set_coef(&deriv, i, coef);
        }

        return deriv;
    }
}

// families of polynomials

// Standard Basis
Polynomial ply_monomial(int n) {
    Polynomial p = ply_create(n);
    int i;
    for (i=0; i<=n-1; i++) {
        ply_set_coef(&p, i, 0.0);
    }
    ply_set_coef(&p, n, 1.0);
    return p;
}

// Legendre (Recursive)
Polynomial ply_legendre(int n) {
    Polynomial p0 = ply_monomial(0);
    Polynomial p1 = ply_monomial(1);

    Polynomial poly_list[n];
    poly_list[0] = p0;
    poly_list[1] = p1;

    int i;
    for (i=1; i<=n-1; i++) {

        Polynomial prod = ply_product(p1, poly_list[i]);
        Polynomial left = ply_scale(2*i+1, prod);


        Polynomial right = ply_scale(-i, poly_list[i-1]);
        Polynomial sum = ply_sum(right, left);
        Polynomial p = ply_scale(1.0/(i+1), sum);

        poly_list[i+1] = p;
    }

    return poly_list[n];
}
//----------------------------------------------------------------------
//--------------------------------------------------------------------














//Matrix
// data structure
Matrix mat_create(int rows, int cols) {
    assert(rows*cols <= MAX_SIZE);
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;

    return mat;
}


Matrix mat_zero(int rows, int cols) {
    Matrix mat = mat_create(rows, cols);
    int i;
    for (i=0; i<rows*cols; i++) {
        mat.data[i] = 0.0;
    }

    return mat;
}

double mat_get_element(Matrix mat, int row, int col) {
    return mat.data[row*(mat.cols) + col];
}

void mat_set_element(Matrix *mat, int row, int col, double element) {
    mat->data[row*(mat->cols) + col] = element;
}

Matrix mat_get_rows(Matrix mat, int rows, int *rows_arr) {
    int cols = mat.cols;
    Matrix row_mat = mat_create(rows, cols);

    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = mat_get_element(mat, rows_arr[i], j);
            mat_set_element(&row_mat, i, j, element);
        }
    }

    return row_mat;
}

Matrix mat_get_cols(Matrix mat, int cols, int *cols_arr) {
    int rows = mat.rows;
    Matrix col_mat = mat_create(rows, cols);

    int i, j;
    for (j=0; j<cols; j++) {
        for (i=0; i<rows; i++) {
            double element = mat_get_element(mat, i, cols_arr[j]);
            mat_set_element(&col_mat, i, j, element);
        }
    }

    return col_mat;
}

Matrix mat_join(Matrix A, Matrix B, int axis) {
    // axis = 0 means vertical join
    // axis = 1 means horizontal join
    int m = A.rows;
    int n = A.cols;
    int r = B.rows;
    int s = B.cols;

    Matrix join;
    if (axis==0) {
        assert(n==s);
        join = mat_create(m+r, n);
        int i, j;
        for (j=0; j<n; j++) {
            for (i=0; i<m; i++) {
                double elementA = mat_get_element(A, i, j);
                mat_set_element(&join, i, j, elementA);
            }

            for (i=0; i<r; i++) {
                double elementB = mat_get_element(B, i, j);
                mat_set_element(&join, i+m, j, elementB);
            }
        }
    }

    else if (axis==1) {
        assert(m==r);
        join = mat_create(m, n+s);
        int i, j;
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                double elementA = mat_get_element(A, i, j);
                mat_set_element(&join, i, j, elementA);
            }

            for (j=0; j<s; j++) {
                double elementB = mat_get_element(B, i, j);
                mat_set_element(&join, i, j+n, elementB);
            }
        }
    }

    return join;
}

void mat_print(Matrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = mat_get_element(mat, i, j);
            printf("%6.2f ", element);
        }
        printf("\n");
    }
}

// create a copy of a matrix in different memory
Matrix mat_copy(Matrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    Matrix cpy = mat_create(rows, cols);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = mat_get_element(mat, i, j);
            mat_set_element(&cpy, i, j, element);
        }
    }

    return cpy;
}
//------------------------------------------------------------------------

// mathematics

// elementary row operations:

// (type 1) swaps rows i,j
void mat_row_op1(Matrix *mat, int i, int j) {
    int cols = mat->cols;
    int rows_arr[1] = {i};
    Matrix temp_i = mat_get_rows(*mat, 1, rows_arr);
    int col;
    for (col=0; col<cols; col=col+1) {
        double element_j = mat_get_element(*mat, j, col);
        mat_set_element(mat, i, col, element_j);
        double element_i = mat_get_element(temp_i, 0, col);
        mat_set_element(mat, j, col, element_i);
    }
}

// (type 2) multiplies row i by a constant k
void mat_row_op2(Matrix *mat, int i, double k) {
    int cols = mat->cols;
    int j;
    for (j=0; j<cols; j=j+1) {
        double element = k*mat_get_element(*mat, i, j);
        mat_set_element(mat, i, j, element);
    }
}

// (type 3) multiplies row j by constant k and adds it to row i
void mat_row_op3(Matrix *mat, int i, int j, double k) {
    int cols = mat->cols;
    int col;
    for (col=0; col<cols; col=col+1) {
        double element = mat_get_element(*mat, i, col) + k*mat_get_element(*mat, j, col);
        mat_set_element(mat, i, col, element);
    }
}

int mat_equal(Matrix A, Matrix B) {
    int same_rows = (A.rows == B.rows);
    int same_cols = (A.cols == B.cols);
    if (!same_rows || !same_cols) {
        return FALSE;
    }

    int rows = A.rows;
    int cols = A.cols;
    int i;
    for (i=0; i<rows*cols; i=i+1) {
        if (A.data[i] != B.data[i])
            return FALSE;
    }

    return TRUE;
}

// standard matrix product
Matrix mat_product(Matrix A, Matrix B) {
    assert(A.cols == B.rows);
    int m,n,p;
    m = A.rows;
    n = A.cols;
    p = B.cols;

    Matrix mat = mat_create(m, p);
    int i,j;
    for (i=0; i<m; i=i+1) {
        for (j=0; j<p; j=j+1) {

            double kron_prod[n];
            int k;
            for (k=0; k<n; k=k+1) {
                kron_prod[k] = mat_get_element(A, i, k)*mat_get_element(B, k, j);
            }

            double element = arr_sum(kron_prod, n);
            mat_set_element(&mat, i, j, element);
        }
    }

    return mat;
}

// Hadamard product
Matrix mat_had_product(Matrix A, Matrix B) {
    assert(A.rows == B.rows);
    assert(A.cols == B.cols);

    int rows = A.rows;
    int cols = A.cols;

    Matrix prod = mat_create(rows, cols);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = mat_get_element(A, i, j)*mat_get_element(B, i, j);
            mat_set_element(&prod, i, j, element);
        }
    }

    return prod;
}

Matrix mat_scale(double c, Matrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    Matrix prod = mat_create(rows, cols);

    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = c*mat_get_element(mat, i, j);
            mat_set_element(&prod, i, j, element);
        }
    }

    return prod;
}

Matrix mat_sum(Matrix A, Matrix B) {
    assert(A.rows == B.rows);
    assert(A.cols == B.cols);
    int m = A.rows;
    int n = A.cols;

    Matrix mat = mat_create(m, n);
    int i,j;
    for (i=0; i<m; i=i+1) {
        for (j=0; j<n; j=j+1) {
            double element = mat_get_element(A, i, j) + mat_get_element(B, i, j);
            mat_set_element(&mat, i, j, element);
        }
    }

    return mat;
}

Matrix mat_transpose(Matrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;

    Matrix trans = mat_create(cols, rows);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            double element = mat_get_element(mat, i, j);
            mat_set_element(&trans, j, i, element);
        }
    }

    return trans;
}

//----------------------------------------------------------------------

// algorithms

// row echelon form
static void sub_ref(Matrix *mat, int start_row, int start_col) {

    // check if mat is zero matrix
    int rows = mat->rows;
    int cols = mat->cols;
    int val = TRUE;
    double element;
    int row = start_row;
    int col = start_col;
    while (col<cols) {
        row = start_row;
        while (row<rows) {
            element = mat_get_element(*mat, row, col);
            if (element != 0) {
                val = FALSE;
                break;
            }
            row = row+1;
        }
        if (!val)
            break;
        col = col+1;
    }
    if (val)
        return;

    // (row, col) is now pivot;
    // row operations to create a pivot of 1 and zeros below pivot
    mat_row_op2(mat, row, 1.0/element);
    mat_row_op1(mat, start_row, row);

    row = row+1;
    while (row<rows) {
        double k = -mat_get_element(*mat, row, col);
        mat_row_op3(mat, row, start_row, k);
        row = row+1;
    }
    // recursive call on submatrix down-right from pivot
    sub_ref(mat, start_row+1, col+1);
}

void mat_ref(Matrix *mat) {
    sub_ref(mat, 0, 0);
}

// reduced row echelon form
static void sub_rref(Matrix *mat, int start_row, int start_col) {
    // assume matrix is in row echelon form
    // check if mat is zero matrix
    int rows = mat->rows;
    int cols = mat->cols;
    int val = TRUE;
    double element;
    int row = start_row;
    int col = start_col;
    while (col<cols) {
        row = start_row;
        while (row<rows) {
            element = mat_get_element(*mat, row, col);
            if (element != 0) {
                val = FALSE;
                break;
            }
            row = row+1;
        }
        if (!val)
            break;
        col = col+1;
    }
    if (val)
        return;

    // now (row, col) is the pivot
    int prev_row = row-1;
    while (0<=prev_row) {
        element = -mat_get_element(*mat, prev_row, col);
        if (element != 0)
            mat_row_op3(mat, prev_row, row, element);
        prev_row = prev_row-1;
    }

    row = row+1;
    col = col+1;
    sub_rref(mat, row, col);
}

void mat_rref(Matrix *mat) {
    mat_ref(mat);
    sub_rref(mat, 0, 0);
}

// solve system Ax=b
Matrix mat_solve_system(Matrix A, Matrix b) {
    int m = A.rows;
    int n = A.cols;
    int l = b.rows;
    assert(m==n);
    assert(m==l);

    Matrix xsol = mat_join(A, b, 1);
    mat_rref(&xsol);
    int cols_arr[1] = {m};
    Matrix x = mat_get_cols(xsol, 1, cols_arr);

    return x;
}

//-------------------------------------------------------------------









//PolyMatrix
PolyMatrix pymat_create(int rows, int cols) {
    assert(rows*cols <= MAX_SIZE);
    PolyMatrix mat;
    mat.rows = rows;
    mat.cols = cols;

    return mat;
}

PolyMatrix pymat_zero(int rows, int cols) {
    PolyMatrix mat = pymat_create(rows, cols);
    int i;
    for (i=0; i<rows*cols; i++) {
        mat.data[i] = ply_zero();
    }

    return mat;
}

Polynomial pymat_get_element(PolyMatrix mat, int row, int col) {
    return mat.data[row*(mat.cols) + col];
}

void pymat_set_element(PolyMatrix *mat, int row, int col, Polynomial element) {
    mat->data[row*(mat->cols) + col] = element;
}

PolyMatrix pymat_get_rows(PolyMatrix mat, int rows, int *rows_arr) {
    int cols = mat.cols;
    PolyMatrix row_mat = pymat_create(rows, cols);

    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial element = pymat_get_element(mat, rows_arr[i], j);
            pymat_set_element(&row_mat, i, j, element);
        }
    }

    return row_mat;
}

PolyMatrix pymat_get_cols(PolyMatrix mat, int cols, int *cols_arr) {
    int rows = mat.rows;
    PolyMatrix col_mat = pymat_create(rows, cols);

    int i, j;
    for (j=0; j<cols; j++) {
        for (i=0; i<rows; i++) {
            Polynomial element = pymat_get_element(mat, i, cols_arr[j]);
            pymat_set_element(&col_mat, i, j, element);
        }
    }

    return col_mat;
}

PolyMatrix pymat_join(PolyMatrix A, PolyMatrix B, int axis) {
    // axis = 0 means vertical join
    // axis = 1 means horizontal join
    int m = A.rows;
    int n = A.cols;
    int r = B.rows;
    int s = B.cols;

    PolyMatrix join;
    if (axis==0) {
        assert(n==s);
        join = pymat_create(m+r, n);
        int i, j;
        for (j=0; j<n; j++) {
            for (i=0; i<m; i++) {
                Polynomial elementA = pymat_get_element(A, i, j);
                pymat_set_element(&join, i, j, elementA);
            }

            for (i=0; i<r; i++) {
                Polynomial elementB = pymat_get_element(B, i, j);
                pymat_set_element(&join, i+m, j, elementB);
            }
        }
    }

    else if (axis==1) {
        assert(m==r);
        join = pymat_create(m, n+s);
        int i, j;
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                Polynomial elementA = pymat_get_element(A, i, j);
                pymat_set_element(&join, i, j, elementA);
            }

            for (j=0; j<s; j++) {
                Polynomial elementB = pymat_get_element(B, i, j);
                pymat_set_element(&join, i, j+n, elementB);
            }
        }
    }

    return join;
}

PolyMatrix pymat_copy(PolyMatrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    PolyMatrix cpy = pymat_create(rows, cols);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial element = pymat_get_element(mat, i, j);
            pymat_set_element(&cpy, i, j, element);
        }
    }

    return cpy;
}

int pymat_equal(PolyMatrix A, PolyMatrix B) {
    int same_rows = (A.rows == B.rows);
    int same_cols = (A.cols == B.cols);
    if (!same_rows || !same_cols) {
        return FALSE;
    }

    int rows = A.rows;
    int cols = A.cols;
    int i;
    for (i=0; i<rows*cols; i=i+1) {
        if (ply_equal(A.data[i], B.data[i])) {
            return FALSE;
        }
    }

    return TRUE;
}

// standard matrix product
PolyMatrix pymat_product(PolyMatrix A, PolyMatrix B) {
    assert(A.cols == B.rows);
    int m,n,p;
    m = A.rows;
    n = A.cols;
    p = B.cols;

    PolyMatrix mat = pymat_create(m, p);
    int i,j;
    for (i=0; i<m; i=i+1) {
        for (j=0; j<p; j=j+1) {

            Polynomial kron_prod[n];
            int k;
            for (k=0; k<n; k=k+1) {
                Polynomial Aik = pymat_get_element(A, i, k);
                Polynomial Bkj = pymat_get_element(B, k, j);
                kron_prod[k] = ply_product(Aik, Bkj);
            }

            Polynomial sum[n-1];
            sum[0] = ply_copy(kron_prod[0]);
            int l;
            for (l=1; l<n-1;l++) {
                sum[i] = ply_sum(kron_prod[i-1], kron_prod[i]);
            }

            Polynomial element = sum[n-2];
            pymat_set_element(&mat, i, j, element);
        }
    }

    return mat;
}

// Hadamard product
PolyMatrix pymat_had_product(PolyMatrix A, PolyMatrix B) {
    assert(A.rows == B.rows);
    assert(A.cols == B.cols);

    int rows = A.rows;
    int cols = A.cols;

    PolyMatrix prod = pymat_create(rows, cols);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial Aij = pymat_get_element(A, i, j);
            Polynomial Bij = pymat_get_element(B, i, j);
            Polynomial element = ply_product(Aij, Bij);
            pymat_set_element(&prod, i, j, element);
        }
    }

    return prod;
}

PolyMatrix pymat_scale(double c, PolyMatrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    PolyMatrix prod = pymat_create(rows, cols);

    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial element = ply_scale(c, pymat_get_element(mat, i, j));
            pymat_set_element(&prod, i, j, element);
        }
    }

    return prod;
}

PolyMatrix pymat_poly_scale(Polynomial p, PolyMatrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    PolyMatrix prod = pymat_create(rows, cols);

    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial element = ply_product(p, pymat_get_element(mat, i, j));
            pymat_set_element(&prod, i, j, element);
        }
    }

    return prod;
}

PolyMatrix pymat_sum(PolyMatrix A, PolyMatrix B) {
    assert(A.rows == B.rows);
    assert(A.cols == B.cols);
    int m = A.rows;
    int n = A.cols;

    PolyMatrix mat = pymat_create(m, n);
    int i,j;
    for (i=0; i<m; i=i+1) {
        for (j=0; j<n; j=j+1) {
            Polynomial Aij = pymat_get_element(A, i, j);
            Polynomial Bij = pymat_get_element(B, i, j);
            Polynomial element = ply_sum(Aij, Bij);
            pymat_set_element(&mat, i, j, element);
        }
    }

    return mat;
}

PolyMatrix pymat_transpose(PolyMatrix mat) {
    int rows = mat.rows;
    int cols = mat.cols;

    PolyMatrix trans = pymat_create(cols, rows);
    int i, j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            Polynomial element = pymat_get_element(mat, i, j);
            pymat_set_element(&trans, j, i, element);
        }
    }

    return trans;
}
//-----------------------------------------------------------
//------------------------------------------------------------












// number fields

// structure

NFNumber nf_create(Polynomial min_poly, Polynomial number) {
    NFNumber x;
    x.min_poly = min_poly;
    if (number.deg > min_poly.deg) {
        x.number = pymod_reduce(number, min_poly);
    }
    else {
        x.number = number;
    }


    return x;
}

void nf_print(NFNumber x) {
    Polynomial poly = x.number;
    int deg = poly.deg;
    double coef0 = ply_get_coef(poly, 0);
    double coefn = ply_get_coef(poly, deg); //PROBLEM
    if (deg!=0) {
        assert(coefn!=0);
        int k;
        int i;
        for (i=0; i<=deg; i++) {
            double coef = ply_get_coef(poly, i);
            if (coef != 0) {
                k = i;
                break;
            }
        }
        if (k<deg) {
            double coef = ply_get_coef(poly, k);
            printf("%.2fa^%d ", coef, k);
            for (i=k+1; i<=deg-1; i++) {
                double coef = ply_get_coef(poly, i);
                if (coef>0) {
                    printf("+ %.2fa^%d ", coef, i);
                }
                else if (coef<0)
                    printf("- %.2fa^%d ", -coef, i);
                }
            if (coefn>0)
                printf("+ %.2fa^%d\n", coefn, deg);
            else
                printf("- %.2fa^%d\n", -coefn, deg);
        }

        else if (k==deg) {
            printf("%.3fa^%d\n", coefn, deg);
        }
    }

    else if (deg==0) {
        printf("%.2f\n", coef0);
    }
}
//-----------------------------------------------

// algebra

NFNumber nf_sum(NFNumber x, NFNumber y) {
    Polynomial z_num = pymod_sum(x.number, y.number, x.min_poly);
    NFNumber z = nf_create(x.min_poly, z_num);

    return z;
}

NFNumber nf_neg(NFNumber x) {
    Polynomial min_poly = x.min_poly;
    Polynomial neg_num = ply_neg(x.number);

    NFNumber neg_x = nf_create(min_poly, neg_num);

    return neg_x;
}

NFNumber nf_product(NFNumber x, NFNumber y) {
    Polynomial z_num = pymod_product(x.number, y.number, x.min_poly);
    NFNumber z = nf_create(x.min_poly, z_num);

    return z;
}

NFNumber nf_inv(NFNumber x) {
    assert(ply_is_zero(x.number)==FALSE);
    Polynomial inv_num = pymod_inv(x.number, x.min_poly);

    NFNumber inv = nf_create(x.min_poly, inv_num);

    return inv;
}
