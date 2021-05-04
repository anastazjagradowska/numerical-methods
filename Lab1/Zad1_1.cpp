#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const int N = 200;
void set_a(double a[][N], double, double);
void solve_x_Gauss(double*, double *, double a[][N]);
void solve_x_GJ(double*, double*, double a[][N]);
void print_matrix(double a[][N]);
int main()
{
    int i;
    double A = 1;
    double v_0 = 0;
    double omega = 1;
    double h = 0.1;
    double x[N] = {0};
    double b[N] = {0}; 
    b[0] = A; 
    b[1] = v_0 * h;
    double a[N][N] = {0};
    set_a(a, omega, h);
    
    solve_x_GJ(x, b, a);
    FILE* fp = fopen("lab1_x.txt", "w");
    for (i = 0; i < N; i++) {
        fprintf(fp, "%f %f\n", h * i, x[i]);
    }
    
    return 0;
}
void set_a(double a[][N], double omega, double h) {
    a[0][0] = 1;
    a[1][0] = -1;
    for(int i = 0; i < N; i++) {
        a[i][i] = 1;
        if (i >= 2) {
          a[i][i - 2] = 1;
          a[i][i - 1] = omega * omega * h * h - 2;
        }
    }
    for(int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        printf("%f ", a[i][j]);
      }
      printf("\n");
    }
    printf("\n");
    
}
void solve_x(double *x, double *b, double a[][N]) {
    for(int i = 0; i < N; i++) {
      x[i] = b[i] / a[i][i];
      for (int j = N - i; j < N; j++) {
        x[i] = x[i] - (a[i][N - j - 1] * x[N - j - 1]) / a[i][i];
      }
    }
    for (int i = 0; i < N; i++) {
      printf("%f\n", x[i]);
    }
}
void solve_x_GJ(double* x, double* b, double a[][N]) {
    double w, w2;
    for(int i = 0; i < N; i++) {
      w = a[i][i];
      for (int j = 0; j < N; j++) {
        a[i][j] = a[i][j] / w;
        b[i] = b[i] / w;
        }
      for(int k = 0; k < N; k++) {
        w2 = a[k][i];
        if (k!=i) {
          b[k] = b[k] - b[i] * w2;
          for (int j = 0; j < N; j++) {
            a[k][j] = a[k][j] - w2 * a[i][j];
          }
        }
      }
        //print_matrix(a);
    }
    for(int i = 0; i < N; i++) {
      x[i] = b[i];
    }
}
void print_matrix(double a[][N]) {
  for(int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%f ", a[i][j]);
    }
    printf("\n");
    }
    printf("\n");
}