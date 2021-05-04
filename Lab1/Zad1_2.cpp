#include <iostream>
#include <cmath>
#define N 5
using namespace std;

void solution(double*, double*, double a[][N]);
void multiple(double a[N][N], double x[N], double c[N]);
double deviation(double*, double*);

int main()
{ 
	FILE* file = fopen("Zad1_2.txt", "w");
	for(double q = 0.2; q <= 5; q += 0.200001)
	{
    double b[N] = {10, 2, 9, 9, 3};
    double bCopy[N] = {10, 2, 9, 9, 3};
    double a[N][N] = { 
      {q * 0.0002, 1, 6, 9, 10},
      {0.0002, 1, 6, 9, 10} ,
      {1, 6, 6, 8, 6},
      {5, 9, 10, 7, 10},
      {3, 4, 9, 7, 9}	
      };
    double aCopy[N][N] = { 
      {q * 0.0002, 1, 6, 9, 10},
      {0.0002, 1, 6, 9, 10} ,
      {1, 6, 6, 8, 6},
      {5, 9, 10, 7, 10},
      {3, 4, 9, 7, 9}
      };
    double c[N] = {0, 0, 0, 0, 0};
		double x[N] = {1, 1, 1 ,1 ,1};
		solution(x, bCopy, aCopy);
    multiple(a, x, c);
    double result = deviation(c,b);
    cout << result << endl;
		fprintf(file, "%2.12f	%2.12f\n", q, result);
	}
  return 0;
}

void multiple(double a[][N], double x[N], double c[N])
{
	 for (int i = 0; i < N; i++) 
	 {
    for (int j = 0; j< N; j++) 
		{
      c[i] += a[i][j] * x[j];
    }
	}
}

void solution(double* x, double* b, double a[][N]) {
  for(int i = 0; i < N; i++) {
    double divider = a[i][i];
    b[i] = b[i] / divider;
    for(int j = 0; j < N; j++) {
      a[i][j] = a[i][j] / divider;
    }
    for(int k = 0; k < N; k++) {
      double factor = a[k][i];
      if (k!=i) {
        b[k] = b[k] - b[i] * factor;
        for (int j = 0; j < N; j++) {
          a[k][j] = a[k][j] - factor * a[i][j];
        }
      }
    }
  }
  for(int i = 0; i < N; i++) {
    x[i] = b[i];
  }
}

double deviation(double* c, double* b)
{
	double deviation = 0;
	for(int i = 0; i < N; i++)
	{
		deviation = deviation + pow((c[i]-b[i]),2);
	}
  return (sqrt(deviation)) / 5.0;
}

