#include <iostream>
#include <iomanip>
#include <cmath>
#define N 6

using namespace std;

const double eps = 1e-12;
double detLU(double LU[N][N]);

double* calculateY(double* x, double* c)
{
    double* y = new double[N];
    for (int i = 0; i < N; i++)
    {
        y[i] = 0;
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            y[i] += pow(x[i], j) * c[j];
        }
    }
    return y;
}

double** calculateA(double* x)
{
    double** a = new double*[N];
    for (int i = 0; i < N; i++) {
        a[i] = new double[N];
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            a[i][j] = pow(x[i], j);
        }
    }
    return a;
}

/**
 * \brief 
 * \param matrix  
 * \return det(matrix)
 */
double luDecomposition(double** matrix)
{
    double lower[N][N], upper[N][N];
    memset(lower, 0, sizeof(lower));
    memset(upper, 0, sizeof(upper));
	
    for (int i = 0; i < N; i++)
    {
        for (int k = i; k < N; k++)
        {
            int sum = 0;
            for (int j = 0; j < i; j++)
            {
                sum += (lower[i][j] * upper[j][k]);
            }
            upper[i][k] = matrix[i][k] - sum;
        }

        for (int k = i; k < N; k++)
        {
            if (i == k)
            {
                lower[i][i] = 1;
            }
            else
            {
                int sum = 0;
                for (int j = 0; j < i; j++)
                {
                    sum += (lower[k][j] * upper[j][i]);
                }
                lower[k][i] = (matrix[k][i] - sum) / upper[i][i];
            }
        }
    }
	cout << setw(6)<< "      Lower Triangular"<< setw(32)<< "Upper Triangular" << endl;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << setw(6) << lower[i][j] << "\t";
        }
        cout << "\t";
    	
        for (int j = 0; j < N; j++)
        {
            cout << setw(6) << upper[i][j] << "\t";
        }
        cout << endl;
    }

    double determinant = detLU(lower) * detLU(upper);
    cout << "Determinant = " << determinant << endl;

    return determinant;
}

bool luDist(double** A)
{

    for (int k = 0; k < N - 1; k++)
    {
        if (fabs(A[k][k]) < eps) return false;

        for (int i = k + 1; i < N; i++)
            A[i][k] /= A[k][k];

        for (int i = k + 1; i < N; i++)
            for (int j = k + 1; j < N; j++)
                A[i][j] -= A[i][k] * A[k][j];
    }

    return true;
}

bool luSolve(int k, double** A, double** X)
{
    double s;

    for (int i = 1; i < N; i++)
    {
        s = 0;

        for (int j = 0; j < i; j++) s += A[i][j] * X[j][k];
        {
            X[i][k] -= s;
        }
    }

    if (fabs(A[N - 1][N - 1]) < eps)
    {
        return false;
    }

    X[N - 1][k] /= A[N - 1][N - 1];

    for (int i = N - 2; i >= 0; i--)
    {
        s = 0;

        for (int j = i + 1; j < N; j++)
        {
            s += A[i][j] * X[j][k];
        }

        if (fabs(A[i][i]) < eps)
        {
            return false;
        }

        X[i][k] = (X[i][k] - s) / A[i][i];
    }
    return true;
}

double detLU(double LU[N][N])
{
    double a = 1;
    for (int i = 0, j = 0; i < N && j < N; i++, j++)
    {
        a *= LU[i][j];
    }
    return a;
}

double** calculateInverse(double** A)
{

    bool ok;
    double** X = new double* [N];

    for (int i = 0; i < N; i++)
	{
		 X[i] = new double[N];
	}
	
    if (luDist(A))
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                X[i][j] = 0;
            }
            X[i][i] = 1;
        }
       
        ok = true;
        for (int i = 0; i < N; i++)
            if (!luSolve(i, A, X))
            {
                ok = false;
                break;
            }
    }
    else ok = false;
    
    cout << "Inverse matrix: " << endl;
    if (ok)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                cout << setw(10) << X[i][j] << " ";
            }
            cout << endl;
        }
    }
    else cout << "DZIELNIK ZERO\N";

    return X;
}

double calculateConditionNumber(double** A, double** inversedA)
{
    double maxA = A[0][0];
    double maxIinversed = inversedA[0][0];

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (A[i][j] > maxA)
			{
                maxA = A[i][j];
			}
			if (inversedA[i][j] > maxIinversed)
			{
                maxIinversed = inversedA[i][j];
			}
		}
	}

    return maxA * maxIinversed;
}

int  main()
{
    double** A, ** X;
    bool ok;
    double x[N] = { 1, 2, 3, 4, 5, 6 };
    double c[N] = { 1, 5, -3, 7, 2, 5 };
    cout << setprecision(5) << fixed;

    A = new double* [N];
    X = new double* [N];
    A = calculateA(x);
    double* y = new double[N];
    y = calculateY(x, c);
    cout << "A:" << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << "Y:" << endl;
    for (int i = 0; i < N; i++)
    {
        cout << y[i] << " ";
        
    }
    cout << endl;

    double detA = luDecomposition(A);
    double ** inverseA = calculateInverse(A);
    cout << "Calculating LU for inversed matrix A..." << endl;
    double detInversedA = luDecomposition(inverseA);
    double conditionNumber = calculateConditionNumber(A, inverseA);
    cout << "Condition number = " << conditionNumber << endl;
	
    return 0;
}
