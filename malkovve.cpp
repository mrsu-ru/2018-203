#include "malkovve.h"

/**
 * Введение в дисциплину
 */
void malkovve::lab1()
{
std::cout<<"hello world";
}



/**
 * Метод Гаусса с выбором главного элемента
 */
void malkovve::lab2()
{
	double Q = 0;

    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            Q = A[j][i] / A[i][i];

            for (int k = i; k < N; k++)
            {
                A[j][k] -= Q * A[i][k];
            }
            b[j] -= Q * b[i];
        }
    }

    for(int i = 0; i < N; i++)
	{
        x[i] = b[i];
	}

    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
		{
			x[i] -= A[i][j] * x[j];
		}

        x[i] /= A[i][i];
	}
}



/**
 * Метод прогонки
 */
void malkovve::lab3()
{
double Q = 0;
	int max;

    for (int i = 0; i < N - 1; i++)
    {
				max = i;

		for (int j = i + 1; j < N; j++)
		{
			if(abs(A[j][i]) > abs(A[max][i]))
			{
				max = j;
			}
		}


				std::swap(A[max], A[i]);
		    std::swap(b[max], b[i]);

        for (int j = i + 1; j < N; j++)
        {
            Q = A[j][i] / A[i][i];
            for (int k = i; k < N; k++)
            {
                A[j][k] -= Q * A[i][k];
            }
            b[j] -= Q * b[i];
        }
    }

    for(int i = 0; i < N; i++)
	{
        x[i] = b[i];
	}

    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
		{
			x[i] -= A[i][j] * x[j];
		}

        x[i] /= A[i][i];
	}
}



/**
 * Метод простых итераций
 */
void malkovve::lab4()
{
    double eps = 1e-13;
    double tau = 1e-5;
    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }
    double x1;
    double *xn = new double[N];
    int step = 0;
    do {
        step++;
        for (int i = 0; i < N; i++) {
            xn[i] = x[i];
            for (int k = 0; k < N; k++)
                xn[i] -= tau*A[i][k] * x[k];
            xn[i] += tau * b[i];
        }
        x1 = 0.;
        for (int i = 0; i < N; i++) {
            x1 += (xn[i]-x[i])*(xn[i]-x[i]);
        }
        for (int i = 0; i < N; i++) {
            x[i] = xn[i];
        }
    } while (sqrt(x1)>eps);
}



/**
 * Метод Якоби или Зейделя
 */
void malkovve::lab5()
{
    //Якоби
    const double eps = 10E-20;
    double* y = new double[N];
    double r = 0;
    for(int i=0; i<N; i++){
        x[i] = 0; }
    do {
        for(int i=0; i<N; i++){
            y[i] = b[i];
            for(int j=0; j<N; j++){
                if(i != j){
                    y[i] -= A[i][j]*x[j];
                }
            }
            y[i] /= A[i][i];
        }
        r = abs(x[0] - y[0]);
        for(int i=0; i<N; i++){
            if(abs(x[i]-y[i]) > r)
                r = sqrt((x[i]-y[i])*(x[i]-y[i]));
            x[i] = y[i];
        }
    } while(r >= eps);
    delete[] y;
}




/**
 * Метод минимальных невязок
 */
void malkovve::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void malkovve::lab7()
{

}


void malkovve::lab8()
{

}


void malkovve::lab9()
{

}

/** 
* Решение нелинейных уравнений 
*/ 

void malkovve::lab10() 
{ 
double y; 
for(int i=0; i < 100; i++) { 
double xd; 
double eps = 1e-5; 
y = i; 
int ind = 0; 
do { 
ind++; 
xd = y; 
y = exp((-y)); 
} while (abs(xd - y) > eps || ind > 1000); 
if (y == y) break; 
} 
cout << y << endl; 
}


std::string malkovve::get_name()
{
  return "Malkov V.E.";
}
