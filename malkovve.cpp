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
double *R = new double [N]; 
double *Delta = new double [N]; 
double *TempX = new double[N]; 
double maxi=0.0, Tau=0.0, TempTau=0.0; 
double eps = 0.0000001; 
for (int i=0; i<N; i++) 
TempX[i]=0; 
do 
{ 
MatrVekt(N, A, TempX, R); 
for(int i=0; i<N; i++) 
{ 
Delta[i]=R[i]-b[i]; 
} 
MatrVekt(N, A, Delta, R); 
Tau=0.0; 
TempTau=0.0; 
for(int i=0; i<N; i++) 
{ 
Tau+=R[i]*Delta[i]; 
TempTau+=R[i]*R[i]; 
} 
Tau=Tau/TempTau; 
for(int i=0; i<N; i++) 
x[i]=TempX[i]-Tau*Delta[i]; 
maxi = fabs(x[0] - TempX[0]); 
for(int i=0; i<N; i++) 
{ 
if(fabs(x[i]-TempX[i])>maxi) 
maxi=fabs(x[i]-TempX[i]); 
TempX[i]=x[i]; 
} 
} 
while (maxi>=eps); 

delete[] R; 
delete[] Delta; 
delete[] TempX; 

}



/**
 * Метод сопряженных градиентов
 */
void malkovve::lab7()
{
double eps = 1e-15;
	double delta, r, rModul;
	double *w = new double[N];
	double *wp = new double[N];
	double *v = new double[N];
	double *result = new double[N];

	for (int i = 0; i<N; i++)
		result[i] = 0;

	do {
		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++)
				w[i] += A[i][j] * result[j];
		}

		for (int i = 0; i < N; i++) {
			v[i] = w[i] - b[i];
		}

		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++)
				w[i] += A[i][j] * v[j];
		}
		
		for (int i = 0; i < N; i++) {
			wp[i] = 0;
			for (int j = 0; j < N; j++) {
				wp[i] += A[i][j] * w[j];
			}
		}
		
		r = 0;
		rModul = 0;
		
		for (int i = 0; i < N; i++) {
			r += w[i] * v[i];
			rModul += wp[i] * w[i];
		}
		
		if (r == rModul)r = 1;
		else r = r / rModul;
		for (int i = 0; i < N; i++)
			x[i] = result[i] - r*v[i];
		delta = abs(x[0] - result[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - result[i])>delta)
				delta = abs(x[i] - result[i]);
			result[i] = x[i];
		}
	} 
	while (eps < delta);
}


void malkovve::lab8()
{
double * Y = new double[N];//предыдущее приближение
	double * y = new double[N];//последующее приближение
	double maxSob,sob,sum;
	double eps = 1e-9;
	for (int i = 0; i < N; i++)
		Y[i] = 0;
	Y[0] = 1;
	do{
		sum = 0;
		for (int i = 0; i < N; i++)
			sum += Y[i] * Y[i];
		sob = sqrt(sum);
		for (int i = 0; i < N; i++)
		{
			y[i] = 0;
			for (int j = 0; j < N; j++)
				y[i] += A[i][j] * Y[j] / sob;
		}
		sum = 0;
		for (int i = 0; i < N; i++)
			sum += y[i] * y[i];
		maxSob = sqrt(sum);
		for (int i = 0; i<N; i++)
			Y[i] = y[i];
	} 
	while (abs(maxSob - sob)>eps);
	cout << maxSob << endl;
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
