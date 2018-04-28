#include "biryukovaes.h"

/**
 * Введение в дисциплину
 */
void biryukovaes::lab1()
{
std::cout<<"Hello, world!";
	
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void biryukovaes::lab2()
{
int i, j, m, l;
    for(i=0;i<N;i++)
        {
          m=i;
          for(l=i+1; l<N; l++)
          if (abs(A[l][i])>abs(A[m][i]))
                m=l;
          swap(A[m],A[i]);
          swap(b[m],b[i]);
           b[i]/=A[i][i];
          for(j=N-1; j>i; A[i][j--]/=A[i][i]);
          A[i][i]=1;
          for(int j=i+1; j<N;j++){
          for(l=N-1;l>i;l--)
          A[j][l]-=A[i][l]*A[j][i];
          b[j]-=b[i]*A[j][i];
          A[j][i]=0;
          }
            cout<<endl;
          }
          x[N-1]=b[N-1];
          for ( int i = N - 2; i >= 0; i-- )
         {
           x[i] = b[i];
           for ( int j = i + 1; j < N; j++ ) {
              x[i] -= A[i][j] * x[j];
              }
          }
}



/**
 * Метод прогонки
 */
void biryukovaes::lab3()
{

    double* new_A = new double[N];
    double* new_b = new double[N];

    new_A[0] = A[0][1] / (-A[0][0]);
    new_b[0] = b[0] / A[0][0];

    for(int i = 1; i < N; i++)
    {
        new_A[i] = A[i][i+1] / (-A[i][i-1] * new_A[i-1] - A[i][i]);
        new_b[i] = (-b[i] + A[i][i-1] * new_b[i-1]) / ( -A[i][i-1] * new_A[i-1] - A[i][i]);
    }

    for(int i = N - 1; i >= 0; i--)
    {
        x[i] = new_A[i] * x[i+1] + new_b[i];
    }

    delete[] new_A;
    delete[] new_b;

}

/**
 * Метод простых итераций
 */
void biryukovaes::lab4()
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
void biryukovaes::lab5()
{ //Якоби
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
void biryukovaes::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void biryukovaes::lab7()
{
	double Eps = 1e-20;
	double Del, a, aAbs;
	double *K = new double[N];
	double *L = new double[N];
	double *M = new double[N];
	double *xrez = new double[N]; 
	for (int i = 0; i<N; i++){
		xrez[i] = 0;
	}
	do {
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * xrez[j];
		}
		for (int i = 0; i < N; i++) {
			L[i] = K[i] - b[i];
		}
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * L[j];
		}
		for (int i = 0; i < N; i++) {
			M[i] = 0;
			for (int j = 0; j < N; j++) {
				M[i] += A[i][j] * K[j];
			}
		}
		a = 0;
		aAbs = 0;
		for (int i = 0; i < N; i++) {
			a += K[i] * L[i];
			aAbs += M[i] * K[i];
		}
		if (a == aAbs)
			a = 1;
		else 
			a = a / aAbs;
		for (int i = 0; i < N; i++)
			x[i] = xrez[i] - a*L[i];
		Del = abs(x[0] - xrez[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - xrez[i])>Del)
				Del = abs(x[i] - xrez[i]);
				xrez[i] = x[i];
		}
	} while (Eps < Del);
}


void biryukovaes::lab8()
{

}

void biryukovaes::lab9()
{

}


/**
 * Решение нелинейных уравнений
 */
void biryukovaes::lab10()
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

std::string biryukovaes::get_name()
{
  return "Biryukova E.S.";
}
