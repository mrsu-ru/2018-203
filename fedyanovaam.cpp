#include "fedyanovaam.h"

/**
 * Введение в дисциплину
 */
void fedyanovaam::lab1()
{
cout<<"Hello!!!";
}


/**
 * Метод Гаусса с выбором главного элемента
 */

void fedyanovaam::lab2()
{
  double p;
    int maxn;

      for (int k=0; k<N-1; k++){
          maxn = k;
          for (int i=k+1; i<N; i++)
        if(abs(A[i][k]) > abs(A[maxn][k])) maxn = i; ///Выбор главного элемента
          std::swap(A[maxn], A[k]);///Меняем строки местами
          std::swap(b[maxn], b[k]);

          for (int i=k+1; i<N; i++){
              p = A[i][k]/A[k][k];
              for (int j=k; j<N; j++)
                  A[i][j] -= p*A[k][j];
              b[i] -= p*b[k];
          }
      }

      for(int i = 0; i<N; i++){
          x[i]=b[i];
      }

      for (int i=N-1; i>=0; i--){//Перебираем массив
          for (int j=i+1;j<N;j++)
              x[i]-=A[i][j]*x[j];
          x[i] /= A[i][i];

      }  

      }

//}




/**
 * Метод прогонки
 */

void fedyanovaam::lab3()
{
 double *P = new double [N]; ///Коэффициенты alfa
  double *Q = new double [N]; ///Коэффициенты betta

    P[0] = -A[0][1]/A[0][0];
    Q[0] = b[0]/A[0][0];

    for(int i=1; i<N; i++){ ///Определим прогоночные коэффициенты
    
      P[i] = A[i][i+1]/(-A[i][i] - A[i][i-1]*P[i-1]);
      Q[i] = (-b[i] + A[i][i-1]*Q[i-1])/(-A[i][i] - A[i][i-1]*P[i-1]);
    }


    x[N-1] = Q[N-1];
    for(int i=N-2; i>=0; i--) ///Определим решение
      x[i] = P[i]*x[i+1] + Q[i];



    x[N-1] = Q[N-1];
    for(int i=N-2; i>=0; i--) ///Определим решение
      x[i] = P[i]*x[i+1] + Q[i];


    delete [] P;
    delete [] Q;
}




/**
 * Метод простых итераций
 */
void fedyanovaam::lab4()
{
double eps = 1e-13;
double tauh = 1e-5;
for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double x1;
	double *xr = new double[N];
	int step = 0;

	do {
		step++;
		for (int i = 0; i < N; i++) {
			xr[i] = x[i];
			for (int k = 0; k < N; k++)
				xr[i] -= tauh*A[i][k] * x[k];
			xr[i] += tauh * b[i];

		}
		x1 = 0.;
		for (int i = 0; i < N; i++) {
			x1 += (xr[i]-x[i])*(xr[i]-x[i]);
		}

		for (int i = 0; i < N; i++) {
			x[i] = xr[i];
		}
	} while (sqrt(x1)>eps);
}



/**
 * Метод Якоби или Зейделя
 */
void fedyanovaam::lab5()
{ //Якоби
    const double eps = 10E-20;

	double* y = new double[N];
	double r = 0; /// норма, определяемая как наибольшая разность компонент столбца иксов соседних итераций.

	for(int i=0; i<N; i++)
	{
		x[i] = 0;
	}

	do
	{
		for(int i=0; i<N; i++)
		{
			y[i] = b[i];
			for(int j=0; j<N; j++)
			{
				if(i != j)
				{
					y[i] -= A[i][j]*x[j];
				}
			}
			y[i] /= A[i][i];
		}

		r = abs(x[0] - y[0]);

		for(int i=0; i<N; i++)
		{
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
void fedyanovaam::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void fedyanovaam::lab7()
{

}


void fedyanovaam::lab8()
{

}
void fedyanovaam::lab9()
{

}


static double f(double x)
{
return ((x*x)/10-exp(2*(x)));
}


void fedyanovaam::lab10()
{ //метод половинного деления
double a=1;
double b=2;
double e=0.001;
double c;
int i=0;
if (f(a)*f(b)>0) cout<<"На этом промежутке корней нет"<<endl;
else {
do
{
{c=(a+b)/2 ;
if (f(a)*f(c)<0)
b=c;
else   a=c;
i++;}}
while((fabs(b-a)>e)&&(f(c)!=0));
cout<<"x="<<c<<"\n";}
}



std::string fedyanovaam::get_name()
{
  return "Fedyanova A.M.";
}
