#include "velmiskinaav.h"

/**
 * Введение в дисциплину
 */
void velmiskinaav::lab1()
{
	std::cout << "Hello world!!!" << std::endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void velmiskinaav::lab2()
{
int i, j, k, l;

    for(i=0;i<N;i++)
        {
          k=i;
          for(l=i+1; l<N; l++)
          if (abs(A[l][i])>abs(A[k][i]))
                k=l;

          std::swap(A[k],A[i]);
          std::swap(b[k],b[i]);
           b[i]/=A[i][i];
          for(j=N-1; j>i; A[i][j--]/=A[i][i]);

          A[i][i]=1;

          for(int j=i+1; j<N;j++)
		  {
          for(l=N-1;l>i;l--)
          A[j][l]-=A[i][l]*A[j][i];
          b[j]-=b[i]*A[j][i];

          A[j][i]=0;
          }
        }
		  
		  
          x[N-1]=b[N-1];
          for ( int i = N - 2; i >= 0; i-- )
         {
           x[i] = b[i];
           for ( int j = i + 1; j < N; j++ ) 
		   {
              x[i] -= A[i][j] * x[j];
           }
         }        
 }



/**
 * Метод прогонки
 */
void velmiskinaav::lab3()
{
 double *P = new double[N];
double *Q = new double[N];
int z;

 P[0]=A[0][1]/(-A[0][0]);
 Q[0]=b[0]/A[0][0];

for(int i=1;i<N;i++)
    {

     P[i] = A[i][i+1]/(-A[i][i-1]*P[i-1]-A[i][i]);
     Q[i] = (-b[i] + A[i][i-1]*Q[i-1])/(-A[i][i-1]*P[i-1]-A[i][i]);

    }
    x[N-1] = (b[N-1] - A[N-1][N - 2] * Q[N - 2]) / (A[N-1][N-1] + A[N-1][N-1] * P[N-1]);
for(int i=N-2;i>=0;i--)
    x[i]=P[i]*x[i+1]+Q[i];
}



/**
 * Метод простых итераций
 */
void velmiskinaav::lab4()
{
	
}



/**
 * Метод Якоби
 */
void velmiskinaav::lab5()
{
double eps =0.00001;

		double* Y = new double[N];
		double norma = 0;

		for(int i=0; i<N; i++){
			x[i] = 0;}
        do
		{
			for(int i=0; i<N; i++)
			{
				Y[i] = b[i];
				for(int j=0; j<N; j++)
				{
					if(i != j)
					{
						Y[i] -= A[i][j]*x[j];
					}
				}
				Y[i] /= A[i][i];
			}

		norma = abs(x[0] - Y[0]);

			for(int i=0; i<N; i++)
			{
				if(abs(x[i]-Y[i]) > norma)
					norma = sqrt((x[i]-Y[i])*(x[i]-Y[i]));
				x[i] = Y[i];
			}
		} while (norma >= eps);
		delete[] Y;
}

void MatrVekt(int N, double **A, double *V, double *R)//перемножения матрицы на вектор
//N- размерность, A- матрица, V- вектор, R- результат
{
for(int i=0; i<N; i++)
        {
        R[i]=0;
        for(int j=0; j<N; j++)
              R[i]+= A[i][j]*V[j];
        }
}

/**
 * Метод минимальных невязок
 */
void velmiskinaav::lab6()
{
//N- размерность, A- матрица, b- вектор свободных членов, x-вектор результат
//int count=0;//  количество итераций
double *R = new double [N];
double *Delta = new double [N];
double *TempX = new double[N];
//double *x = new double[N];
double maxi=0.0, Tau=0.0, TempTau=0.0;
double eps = 0.0000001;
for (int i=0; i<N; i++)
    TempX[i]=0;//первое приближение задаём нулевым
do
{
MatrVekt(N, A, TempX, R);
for(int i=0; i<N; i++)
    {
    Delta[i]=R[i]-b[i];//Вектор невязок
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
void velmiskinaav::lab7()
{

}


void velmiskinaav::lab8()
{

}

void velmiskinaav::lab9()
{

}

//метод касательных
void velmiskinaav::lab10()
{
   double f, x, df;
double  EPS = 0.000001;
  do {
    f = pow(x,3) + 3*x + 1;
    df = 3*pow(x, 2)+3;
    x = x - f / df;
  } while (fabs(f) > EPS);

  cout << "x0= "<< x;
}


std::string velmiskinaav::get_name()
{
  return "Velmiskina A. V.";
}
