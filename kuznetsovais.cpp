#include "kuznetsovais.h"

/**
 * Введение в дисциплину
 */
void kuznetsovais::lab1()
{

}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kuznetsovais::lab2()
{
int i, j, k, l;

    for(i=0;i<N;i++)
        {
          k=i;
          for(l=i+1; l<N; l++)
          if (abs(A[l][i])>abs(A[k][i]))
                k=l;

          swap(A[k],A[i]);
          swap(b[k],b[i]);
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

//}



/**
 * Метод прогонки
 */
void kuznetsovais::lab3()
{
double *C = new double[N];
double *B = new double[N];
int z;

 C[0]=A[0][1]/(-A[0][0]);
 B[0]=b[0]/A[0][0];

for(int i=1;i<N;i++)
    {

     C[i] = A[i][i+1]/(-A[i][i-1]*C[i-1]-A[i][i]);
     B[i] = (-b[i] + A[i][i-1]*B[i-1])/(-A[i][i-1]*C[i-1]-A[i][i]);

    }
    x[N-1] = (b[N-1] - A[N-1][N - 2] * B[N - 2]) / (A[N-1][N-1] + A[N-1][N-1] * C[N-1]);
for(int i=N-2;i>=0;i--)
    x[i]=C[i]*x[i+1]+B[i];
}



/**
 * Метод простых итераций
 */
void kuznetsovais::lab4()
{
    double eps=0.0000001;
    double summ,L;
    int k;
    double x1[N];
    for(int i=0; i<N; i++) x1[i] = 0;
do
{
for(int i=0; i<N; i++)
{
summ = 0;k=0;
for(int j=0; j<N; j++)
if(i!=j)
summ += A[i][j] * x1[j];
x[i] = (1/A[i][i]) * (x1[i] - summ);
L+=(x1[i]-x[i])*(x1[i]-x[i]);
}
k++;
}while(sqrt(L)>eps);

}



/**
 * Метод Якоби или Зейделя
 */
void kuznetsovais::lab5()
{
double eps =0.001;

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
//}



/**
 * Метод минимальных невязок
 */
void kuznetsovais::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void kuznetsovais::lab7()
{

}


void kuznetsovais::lab8()
{

}

void kuznetsovais::lab9()
{

}


std::string kuznetsovais::get_name()
{
  return "Kuznetsova I.S.";
}
