#include "syusinaev.h"





/**
 * Введение в дисциплину
 */
void syusinaev::lab1()
{
	cout<<"It's working!!!";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void syusinaev::lab2()
{
	   int i, l, k;

 for(i=0;i<N;i++)
        {
          for(k=i+1,l=i ;k<N; k++)
              if (abs(A[k][i])>abs(A[l][i])) l=k;

          swap(A[l],A[i]);
          swap(b[l],b[i]);

           b[i]/=A[i][i];

          for(  int j=N-1 ;j>i;A[i][j--]/=A[i][i]);
          A[i][i]=1;

          for(int j=i+1; j<N;j++){
            for(k=N-1;k>i;k--)
                A[j][k]-=A[i][k]*A[j][i];
              b[j]-=b[i]*A[j][i];
            A[j][i]=0;
          }
       cout<<endl;

        }//прямой ход
       x[N-1]=b[N-1];
       for ( int i = N - 2; i >= 0; i-- )
         {
           x[i] = b[i];
           for ( int j = i + 1; j < N; j++ ) {
              x[i] -= A[i][j] * x[j];
              }
                  cout<<endl;

         }//обратный ход

}



/**
 * Метод прогонки
 */
void syusinaev::lab3()
{
double y;
double* Al= new double[N];
double* Be= new double[N];
 y = A[0][0];
  Al[0] = -A[0][1] / y;
  BE[0] = b[0] / y  ;
  for (int i = 1; i < N; i++) {
    y = A[i][i] + A[i][i - 1] * Al[i - 1];
    Al[i] = -A[i][i + 1] / y;
    Be[i] = (b[i] - A[i][i - 1] * Be[i - 1]) / y;
  }
  x[N] = (b[N] - A[N][N - 1] * BE[N - 1]) / (A[N][N] + A[N][N - 1] * Al[N - 1]);
  for (int i = N - 1; i >= 0; i--) {
    x[i] = Al[i] * x[i + 1] + Be[i];
  }

}

/**
 * Метод простых итераций
 */
void syusinaev::lab4()
{
double eps = 1e-9;
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
		printf("err = %f, step = %d\n", x1, step);
	} while (sqrt(x1)>eps);
}


/**
 * Метод Якоби или Зейделя
 */
void syusinaev::lab5()
{//Якоби
double eps=0.001;

 double* temp = new double[N];
	double norm; 

	do {
		for (int i = 0; i < N; i++) {
			temp[i] = b[i];
			for (int g = 0; g < N; g++) {
				if (i != g)
					temp[i] -= A[i][g] * x[g];
			}
			temp[i] /= A[i][i];
		}
        norm = abs(x[0] - temp[0]);
		for (int h = 0; h < N; h++) {
			if (fabs(x[h] - temp[h]) > norm)
				norm = abs(x[h] - temp[h]);
			x[h] = temp[h];
		}
	} while (norm > eps);
	delete[] temp;
}



/**
 * Метод минимальных невязок
 */
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
void syusinaev::lab6()
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
void syusinaev::lab7()
{

}

/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void syusinaev::lab8()
{

}
/**
 * Нахождение наибольшего по модолю собственного значения матрицы
 */

void syusinaev::lab9()
{

}
static double f(double x) 
{ 
return (2*(log(x))-(x/2)); 
} 
static double f1(double x) 
{ 
return ((2.0/x) - (1.0/2.0)); 
} 
static double f2(double x) 
{ 
return ((-2.0)/(x*x)); 
} 

void syusinaev::lab10()
{
double a=1,b=2 c,e=0.001; 
if( (2*(log(a))-(a/2))*((-2.0)/(a*a))>0 ) c=a; 
else c=b; 
do { 
c=c-(2*(log(c))-(c/2))/((2.0/c) - (1.0/2.0));  
} 
while (abs((2*(log(c))-(c/2)))>=e); 
}

std::string syusinaev::get_name()
{
  return "Syusina E.V.";
}















