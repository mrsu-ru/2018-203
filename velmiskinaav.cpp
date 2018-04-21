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
        //прямой ход
        for(i = 0; i < N; i++)
        {
         for(k = i + 1, l = i; k < N; k++)
          if (std::abs(A[k][i]) > std::abs(A[l][i]))
                l = k;

          std::swap(A[l],A[i]);
          std::swap(b[l],b[i]);

          b[i]/=A[i][i];
          for(j=N-1; j>i; A[i][j--]/=A[i][i]);

          A[i][i]=1;

          for(int j=i+1; j<N;j++)
            {
            for(k=N-1;k>i;k--)
            A[j][k]-=A[i][k]*A[j][i];
            b[j]-=b[i]*A[j][i];
            A[j][i]=0;
            }

        }
          //обратный ход
          double *x=new double[N];

          x[N-1]=b[N-1];
          for ( int i = N - 2; i >= 0; i--)
         {
           x[i] = b[i];
           for ( int j = i + 1; j < N; j++)
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
 int i;
 double znam;
 //double *x=new double[N];

    b[0] /= A[0][0];//Q[1]
    A[0][1] /= -A[0][0];//P[1]

    for(i = 1;i < N-1;i++)
    {
        znam = - A[i][i] - A[i][i - 1] * A[i - 1][i]; //общий знаменатель для формул нахождения P[i], Q[i]
        A[i][i + 1] /= znam; //P[i]
        b[i] = (A[i][i - 1] * b[i - 1] - b[i]) / znam; //Q[i]
    }
        //строка ниже для вычисления Q[N]
    b[N - 1] = (A[N - 1][N - 2] * b[N - 2] - b[N - 1]) / (-A[N - 1][N - 1] -A[N - 1][N - 2] * A[N - 2][N - 1]);


        //обратный ход
    for(i = N - 2; i > -1; i--)
    {
       x[i] = b[i] + b[i + 1] * A[i][i + 1];
    }
}



/**
 * Метод простых итераций
 */
void velmiskinaav::lab4()
{
double norma; //чебышевская норма вектора
 double xn[N]={0};//вектор для текущей итерации, начальное значение
       //должно быть равно начальному приближению
//double *x=new double[N];


 for(int i=0; i < N;i++)
  {
   x[i]=-b[i];

   for(int j=0;j < N;j++)
   {
    if(i!=j)
     x[i]+=A[i][j]*x[j];
   }

   x[i]/=-A[i][i];
  }

  for(int i=0;i < N;i++)
  {
   if(fabs(x[i]-xn[i]) > norma)
    norma=fabs(x[i]-xn[i]);
   xn[i]=x[i];
  }
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
int count=0;//  количество итераций
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

double f(double x) {

    return   pow(x,3) + 3*x + 1;
}

double f1(double x) {

    return   3*pow(x, 2)+3;
}

double f2(double x) {

    return  6*x;
}
void velmiskinaav::lab10()
{
   int a = -3;
    int b = -2;
    double c;
    double eps = 0.0001;

    if(f(a)*f2(a) > 0) c = a;
    else c = b;
    do {
        c=c-f(c)/f1(c);

    }
    while (fabs(f(c))>=eps);
        std::cout<<"c="<<c<<"\n";
}


std::string velmiskinaav::get_name()
{
  return "Velmiskina A. V.";
}
