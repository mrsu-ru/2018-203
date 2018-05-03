#include "serguninaes.h"

/**
 * Введение в дисциплину
 */
void serguninaes::lab1()
{

}


/**
 * Метод Гаусса с выбором главного элемента
 */
void serguninaes::lab2()
{
    int i, j, maximal;

 for(i=0;i<N;i++)
        {
          for(j=i+1,maximal=i ;j<N; j++)
              if (abs(A[j][i])>abs(A[j][i])) maximal=j;
          if (A[maximal][i]==0) return;

          swap(A[maximal],A[i]);
          swap(b[maximal],b[i]);

           b[i]=b[i]/A[i][i];

          for(  int k=N-1 ;k>i; A[i][k--]/=A[i][i]);
          A[i][i]=1;

          for(int k=i+1; k<N;k++){
            for(j=N-1;j>i;j--)
                A[k][j]-=A[i][j]*A[k][i];
              b[k]-=b[i]*A[k][i];
            A[k][i]=0;
          }
       cout<<endl;

        }

       x[N]=b[N];

       for ( int i = N - 3; i >= 0; i-- )
         {
           x[i] = b[i];
           for ( int k = i + 1; k < N; k++ ) {
              x[i] -= A[i][k] * x[k];
              }
                  cout<<endl;

         }
}



/**
 * Метод прогонки
 */
void serguninaes::lab3()
{
    double* AA = new double [N]; //Коэффициенты "альфа"
double* BB = new double [N]; //Коэффициенты "бетта"

AA[0] = -A[0][1]/A[0][0];
BB[0] = b[0]/A[0][0];

for(int i=1; i<N; i++) //Определяем прогоночные коэффициенты
{
AA[i] = A[i][i+1]/(-A[i][i] - A[i][i-1]*AA[i-1]);
BB[i] = (-b[i] + A[i][i-1]*BB[i-1])/(-A[i][i] - A[i][i-1]*AA[i-1]);
}

x[N-1] = BB[N-1];
for(int i=N-2; i>=0; i--) //Определяем решение
x[i] = AA[i]*x[i+1] + BB[i];

}



/**
 * Метод простых итераций
 */
void serguninaes::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void serguninaes::lab5()
{
    //Метод Якоби
double Eps = 0.0000001;
double* a= new double[N];
double n=0; // норма, определяемая как наибольшая разность компонент столбца иксов соседних итераций

for(int i=0; i<N; i++) //начальное приближение
{
x[i] = 0;
}

do {
for (int i = 0; i < N; i++) {
a[i] = b[i];

for (int k = 0; k < N; k++) {
if (i != k) {a[i] -= A[i][k] * x[k];}
}
a[i] =a[i]/ A[i][i];
}

n = x[0] - a[0];
n=abs(n);

for (int j = 0; j < N; j++) {
if (abs(x[j] - a[j]) > n)
n = x[j] - b[j];
n=abs(n);
x[j] = a[j];
}
} while (n >= Eps);

delete[] a;

}



/**
 * Метод минимальных невязок
 */
void serguninaes::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void serguninaes::lab7()
{

}


void serguninaes::lab8()
{

}


void serguninaes::lab9()
{

}


static double f(double x) {
  return x*x/(10.-exp(2*x));
}

void serguninaes::lab10()
{ //vtnjl gjkjdby ltktybz

double a=1;
double b=2;
double Eps = 0.0000001;
double c;
int i=0;

do
{ c=(a+b)/2 ;
if (f(a)*f(c)<0)
b=c;
else
a=c;

i++;} while((abs(b-a)>Eps)&&(f(c)!=0));

cout<<"x = "<<c<<"\n";

}


std::string serguninaes::get_name()
{
  return "Sergunina E.S.";
}
