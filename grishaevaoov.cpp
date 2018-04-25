#include "grishaevaoov.h"

/**
 * Введение в дисциплину
 */
void grishaevaov::lab1()
{
std::cout<<"Hello, world!!!";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void grishaevaov::lab2()
{
double maxi;
  int k, ii;
 const double eps = 0.000000000000000000000000000000001;

  k = 0;
  while (k < N)
    {
    maxi = abs(A[k][k]);
    ii = k;
    for (int i = k + 1; i < N; i++) {
      if (abs(A[i][k]) > maxi) {
        maxi = abs(A[i][k]);
        ii = i;
      }
    }
  	for (int j = 0; j < N; j++) {
       swap (A[k][j], A[ii][j]);
    }
    swap(b[k], b[ii]);
    for (int i = k; i < N; i++) {
      double c = A[i][k];
      if (abs(c) < eps)
      continue;
      for (int j = 0; j < N; j++) {
        A[i][j] = A[i][j] / c;
      }
      b[i] = b[i] / c;
      if (i == k)
       continue; 
      for (int j = 0; j < N; j++) {
        A[i][j] = A[i][j] - A[k][j];
      }
      b[i] = b[i] - b[k];
    }
    k++;
  }
  for (k = N - 1; k >=0; k--) {
  x[k] = b[k];
    for (int i = k-1; i >=0; i--) {
       double c=A[i][k];
            for (int j=0; j<N; j++)
           A[i][j]=A[k][j]*c+A[i][j];
      b[i]=-b[k]*c+b[i];
}
  }
}



/**
 * Метод прогонки
 */
void grishaevaov::lab3()
{
int n=N-1;
  double AA[N];
    double B[N];
int z;
 AA[0]=-A[0][1]/A[0][0];
 B[0]=b[0]/A[0][0];

for(int i=1;i<N;i++)
    {

     AA[i] = -A[i][i+1]/(A[i][i-1]*AA[i-1]+A[i][i]);
     B[i] = (b[i] - A[i][i-1]*B[i-1])/(A[i][i-1]*AA[i-1]+A[i][i]);

    }
    x[N]=b[N];
for(int i=N-1;i>=0;i--)
    x[i]=AA[i]*x[i+1]+B[i];

}



/**
 * Метод простых итераций
 */
void grishaevaov::lab4()
{
double Eps=0.00000001;
double Err;
double *nx = new double[N];
for (int i=0;i<N;i++)
	x[i]=b[i];
int step=0;
	
do{
step++;
  for(int i=0;i < N;i++)
  {
   nx[i]=-b[i];
 
   for(int j=0;j < N;j++)
   {
    if(i!=j)
     nx[i]+=A[i][j]*x[j];
   }
 
   nx[i]/=-A[i][i];
  }
  Err=0;
for(int i=0; i<N; i++) { 
if(std::abs(x[i]-nx[i]) > Err)
Err = std::abs(x[i]-nx[i]);
}
for(int i=0; i<N; i++) 
	x[i]=nx[i];
std::cout<<step<<"     "<<Err<<endl;
}while (Err>Eps);

delete[] nx;
}



/**
 * Метод Якоби или Зейделя
 */
void grishaevaov::lab5()
{
double eps = 0.00000001; 
double px[N]; 
double sum = 0; 
int flag = 1; 
for(int i = 0; i < N; i++) 
{ 
x[i] = 0; 
} 
do 
{ 
for(int i = 0; i < N; i++) 
{ 
px[i] = x[i]; 
} 
for(int i = 0; i < N; i++) 
{ 
sum = 0; 
for(int j = 0; j < i; j++) 
{ 
sum =sum + A[i][j] * x[j]; 
} 
for(int j = i+1; j < N; j++) 
sum =sum + A[i][j] * px[j]; 
x[i] = (b[i] - sum) / A[i][i]; 
} 
flag = 1; 
for(int i = 0; i < N; i++) 
{ 
if(fabs(x[i] - px[i]) < eps) 
{ 
flag = 0; 
break; 
} 
} 
} while(flag == 1); 
}



/**
 * Метод минимальных невязок
 */
void grishaevaov::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void grishaevaov::lab7()
{

}


void grishaevaov::lab8()
{

}

void grishaevaov::lab9()
{

}
void grishaevaov::lab10()
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

std::string grishaevaov::get_name()
{
  return "Grishaeva O.V.";
}
