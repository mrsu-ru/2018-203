#include "borisovrs.h"

/**
 * Введение в дисциплину
 */
void borisovrs::lab1()
{
std::cout<<"hello world";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void borisovrs::lab2()
{
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

double determ(double** a, int n)
{
int i, j;
double det = 0;
double** matr;
if (n == 1)
{
det = a[0][0];
}
else if (n == 2)
{
det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
}
else
{
matr = new double*[n - 1];
for (i = 0; i<n; ++i)
{
for (j = 0; j<n - 1; ++j)
{
if (j<i)
matr[j] = a[j];
else
matr[j] = a[j + 1];
}
det += pow((double)-1, (i + j));
matr, n - 1*a[i][n - 1];
}
delete[] matr;
}
return det;
}

double * gauss(double **a, double *y, int n)
{
double *x, max;
int k, index;
const double eps = 0.00001; // точность
x = new double[n];
k = 0;
while (k < n)
{
// Поиск строки с максимальным a[i][k]
max = abs(a[k][k]);
index = k;
for (int i = k + 1; i < n; i++)
{
if (abs(a[i][k]) > max)
{
max = abs(a[i][k]);
index = i;
}
}

if (max < eps)
{
cout << "Решения нет " << endl;

return 0;
}
for (int j = 0; j < n; j++)
{
double temp = a[k][j];
a[k][j] = a[index][j];
a[index][j] = temp;
}
double temp = y[k];
y[k] = y[index];
y[index] = temp;
// Нормализация уравнений
for (int i = k; i < n; i++)
{
double temp = a[i][k];
if (abs(temp) < eps) continue;
for (int j = 0; j < n; j++)
a[i][j] = a[i][j] / temp;
y[i] = y[i] / temp;
if (i == k) continue;
for (int j = 0; j < n; j++)
a[i][j] = a[i][j] - a[k][j];
y[i] = y[i] - y[k];
}
k++;
}

// обратная подстановка
for (k = n - 1; k >= 0; k--)
{
x[k] = y[k];
for (int i = 0; i < k; i++)
y[i] = y[i] - a[i][k] * x[k];
}
return x;
}

int main()
{
setlocale(LC_ALL, "rus");
double **a, *y, *x;
int n;
cout << "Razmernost': ";
cin >> n;
a = new double*[n];
y = new double[n];
for (int i = 0; i < n; i++)
{
a[i] = new double[n];
for (int j = 0; j < n; j++)
{
a[i][j] = rand() % 10;
}
}
for (int i = 0; i < n; i++)
{
for (int j = 0; j < n; j++) {
cout.width(3);
cout << a[i][j];
}
cout <<endl;
}
for (int i = 0; i < n; i++)
{

y[i] = rand() % 10;

}
cout << endl;
for (int i = 0; i < n; i++)
{
cout.width(3);
cout << y[i] << endl;

}
cout << endl;
cout << determ(a, n) << endl;
x = gauss(a, y, n);
cout << endl;
if (x != 0)
{

for (int i = 0; i < n; i++)
cout << x[i] << endl;
}


for (int i = 0; i < n; i++)
delete[] a[i];
delete[] a;
delete[] x;
delete[] y;
}


}



/**
 * Метод прогонки
 */
void borisovrs::lab3()
{

}



/**
 * Метод простых итераций
 */
void borisovrs::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void borisovrs::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void borisovrs::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void borisovrs::lab7()
{

}


void borisovrs::lab8()
{

}


std::string borisovrs::get_name()
{
  return "Borisov R.S.";
}
