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
    int t;
for (int k = 0; k < N; k++)
{
//максимальный элемент по модулю в столбце и меняем строки местами
int max=k;
for(int i=k+1;i<N;i++)
if(abs(A[i][k]) > abs(A[max][k]))
max=i;
for(int i=0;i<N;i++)
std::swap(A[k][i],A[max][i]);
std::swap(b[k],b[max]);

t = A[k][k]; //получаем единицы на [k][k]
for (int j = 0; j < N; j++)
A[k][j]= A[k][j]/t;
b[k]=b[k]/t;

for (int i = k + 1; i < N; i++)
{
t = A[i][k];
//получаем нули под единицами
for (int j = 0; j < N; j++)
{
A[i][j] -= A[k][j] * t;
}
b[i] -= b[k] * t;
}
}

//обратный ход
for (int i = N - 1; i > 0; i--)
{
for (int j= i - 1; j >= 0; j--)
{
t = A[j][i];
//получаем нули над единицами
for (int k = 0; k < N; k++)
A[j][k] -= A[i][k] * t;
// b[i] -= b[k] * t;
b[i] -= b[i] * t; // ZHRV:  'k' was changed to 'i'
}
}

for(int i=0; i<N; i++) x[i]=b[i];

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


std::string serguninaes::get_name()
{
  return "Sergunina E.S.";
}
