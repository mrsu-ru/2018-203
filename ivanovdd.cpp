#include "ivanovdd.h"
using namespace std;
/**
 * Введение в дисциплину
 */
void ivanovdd::lab1()
{
std::cout << "Hello World!!!";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void ivanovdd::lab2()
{

for (int k = 0; k < N; k++) {
	int idmax = -1;
	double max = 0;
for (int l=0; l<N; l++) {
	if (abs(A[k][l]) >= max) {
		max = abs(A[k][l]);
		idmax = l;
}
}

if(idmax != -1){
for (int j = 0; j < N; j++)
{
 swap(A[j][idmax], A[j][k]);
}
swap(B[idmax], B[k]);


for (int i = k + 1; i < N; i++)
{
double c_ki = A[i][k]/A[k][k];
for (int j = k; j < N; j++) {
A[i][j] -= A[k][j]*c_ki;
}
b[i] -= b[k]*c_ki;
}
}
}

for (int k = N-1; k >= 0; k--) {
x[k] = b[k];
for (int i = k+1; i < N; i++) {
x[k] -= A[k][i]*x[i];
}
x[k] /= A[k][k];
}
}



/**
 * Метод прогонки
 */
void ivanovdd::lab3()
{

}



/**
 * Метод простых итераций
 */
void ivanovdd::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void ivanovdd::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void ivanovdd::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void ivanovdd::lab7()
{

}


void ivanovdd::lab8()
{

}


std::string ivanovdd::get_name()
{
  return "Ivanov D.D.";
}
