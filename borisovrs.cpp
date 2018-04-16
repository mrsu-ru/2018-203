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
 double t;

for (int k = 0; k < N; k++) { 
	int maxel=k;
		for(int i=k+1;i<N;i++)
			if(abs(A[i][k]) > abs(A[maxel][k]))
				maxel=i;//ищем наибольший элемент и запоминаем номер строки
	for(int i=0;i<N;i++)
	std::swap(A[k][i],A[maxel][i]);//меняем Kю строку со строкой на которой максимальный элемент
	std::swap(b[k],b[maxel]);
	t = A[k][k];
	for (int j = 0; j < N; j++)
		A[k][j] = A[k][j] / t;//получаем 1 на главной диагонали
	b[k] =b[k]/t;

	for (int i = k + 1; i < N; i++) {
		t = A[i][k];
		for (int j = 0; j < N; j++) {
			A[i][j] =A[i][j]- A[k][j] * t;//Придаем матрице треугольный
		}
	b[i] =b[i]- b[k] * t;
	}
}

//получаем матрицу треугольниго вида и осуществляем обратный ход 
for (int k = N - 1; k > 0; k--){
	for (int i = k - 1; i >= 0; i--){
		t = A[i][k];
		for (int j = 0; j < N; j++)
			A[i][j] =A[i][j]- A[k][j] * t;
		b[i] =b[i] - b[k] * t;
	}
}

for(int i=0; i<N; i++)
	x[i]=b[i];


}



/**
 * Метод прогонки
 */
void borisovrs::lab3()
{
double *a = new double[N];
double *beta = new double[N];
a[0] = -A[0][1] / A[0][0];
beta[0] = b[0] /A[0][0];
for (int i = 1; i < N-1; i++) {
	y = A[i][i] + A[i][i - 1] * a[i - 1];
	a[i] = -A[i][i + 1] / y;
	beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / A[0][0];
}
x[N-1] = (b[N-1] - A[N-1][N-2] * beta[N-2]) / (A[N-1][N-1] + A[N-1][N-2] * a[N-2]);
for (int i = N-1 - 1; i >= 0; i--) {
		x[i] = a[i] * x[i + 1] + beta[i];
}


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


void borisovrs::lab9()
{

}


std::string borisovrs::get_name()
{
  return "Borisov R.S.";
}
