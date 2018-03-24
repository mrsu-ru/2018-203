#include "bulychevaoa.h"

/**
 * Введение в дисциплину
 */
void bulychevaoa::lab1()
{
std::cout<<"hello world";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void bulychevaoa::lab2()
{
//поиск максимального
	for (int i = 0; i < N; i++) {
		double maxEl = 0;
		int indRow = i;
		for (int j = i + 1; j<N; j++)
			if (maxEl<abs(A[j][i])) {
				indRow = j;
				maxEl = abs(A[j][i]);
			}

//меняем столбцы
		if (indRow != i) {
			for (int j = i; j < N; j++) {
				swap(A[i][j], A[indRow][j]);
			}
			swap(b[i], b[indRow]);
		}

//проходим Гаусом по матрице
		maxEl = A[i][i];
		b[i] /= A[i][i];
		A[i][i] = 1;
		for (int j = i + 1; j<N; j++)
			A[i][j] /= maxEl;

		for (int k = i + 1; k < N; k++) {
			double multiplier = A[k][i];
			A[k][i] = 0;
			if (multiplier != 0) {
				for (int j = i + 1; j < N; j++)
					A[k][j] -= multiplier*A[i][j];
				b[k] -= multiplier*b[i];
			}
		}
	}
//обратное вычисление корнеЙ
	for (int k = N - 1; k >= 0; k--) {
		x[k] = 0;
		double sum = b[k];
		for (int j = N - 1; j>k; j--)
			sum -= A[k][j] * x[j];
		sum -= b[k] * x[k];
		x[k] = sum;
	}

}



/**
 * Метод прогонки
 */
void bulychevaoa::lab3()
{

}



/**
 * Метод простых итераций
 */
void bulychevaoa::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void bulychevaoa::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void bulychevaoa::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void bulychevaoa::lab7()
{

}


void bulychevaoa::lab8()
{

}


void bulychevaoa::lab9()
{

}


std::string bulychevaoa::get_name()
{
  return "Bulycheva O.A.";
}
