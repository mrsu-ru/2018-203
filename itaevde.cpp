#include "itaevde.h"

/**
 * Введение в дисциплину
 */
void itaevde::lab1()
{
std::cout << "hello world";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void itaevde::lab2()
{
double max;
	int k, index;
	const double eps = 0.00001;
	k = 0;
	while (k < N){
		max = abs(A[k][k]);
		index = k;
		for (int i = k + 1; i < N; i++){
			if (abs(A[i][k]) > max){
				max = abs(A[i][k]);
				index = i;
			}
		}
		
		if (max < eps){
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return;
		}
		for (int j = 0; j < N; j++)
		    swap(A[k][j], A[index][j]);
		swap(b[k],b[index]);

		for (int i = k; i < N; i++){
			double temp = A[i][k];
			if (abs(temp) < eps) 
				continue; 
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] / temp;
			b[i] = b[i] / temp;
			if (i == k)
				continue; 
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] - A[k][j];
			b[i] = b[i] - b[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = N - 1; k >= 0; k--){
		x[k] = b[k];
		for (int i = 0; i < k; i++)
			b[i] = b[i] - A[i][k] * x[k];
	}
}



/**
 * Метод прогонки
 */
void itaevde::lab3()
{

}



/**
 * Метод простых итераций
 */
void itaevde::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void itaevde::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void itaevde::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void itaevde::lab7()
{

}


void itaevde::lab8()
{

}


void itaevde::lab9()
{

}


std::string itaevde::get_name()
{
  return "Itaev D.E.";
}
