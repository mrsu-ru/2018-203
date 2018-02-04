#include "polyakovda.h"

/**
 * Введение в дисциплину
 */
void polyakovda::lab1()
{
	for (int i = 0; i < N; i++){
		double maxEl = 0;
		int indRow = i;
		for (int j = i + 1; j<N; j++)
			if (maxEl<abs(A[j][i])){
				indRow = j;
				maxEl = abs(A[j][i]);
			}


		if (indRow != i){
			for (int j = i; j < N; j++){
				swap (A[i][j], A[indRow][j]);
			}
			swap(b[i], b[indRow]);
		}

		maxEl = A[i][i];
		b[i] /= A[i][i];
		A[i][i] = 1;
		for (int j = i + 1; j<N; j++)
			A[i][j] /= maxEl;

		for (int k = i + 1; k < N; k++){
			double multiplier = A[k][i];
			A[k][i] = 0;
			if (multiplier != 0) {
				for (int j = i + 1; j < N; j++)
					A[k][j] -= multiplier*A[i][j];
				b[k] -= multiplier*b[i];
			}
		}
	}

	for (int k = N - 1; k >= 0; k--){
		x[k] = 0;
		double sum = b[k];
		for (int j = N-1; j>k; j--)
			sum -= A[k][j] * x[j];
		sum -= b[k] * x[k];
		x[k] = sum;
	}

}


/**
 * Метод Гаусса с выбором главного элемента
 */
void polyakovda::lab2()
{
	double *apr = new double[N];
	double *bpr = new double[N];
	double y = A[0][0];
	apr[0] = -A[0][1] / y;
	bpr[0] = b[0] / y;
	for (int i = 1; i < N-1; i++) {
		y = A[i][i] + A[i][i - 1] * apr[i - 1];
		apr[i] = -A[i][i + 1] / y;
		bpr[i] = (b[i] - A[i][i - 1] * bpr[i - 1]) / y;
	}

	x[N-1] = (b[N-1] - A[N-1][N-2] * bpr[N-2]) / (A[N-1][N-1] + A[N-1][N-2] * apr[N-2]);
	for (int i = N-1 - 1; i >= 0; i--) {
		x[i] = apr[i] * x[i + 1] + bpr[i];
	}
}



/**
 * Метод прогонки
 */
void polyakovda::lab3()
{

}



/**
 * Метод простых итераций
 */
void polyakovda::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void polyakovda::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void polyakovda::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void polyakovda::lab7()
{


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
}


void polyakovda::lab8()
{
/**
 * Нахождение наибольшего по модолю собственного значения матрицы
 */
}


std::string polyakovda::get_name()
{
  return "Polyakov D.A.";
}
