#include "bulychevaoa.h"
#include <stdio.h>
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
 * Метод простых итераций
 */
void bulychevaoa::lab4()
{
double eps = 1e-9;
double tauh = 1e-5;	

	/*for (int i = 0; i < N; i++) {
		double maxEl = A[i][i];
		int indRow = i;
		for (int j = i + 1; j < N; j++)
			if (maxEl < abs(A[j][i])) {
				indRow = j;
				maxEl = abs(A[j][i]);
			}


		if (indRow != i) {
			for (int j = i; j < N; j++) {
				swap(A[i][j], A[indRow][j]);
			}
			swap(b[i], b[indRow]);
		}

		double summ = 0;
		for (int j = 0; j < N; j++)
			summ += abs(A[i][j]);
		if (2 * abs(A[i][i]) < summ) { cout << "Error" << std::endl; }

		maxEl = A[i][i];
		b[i] /= A[i][i];
		A[i][i] = 0;
		for (int j = 0; j<N; j++)
			if (j != i) A[i][j] /= maxEl;

	}
	
	*/
	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double x1;
	double *xr = new double[N];
	int step = 0;
	
	do {
		step++;
		for (int i = 0; i < N; i++) {
			xr[i] = x[i];
			for (int k = 0; k < N; k++)
				xr[i] -= tauh*A[i][k] * x[k];
			xr[i] += tauh * b[i];
			
		}
		x1 = 0.;
		for (int i = 0; i < N; i++) {
			x1 += (xr[i]-x[i])*(xr[i]-x[i]);
		}
		
		for (int i = 0; i < N; i++) {
			x[i] = xr[i];
		}
		printf("err = %f, step = %d\n", x1, step);
	} while (sqrt(x1)>eps);
}



/**
 * Метод Якоби или Зейделя
 */
void bulychevaoa::lab5()
{
double eps = 1e-9;

	for (int i = 0; i < N; i++) {
		double maxEl = A[i][i];
		int indRow = i;
		for (int j = i + 1; j < N; j++)
			if (maxEl < abs(A[j][i])) {
				indRow = j;
				maxEl = abs(A[j][i]);
			}


		if (indRow != i) {
			for (int j = i; j < N; j++) {
				swap(A[i][j], A[indRow][j]);
			}
			swap(b[i], b[indRow]);
		}

		double summ = 0;
		for (int j = 0; j < N; j++)
			summ += abs(A[i][j]);
		if (2 * abs(A[i][i]) < summ) { cout << "Error" << std::endl; }

		maxEl = A[i][i];
		b[i] /= A[i][i];
		A[i][i] = 0;
		for (int j = 0; j<N; j++)
			if (j != i) A[i][j] /= maxEl;

	}


	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double x1 = b[0];
	double *xr = new double[N];

	do {
		for (int i = 0; i < N; i++) {
			xr[i] = 0;
			for (int k = 0; k < N; k++){
                if (i != k)
				xr[i] -= A[i][k] * x[k];}
			xr[i] += b[i];
		}
		x1 = x[0];

		for (int i = 0; i < N; i++) {
			x[i] = xr[i];
		}
	} while (abs(x[0] - x1)>eps);
	
	//while ((abs(x[0] - x1)/abs(x[0])>eps); - для Зейделя
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
