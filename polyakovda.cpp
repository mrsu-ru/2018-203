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
	double eps = 0.001;

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
		if (2 * abs(A[i][i]) < summ) { cout << "Error" << endl; system("pause"); return 0; }

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
			for (int k = 0; k < N; k++)
				xr[i] -= A[i][k] * x[k];
			xr[i] += b[i];
		}
		x1 = x[0];
		for (int i = 0; i < N; i++) {
			x[i] = xr[i];
		}
	} while (abs(x[0] - x1)>eps);
}



/**
 * Метод простых итераций
 */
void polyakovda::lab4()
{
	double eps = 0.001;
	double *xResult = new double[N];
	double delta;
	
	for (int i = 0; i < N; i++) {
		x[i] = b[i];
	}

	do {
		for (int i = 0; i < N; i++) {
			xResult[i] = b[i];
			for (int j = 0; j < N; j++) {
				if (i != j)
					xResult[i] -= A[i][j] * x[j];
			}
			xResult[i] /= A[i][i];
		}
		delta = abs(x[0] - xResult[0]);
		for (int i = 0; i < N; i++) {
			if (delta < abs(x[i] - xResult[i]))
				delta = abs(x[i] - xResult[i]);
			x[i] = xResult[i];
		}
	} while (eps < delta);
}



/**
 * Метод Якоби или Зейделя
 */
void polyakovda::lab5()
{
	double eps = 0.001;
	double delta, r, rModul;

	double *w = new double[N];
	double *v = new double[N];
	double *result = new double[N];

	for (int i = 0; i<N; i++)
		result[i] = 0;

	do{
		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++)
				w[i] += A[i][j] * result[j];
		}

		for (int i = 0; i < N; i++) {
			v[i] = w[i] - b[i];
		}

		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++)
				w[i] += A[i][j] * v[j];
		}

		r = 0.0;
		rModul = 0.0;
		for (int i = 0; i < N; i++) {
			r += w[i] * v[i];
			rModul += w[i] * w[i];
		}
		if (r==rModul) {r=1;}
		else {r = r / rModul;}
		for (int i = 0; i < N; i++)
			x[i] = result[i] - r*v[i];
		delta = abs(x[0] - result[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - result[i])>delta)
				delta = abs(x[i] - result[i]);
			result[i] = x[i];
		}
	} while (eps < delta);
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
