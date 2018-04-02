#include "salinaa.h"

/**
 * Введение в дисциплину
 */
void salinaa::lab1()
{
std::cout<<"Hello world";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void salinaa::lab2()
{
double max;
	int k, index;
	const double eps = 0.00001;
	k = 0;
	while (k < N)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(A[k][k]);
		index = k;
		for (int i = k + 1; i < N; i++)
		{
			if (abs(A[i][k]) > max)
			{
				max = abs(A[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return;
		}
		for (int j = 0; j < N; j++)
		{
		    swap(A[k][j], A[index][j]);

		}
		swap(b[k],b[index]);


		for (int i = k; i < N; i++)
		{
			double temp = A[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] / temp;
			b[i] = b[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] - A[k][j];
			b[i] = b[i] - b[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = N - 1; k >= 0; k--)
	{
		x[k] = b[k];
		for (int i = 0; i < k; i++)
			b[i] = b[i] - A[i][k] * x[k];
	}
}



/**
 * Метод прогонки
 */
void salinaa::lab3()
{
	int N1=N-1;
	double *alfa = new double[N];
	double *beta = new double[N];
	double y = A[0][0];
	alfa[0] = -A[0][1] / y;
	beta[0] = b[0] / y;
	for (int i = 1; i < N1; i++) {
		y = A[i][i] + A[i][i - 1] * alfa[i - 1];
		alfa[i] = -A[i][i + 1] / y;
		beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
	}

	x[N1] = (b[N1] - A[N1][N1-1] * beta[N1-1]) / (A[N1][N1] + A[N1][N1-1] * alfa[N1-1]);
	for (int i = N1-1; i >= 0; i--) {
		x[i] = alfa[i] * x[i + 1] + beta[i];
	}
}



/**
 * Метод простых итераций
 */
void salinaa::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void salinaa::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void salinaa::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void salinaa::lab7()
{

}


void salinaa::lab8()
{

}


void salinaa::lab9()
{

}


std::string salinaa::get_name()
{
  return "Salin A.A.";
}
