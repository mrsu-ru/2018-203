#include "scherbakovdv.h"

/**
 * Введение в дисциплину
 */
void scherbakovdv::lab1()
{
	
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void scherbakovdv::lab2()
{
	//A[i][i] - матрица, b[i] - столбец свободных членов
	memcpy(x,b)
	for (int i = 0; i < N; i++)
		{
			if (!A[i][i])
			{
				double* MaxRow = A[i];
				for (int j = i + 1; j < N; j++)
					if (A[j][i] > MaxRow[i])
					{
						double* tmp = MaxRow;
						MaxRow = A[j];
						A[j] = tmp;
					}
				if (!MaxRow[i])
					return;
			}
			for (int j = N - 1; j > i; j--)
			{
				A[i][j] /= A[i][i];
				for (int k = i + 1; k < N; k++)
					A[k][j] -= A[k][i] * A[i][j];
			}
			b[i] /= A[i][i];
			A[i][i] = 1;
			for (int j = N - 1; j > i; j--)
			{
				b[j] -= b[i] * A[j][i];
				A[j][i] = 0;
			}
		}
		for (int i = N - 1; i >= 0; i--)
		{
			for (int j = i - 1; j >= 0; j--)
			{
				b[j] -= b[i] * A[j][i];
				A[j][i] = 0;
			}
		}
		memcpy(x,b,sizeof(double)*N);
}



/**
 * Метод прогонки
 */
void scherbakovdv::lab3()
{

}



/**
 * Метод простых итераций
 */
void scherbakovdv::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void scherbakovdv::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void scherbakovdv::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void scherbakovdv::lab7()
{

}


void scherbakovdv::lab8()
{

}


std::string scherbakovdv::get_name()
{
  return "Scherbakov D.V.";
}
