#include "nefedovms.h"

/**
 * Введение в дисциплину
 */
void nefedovms::lab1()
{
std::cout<<"hello world";
std::cout<<"hello world";
std::cout<<"hello world";
}



/**
 * Метод Гаусса с выбором главного элемента
 */
void nefedovms::lab2()
{
	double Q = 0;

    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            Q = A[j][i] / A[i][i];

            for (int k = i; k < N; k++)
            {
                A[j][k] -= Q * A[i][k];
            }
            b[j] -= Q * b[i];
        }
    }

    for(int i = 0; i < N; i++)
	{
        x[i] = b[i];
	}

    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
		{
			x[i] -= A[i][j] * x[j];
		}

        x[i] /= A[i][i];
	}
}



/**
 * Метод прогонки
 */
void nefedovms::lab3()
{
double Q = 0;
	int max;

    for (int i = 0; i < N - 1; i++)
    {
				max = i;

		for (int j = i + 1; j < N; j++)
		{
			if(abs(A[j][i]) > abs(A[max][i]))
			{
				max = j;
			}
		}


				std::swap(A[max], A[i]);
		    std::swap(b[max], b[i]);

        for (int j = i + 1; j < N; j++)
        {
            Q = A[j][i] / A[i][i];
            for (int k = i; k < N; k++)
            {
                A[j][k] -= Q * A[i][k];
            }
            b[j] -= Q * b[i];
        }
    }

    for(int i = 0; i < N; i++)
	{
        x[i] = b[i];
	}

    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
		{
			x[i] -= A[i][j] * x[j];
		}

        x[i] /= A[i][i];
	}
}



/**
 * Метод простых итераций
 */
void nefedovms::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void nefedovms::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void nefedovms::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void nefedovms::lab7()
{

}


void nefedovms::lab8()
{

}


void nefedovms::lab9()
{

}


std::string nefedovms::get_name()
{
  return "Nefedov M.S.";
}
