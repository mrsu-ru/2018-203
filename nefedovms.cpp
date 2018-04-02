#include "nefedovms.h"

/**
 * Введение в дисциплину
 */
void nefedovms::lab1()
{
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

    double* new_A = new double[N];
    double* new_b = new double[N];

    new_A[0] = A[0][1] / (-A[0][0]);
    new_b[0] = b[0] / A[0][0];

    for(int i = 1; i < N; i++)
    {
        new_A[i] = A[i][i+1] / (-A[i][i-1] * new_A[i-1] - A[i][i]);
        new_b[i] = (-b[i] + A[i][i-1] * new_b[i-1]) / ( -A[i][i-1] * new_A[i-1] - A[i][i]);
    }

    for(int i = N - 1; i >= 0; i--)
    {
        x[i] = new_A[i] * x[i+1] + new_b[i];
    }

    delete[] new_A;
    delete[] new_b;

}



/**
 * Метод простых итераций
 */
void nefedovms::lab4()
{

}



/**
 * Метод Якоби
 */
void nefedovms::lab5()
{
double eps = 1.e-10;
double nx[N];
double razn = 0;
int flag=1;
	for(int i = 0; i < N; i++)
	{
		x[i] = 0;
	}
	do
	{
		for(int i = 0; i < N; i++)
		{
			nx[i] = b[i];
			for(int j = 0; j < N; j++)
			{
				if(i != j)
				{
					nx[i] = nx[i] - A[i][j] * x[j];
				}
			}
			nx[i] = nx[i]/(A[i][i]);
		}
		razn = fabs( nx[0] - x[0] );
		for(int i = 0; i < N; i++)
		{
			if(fabs(nx[i] - x[i]) > razn)
			{
				razn = fabs(nx[i] - x[i]);
			}
			x[i] = nx[i];
		}
		if(razn < eps)
            flag=0;
	} while(flag==1);
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
