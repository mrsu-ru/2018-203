#include "velmiskinaav.h"

/**
 * Введение в дисциплину
 */
void velmiskinaav::lab1()
{
	std::cout << "Hello world!!!" << std::endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void velmiskinaav::lab2()
{
int i, j, k, l;
        //прямой ход
        for(i = 0; i < N; i++)
        {
         for(k = i + 1, l = i; k < N; k++)
          if (std::abs(A[k][i]) > std::abs(A[l][i]))
                l = k;

          std::swap(A[l],A[i]);
          std::swap(b[l],b[i]);

          b[i]/=A[i][i];
          for(j=N-1; j>i; A[i][j--]/=A[i][i]);

          A[i][i]=1;

          for(int j=i+1; j<N;j++)
            {
            for(k=N-1;k>i;k--)
            A[j][k]-=A[i][k]*A[j][i];
            b[j]-=b[i]*A[j][i];
            A[j][i]=0;
            }

        }
          //обратный ход
          double *x=new double[N];

          x[N-1]=b[N-1];
          for ( int i = N - 2; i >= 0; i--)
         {
           x[i] = b[i];
           for ( int j = i + 1; j < N; j++)
            {
              x[i] -= A[i][j] * x[j];
            }
         }
}



/**
 * Метод прогонки
 */
void velmiskinaav::lab3()
{
 int i;
 double znam;

    b[0] /= A[0][0];//Q[1]
    A[0][1] /= -A[0][0];//P[1]

    for(i = 1;i < N-1;i++)
    {
        znam = - A[i][i] - A[i][i - 1] * A[i - 1][i]; //общий знаменатель для формул нахождения P[i], Q[i]
        A[i][i + 1] /= znam; //P[i]
        b[i] = (A[i][i - 1] * b[i - 1] - b[i]) / znam; //Q[i]
    }
        //строка ниже для вычисления Q[N]
    b[N - 1] = (A[N - 1][N - 2] * b[N - 2] - b[N - 1]) / (-A[N - 1][N - 1] -A[N - 1][N - 2] * A[N - 2][N - 1]);


        //обратный ход
    for(i = N - 2; i > -1; i--)
    {
        b[i] += b[i + 1] * A[i][i + 1];
    }
}



/**
 * Метод простых итераций
 */
void velmiskinaav::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void velmiskinaav::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void velmiskinaav::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void velmiskinaav::lab7()
{

}


void velmiskinaav::lab8()
{

}


std::string velmiskinaav::get_name()
{
  return "Velmiskina A. V.";
}
