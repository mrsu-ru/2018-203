#include "tarasovams.h"

/**
 * Введение в дисциплину
 */
void tarasovams::lab1()
{

}


/**
 * Метод Гаусса с выбором главного элемента
 */
void tarasovams::lab2()
{
    int i, j, m, c;
    int n = N;
    for(i=0;i<n;i++)
        {
          m=i;
          for(c=i+1; c<n; c++)
          if (abs(A[c][i])>abs(A[m][i]))
                m=c;



          swap(A[m],A[i]);
          swap(b[m],b[i]);
           b[i]/=A[i][i];
          for(j=n-1; j>i; A[i][j--]/=A[i][i]);

          A[i][i]=1;

          for(int j=i+1; j<n;j++){
          for(c=n-1;c>i;c--)
          A[j][c]-=A[i][c]*A[j][i];
          b[j]-=b[i]*A[j][i];

          A[j][i]=0;
          }

            cout<<endl;

          }
          x[n-1]=b[n-1];
          for ( int i = n - 2; i >= 0; i-- )
         {
           x[i] = b[i];
           for ( int j = i + 1; j < n; j++ ) {
              x[i] -= A[i][j] * x[j];
              }
              }

        return;
        }

/**
 * Метод прогонки
 */
void tarasovams::lab3()
{

}



/**
 * Метод простых итераций
 */
void tarasovams::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void tarasovams::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void tarasovams::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void tarasovams::lab7()
{

}


void tarasovams::lab8()
{

}

void tarasovams::lab9()
{

}


std::string tarasovams::get_name()
{
  return "Tarasova M. S.";
}
