#include "tarasovams.h"

/**
 * ¬ведение в дисциплину
 */
void tarasovams::lab1()
{

}


/**
 * ћетод √аусса с выбором главного элемента
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

        }

/**
 * ћетод прогонки
 */
void tarasovams::lab3()
{
double* Alpha = new double[N];
double* Betta = new double[N];

    Alpha[0] = -A[0][1]/A[0][0];
    Betta[0] = b[0]/A[0][0];

    for(int i=1; i<N; i++)
    {
        Alpha[i] = -A[i][i+1]/(A[i][i-1]*Alpha[i-1]-A[i][i]);
        Betta[i] = (-b[i] + A[i][i-1]*Betta[i-1])/(-A[i][i-1]*Alpha[i-1]-A[i][i]);
    }

    for(int i=N-1; i>=0; i--)

    x[i] = Alpha[i]*x[i+1]+Betta[i];
    delete[] Alpha;
    delete[] Betta;
}





/**
 * ћетод простых итераций
 */
void tarasovams::lab4()
{

}



/**
 * ћетод якоби или «ейдел€
 */
void tarasovams::lab5()
{

}



/**
 * ћетод минимальных нев€зок
 */
void tarasovams::lab6()
{

}



/**
 * ћетод сопр€женных градиентов
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
