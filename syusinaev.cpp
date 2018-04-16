#include "syusinaev.h"

/**
 * Введение в дисциплину
 */
void syusinaev::lab1()
{
	lab2();

}


/**
 * Метод Гаусса с выбором главного элемента
 */
void syusinaev::lab2()
{
	   int i, l, k;

 for(i=0;i<N;i++)
        {
          for(k=i+1,l=i ;k<N; k++)
              if (abs(A[k][i])>abs(A[l][i])) l=k;
          if (A[l][i]==0) return;

          swap(A[l],A[i]);
          swap(b[l],b[i]);

           b[i]/=A[i][i];

          for(  int j=N-1 ;j>i;A[i][j--]/=A[i][i]);
          A[i][i]=1;

          for(int j=i+1; j<N;j++){
            for(k=N-1;k>i;k--)
                A[j][k]-=A[i][k]*A[j][i];
              b[j]-=b[i]*A[j][i];
            A[j][i]=0;
          }
       cout<<endl;

        }//прямой ход
       x[N-1]=b[N-1];
       for ( int i = N - 2; i >= 0; i-- )
         {
           x[i] = b[i];
           for ( int j = i + 1; j < N; j++ ) {
              x[i] -= A[i][j] * x[j];
              }
                  cout<<endl;

         }//обратный ход



}



/**
 * Метод прогонки
 */
void syusinaev::lab3()
{
	
}








/**
 * Метод простых итераций
 */
void syusinaev::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void syusinaev::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void syusinaev::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void syusinaev::lab7()
{

}


void syusinaev::lab8()
{

}
<<<<<<< HEAD
<<<<<<< HEAD
/*
void syusinaev::lab9()
{

}
*/


void syusinaev::lab9()
{

}


std::string syusinaev::get_name()
{
  return "Syusina E.V.";
}
