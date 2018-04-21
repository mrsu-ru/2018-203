#include "fedyanovaam.h"

/**
 * Введение в дисциплину
 */
void fedyanovaam::lab1()
{
cout<<"Hello!!!";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
<<<<<<< HEAD
void fedyanovaam::lab2()
{
 double p;
=======
void fedyanovaam::lab2(){
  double p;
>>>>>>> fa39b47249585ff972fc100c09bc1da99c5ec339
    int maxn;

      for (int k=0; k<N-1; k++){
          maxn = k;
          for (int i=k+1; i<N; i++)
        if(abs(A[i][k]) > abs(A[maxn][k])) maxn = i; ///Выбор главного элемента
          std::swap(A[maxn], A[k]);///Меняем строки местами
          std::swap(b[maxn], b[k]);

          for (int i=k+1; i<N; i++){
              p = A[i][k]/A[k][k];
              for (int j=k; j<N; j++)
                  A[i][j] -= p*A[k][j];
              b[i] -= p*b[k];
          }
      }

      for(int i = 0; i<N; i++){
          x[i]=b[i];
      }

      for (int i=N-1; i>=0; i--){//Перебираем массив
          for (int j=i+1;j<N;j++)
              x[i]-=A[i][j]*x[j];
          x[i] /= A[i][i];
<<<<<<< HEAD
      }  
=======
      }
>>>>>>> fa39b47249585ff972fc100c09bc1da99c5ec339
}




/**
 * Метод прогонки
 */
<<<<<<< HEAD
void fedyanovaam::lab3()
{
 double *P = new double [N]; ///Коэффициенты alfa
=======
void fedyanovaam::lab3(){
  double *P = new double [N]; ///Коэффициенты alfa
>>>>>>> fa39b47249585ff972fc100c09bc1da99c5ec339
  double *Q = new double [N]; ///Коэффициенты betta

    P[0] = -A[0][1]/A[0][0];
    Q[0] = b[0]/A[0][0];

    for(int i=1; i<N; i++){ ///Определим прогоночные коэффициенты
    
      P[i] = A[i][i+1]/(-A[i][i] - A[i][i-1]*P[i-1]);
      Q[i] = (-b[i] + A[i][i-1]*Q[i-1])/(-A[i][i] - A[i][i-1]*P[i-1]);
    }
<<<<<<< HEAD

    x[N-1] = Q[N-1];
    for(int i=N-2; i>=0; i--) ///Определим решение
      x[i] = P[i]*x[i+1] + Q[i];

=======

    x[N-1] = Q[N-1];
    for(int i=N-2; i>=0; i--) ///Определим решение
      x[i] = P[i]*x[i+1] + Q[i];

>>>>>>> fa39b47249585ff972fc100c09bc1da99c5ec339
    delete [] P;
    delete [] Q;
}




/**
 * Метод простых итераций
 */
void fedyanovaam::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void fedyanovaam::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void fedyanovaam::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void fedyanovaam::lab7()
{

}


void fedyanovaam::lab8()
{

}
void fedyanovaam::lab9()
{

}

std::string fedyanovaam::get_name()
{
  return "Fedyanova A.M.";
}
