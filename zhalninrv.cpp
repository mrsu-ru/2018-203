#include "zhalninrv.h"
#include <cstdio>
/**
 * Введение в дисциплину
 */
void zhalninrv::lab1()
{
  std::cout << "Hello World!!!" << std::endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void zhalninrv::lab2()
{
    for (int k = 0; k < N-1; k++) {
        for (int i = k+1; i < N; i++) {
		  // choose of main element
		  // ...
		  //
          double c_ki = A[i][k]/A[k][k];
          for (int j = k; j < N; j++) {
            A[i][j] -= A[k][j]*c_ki;
          }
          b[i] -= b[k]*c_ki;
        }
    }

    for (int k = N-1; k >= 0; k--) {
      x[k] = b[k];
      for (int i = k+1; i < N; i++) {
        x[k] -= A[k][i]*x[i];
      }
      x[k] /= A[k][k];
    }

}



/**
 * Метод прогонки
 */
void zhalninrv::lab3()
{

}



/**
 * Метод простых итераций
 */
void zhalninrv::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void zhalninrv::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void zhalninrv::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void zhalninrv::lab7()
{

}


void zhalninrv::lab8()
{

}

void zhalninrv::lab9()
{

}


std::string zhalninrv::get_name()
{
  return "Zhalnin R.V.";
}
