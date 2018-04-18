#include "biryukovaes.h"

/**
 * Введение в дисциплину
 */
void biryukovaes::lab1()
{
std::cout<<"Hello, world!";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void biryukovaes::lab2()
{
int i, j, m, l;
    for(i=0;i<N;i++)
        {
          m=i;
          for(l=i+1; l<N; l++)
          if (abs(A[l][i])>abs(A[m][i]))
                m=l;
          swap(A[m],A[i]);
          swap(b[m],b[i]);
           b[i]/=A[i][i];
          for(j=N-1; j>i; A[i][j--]/=A[i][i]);
          A[i][i]=1;
          for(int j=i+1; j<N;j++){
          for(l=N-1;l>i;l--)
          A[j][l]-=A[i][l]*A[j][i];
          b[j]-=b[i]*A[j][i];
          A[j][i]=0;
          }
            cout<<endl;
          }
          x[N-1]=b[N-1];
          for ( int i = N - 2; i >= 0; i-- )
         {
           x[i] = b[i];
           for ( int j = i + 1; j < N; j++ ) {
              x[i] -= A[i][j] * x[j];
              }
          }
}



/**
 * Метод прогонки
 */
void biryukovaes::lab3()
{
 double* a=new double[N];//Прогоночные коэффициеты
 double* b=new double[N];//
    //Прямая прогонка
    a[0]=-A[0][1]/A[0][0];
    b[0]=b[0]/A[0][0];
   for(int i=1;i<N;i++)
    {
     a[i]=-A[i][i+1]/(A[i][i-1]*a[i-1]+A[i][i]);
     b[i]=(b[i]-A[i][i-1]*b[i-1])/(A[i][i-1]*a[i-1]+A[i][i]);
    }
	//Обратная прогонка
    x[N-1] = b[N-1];
    for(int i=N-2; i>=0; i--)
     x[i] = a[i]*x[i+1]+b[i];
      delete[] a;
      delete[] b;
}



/**
 * Метод простых итераций
 */
void biryukovaes::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void biryukovaes::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void biryukovaes::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void biryukovaes::lab7()
{

}


void biryukovaes::lab8()
{

}

void biryukovaes::lab9()
{

}


std::string biryukovaes::get_name()
{
  return "Biryukova E.S.";
}
