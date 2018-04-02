#include "mescheryakovam.h"

/**
 * Введение в дисциплину
 */

 void mescheryakovam::lab1()
{
std::cout<<"Hello world";
}


/**
* Метод Гаусса с выбором главного элемента
*/
void mescheryakovam::lab2()
{
double t;

for (int k = 0; k < N; k++) { //Ищем самый большой элемент стоящий на к'атом месте по всем строкам, по умолчанию он k'ый
	int maximal=k;
		for(int i=k+1;i<N;i++)
			if(abs(A[i][k]) > abs(A[maximal][k]))
				maximal=i;
//теперь нужно поменять k'ую строчку и строчку с максимальным элементом местами
	for(int i=0;i<N;i++)
	std::swap(A[k][i],A[maximal][i]);
	std::swap(b[k],b[maximal]);

//Прямой ход
t = A[k][k];
for (int j = 0; j < N; j++)
	A[k][j] = A[k][j] / t;
b[k] =b[k]/t;

for (int i = k + 1; i < N; i++) {
	t = A[i][k];
//Вычитаем из всех строк лежащих ниже k'ой к'ую строку помноженную на k'ый элемент строки,
// из которой вычитаем, что даёт нам ноль в этом элементе после вычитания и матрица постепенно приобретает треугольный вид

	for (int j = 0; j < N; j++) {
		A[i][j] =A[i][j]- A[k][j] * t;
	}
b[i] =b[i]- b[k] * t;
}
}

//Матрица треугольного вида найдена, можно вычислять искомые значения элементов матрицы x

//Обратный ход
for (int k = N - 1; k > 0; k--)
{
for (int i = k - 1; i >= 0; i--)
{
t = A[i][k];

for (int j = 0; j < N; j++)
A[i][j] =A[i][j]- A[k][j] * t;
b[i] =b[i] - b[k] * t;
}
}

for(int i=0; i<N; i++)
x[i]=b[i];
}



/**
 * Метод прогонки
 */
void mescheryakovam::lab3()
{
	
} 



/**
 * Метод простых итераций
 */
void mescheryakovam::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void mescheryakovam::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void mescheryakovam::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void mescheryakovam::lab7()
{

}


void mescheryakovam::lab8()
{

}


void mescheryakovam::lab9()
{

}


std::string mescheryakovam::get_name()
{
  return "Mescheryakov A.M.";
}
