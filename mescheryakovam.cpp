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
	//Адаптация метода Гаусса,в случае трехдиагональной матрицы
	double *upper, *middle, *lower; // верхняя, главная и нижняя диагонали 
	upper = new double[N]; 
	middle = new double[N]; 
	lower = new double[N]; 
	double k; 

	lower[0] = 0; 
	upper[N-1] = 0; 

	for (int i = 0; i < N; i++)//Заполнение"диагональных" массивов 
	{ 
		if (i - 1 >= 0 && i - 1 < N) 
		upper[i] = A[i-1][i];  
		middle[i] = A[i][i];  
		if (i + 1 >= 0 && i + 1 < N) 
		lower[i] = A[i+1][i]; //нижняя 
	} 
	//Прямой ход
	for (int i = 1; i < N; i++) //Вычисляем коэффициенты прогонки 
	{ 
		k = lower[i]/middle[i-1]; 
		middle[i] = middle[i] - k*upper[i-1]; 
		b[i] = b[i] - k*b[i-1]; 
	}	 

	//Обратный ход
	x[N-1] = b[N-1]/middle[N-1]; //Вычисляется решение 

	for (int i = N - 2; i >= 0; i--) //рекуррентная формула для вычисления остальных неизвестных
		x[i]=(b[i]-upper[i]*x[i+1])/middle[i]; 

	delete[] upper, middle, lower;
} 



/**
 * Метод простых итераций
 */
void mescheryakovam::lab4()
{
	double mis=1e-15;//порядок ошибки
	double error;
	double *nx = new double[N];//для хранения промежуточных значений

	for (int i=0;i<N;i++)//для первичного приближения возьмём столбец свободных членов
		x[i]=b[i];
	int step=0;
	
	do
	{
		step++;
	for(int i=0;i < N;i++)
	{
		nx[i]=-b[i];		
		for(int j=0;j < N;j++)
		{
			if(i!=j)
			nx[i]+=A[i][j]*x[j];
		}
		nx[i]/=-A[i][i];
	}
	error=0;
	for(int i=0; i<N; i++) 
	{ 
		if(std::abs(x[i]-nx[i]) > error)//Максимальная разница между элементами решения 
			error = std::abs(x[i]-nx[i]);
	}
	for(int i=0; i<N; i++) 
		x[i]=nx[i];
	std::cout<<step<<"    "<<error<<endl;
	}
	while (error>mis);
	delete[] nx;
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
