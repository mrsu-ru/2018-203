#include "bagapovar.h"

/**
 * Введение в дисциплину
 */

 void bagapovar::lab1()
{
std::cout<<"Hello world";
}


/**
* Метод Гаусса с выбором главного элемента
*/
void bagapovar::lab2()
{
double t;

for (int k = 0; k < N; k++) { //Ищем самый большой элемент стоящий на к'атом месте по всем строкам, по умолчанию он k'ый
	int maxel=k;
		for(int i=k+1;i<N;i++)
			if(abs(A[i][k]) > abs(A[maxel][k]))
				maxel=i;
//Нашли, теперь надо поменять k'ую строчку и строчку с максимальным элементом местами
	for(int i=0;i<N;i++)
	std::swap(A[k][i],A[maxel][i]);
	std::swap(b[k],b[maxel]);

//Прямой ход
t = A[k][k];
//Делим все элементы k'ой строки на a[k][k] элемент, потому что он диагональный и мы хотим на нём получить значение 1
//для удобства дальнейших вычислений
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

//Матрица треугольного вида готова, теперь можем последовательно, снизу вверх вычислять искомые значения элементов матрицы x
//Осуществляем обратный ход
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
void bagapovar::lab3()
{//Это частный случай для метода Гаусса, используется, когда матрица трёхдиагональная
double *up, *mid, *low; // верхняя, главная, нижняя диагонали 
up = new double[N]; 
mid = new double[N]; 
low = new double[N]; 
double k; 

low[0] = 0; 
up[N-1] = 0; 

for (int i = 0; i < N; i++)//Заполняем "диагональные" массивы 
{ 
	if (i - 1 >= 0 && i - 1 < N) 
	up[i] = A[i-1][i]; //верхняя 
	mid[i] = A[i][i]; //главная 
	if (i + 1 >= 0 && i + 1 < N) 
	low[i] = A[i+1][i]; //нижняя 
} 
//Прямая прогонка 
for (int i = 1; i < N; i++) //Вычисляем коэффициенты прогонки 
{ 
	k = low[i]/mid[i-1]; 
	mid[i] = mid[i] - k*up[i-1]; 
	b[i] = b[i] - k*b[i-1]; 
} 

//Обратная прогонка 
x[N-1] = b[N-1]/mid[N-1]; //Вычисляется решение 

for (int i = N - 2; i >= 0; i--) //Решение для следующих частей 
	x[i]=(b[i]-up[i]*x[i+1])/mid[i]; 

delete[] up, mid, low; 
} 



/**
 * Метод простых итераций
 */
void bagapovar::lab4()
{
double Eps=1e-15;//порядок ошибки
double Err;
double *nx = new double[N];//для хранения промежуточных значений

for (int i=0;i<N;i++)//для первичного приближения возьмём столбец свободных членов
	x[i]=b[i];
int step=0;
	
do{
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
  Err=0;
for(int i=0; i<N; i++) { 
if(std::abs(x[i]-nx[i]) > Err)//Максимальная разница между элементами решения 
Err = std::abs(x[i]-nx[i]);
}
for(int i=0; i<N; i++) 
	x[i]=nx[i];
std::cout<<step<<"    "<<Err<<endl;
}while (Err>Eps);

delete[] nx;
}


/**
 * Метод Якоби или Зейделя
 */
void bagapovar::lab5()
{ 
//Метод Якоби
double *oldx = new double[N]; 

for (int i=0; i<N; i++) { 
x[i]=0; // первоначальное новое решение 
} 

double Err=0.0; 
double eps=1e-20; 
int k=0; 

do { 
k++; 
Err=0.0; 
for(int i=0; i<N; i++) 
oldx[i]=x[i]; // здесь записывается предыдущее решение 
for(int i=0; i<N; i++) 
{ 
double s=0; //вычисляем s, но мы не берём диагональные элементы 
for(int j=0; j<i; j++) 
s += A[i][j] * oldx[j]; 
for(int j=i+1; j<N; j++) 
s += A[i][j] * oldx[j]; 
x[i]=(b[i] - s)/A[i][i]; // вычисляется новое решение 
} 
Err= std::abs(oldx[0]-x[0]); 
for(int i=0; i<N; i++) 
{ 
if(std::abs(oldx[i]-x[i]) > Err) 
Err = std::abs(oldx[i]-x[i]);//максимальная разница между предыдущим решением и текущим. 
} 
} while(Err >= eps); 
delete [] oldx; 
} 



/**
 * Метод минимальных невязок
 */
void bagapovar::lab6()
{
	
	double Eps = 1e-18;//заданная погрешность
	double Del, Res, Abs;//погрешность, невязка, модуль

	double *K = new double[N];//w
	double *L = new double[N];//v
	double *xrez = new double[N];
	
	
	//задаём первоначальное приближение
	for (int i = 0; i<N; i++)
		xrez[i] = 0;

	//цикл для нахождения корней
	do{
		//находим редуцированную систему(одна часть)
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * xrez[j];
		}

		//находим редуцированную систему(вторая часть)
		for (int i = 0; i < N; i++) {
			L[i] = K[i] - b[i];//нахождение вектора невязки
		}

		
		//нахождение скалярного произведения матрицы системы и вектора невязки
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * L[j];
		}

		Res = 0;
		Abs = 0;
		
		//нахождение значения итерационного параметра
		for (int i = 0; i < N; i++) {
			Res += K[i] * L[i];
			Abs += K[i] * K[i];
		}
		
		if (Res==Abs) Res=1;
		else {
		Res = Res / Abs;
		}
		//получение приближения решения
		for (int i = 0; i < N; i++)
			x[i] = xrez[i] - Res*L[i];
		
		//Проверка на уменьшение погрешности
		Del = abs(x[0] - xrez[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - xrez[i])>Del)
				Del = abs(x[i] - xrez[i]);
			xrez[i] = x[i];
		}
	} while (Eps < Del);
}



/**
 * Метод сопряженных градиентов
 */
void bagapovar::lab7()
{
	double Eps = 1e-20;//заданная погрешность
	double Del, s, sAbs;//погрешность итерации, скалярный шаг, модуль шага


	double *K = new double[N];
	double *L = new double[N];
	double *M = new double[N];
	double *xrez = new double[N];//итерационные решения
	
	
	//задание начального приближения 
	for (int i = 0; i<N; i++){
		xrez[i] = 0;
	}
	
	
	do {
		//нахождение скалярного произведения матрицы системы и вектора приближенного решения
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * xrez[j];
		}
		
		//нахождение градиента
		for (int i = 0; i < N; i++) {
			L[i] = K[i] - b[i];
		}
		
		//скалярное произведение матрицы системы и градиента
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * L[j];
		}
		
		
		for (int i = 0; i < N; i++) {
			M[i] = 0;
			for (int j = 0; j < N; j++) {
				M[i] += A[i][j] * K[j];
			}
		}
		
		s = 0;
		sAbs = 0;
		
		//нахождение величины смещения по направлению градиента(скалярного шага)
		for (int i = 0; i < N; i++) {
			s += K[i] * L[i];
			sAbs += M[i] * K[i];
		}
		if (s == sAbs)
			s = 1;
		else 
			s = s / sAbs;
		//записываем новое приближенное решение
		for (int i = 0; i < N; i++)
			x[i] = xrez[i] - s*L[i];
		
		//проверка на уменьшение погрешности
		Del = abs(x[0] - xrez[0]);
		
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - xrez[i])>Del)
				Del = abs(x[i] - xrez[i]);
				xrez[i] = x[i];
		}
	} while (Eps < Del);
}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */

void bagapovar::lab8()
{

}

/**
 * Нахождение наибольшего по модолю собственного значения матрицы
 */
void bagapovar::lab9()
{

}


std::string bagapovar::get_name()
{
  return "Bagapov A.R.";
}
