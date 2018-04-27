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
//Метод Якоби
double *oldx = new double[N]; 

for (int i=0; i<N; i++) { 
x[i]=0; // первоначальное решение 
} 

double Err=0.0; 
double eps=1e-20; 
int k=0; 

do { 
k++; 
Err=0.0; 
for(int i=0; i<N; i++) 
oldx[i]=x[i]; //предыдущее решение сюда
for(int i=0; i<N; i++) 
{ 
double s=0; //вычисляем s, но не берём диагональные элементы 
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
Err = std::abs(oldx[i]-x[i]);//максимальная допустимая разница между предыдущим решением и текущим. 
} 
} while(Err >= eps); 
delete [] oldx; 
} 



/**
 * Метод минимальных невязок
 */
void mescheryakovam::lab6()
{
	
	double Eps = 1e-18;//погрешность
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
void mescheryakovam::lab7()
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

void mescheryakovam::lab8()
{

}

/**
 * Нахождение наибольшего по модолю собственного значения матрицы
 */
void mescheryakovam::lab9()
{
	double * Y = new double[N];//первый вектор приближения
	double * M = new double[N];//второй вектор приближения
	double maxL, L, sum;
	double EPS = 1e-15;
	
	//первичное приближение начального вектора
	for (int i = 0; i < N; i++)
		Y[i] = 0;
	Y[0] = 1;
	
	do{
		sum = 0;
		//нахождение скалярного произведения векторов приближения 
		for (int i = 0; i < N; i++)
			sum += Y[i] * Y[i];
		
		L = sqrt(sum);//норма вектора приближения
		
		//построение последовательности векторов
		for (int i = 0; i < N; i++)
		{
			M[i] = 0;
			for (int j = 0; j < N; j++)
				M[i] += A[i][j] * Y[j] / L;
		}
		sum = 0;
		
		//сравнение нормы полученного вектора с заданной погрешностью
		for (int i = 0; i < N; i++)
			sum += M[i] * M[i];
		maxL = sqrt(sum);
		
		for (int i = 0; i<N; i++)
			Y[i] = M[i];
	} while (abs(maxL - L)>EPS);
}

double static f(double x)
{
    double f=x*x*x-2*x*x-5*x+6;
	return f;
}
//производная этой функции 
double static df(double x)
{
    double df=3*x*x-4*x-5;
	return df;
}

//сжимающее отображение 
double static g(double x)
{
    return x - f(x)/df(x);
}
 



void mescheryakovam::lab10()
{
    double x;
    double eps=1e-30;
    cout<<"Enter initial root value   : ";//начальное приближение корня
	cin>>x;
    //итерации приближения
	for(double iter = 1; eps < fabs(f(x)); iter++){//double,потому что может быть очень много итераций
        system("cls");
        
		cout<<"Iteration : "<<iter<<endl;
        
		if(df(x) == 0)//если df(x) == 0, то будет деление на ноль x - f(x)/df(x)
            break;
        
		cout<<"x = "<<x<<endl;
        cout<<"g(x) = "<<g(x) <<endl;
        cout<<"df(x)= "<<df(x)<<endl;
        cout<<"f(x) = "<<f(x) <<endl;
        
		x = g(x);//новое приближение
    }
    system("pause");
}

std::string mescheryakovam::get_name()
{
  return "Mescheryakov A.M.";
}
