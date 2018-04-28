#include "kozlovdn.h"

/**
 * Введение в дисциплину
 */
void kozlovdn::lab1()
{
std::cout<<"hello world";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kozlovdn::lab2()
{
double max;
	int k, index;
	const double eps = 0.00001;
	k = 0;
	while (k < N)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(A[k][k]);
		index = k;
		for (int i = k + 1; i < N; i++)
		{
			if (abs(A[i][k]) > max)
			{
				max = abs(A[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return;
		}
		for (int j = 0; j < N; j++)
		{
		    swap(A[k][j], A[index][j]);

		}
		swap(b[k],b[index]);


		for (int i = k; i < N; i++)
		{
			double temp = A[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] / temp;
			b[i] = b[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] - A[k][j];
			b[i] = b[i] - b[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = N - 1; k >= 0; k--)
	{
		x[k] = b[k];
		for (int i = 0; i < k; i++)
			b[i] = b[i] - A[i][k] * x[k];
	}
}



/**
 * Метод прогонки
 */
void kozlovdn::lab3()
{
	int N1=N-1;
	double *alfa = new double[N];
	double *beta = new double[N];
	double y = A[0][0];
	alfa[0] = -A[0][1] / y;
	beta[0] = b[0] / y;
	for (int i = 1; i < N1; i++) {
		y = A[i][i] + A[i][i - 1] * alfa[i - 1];
		alfa[i] = -A[i][i + 1] / y;
		beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
	}

	x[N1] = (b[N1] - A[N1][N1-1] * beta[N1-1]) / (A[N1][N1] + A[N1][N1-1] * alfa[N1-1]);
	for (int i = N1-1; i >= 0; i--) {
		x[i] = alfa[i] * x[i + 1] + beta[i];
	}
}





/**
 * Метод простых итераций
 */
void kozlovdn::lab4()
{
double eps = 1e-9;
    
	for (int i = 0; i < N; i++) {
		double maxel = A[i][i];
		int indRow = i;
		for (int j = i + 1; j < N; j++)
			if (maxel < abs(A[j][i])) {
				indRow = j;
				maxel = abs(A[j][i]);
			}


		if (indRow != i) {
			for (int j = i; j < N; j++) {
				swap(A[i][j], A[indRow][j]);
			}
			swap(b[i], b[indRow]);
		}
		
		double summ = 0;
		for (int j = 0; j < N; j++)
			summ += abs(A[i][j]);
		if (2 * abs(A[i][i]) < summ) {
                cout << "Error" << std::endl;
        }

		maxel = A[i][i];
		b[i] /= A[i][i];
		A[i][i] = 0;
		for (int j = 0; j<N; j++)
			if (j != i)
                A[i][j] /= maxel;

	}
	
	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	
	double xx = b[0];
	double *results = new double[N];
	do {
		for (int i = 0; i < N; i++) {
			results[i] = 0;
			for (int k = 0; k < N; k++)
				results[i] -= A[i][k] * x[k];
			results[i] += b[i];
		}
		xx = x[0];
		for (int i = 0; i < N; i++) {
			x[i] = results[i];
		}
	} 
	while (abs(x[0] - xx)>eps);
}



/**
 * Метод Якоби или Зейделя
 */
void kozlovdn::lab5()
{

	double eps = 1e-8;
	double* results = new double[N];
	double norm;

	for (int i = 0; i < N; i++) {
		results[i] = b[i];
	}
	do {
		for (int i = 0; i < N; i++) {
			x[i] = b[i];
			for (int j = 0; j < N; j++) {
				if (i != j)
					x[i] -= A[i][j] * results[j];
			}
			x[i] /= A[i][i];
		}
        norm = fabs(results[0] - x[0]);
		for (int i = 0; i < N; i++) {
			if (fabs(results[i] - x[i]) > norm)
				norm = fabs(results[i] - x[i]);
			results[i] = x[i];
		}
	} 
	while (norm > eps);
	delete[] results;
}



/**
 * Метод минимальных невязок
 */
void kozlovdn::lab6()
{
	double eps = 1e-5;
	double delta, r, rModul; //погрешность, невязка, модуль

	double *w = new double[N];
	double *v = new double[N];
	double *result = new double[N];

	//задаём первоначальное приближение
	for (int i = 0; i<N; i++)
		result[i] = 0;
	//цикл для нахождения корней
	do{
	//находим редуцированную систему(одна часть)
		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++)
				w[i] += A[i][j] * result[j];
		}
		//находим редуцированную систему(вторая часть)
		for (int i = 0; i < N; i++) {
			v[i] = w[i] - b[i];//нахождение вектора невязки
		}
		//нахождение скалярного произведения матрицы системы и вектора невязки
		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++)
				w[i] += A[i][j] * v[j];
		}

		r = 0;
		rModul = 0;
		//нахождение значения итерационного параметра
		for (int i = 0; i < N; i++) {
			r += w[i] * v[i];
			rModul += w[i] * w[i];
		}
		if (r==rModul) {r=1;
		}
		else {r = r / rModul;
		}
		//получение приближения решения
		for (int i = 0; i < N; i++)
			x[i] = result[i] - r*v[i];
			//Проверка на уменьшение погрешности
		delta = abs(x[0] - result[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - result[i])>delta)
				delta = abs(x[i] - result[i]);
			result[i] = x[i];
		}
	} while (eps < delta);
}



/**
 * Метод сопряженных градиентов
 */
void kozlovdn::lab7()
{

}


void kozlovdn::lab8()
{

}


void kozlovdn::lab9()
{
	
}
static double f(double x)
{
 return (x*x-10);
}

static double f2(double x)
{
 return (2*x);
}

void kozlovdn::lab10()
{
	double a=1;
	double b=2;
	double eps=0.001;
	double z;
	do {
		if( f(a)*f2(a)>0){
			b = b - ((a-b) * f(b))/(f(a) - f(b));
			z=b;
		}
		else if (f(b)*f2(b)>0) {
			a = a - ((b-a) * f(a))/(f(b) - f(a));
			z=a;
		}
	} 
	while (fabs(f(z))>=eps);
	cout<<"x = "<<z<<"\n";
}


std::string kozlovdn::get_name()
{
  return "Kozlov D.N.";
}
