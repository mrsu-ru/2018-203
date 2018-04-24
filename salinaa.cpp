#include "salinaa.h"

/**
 * Введение в дисциплину
 */
 
 double static f(double x)
{
    double f=x*x-2*x+1;
	return f;
}

double static df(double x)
{
    double df=2*x-2;
	return df;
}
 
void salinaa::lab1()
{
int iter = 0;
  double x=0, eps=1e-5;
  do {
  x = x - f(x) / df(x);
  iter++;
  } while (fabs(f(x)) > eps && iter<20000);
  cout << iter << " iterations " << endl;
  cout << x;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void salinaa::lab2()
{
double *x, max;
	int k, index;
	const double eps = 0.00001;
	x = new double[N];
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
			return 0;
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
	return x;
}



/**
 * Метод прогонки
 */
void salinaa::lab3()
{
double *x;
	x = new double[N];
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
	return x;
}



/**
 * Метод простых итераций
 */
void salinaa::lab4()
{
double eps = 1e-5;
    
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
	} while (abs(x[0] - xx)>eps);
	
}



/**
 * Метод Якоби или Зейделя
 */
void salinaa::lab5()
{
double eps = 1e-5;
	
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
	} while (norm > eps);
	delete[] results;
}



/**
 * Метод минимальных невязок
 */
void salinaa::lab6()
{
double eps = 1e-5;
    double norm, tao, taoMod;
	double* result = new double[N];
	double *Ark = new double[N];
	double *rk = new double[N];

	for (int i = 0; i < N; i++) {
		result[i] = 0;
	}
	do {
		for (int i = 0; i < N; i++) {
			Ark[i] = 0;
			for (int j = 0; j < N; j++)
					Ark[i] += A[i][j] * result[j];
		}

		for (int i = 0; i < N; i++) {
					rk[i] = Ark[i]-b[i];
		}

		for (int i = 0; i < N; i++) {
			Ark[i] = 0;
			for (int j = 0; j < N; j++)
				Ark[i] += A[i][j] * rk[j];
		}

        tao = 0;
		taoMod = 0;
		for (int i = 0; i < N; i++) {
			tao += Ark[i] * rk[i];
			taoMod += Ark[i] * Ark[i];
		}
		if (tao==taoMod) {
                tao=1;
                }
		else {
		    tao = tao / taoMod;
        }

		for (int i = 0; i < N; i++)
			x[i] = result[i] - tao*rk[i];
		norm = abs(x[0] - result[0]);

		for (int i = 0; i < N; i++) {
			if (abs(x[i] - result[i])>norm)
				norm = abs(x[i] - result[i]);
			result[i] = x[i];
		}
	} while (eps < norm);
}



/**
 * Метод сопряженных градиентов
 */
void salinaa::lab7()
{
double eps = 1e-15;
	double norm, s, sMod;
   
	double *d = new double[N];
	double *g = new double[N];
	double *dd = new double[N];
	double *result = new double[N];



	for (int i = 0; i<N; i++){
		result[i] = 0;
	}


	do {
		for (int i = 0; i < N; i++) {
			d[i] = 0;
			for (int j = 0; j < N; j++)
				d[i] += A[i][j] * result[j];
		}

		for (int i = 0; i < N; i++) {
			g[i] = d[i] - b[i];
		}

		for (int i = 0; i < N; i++) {
			d[i] = 0;
			for (int j = 0; j < N; j++)
				d[i] += A[i][j] * g[j];
		}

		for (int i = 0; i < N; i++) {
			dd[i] = 0;
			for (int j = 0; j < N; j++) {
				dd[i] += A[i][j] * d[j];
			}
		}

		s = 0;
		sMod = 0;
		for (int i = 0; i < N; i++) {
			s += d[i] * g[i];
			sMod += dd[i] * d[i];
		}
		if (s == sMod)
			s = 1;
		else
			s = s / sMod;

		for (int i = 0; i < N; i++)
			x[i] = result[i] - s*g[i];


		norm = abs(x[0] - result[0]);

		for (int i = 0; i < N; i++) {
			if (abs(x[i] - result[i])>norm)
				norm = abs(x[i] - result[i]);
				result[i] = x[i];
		}
	} while (eps < norm);
	
}


void salinaa::lab8()
{

}

void salinaa::lab9()
{

}


std::string salinaa::get_name()
{
  return "Salin A.A.";
}
