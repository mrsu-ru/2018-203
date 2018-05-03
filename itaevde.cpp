#include "itaevde.h"

/**
 * Введение в дисциплину
 */
void itaevde::lab1()
{
double x0;
	for(int i=0;i<100;i++){
	double xd;
	double eps = 1e-5;
	x0 = i;
	int ind = 0;
	do {
		ind++;
		xd = x0;
		x0 = exp((-x0));
	} while (abs(xd - x0) > eps || ind>1000);
	if (x0==x0) break;
	}
	cout << x0 << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void itaevde::lab2()
{
double max;
	int k, index;
	const double eps = 0.00001;
	k = 0;
	while (k < N){
		max = abs(A[k][k]);
		index = k;
		for (int i = k + 1; i < N; i++){
			if (abs(A[i][k]) > max){
				max = abs(A[i][k]);
				index = i;
			}
		}
		
		if (max < eps){
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return;
		}
		for (int j = 0; j < N; j++)
		    swap(A[k][j], A[index][j]);
		swap(b[k],b[index]);

		for (int i = k; i < N; i++){
			double temp = A[i][k];
			if (abs(temp) < eps) 
				continue; 
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] / temp;
			b[i] = b[i] / temp;
			if (i == k)
				continue; 
			for (int j = 0; j < N; j++)
				A[i][j] = A[i][j] - A[k][j];
			b[i] = b[i] - b[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = N - 1; k >= 0; k--){
		x[k] = b[k];
		for (int i = 0; i < k; i++)
			b[i] = b[i] - A[i][k] * x[k];
	}
}



/**
 * Метод прогонки
 */
void itaevde::lab3()
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
void itaevde::lab4()
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
void itaevde::lab5()
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
		for (int j = 0; j < N; j++)
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
	} while (abs(x[0] - xx) > eps);
}



/**
 * Метод минимальных невязок
 */
void itaevde::lab6()
{
	double eps = 1e-5;
	double delta, r, rModul;

	double *w = new double[N];
	double *v = new double[N];
	double *result = new double[N];

	for (int i = 0; i < N; i++)
		result[i] = 0;

	do{
		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++)
				w[i] += A[i][j] * result[j];
		}

		for (int i = 0; i < N; i++) {
			v[i] = w[i] - b[i];
		}

		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++)
				w[i] += A[i][j] * v[j];
		}

		r = 0.0;
		rModul = 0.0;
		for (int i = 0; i < N; i++) {
			r += w[i] * v[i];
			rModul += w[i] * w[i];
		}
		if (r == rModul) {r = 1;}
		else {r = r / rModul;}
		for (int i = 0; i < N; i++)
			x[i] = result[i] - r*v[i];
		delta = abs(x[0] - result[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - result[i]) > delta)
				delta = abs(x[i] - result[i]);
			result[i] = x[i];
		}
	}  while (eps < delta);
}



/**
 * Метод сопряженных градиентов
 */
void itaevde::lab7()
{
	double eps = 1e-5;
	double delta, r, rModul;
	
	double *w = new double[N];
	double *wp = new double[N];
	double *v = new double[N];
	double *result = new double[N];

	for (int i = 0; i < N; i++)
		result[i] = 0;

	do {
		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++)
				w[i] += A[i][j] * result[j];
		}

		for (int i = 0; i < N; i++) {
			v[i] = w[i] - b[i];
		}

		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++)
				w[i] += A[i][j] * v[j];
		}
		
		for (int i = 0; i < N; i++) {
			wp[i] = 0;
			for (int j = 0; j < N; j++) {
				wp[i] += A[i][j] * w[j];
			}
		}
		
		r = 0.0;
		rModul = 0.0;
		
		for (int i = 0; i < N; i++) {
			r += w[i] * v[i];
			rModul += wp[i] * w[i];
		}
		
		if (r == rModul) r = 1;
		else r = r / rModul;
		for (int i = 0; i < N; i++)
			x[i] = result[i] - r*v[i];
		delta = abs(x[0] - result[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - result[i]) > delta)
				delta = abs(x[i] - result[i]);
			result[i] = x[i];
		}
	} while (eps < delta);
}


void itaevde::lab8()
{
	double *sob = new double[N];
	double * b = new double[N];
	double * z = new double[N];
	
	for (int i = 0; i < N; i++){
		z[i] = 0.0;
		b[i] = A[i][i];
		sob[i] = A[i][i];
	}
	
	for (int i = 0; i < 100; i++){
		double sm = 0.0;
		for (int p = 0; p < N - 1; p++){
			for (int q = p + 1; q < N; q++){
				sm += abs(A[p][q]);
			}
		}
		if (sm == 0) break;
		double tresh = sm / (5*N*N);
		for (int p = 0; p < N - 1; p++){
			for (int q = p + 1; q < N; q++){
				double g = 1e12 * abs(A[p][q]);
				if (i >= 3 && abs(sob[p]) > g && abs(sob[q]) > g) A[p][q] = 0.;
				else
					if (abs(A[p][q]) > tresh){
						double theta = (sob[q] - sob[p]) / (2.0 * A[p][q]);
						double t = 1.0 / (abs(theta) + sqrt(1.0 + theta*theta));
						if (theta < 0) t = -t;
					    double s = t / sqrt(1.0 + t*t);
						double tau = s / (1.0 + 1.0 / sqrt(1.0 + t*t));
						z[p] -= t * A[p][q];
						z[q] += t * A[p][q];
						sob[p] -= t * A[p][q];
						sob[q] += t * A[p][q];
						A[p][q] = 0.0;
						for (int j = 0; j < p; j++){
							A[j][p] = A[j][p] - s * (A[j][q] + A[j][p] * tau);
							A[j][q] = A[j][q] + s * (A[j][p] - A[j][q] * tau);
						}
						for (int j = p + 1; j < q; j++){
							A[p][j] = A[p][j] - s * (A[j][q] + A[p][j] * tau);
							A[j][q] = A[j][q] + s * (A[p][j] - A[j][q] * tau);
						}
						for (int j = q + 1; j < N; j++){
							A[p][j] = A[p][j] - s * (A[q][j] + A[p][j] * tau);
							A[q][j] = A[q][j] + s * (A[p][j] - A[q][j] * tau);
						}
					}
			}
		}
		for (int p = 0; p < N; p++){
			sob[p] = b[p] + z[p];
		}
	}

	for (int p = 0; p < N; ++p) {
		cout << sob[p] << endl;
	}

}


void itaevde::lab9()
{
	double * Y = new double[N];
	double * y = new double[N];
	double maxSob,sob,sum;
	double eps = 1e-5;
	
	for (int i = 0; i < N; i++)
		Y[i] = 0;
	
	Y[0] = 1;
	
	do{
		sum = 0;
		for (int i = 0; i < N; i++)
			sum += Y[i] * Y[i];
		sob = sqrt(sum);
		for (int i = 0; i < N; i++)
		{
			y[i] = 0;
			for (int j = 0; j < N; j++)
				y[i] += A[i][j] * Y[j] / sob;
		}
		sum = 0;
		for (int i = 0; i < N; i++)
			sum += y[i] * y[i];
		maxSob = sqrt(sum);
		for (int i = 0; i < N; i++)
			Y[i] = y[i];
	} while (abs(maxSob - sob) > eps);

	cout << maxSob << endl;
}


std::string ivanovii::get_name()
{
  return "Itaev D.E.";
}
