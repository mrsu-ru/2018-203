#include "borisovdv.h" 

/** 
* Введение в дисциплину 
*/ 
void borisovdv::lab1() 
{ 
std::cout<<"Hello, world!"; 
} 


/** 
* Метод Гаусса с выбором главного элемента 
*/ 
void borisovdv::lab2() 
{ 
double maxi; 
int k, ii; 
const double eps = 0.000000000000000000000000000000001; 

k = 0; 
while (k < N) 
{ 
maxi = abs(A[k][k]); 
ii = k; 
for (int i = k + 1; i < N; i++) { 
if (abs(A[i][k]) > maxi) { 
maxi = abs(A[i][k]); 
ii = i; 
} 
} 
for (int j = 0; j < N; j++) { 
swap (A[k][j], A[ii][j]); 
} 
swap(b[k], b[ii]); 
for (int i = k; i < N; i++) { 
double c = A[i][k]; 
if (abs(c) < eps) 
continue; 
for (int j = 0; j < N; j++) { 
A[i][j] = A[i][j] / c; 
} 
b[i] = b[i] / c; 
if (i == k) 
continue; 
for (int j = 0; j < N; j++) { 
A[i][j] = A[i][j] - A[k][j]; 
} 
b[i] = b[i] - b[k]; 
} 
k++; 
} 
for (k = N - 1; k >=0; k--) { 
x[k] = b[k]; 
for (int i = k-1; i >=0; i--) { 
double c=A[i][k]; 
for (int j=0; j<N; j++) 
A[i][j]=A[k][j]*c+A[i][j]; 
b[i]=-b[k]*c+b[i]; 
} 
} 
} 



/** 
* Метод прогонки 
*/ 
void borisovdv::lab3() 
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
void borisovdv::lab4() 
{ 
double eps = 1e-13;
double tau = 1e-5;
for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double x1;
	double *xn = new double[N];
	int step = 0;
	do {
		step++;
		for (int i = 0; i < N; i++) {
			xn[i] = x[i];
			for (int k = 0; k < N; k++)
				xn[i] -= tau*A[i][k] * x[k];
			xn[i] += tau * b[i];
		}
		x1 = 0.;
		for (int i = 0; i < N; i++) {
			x1 += (xn[i]-x[i])*(xn[i]-x[i]);
		}
		for (int i = 0; i < N; i++) {
			x[i] = xn[i];
		}
	} while (sqrt(x1)>eps);
} 



/** 
* Метод Якоби или Зейделя 
*/ 
void borisovdv::lab5() 
{ 
  const double eps = 10E-20;
	double* y = new double[N];
	double r = 0; 
	for(int i=0; i<N; i++){
		x[i] = 0; }
	do {
		for(int i=0; i<N; i++){
			y[i] = b[i];
			for(int j=0; j<N; j++){
				if(i != j){
					y[i] -= A[i][j]*x[j];
				}
			}
			y[i] /= A[i][i];
		}
		r = abs(x[0] - y[0]);
		for(int i=0; i<N; i++){
			if(abs(x[i]-y[i]) > r)
				r = sqrt((x[i]-y[i])*(x[i]-y[i]));
			x[i] = y[i];
		}
	} while(r >= eps);
	delete[] y;
} 

/** 
* Метод минимальных невязок 
*/ 
void borisovdv::lab6() 
{ 
double Eps = 1e-18;
	double Del, Res, Abs;

	double *K = new double[N];
	double *L = new double[N];
	double *xrez = new double[N];
	
	
	
	for (int i = 0; i<N; i++)
		xrez[i] = 0;

	
	do{
		
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * xrez[j];
		}

		
		for (int i = 0; i < N; i++) {
			L[i] = K[i] - b[i];
		}

		
	
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * L[j];
		}

		Res = 0;
		Abs = 0;
		
		
		for (int i = 0; i < N; i++) {
			Res += K[i] * L[i];
			Abs += K[i] * K[i];
		}
		
		if (Res==Abs) Res=1;
		else {
		Res = Res / Abs;
		}
		
		for (int i = 0; i < N; i++)
			x[i] = xrez[i] - Res*L[i];
		
		
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
void borisovdv::lab7() 
{ 
double eps = 1e-15;
	double delta, r, rModul;
	double *w = new double[N];
	double *wp = new double[N];
	double *v = new double[N];
	double *result = new double[N];

	for (int i = 0; i<N; i++)
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
		
		r = 0;
		rModul = 0;
		
		for (int i = 0; i < N; i++) {
			r += w[i] * v[i];
			rModul += wp[i] * w[i];
		}
		
		if (r == rModul)r = 1;
		else r = r / rModul;
		for (int i = 0; i < N; i++)
			x[i] = result[i] - r*v[i];
		delta = abs(x[0] - result[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - result[i])>delta)
				delta = abs(x[i] - result[i]);
			result[i] = x[i];
		}
	} 
	while (eps < delta);
} 


void borisovdv::lab8() 
{ 
double * Y = new double[N];//предыдущее приближение
	double * y = new double[N];//последующее приближение
	double maxSob,sob,sum;
	double eps = 1e-9;
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
		for (int i = 0; i<N; i++)
			Y[i] = y[i];
	} 
	while (abs(maxSob - sob)>eps);
	cout << maxSob << endl;
} 

void borisovdv::lab9() 

{
	
}

static double f(double x)
{
 return (x*x-2);
}

static double f2(double x)
{
 return (2*x);

} 
void borisovdv::lab10() 
{ 
double y; 
for(int i=0; i < 100; i++) { 
double xd; 
double eps = 1e-5; 
y = i; 
int ind = 0; 
do { 
ind++; 
xd = y; 
y = exp((-y)); 
} while (abs(xd - y) > eps || ind > 1000); 
if (y == y) break; 
} 
cout << y << endl; 
} 

std::string borisovdv::get_name() 
{ 
return "Borisov D.V."; 
}