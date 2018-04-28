#include "borisovrs.h"

/**
 * Введение в дисциплину
 */
void borisovrs::lab1()
{
    std::cout<<"hello world";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void borisovrs::lab2()
{
 double t;

for (int k = 0; k < N; k++) { 
	int maxel=k;
		for(int i=k+1;i<N;i++)
			if(abs(A[i][k]) > abs(A[maxel][k]))
				maxel=i;//ищем наибольший элемент и запоминаем номер строки
	for(int i=0;i<N;i++)
	std::swap(A[k][i],A[maxel][i]);//меняем Kю строку со строкой на которой максимальный элемент
	std::swap(b[k],b[maxel]);
	t = A[k][k];
	for (int j = 0; j < N; j++)
		A[k][j] = A[k][j] / t;//получаем 1 на главной диагонали
	b[k] =b[k]/t;

	for (int i = k + 1; i < N; i++) {
		t = A[i][k];
		for (int j = 0; j < N; j++) {
			A[i][j] =A[i][j]- A[k][j] * t;//Придаем матрице треугольный
		}
	b[i] =b[i]- b[k] * t;
	}
}

//получаем матрицу треугольниго вида и осуществляем обратный ход 
for (int k = N - 1; k > 0; k--){
	for (int i = k - 1; i >= 0; i--){
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
void borisovrs::lab3()
{
double *a = new double[N];
double *beta = new double[N];
a[0] = -A[0][1] / A[0][0];
beta[0] = b[0] /A[0][0];
for (int i = 1; i < N; i++) {
	double y = A[i][i] + A[i][i - 1] * a[i - 1];
	a[i] = -A[i][i + 1] / y;
	beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
}
x[N-1] = (b[N-1] - A[N-1][N-2] * beta[N-2]) / (A[N-1][N-1] + A[N-1][N-1] * a[N-1]);
for (int i = N-1 - 1; i >= 0; i--) {
		x[i] = a[i] * x[i + 1] + beta[i];
}

/*
double* Alpha = new double[N];
double* Betta = new double[N];

    Alpha[0] = -A[0][1]/A[0][0];
    Betta[0] = b[0]/A[0][0];

    for(int i=1; i<N; i++)
    {
        Alpha[i] = -A[i][i+1]/(A[i][i] + A[i][i - 1] * Alpha[i - 1]);
        Betta[i] = (b[i] - A[i][i-1]*Betta[i-1])/(A[i][i] + A[i][i - 1] * Alpha[i - 1]);
    }
    x[N-1] = (b[N-1] - A[N-1][N - 2] * Betta[N - 2]) / (A[N-1][N-1] + A[N-1][N-1] * Alpha[N-1]);
    for(int i=N-2; i>=0; i--)

    x[i] = Alpha[i]*x[i+1]+Betta[i];
    delete[] Alpha;
    delete[] Betta;*/
}



/**
 * Метод простых итераций
 */
void borisovrs::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void borisovrs::lab5()
{
//Метод Якоби
	double *oldx = new double[N]; 
	for (int i=0; i<N; i++) { 
		x[i]=0;} // заполняем решение нулями 
	double Err=0.0; 
	double eps=1e-20; // погрешность
	int k=0; 
	do { 
		k++; 
		Err=0.0; 
		for(int i=0; i<N; i++) 
			oldx[i]=x[i]; // предыдущее решение 
			for(int i=0; i<N; i++) 
			{ 
				double s=0;
				for(int j=0; j<i; j++) 
					s += A[i][j] * oldx[j]; //под главной диагональю
				for(int j=i+1; j<N; j++) 
					s += A[i][j] * oldx[j]; //над главной диагональю
				x[i]=(b[i] - s)/A[i][i]; // вычисляется новое решение 
			}			 
			Err=std::abs(oldx[0]-x[0]); 
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
void borisovrs::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void borisovrs::lab7()
{

}


void borisovrs::lab8()
{

}


void borisovrs::lab9()
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
void borisovrs::lab10()
{
	double A=1;
	double B=2;
	double eps=0.001,t;
	do {
		if( f(A)*f2(A)>0){
			B = B - ((A-B) * f(B))/(f(A) - f(B));
			t=B;
		}
		else if (f(B)*f2(B)>0) {
			A = A - ((B-A) * f(A))/(f(B) - f(A));
			t=A;
		}
	} while (fabs(f(t))>=eps);
	cout<<"x = "<<t<<"\n";
}

std::string borisovrs::get_name()
{
  return "Borisov R.S.";
}
