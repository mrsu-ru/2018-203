#include "scherbakovdv.h"
#define MSIZE sizeof(double)*N
#include <cstdio>
/**
 * Поиск скалярного произведения
 */
double scherbakovdv::scala(double* l, double* r){
	double res=0;
	for (int i=0;i<N;i++)
		res+=l[i]*r[i];
	return res;
}

/**
 * Парочка произведений матриц и векторов
 */
double* scherbakovdv::mul(double** mat, double* vec){
	double* res = new double[N];
	memset(res,0,sizeof(double)*N);
	for (int i=0;i<N;i++)
		for (int j=0;j<N;j++)
			res[i]+=mat[i][j]*vec[j];
	return res;
}
double** scherbakovdv::mul(double** mat1, double** mat2){
	double** res = new double*[N];
	for (int i=0;i<N;i++)
	{
		res[i] = new double[N];
		for (int j=0;j<N;j++)
		{
			res[i][j]=0;
			for (int k=0;k<N;k++)
				res[i][j]+=mat1[i][k]*mat2[k][j];
		}
	}
	return res;
}

/**
 * Введение в дисциплину
 */
void scherbakovdv::lab1()
{
	printf("Hello world");
}


/**
 * Метод Гаусса с выбором главного элемента - done and passed
 */
void scherbakovdv::lab2()
{
	//A[i][i] - матрица, b[i] - столбец свободных членов
	for (int i = 0; i < N; i++)
		{
			if (!A[i][i])
			{
				double* MaxRow = A[i];
				for (int j = i + 1; j < N; j++)
					if (A[j][i] > MaxRow[i])
					{
						double* tmp = MaxRow;
						MaxRow = A[j];
						A[j] = tmp;
					}
				if (!MaxRow[i])
					return;
			}
			for (int j = N - 1; j > i; j--)
			{
				A[i][j] /= A[i][i];
				for (int k = i + 1; k < N; k++)
					A[k][j] -= A[k][i] * A[i][j];
			}
			b[i] /= A[i][i];
			A[i][i] = 1;
			for (int j = N - 1; j > i; j--)
			{
				b[j] -= b[i] * A[j][i];
				A[j][i] = 0;
			}
		}
		for (int i = N - 1; i >= 0; i--)
		{
			for (int j = i - 1; j >= 0; j--)
			{
				b[j] -= b[i] * A[j][i];
				A[j][i] = 0;
			}
		}
		memcpy(x,b,sizeof(double)*N);
}



/**
 * Метод прогонки - done and passed
 */
void scherbakovdv::lab3()
{
	double *al = new double[N], *bt = new double[N];
	//Прямой ход
	al[0] = -A[0][1]/A[0][0];
	bt[0] = b[0]/A[0][0];
	for (int i=1;i<N-1;i++){
		al[i] = -A[i][i+1]/(A[i][i] + A[i][i-1]*al[i-1]);
		bt[i] = (b[i] - A[i][i-1]*bt[i-1])/(A[i][i] + A[i][i-1]*al[i-1]);
	}
	x[N-1] = (b[N-1]-A[N-1][N-2]*bt[N-1])/(A[N-1][N-1]+A[N-1][N-2]*al[N-2]);
	//Обратный ход
	for (int i=N-2;i!=-1;i--)
		x[i] = al[i]*x[i+1]+bt[i];
	delete[] al;
	delete[] bt;
}



/**
 * Метод простых итераций - working
 */
void scherbakovdv::lab4()
{
	//Условие работоспособности метода: корень суммы квадратов всех элементов матрицы меньше единицы
	double diff;
	//Допустимая погрешность
	const double Eps=0.1E-10;
	//Аварийный счётчик
	int counter=0;
	//Переходной массив иксов
	double* xOld = new double[N];
	// Матрица итерационной формы системы
	// Для неё условием сходимости является не прывышение суммой всех элементов строки единицы
	double** AA = new double*[N];
	double* bb = new double[N];
	for (int i=0;i<N;i++){
		AA[i]=new double[N];
		for (int j=0;j<N;j++)
			AA[i][j]=A[i][j];
		bb[i]=b[i];
	}
	//Первый этап - приведение матрицы к итерационной форме
	for (int i=0;i<N;i++){
		double maxEl = fabs(AA[i][i]);
		int row = i;
		for (int j=i+1;j<N;j++)
			if (maxEl<fabs(AA[j][i])){
				row = j;
				maxEl=fabs(AA[j][i]);
			}
		if (row!=i){
			for (int j = i; j < N; j++)
				swap(AA[i][j], AA[row][j]);
			swap(bb[i], bb[row]);
		}
		maxEl = AA[i][i];
		bb[i] /= AA[i][i];
		AA[i][i] = 0;
		for (int j = 0; j<N; j++)
			AA[i][j] /= maxEl;
	}
	//Второй этап - получение x
	memset(xOld,0,MSIZE);
	do {
		for (int i=0;i<N;i++){
			x[i]=bb[i];
			for (int j=0;j<N;j++)
				x[i]-=AA[i][j]*xOld[j];
		}
		diff=fabs(x[0]-xOld[0]);
		for (int i=1;i<N;i++)
			if (fabs(x[i]-xOld[i])>diff)
				diff=fabs(x[i]-xOld[i]);
		counter++;
		memcpy(xOld,x,sizeof(double)*N);
	} while ((diff>Eps)&&(counter<1000));
	printf("Counter=%d, diff=%.2f",counter,diff);
	delete[] xOld;
}



/**
 * Метод Якоби или Зейделя - done
 */
void scherbakovdv::lab5()
{
	//Достаточное признак сходимости метода Якоби - диагональное преобладание
	//Средняя разность между новыми и старыми значениями
	double diff;
	//Допустимая погрешность
	const double Eps=0.1E-10;
	//Аварийный счётчик
	int counter=0;
	//Переходной массив иксов
	double* xOld = new double[N];
	memcpy(xOld,b,sizeof(double)*N);
	do {
		for (int i=0;i<N;i++) {
			x[i]=b[i];
			for (int j=0;j<N;j++)
				if (i!=j)
					x[i]-=A[i][j]*xOld[j];
			x[i]/=A[i][i];
		}
		diff=fabs(x[0]-xOld[0]);
		counter++;
		memcpy(xOld,x,sizeof(double)*N);
	} while ((diff>Eps)&&(counter<100));
	delete[] xOld;
	if (diff>Eps) {
		perror("diff>Eps");
	} else if (counter==100) {
		perror("Out of cycle");
	}
}



/**
 * Метод минимальных невязок - done
 */
void scherbakovdv::lab6()
{
	//Допустимая погрешность
	const double Eps=0.1E-10;
	printf("LAB 6: Eps is %.2e",Eps);
	//Аварийный счётчик
	int counter=0;
	//vec - вектор невязок (интересно, почему так называется)
	double *vec = new double[N], diff=0;
	//начальное приближение будет b
	memcpy(x,b,sizeof(double)*N);
	do {
		memcpy(vec,b,sizeof(double)*N);
		for (int i=0;i<N;i++)
			for (int j=0;j<N;j++)
				vec[i]-=A[i][j]*x[j];
		//Поиск параметра тау в формуле X^(k+1)=X^(k)-tau*vec
		//Из-за обилия A*vec я счёл разумным добавить вектор Avec
		double tau=0, *Avec = mul(A,vec);
		tau=-scala(Avec,vec)/scala(Avec,Avec);
		delete[] Avec;
		for (int i=0;i<N;i++)
		{
			x[i]-=tau*vec[i];
			//diff=x[i]-xOld[i], но так экономнее
			diff+=fabs(tau*vec[i]);
		}
		//Для чистоты эксперимента будем проверять среднее арифметическое разности
		diff/=N;
		counter++;
	} while ((diff>Eps)&&(counter<1000));
	delete[] vec;
}



/**
 * Метод сопряженных градиентов - done
 */
void scherbakovdv::lab7()
{
	//Допустимая погрешность
	const double Eps=0.1E-10;
	printf("LAB 7: Eps is %.2e",Eps);
	//Аварийный счётчик
	int counter=0;
	double diff;
	double *AGrad = new double[N];
	double *SDir = new double[N];
	double *xOld = new double[N];
	memcpy(AGrad,b,MSIZE);
		for (int i=0;i<N;i++)
			for (int j=0;j<N;j++)
				AGrad[i]-=A[i][j]*x[j];
	memcpy(SDir,AGrad,MSIZE);
	do {
		double *AGrOld = new double[N];
		double* A_SDir = mul(A,SDir);
		double koeff = scala(AGrad,AGrad)/scala(A_SDir,SDir);
		memcpy(xOld,x,MSIZE);
		memcpy(AGrOld,AGrad,MSIZE);
		for (int i=0;i<N;i++){
			x[i]+=koeff*SDir[i];
			AGrad[i]-=koeff*A_SDir[i];
		}
		delete[] A_SDir;
		koeff = scala(AGrad,AGrad)/scala(AGrOld,AGrOld);
		delete[] AGrOld;
		for (int i=0;i<N;i++)
			SDir[i]=AGrad[i]+koeff*SDir[i];
		diff=fabs(x[0]-xOld[0]);
		for (int i=1;i<N;i++)
			if (fabs(x[i]-xOld[i])>diff)
				diff=fabs(x[i]-xOld[i]);
		counter++;
	} while ((diff>Eps)&&(counter<100));
	// printf("Counter equals:%d",counter);
	delete[] AGrad;
	delete[] SDir;
	delete[] xOld;
}



/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void scherbakovdv::lab8()
{

}



/**
 * Нахождение наибольшего по модолю собственного значения матрицы
 */
void scherbakovdv::lab9() 
{
	
}



std::string scherbakovdv::get_name()
{
  return "Scherbakov D.V.";
}
