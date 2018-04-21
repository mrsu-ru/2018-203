#include "tarasovams.h"

/**
 * ¬ведение в дисциплину
 */
void tarasovams::lab1()
{

}


/**
 * ћетод √аусса с выбором главного элемента
 */
void tarasovams::lab2()
{
    int i, j, m, c;
    int n = N;
    for(i=0;i<n;i++)
        {
          m=i;
          for(c=i+1; c<n; c++)
          if (abs(A[c][i])>abs(A[m][i]))
                m=c;


          swap(A[m],A[i]);
          swap(b[m],b[i]);
           b[i]/=A[i][i];
          for(j=n-1; j>i; A[i][j--]/=A[i][i]);

          A[i][i]=1;

          for(int j=i+1; j<n;j++){
          for(c=n-1;c>i;c--)
          A[j][c]-=A[i][c]*A[j][i];
          b[j]-=b[i]*A[j][i];

          A[j][i]=0;
          }

            cout<<endl;

          }
          x[n-1]=b[n-1];
          for ( int i = n - 2; i >= 0; i-- )
         {
           x[i] = b[i];
           for ( int j = i + 1; j < n; j++ ) {
              x[i] -= A[i][j] * x[j];
              }
              }

        }


/**
 * ћетод прогонки
 */
void tarasovams::lab3()
{
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
    delete[] Betta;
}





/**
 * ћетод простых итераций
 */
void tarasovams::lab4()
{
double eps = 1e-9;
double tauh = 1e-5;
for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double x1;
	double *xr = new double[N];
	int step = 0;

	do {
		step++;
		for (int i = 0; i < N; i++) {
			xr[i] = x[i];
			for (int k = 0; k < N; k++)
				xr[i] -= tauh*A[i][k] * x[k];
			xr[i] += tauh * b[i];

		}
		x1 = 0.;
		for (int i = 0; i < N; i++) {
			x1 += (xr[i]-x[i])*(xr[i]-x[i]);
		}

		for (int i = 0; i < N; i++) {
			x[i] = xr[i];
		}
		printf(x1, step);
	} while (sqrt(x1)>eps);
}




/**
 * ћетод якоби или «ейдел€
 */
void tarasovams::lab5()
{

double e =0.001;
double* a = new double[N];
double Norm = 0;

        for(int i=0; i<N; i++)
        {
        x[i] = 0;
        }
        do
		{
			for(int i=0; i<N; i++)
			{
				a[i] = b[i];
				for(int j=0; j<N; j++)
				{
                if(i != j)
					{
                    a[i] -= A[i][j]*x[j];
					}
				}
				a[i] /= A[i][i];
			}

		Norm = abs(x[0] - a[0]);

			for(int i=0; i<N; i++)
			{
				if(abs(x[i]-a[i]) > Norm)
                Norm = sqrt((x[i]-a[i])*(x[i]-a[i]));
				x[i] = a[i];
			}
		}
		while (Norm >= e);
		delete[] a;

}



/**
 * ћетод минимальных нев€зок
 */
void tarasovams::lab6()
{

}



/**
 * ћетод сопр€женных градиентов
 */
void tarasovams::lab7()
{

}


void tarasovams::lab8()
{

}

void tarasovams::lab9()
{

}
double f(double x)
{
 return (5-x*exp(x));
}

double f2(double x)
{
 return (-(x+2)*exp(x));
}
void tarasovams::lab10()
{

double c;
do {
if(f(a)*f2(a)>0){
b = b - ((a-b) * f(b))/(f(a) - f(b));

c=b;
}
else if (f(b)*f2(b)>0) {
a = a - ((b-a) * f(a))/(f(b) - f(a));

c=a;
}
} while (fabs(f(c))>=e);

cout<<"x = "<<c<<"\n";
}


std::string tarasovams::get_name()
{
  return "Tarasova M. S.";
}
