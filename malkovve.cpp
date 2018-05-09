#include "malkovve.h"

/**
 * Введение в дисциплину
 */
void malkovve::lab1()
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



/**
 * Метод Гаусса с выбором главного элемента
 */
void malkovve::lab2()
{
	double Q = 0;

    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            Q = A[j][i] / A[i][i];

            for (int k = i; k < N; k++)
            {
                A[j][k] -= Q * A[i][k];
            }
            b[j] -= Q * b[i];
        }
    }

    for(int i = 0; i < N; i++)
	{
        x[i] = b[i];
	}

    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
		{
			x[i] -= A[i][j] * x[j];
		}

        x[i] /= A[i][i];
	}
}



/**
 * Метод прогонки
 */
void malkovve::lab3()
{
double Q = 0;
	int max;

    for (int i = 0; i < N - 1; i++)
    {
				max = i;

		for (int j = i + 1; j < N; j++)
		{
			if(abs(A[j][i]) > abs(A[max][i]))
			{
				max = j;
			}
		}


				std::swap(A[max], A[i]);
		    std::swap(b[max], b[i]);

        for (int j = i + 1; j < N; j++)
        {
            Q = A[j][i] / A[i][i];
            for (int k = i; k < N; k++)
            {
                A[j][k] -= Q * A[i][k];
            }
            b[j] -= Q * b[i];
        }
    }

    for(int i = 0; i < N; i++)
	{
        x[i] = b[i];
	}

    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
		{
			x[i] -= A[i][j] * x[j];
		}

        x[i] /= A[i][i];
	}
}



/**
 * Метод простых итераций
 */
void malkovve::lab4()
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
void malkovve::lab5()
{
    //Якоби
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
void malkovve::lab6()
{
    double eps = 1e-5;
    double delta, r, rModul;
    
    double *w = new double[N];
    double *v = new double[N];
    double *result = new double[N];
    
    for (int i = 0; i<N; i++)
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
        if (r==rModul) {r=1;}
        else {r = r / rModul;}
        for (int i = 0; i < N; i++)
        x[i] = result[i] - r*v[i];
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
void malkovve::lab7()
{
    double eps = 1e-5;
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
        r = 0.0;
        rModul = 0.0;
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
    } while (eps < delta);
}


void malkovve::lab8()
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


void malkovve::lab9()
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
        for (int i = 0; i<N; i++)
        Y[i] = y[i];
    } while (abs(maxSob - sob)>eps);
    
    cout << maxSob << endl;
}

/** 
* Решение нелинейных уравнений 
*/ 


std::string malkovve::get_name()
{
  return "Malkov V.E.";
}
