#include "seninvs.h"

void seninvs::MxV(int N, double **A, double *V, double *R)
{ 
	for(int i=0; i<N; i++) 
	{ 
		R[i]=0; 
		for(int j=0; j<N; j++) 
			R[i]+= A[i][j]*V[j]; 
	} 
}

/**
 * Введение в дисциплину
 */
void seninvs::lab1() {
  std::cout << "Hello World";
}

/**
 * Метод Гаусса с выбором главного элемента
 */
void seninvs::lab2() {
  for (int i = 0; i < N; i++)
    x[i] = b[i];
  long double l;
  int ind;
  for (int k = 0; k < N - 1; k++) {
    ind = k;
    for (int i = k + 1; i < N; i++)
      if (abs(A[i][k]) > abs(A[ind][k]))
        ind = i;
    std::swap(A[ind], A[k]);
    std::swap(b[ind], b[k]);
    for (int i = k + 1; i < N; i++) {
      l = A[i][k] / A[k][k];
      for (int j = k; j < N; j++) {
        A[i][j] -= l * A[k][j];
      }
      x[i] -= l * x[k];
    }
  }

  for (int i = N - 1; i >= 0; i--) {
    for (int j = i + 1; j < N; j++)
      x[i] -= A[i][j] * x[j];
    x[i] /= A[i][i];
  }
}

/**
 * Метод прогонки
 */
void seninvs::lab3() {
  long double* P = new long double[N];
  long double* Q = new long double[N];
  for (int i = 0; i < N; i++) {
    P[i] = 0;
    Q[i] = 0;
  }

  P[0] = A[0][1] / (-A[0][0]);
  Q[0] = b[0] / A[0][0];

  for (int i = 1; i < N; i++) {
    P[i] = A[i][i + 1] / (-A[i][i] - A[i][i - 1] * P[i - 1]);
    Q[i] =
        (-b[i] + A[i][i - 1] * Q[i - 1]) / (-A[i][i] - A[i][i - 1] * P[i - 1]);
  }

  x[N - 1] = Q[N - 1];
  for (int i = N - 2; i >= 0; i--)
    x[i] = P[i] * x[i + 1] + Q[i];

  delete[] P;
  delete[] Q;
}
/**
 * Метод простых итераций
 */
void seninvs::lab4() {
  double eps = 1e-13;
  double t = 1e-5;
  for (int i = 0; i < N; i++) {
    x[i] = 0;
  }
  double x1;
  double* xr = new double[N];
  int st = 0;

  do {
    st++;
    for (int i = 0; i < N; i++) {
      xr[i] = x[i];
      for (int k = 0; k < N; k++) xr[i] -= t * A[i][k] * x[k];
      xr[i] += t * b[i];
    }
    x1 = 0.;
    for (int i = 0; i < N; i++) {
      x1 += (xr[i] - x[i]) * (xr[i] - x[i]);
    }

    for (int i = 0; i < N; i++) {
      x[i] = xr[i];
    }
  } while (sqrt(x1) > eps);
}

/**
 * Метод Якоби
 */
void seninvs::lab5() {
  long double eps = 0.0001;
  long double* p = new long double[N];
  long double norl;
  for (int i = 0; i < N; i++)
    x[i] = 0;
  do {
    for (int i = 0; i < N; i++) {
      p[i] = b[i];
      for (int j = 0; j < N; j++)
        if (i != j)
          p[i] -= A[i][j] * x[j];
      p[i] /= A[i][i];
    }
    norl = fabs(x[0] - p[0]);
    for (int h = 0; h < N; h++) {
      if (fabs(x[h] - p[h]) > norl)
        norl = fabs(x[h] - p[h]);
      x[h] = p[h];
    }
  } while (norl > eps);
  delete[] p;
}
/**
 * Метод минимальных невязок
 */
void seninvs::lab6() {
  double* r = new double[N];
  double* delt = new double[N];
  double* temp = new double[N];

  double max = 0.0, tau = 0.0, tempt = 0.0;
  double eps = 0.0000001;
  for (int i = 0; i < N; i++)
    temp[i] = 0;
  do {
    MxV(N, A, temp, r);
    for (int i = 0; i < N; i++) {
      delt[i] = r[i] - b[i];
    }
    MxV(N, A, delt, r);
    tau = 0.0;
    tempt = 0.0;
    for (int i = 0; i < N; i++) {
      tau += r[i] * delt[i];
      tempt += r[i] * r[i];
    }
    tau = tau / tempt;
    for (int i = 0; i < N; i++)
      x[i] = temp[i] - tau * delt[i];
    max = fabs(x[0] - temp[0]);
    for (int i = 0; i < N; i++) {
      if (fabs(x[i] - temp[i]) > max)
        max = fabs(x[i] - temp[i]);
      temp[i] = x[i];
    }
  } while (max >= eps);

  delete[] r;
  delete[] delt;
  delete[] temp;
}

/**
 * Метод сопряженных градиентов
 */
void seninvs::lab7() {}

void seninvs::lab8() {}

/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */

void seninvs::lab9() {
  double* y = new double[N];
  double ms, s, sum;
  double eps = 1e-9;
  for (int i = 0; i < N; i++) y[i] = 0;
  y[0] = 1;
  do {
    sum = 0;
    for (int i = 0; i < N; i++) sum += y[i] * y[i];
    s = sqrt(sum);
    for (int i = 0; i < N; i++) {
      y[i] = 0;
      for (int j = 0; j < N; j++) y[i] += A[i][j] * y[j] / s;
    }
    sum = 0;
    for (int i = 0; i < N; i++) sum += y[i] * y[i];
    ms = sqrt(sum);
    for (int i = 0; i < N; i++) y[i] = y[i];
  } while (abs(ms - s) > eps);

  cout << ms << endl;
}

/**
 * Метод касательных
 */
void seninvs::lab10() {
  double f, x, df;
  double eps = 0.000001;
  do {
    f = pow(x, 3) + 3 * x + 1;
    df = 3 * pow(x, 2) + 3;
    x = x - f / df;
  } while (fabs(f) > eps);

  cout << "x0= " << x;
}

std::string seninvs::get_name() {
  return "Senin V. S.";
}