#include "pankratovdv.h"

/**
 * Решение нелинейных уравнений
 */

void pankratovdv::lab1() {
	double x0;
	for (int i = 0; i < 100; i++) {
        double xd;
        double eps = 1e-5;
        x0 = i;
        int ind = 0;
        do {
            ind++;
            xd = x0;
            x0 = exp(-x0);
        } while (abs(xd - x0) > eps || ind > 1000);
        if (x0 == x0) {
            break;
        }
    }
	cout << x0 << '\n';
}

/**
 * Метод Гаусса с выбором главного элемента
 */

void pankratovdv::lab2() {
	for (int i = 0; i < N; i++) {
		double mEL = 0;
		int indRow = i;
		for (int j = i + 1; j < N; j++)
			if (mEL < abs(A[j][i])) {
				indRow = j;
				mEL = abs(A[j][i]);
			}
        if (indRow != i) {
            for (int j = i; j < N; j++) {
                swap (A[i][j], A[indRow][j]);
            }
            swap(b[i], b[indRow]);
        }
        mEL = A[i][i];
        b[i] /= A[i][i];
        A[i][i] = 1;
        for (int j = i + 1; j < N; j++) {
            A[i][j] /= mEL;
        }
		for (int k = i + 1; k < N; k++) {
			double mult = A[k][i];
			A[k][i] = 0;
			if (mult != 0) {
				for (int j = i + 1; j < N; j++) {
					A[k][j] -= mult * A[i][j];
				}
				b[k] -= mult * b[i];
			}
		}
	}
	for (int k = N - 1; k >= 0; k--) {
		x[k] = 0;
		double sum = b[k];
		for (int j = N-1; j>k; j--) {
			sum -= A[k][j] * x[j];
		}
		sum -= b[k] * x[k];
		x[k] = sum;
	}
}

/**
 * Метод прогонки
 */

void pankratovdv::lab3() {
	double *ap = new double[N];
	double *bp = new double[N];
	double y = A[0][0];
	ap[0] = -A[0][1] / y;
	bp[0] = b[0] / y;
	for (int i = 1; i < N - 1; i++) {
		y = A[i][i] + A[i][i - 1] * ap[i - 1];
		ap[i] = -A[i][i + 1] / y;
		bp[i] = (b[i] - A[i][i - 1] * bp[i - 1]) / y;
	}
	x[N - 1] = (b[N - 1] - A[N - 1][N - 2] * bp[N - 2]) / (A[N - 1][N - 1] + A[N - 1][N - 2] * ap[N - 2]);
	for (int i = N - 2; i >= 0; i--) {
		x[i] = ap[i] * x[i + 1] + bp[i];
	}
}

/**
 * Метод простых итераций
 */

void pankratovdv::lab4() {
	double eps = 1e-5;
	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double x1 = b[0];
	double *xr = new double[N];
	do {
		for (int i = 0; i < N; i++) {
			xr[i] = 0;
			for (int k = 0; k < N; k++) {
				xr[i] -= A[i][k] * x[k];
			}
			xr[i] += b[i];
		}
		x1 = x[0];
		for (int i = 0; i < N; i++) {
			x[i] = xr[i];
		}
	} while (abs(x[0] - x1) > eps);
}

/**
 * Метод Якоби или Зейделя
 */

void pankratovdv::lab5() {
	double eps = 1e-5;
	double *xR = new double[N];
	double dt = 0.0;
	for (int i = 0; i < N; i++) {
		x[i] = b[i];
	}
	do {
		for (int i = 0; i < N; i++) {
			xR[i] = b[i];
			for (int j = 0; j < N; j++) {
				if (i != j) {
					xR[i] -= A[i][j] * x[j];
				}
			}
			xR[i] /= A[i][i];
		}
		dt = abs(x[0] - xR[0]);
		for (int i = 0; i < N; i++) {
			if (dt < abs(x[i] - xR[i])) {
				dt = abs(x[i] - xR[i]);
			}
			x[i] = xR[i];
		}
	} while (eps < dt);
}

/**
 * Метод минимальных невязок
 */

void pankratovdv::lab6() {
	double eps = 1e-5;
	double dt = 0.0, r = 0.0, rM = 0.0;
	double *w = new double[N];
	double *v = new double[N];
	double *res = new double[N];
	for (int i = 0; i<N; i++) {
		res[i] = 0;
	}
	do {
		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++) {
				w[i] += A[i][j] * res[j];
			}
		}
		for (int i = 0; i < N; i++) {
			v[i] = w[i] - b[i];
		}
		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++) {
				w[i] += A[i][j] * v[j];
			}
		}
		r = 0.0;
		rM = 0.0;
		for (int i = 0; i < N; i++) {
			r += w[i] * v[i];
			rM += w[i] * w[i];
		}
		if (r == rM) {
            r = 1;
        } else {
            r = r / rM;
        }
		for (int i = 0; i < N; i++) {
			x[i] = res[i] - r * v[i];
		}
		dt = abs(x[0] - res[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - res[i]) > dt) {
				dt = abs(x[i] - res[i]);
			}
			res[i] = x[i];
		}
	} while (eps < dt);
}

/**
 * Метод сопряженных градиентов
 */

void pankratovdv::lab7() {
	double eps = 1e-5;
	double dt, r = 0.0, rM = 0.0;
	double *w = new double[N];
	double *wp = new double[N];
	double *v = new double[N];
	double *res = new double[N];
	for (int i = 0; i < N; i++) {
		res[i] = 0;
	}
	do {
		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++) {
				w[i] += A[i][j] * res[j];
			}
		}
		for (int i = 0; i < N; i++) {
			v[i] = w[i] - b[i];
		}
		for (int i = 0; i < N; i++) {
			w[i] = 0;
			for (int j = 0; j < N; j++) {
				w[i] += A[i][j] * v[j];
			}
		}
		for (int i = 0; i < N; i++) {
			wp[i] = 0;
			for (int j = 0; j < N; j++) {
				wp[i] += A[i][j] * w[j];
			}
		}
		r = 0.0;
		rM = 0.0;
		for (int i = 0; i < N; i++) {
			r += w[i] * v[i];
			rM += wp[i] * w[i];
		}
		if (r == rM) {
            r = 1;
		} else {
		    r = r / rM;
		}
		for (int i = 0; i < N; i++) {
			x[i] = res[i] - r * v[i];
		}
		dt = abs(x[0] - res[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - res[i]) > dt) {
				dt = abs(x[i] - res[i]);
			}
			res[i] = x[i];
		}
	} while (eps < dt);
}

/**
 * Метод вращения для нахождения собственных значений матрицы
 */

void pankratovdv::lab8() {
    double *sob = new double[N];
	double *b = new double[N];
	double *z = new double[N];
	for (int i = 0; i < N; i++) {
		z[i] = 0.0;
		b[i] = A[i][i];
		sob[i] = A[i][i];
	}
	for (int i = 0; i < 100; i++) {
		double sm = 0.0;
		for (int p = 0; p < N - 1; p++) {
			for (int q = p + 1; q < N; q++) {
				sm += abs(A[p][q]);
			}
		}
		if (sm == 0) {
            break;
		}
		double tr = sm / (5 * N * N);
		for (int p = 0; p < N - 1; p++) {
			for (int q = p + 1; q < N; q++) {
				double g = 1e12 * abs(A[p][q]);
				if (i >= 3 && abs(sob[p]) > g && abs(sob[q]) > g) {
                    A[p][q] = 0.0;
				} else {
					if (abs(A[p][q]) > tr) {
						double th = (sob[q] - sob[p]) / (2.0 * A[p][q]);
						double t = 1.0 / (abs(th) + sqrt(1.0 + th * th));
						if (th < 0) {
                            t *= -1;
						}
					    double s = t / sqrt(1.0 + t * t);
						double tu = s / (1.0 + 1.0 / sqrt(1.0 + t * t));
						z[p] -= t * A[p][q];
						z[q] += t * A[p][q];
						sob[p] -= t * A[p][q];
						sob[q] += t * A[p][q];
						A[p][q] = 0.0;
						for (int j = 0; j < p; j++) {
							A[j][p] = A[j][p] - s * (A[j][q] + A[j][p] * tu);
							A[j][q] = A[j][q] + s * (A[j][p] - A[j][q] * tu);
						t
						for (int j = p + 1; j < q; j++) {
							A[p][j] = A[p][j] - s * (A[j][q] + A[p][j] * tu);
							A[j][q] = A[j][q] + s * (A[p][j] - A[j][q] * tu);
						}
						for (int j = q + 1; j < N; j++) {
							A[p][j] = A[p][j] - s * (A[q][j] + A[p][j] * tu);
							A[q][j] = A[q][j] + s * (A[p][j] - A[q][j] * tu);
						}
					}
				}
			}
		}
		for (int p = 0; p < N; p++) {
			sob[p] = b[p] + z[p];
		}
	}
	for (int p = 0; p < N; ++p) {
        cout << sob[p] << '\n';
	}
}

/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */

void pankratovdv::lab9() {
	double * Y = new double[N];
	double * y = new double[N];
	double mS, sob, sum;
	double eps = 1e-5;
	for (int i = 0; i < N; i++) {
		Y[i] = 0;
	}
	Y[0] = 1;
	do {
		sum = 0;
		for (int i = 0; i < N; i++) {
			sum += Y[i] * Y[i];
		}
		sob = sqrt(sum);
		for (int i = 0; i < N; i++) {
			y[i] = 0;
			for (int j = 0; j < N; j++) {
				y[i] += A[i][j] * Y[j] / sob;
			}
		}
		sum = 0;
		for (int i = 0; i < N; i++) {
			sum += y[i] * y[i];
		}
		mS = sqrt(sum);
		for (int i = 0; i<N; i++) {
			Y[i] = y[i];
		}
	} while (abs(mS - sob) > eps);

	cout << mS << endl;
}

std::string pankratovdv::get_name()
{
  return "Pankratov D.V.";
}
