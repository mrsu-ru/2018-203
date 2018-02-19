#include "lab.h"
#include <cstdlib>
#include <cstring>
#include <iostream>

using std::cout;
using std::endl;


void lab::read_file() {
    N = 100;
    A = new double*[N];
    for (int i = 0;  i < N; i++) {
        A[i] = new double[N];
        memset(A[i], 0, sizeof(double)*N);
    }

    b = new double[N];
    x = new double[N];

    memset(x, 0, sizeof(double)*N);

    double h = 3.1415926535/(N-1);
    double h2 = h*h;

    A[0][0] = 20.0;
    A[0][1] = -1.0;
    b[0]    = -sin(h)*h2;
    for (int i = 1; i < N-1; i++) {
        A[i][i-1]   =  -1.0;
        A[i][i]     =  20.0;
        A[i][i+1]   =  -1.0;

        b[i]        = -sin(h*i)*h2;
    }
    A[N-1][N-2] =  -1.0;
    A[N-1][N-1] =  20.0;
    b[N-1]      = -A[N-1][N-2]*sin(h*(N-2))*h2;

    A_or = new double*[N];
    for (int i = 0;  i < N; i++) {
        A_or[i] = new double[N];
        memcpy(A_or[i], A[i], sizeof(double)*N);
    }
    b_or = new double[N];
    memcpy(b_or, b, sizeof(double)*N);
}


void lab::run(int lab_number)
{
    std::cout << get_name() << " passes lab #" << lab_number << std::endl;

    switch (lab_number) {
        case 1: this->lab1(); break;
        case 2: this->lab2(); break;
        case 3: this->lab3(); break;
        case 4: this->lab4(); break;
        case 5: this->lab5(); break;
        case 6: this->lab6(); break;
        case 7: this->lab7(); break;
        case 8: this->lab8(); break;
        case 9: this->lab9(); break;
    }
}


void lab::write_result()
{
    cout << endl << "RESULT:" << endl << "*************************" << endl << endl;
    for (int i = 0; i < N; i++) {
        cout << "x[";
        cout.width(4);
        cout << i << "] = ";
        cout.precision(5);
        cout.width(16);
        cout << x[i] << endl;
    }
}


void lab::check_result()
{
    double *r = new double[N];

    for (int i = 0; i < N; i++) {
        r[i] = 0.0;
        for (int j = 0; j < N; j++) {
            r[i] += A_or[i][j]*x[j];
        }
        r[i] -= b_or[i];
    }
    cout << endl << "RESIDUAL:" << endl << "*************************" << endl << endl;
    double maxR = 0, sqrtR = 0, sumR = 0;
    for (int i = 0; i < N; i++) {
        cout << "r[";
        cout.width(4);
        cout << i << "] = ";
        cout.precision(5);
        cout.width(16);
        cout << r[i] << endl;

        if (fabs(r[i]) > maxR) maxR = fabs(r[i]);
        sumR += fabs(r[i]);
        sqrtR += r[i]*r[i];
    }
    sqrtR = sqrt(sqrtR);

    cout << "sumR  = " << sumR  << endl;
    cout << "sqrtR = " << sqrtR << endl;
    cout << "maxR  = " << maxR  << endl;

    delete[] r;
}


std::string lab::get_name()
{
  return std::string("Unknown People");
}


/**
 *
 */

lab::~lab()
{
  if (A != NULL) {
    for (int i = 0; i < N; i++) {
      if (A[i] != NULL) {
        delete[] A[i];
      }
    }
    delete[] A;
  }
  if (b != NULL) {
    delete[] b;
  }
  if (x != NULL) {
    delete[] x;
  }


  if (A_or != NULL) {
    for (int i = 0; i < N; i++) {
      if (A_or[i] != NULL) {
        delete[] A[i];
      }
    }
    delete[] A_or;
  }
  if (b_or != NULL) {
    delete[] b_or;
  }

}
