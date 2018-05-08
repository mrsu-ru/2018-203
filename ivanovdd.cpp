#include "ivanovdd.h"
using namespace std;
/**
 * Введение в дисциплину
 */
void ivanovdd::lab1()
{
std::cout << "Hello World!!!" << std::endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void ivanovdd::lab2()
{
for (int k = 0; k < N; k++) {
	int idmax = -1;
	double max = 0;
for (int l=0; l<N; l++) {
	if (abs(A[k][l]) >= max) {
	max = abs(A[k][l]);
	idmax = l;
}
}

if(idmax != -1){
for (int j = 0; j < N; j++)
{
 swap(A[j][idmax], A[j][k]);
}
swap(b[idmax], b[k]);


for (int i = k + 1; i < N; i++)
{
double c_ki = A[i][k]/A[k][k];
for (int j = k; j < N; j++) {
A[i][j] -= A[k][j]*c_ki;
}
b[i] -= b[k]*c_ki;
}
}
}

for (int k = N-1; k >= 0; k--) {
x[k] = b[k];
for (int i = k+1; i < N; i++) {
x[k] -= A[k][i]*x[i];
}
x[k] /= A[k][k];
}
}



/**
 * Метод прогонки
 */
void ivanovdd::lab3()
{
int n1=N-1;
double *alfa = new double[N];
double *beta = new double[N];
double y = A[0][0];
alfa[0] = -A[0][1] / y;
beta[0] = b[0] / y;
for (int i = 1; i < n1; i++) {
y = A[i][i] + A[i][i - 1] * alfa[i - 1];
alfa[i] = -A[i][i + 1] / y;
beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
}
 
x[n1] = (b[n1] - A[n1][n1-1] * beta[n1-1]) / (A[n1][n1] + A[n1][n1-1] * alfa[n1-1]);
for (int i = n1-1; i >= 0; i--) {
x[i] = alfa[i] * x[i + 1] + beta[i];
}
}



/**
 * Метод простых итераций
 */
void ivanovdd::lab4()
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
void ivanovdd::lab5()
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
 * Метод минимальных невязок
 */
void ivanovdd::lab6()
{
// double* r = new double[N];
// double* delt = new double[N];
// double* temp = new double[N];

// double max = 0.0, tau = 0.0, tempt = 0.0;
// double eps = 0.0000001;
// for (int i = 0; i < N; i++)
// temp[i] = 0;
// do {
// MxV(N, A, temp, r);
// for (int i = 0; i < N; i++) {
// delt[i] = r[i] - b[i];
// }
// MxV(N, A, delt, r);
// tau = 0.0;
// tempt = 0.0;
// for (int i = 0; i < N; i++) {
// tau += r[i] * delt[i];
// tempt += r[i] * r[i];
// }
// tau = tau / tempt;
// for (int i = 0; i < N; i++)
// x[i] = temp[i] - tau * delt[i];
// max = fabs(x[0] - temp[0]);
// for (int i = 0; i < N; i++) {
// if (fabs(x[i] - temp[i]) > max)
// max = fabs(x[i] - temp[i]);
// temp[i] = x[i];
// }
//  } while (max >= eps);

// delete[] r;
// delete[] delt;
// delete[] temp;
}



/**
 * Метод сопряженных градиентов
 */
void ivanovdd::lab7()
{

}


void ivanovdd::lab8()
{

}

void ivanovdd::lab9()
{

}

/**
 * Метод решения нелинейных уравнений
 */
 void ivanovdd::lab10()
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
std::string ivanovdd::get_name()
{
  return "Ivanov D.D.";
}
