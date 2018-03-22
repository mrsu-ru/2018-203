#pragma once
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <string>

using namespace std;

class lab
{
public:
  //void init(int N, double **A, double *b, double *x);
  void read_file();
  void run(int lab_number);
  void write_result();
  void check_result();
  virtual ~lab();
protected:
  int N;
  double **A = NULL, *b = NULL, *x = NULL;
  double **A_or = NULL, *b_or = NULL;

  virtual std::string get_name();


  /**
   * Метод Гаусса
   */
  virtual void lab1() = 0;
  /**
   * Метод Гаусса с выбором главного элемента
   */
  virtual void lab2() = 0;
  /**
   * Метод квадратного корня (метод Холецкого)
   */
  virtual void lab3() = 0;
  /**
   * Метод прогонки
   */
  virtual void lab4() = 0;
  /**
   * Метод Якоби
   */
  virtual void lab5() = 0;
  /**
   * Метод Зейделя
   */
  virtual void lab6() = 0;
  /**
   * Один из градиентных методов
   */
   virtual void lab7() = 0;
   virtual void lab8() = 0;
   virtual void lab9() = 0;
};
