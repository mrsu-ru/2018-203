#pragma once
#include "lab.h"

class scherbakovdv : public lab
{
  virtual double scala(double* l, double* r);
  virtual double* mul(double** mat, double* vec);
  virtual double** mul(double** mat1, double** mat2);
  /**
   * Метод Гаусса
   */
  virtual void lab1();
  /**
   * Метод Гаусса с выбором главного элемента
   */
  virtual void lab2();
  /**
   * Метод квадратного корня (метод Холецкого)
   */
  virtual void lab3();
  /**
   * Метод прогонки
   */
  virtual void lab4();
  /**
   * Метод Якоби
   */
  virtual void lab5();
  /**
   * Метод Зейделя
   */
  virtual void lab6();
  /**
   * Один из градиентных методов
   */
  virtual void lab7();
  /**
   * Метод вращения для нахождения собственных значений матрицы
   */
  virtual void lab8();
  /**
   * Нахождение наибольшего по модолю собственного значения матрицы
   */
  virtual void lab9();


  virtual std::string get_name();

};
