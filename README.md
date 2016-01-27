# QODES
Quick Ordinary Differential Equation Solver (QODES) is an ODE solving library with a focus on usability and ease.


##
Background
Each type of solver is a derived Algorithm. The Algorithm base class is abstract and contains virtual methods for solving, displaying, and setting initial conditions. 

## Supported Algorithms

Format is Algorithm : Class Name

Runge Kutta 4 : RK4

Runge Kutta 3/8 : RK38

Runge–Kutta–Fehlberg : RK45

Forward Euler : ForwardEuler



## Usage

First, create a new interface to the desired algorithm, shown below. During the creation, the step size, final x value, and initial condition are supplied.
  ```c++
Algorithm *RK = new RK38(0.10, 10, 2);
```
  where
  
  0.10 is the step size
  
  10 is the final x
  
  2 is the initial condition
  
  Once the algorithm is created, the ordinary differential equation must be configured. The ODE is represented by a function with a double return type. 
  ```c++
  double MyFunction(double x, double y)
{
	double result = x + y;
	return result;
}
```
The next step is to point the Eqn.differential_equation to the appropriate location.
  ```c++
  RK->Eqn.differential_equation = MyFunction;
```
  
  The final step is to initiate the Solve() method. This will solve the ODE pointed to by Eqn.differential equation in conjunction with the parameters in the constructor.
  ```c++
  RK->Solve();
```
    
    
##Example main.cpp
This example solves the following differential equation at x=10 with a step of 0.5 and initial condition y(0)=1.
    

  ```c++
#include "stdafx.h"
#include "QODE.hpp"
double MyFunction(double x, double y);

int main()
{
	std::cout<< "RK4" << std::endl;
	Algorithm *ClassicRK = new RK4(0.50, 10, 1);
	ClassicRK->Eqn.differential_equation = MyFunction;
	ClassicRK->Solve();
	return 0;
}

double MyFunction(double x, double y)
{
	double result = x + y;
	return result;
}

    ```
  
  
