![alt text](https://github.com/ThomasThelen/QODES/raw/master/qodes.png)

[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)  [![license](https://img.shields.io/github/license/mashape/apistatus.svg)]()


# QODES
Quick Ordinary Differential Equation Solver (QODES) is an ODE solving library with a focus on usability and ease.


## Structure

Solver - A type of numerical method. ie Rune Kutta 545

Each solver is a derived class from the abstract Algorithm class. Each solver is initialized with the step size, target x, and the initial condition.

The source can be found in `QODES/src/`

## Supported Algorithms

I probably used some of the algorithms in [this](https://github.com/ThomasThelen/Runge-Kutta) repository.

Format is Algorithm : Class Name

Runge Kutta 4 : RK4

Runge Kutta 3/8 : RK38

Runge–Kutta–Fehlberg : RK45

Runge-Kutta-Dormand-Prince

Forward Euler : ForwardEuler


## Compiling
To compile the project, include the QUODES.hpp header file in your main document.


## Usage

The library works by interfacing the Algorithm base class. See the code below for an example. The constructor takes the step size, final x value, and initial condition. For adaptive methods, supply the initial step size and subsequent sizes will be automatically computed.
 
 ```c++
auto RK = std::make_shared<Algorithm<double>> RK38<double>(0.10, 10, 2);
```
  where
  
  0.10 is the step size
  
  10 is the final x
  
  2 is the initial condition
  
  Once the algorithm is created, the ordinary differential equation must be configured. The ODE is represented by a function with a double return type. 
  ```c++
  template <class T>
  T MyFunction(T x, T y)
{
	T result = x + y;
	return result;
}
```
The next step is to point the Eqn.differential_equation to the appropriate location. This is done for every solution method.
  ```c++
  RK->Eqn.differential_equation = MyFunction;
```
  
  The final step is to initiate the Solve() method. This will solve the ODE pointed to by Eqn.differential equation in conjunction with the parameters in the constructor.
  ```c++
  RK->Solve();
```
    
    
## Example main.cpp
This example solves the following differential equation at x=10 with a step of 0.5 and initial condition y(0)=1.
    

  ```c++
#include "stdafx.h"
#include "QODE.hpp"

int main()
{
	std::cout<< "RK4" << std::endl;
	Algorithm *ClassicRK = new RK4(0.50, 10, 1);
	ClassicRK->Eqn.differential_equation = MyFunction;
	ClassicRK->Solve();
	delete ClassicRK;
	return 0;
}

template <class T>
T MyFunction(T x, T y)
{
	T result = x;
	return result;
}

    ```
  
  
