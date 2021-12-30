![alt text](https://github.com/ThomasThelen/QODES/raw/master/qodes.png)

[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)  [![license](https://img.shields.io/github/license/mashape/apistatus.svg)]()
[![CodeFactor](https://www.codefactor.io/repository/github/thomasthelen/qodes/badge)](https://www.codefactor.io/repository/github/thomasthelen/qodes)

# QODES
Quick Ordinary Differential Equation Solver (QODES) is a header-only ODE solving library with a focus on end user usability. The library supports a small, but standard set of numerical 
methods for solving ODEs. These include the Forward Euler method and the Rune Kutta class.

Note that there are much more robust libraries out there like Odient; this was largely a project that I had worked on for fun and can be used for reference.
# Building
To build the example, run the following from inside the `build/` directory,

```
cmake ..
cmake --build .
```

# Quick Start
To use this, download `QODES.h` and include it in your C++ project. A list of classes are provided below that represent the different solution algorithms.

```
Runge Kutta 4 : RK4
Runge Kutta 3/8 : RK38
Runge–Kutta–Fehlberg : RK45
Runge-Kutta-Dormand-Prince: RKDP
Forward Euler : ForwardEuler
```

The constructor of these classes takes the step size, the target x value in question, and an initial condition. For adaptive methods, supply the initial step size and subsequent sizes will be automatically computed.

The general form of creating an algorithm follows, where 0.10 is the step size, 10 is the x value that we want the slope at, and 2 is the initial condition `y(0)=2`. In this case, the Runge Kutta 3/8 method is chosen.
 ```c++
auto RK = std::make_shared<RK38<double>> RK38<double>(0.10, 10, 2);
```

 Once the algorithm is created, the differential should be created. The following is a good template to use; note that `y` is _not_ neccessary in the ODE but should still be a function parameter.
  ```c++
  template <class T>
  T MyFunction(T x, T y)
{
	T result = x + y;
	return result;
}
```
The final steps are to set the Eqn.differential_equation to the ODE from above and call the `Solve()` method on the solver.
  ```c++
  RK->Eqn.differential_equation = MyFunction;
```
  ```c++
  RK->Solve();
```

## Example main.cpp
This example solves a basic differential equation at x=10 with a step of 0.5 and initial condition y(0)=1.
  ```c++
#include "stdafx.h"
#include "QODE.hpp"

int main()
{
	std::cout<< "Solving dy/dx=x with the RK4 method..." << std::endl;
	Algorithm ClassicRK = RK4(0.50, 10, 1);
	ClassicRK.Eqn.differential_equation = MyFunction;
	ClassicRK.Solve();
	return 0;
}

template <class T>
T MyFunction(T x, T y)
{
	T result = x;
	return result;
}
```
  
# Sources
http://depa.fquim.unam.mx/amyd/archivero/DormandPrince_19856.pdf

http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html