#include "QODES.h"


template <class T>
T MyFunction(T x, T y)
{
	// function of dy/dx=x
	return (x);
}
int main()
{

	std::cout << "RK38" << std::endl;
	auto RungeKutta38 = RK38<double>(0.10, 2.0, -1.0);
	RungeKutta38.m_equation.differential_equation = MyFunction;
	RungeKutta38.Solve();
	
	cout << endl << "RK4" << endl;
	auto RungeKutta4 = RK4<double>(0.10, 2.0, -1.0);
	RungeKutta4.m_equation.differential_equation = MyFunction;
	RungeKutta4.Solve();

	cout << endl << "Forward Euler" << endl;
	auto ForwardEulerMethod = ForwardEuler<double>(0.10, 2.0, -1.0);
	ForwardEulerMethod.m_equation.differential_equation = MyFunction;
	ForwardEulerMethod.Solve();


	cout << endl << "RK45" << endl;
	auto RungeKutta45 = RK45<double>(0.1, 2.0, -1.0);
	RungeKutta45.m_equation.differential_equation = MyFunction;
	RungeKutta45.Solve();

	
	cout << endl << "RKDP" << endl;
	auto DormandPrince = RKDP<double>(0.1, 2.0, -1.0);
	DormandPrince.m_equation.differential_equation = MyFunction;
	DormandPrince.Solve();
	
	return 0;
}


