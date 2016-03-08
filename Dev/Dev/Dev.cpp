#include "stdafx.h"
#include "../../QODES.hpp"


template <class T>
T MyFunction(T x, T y)
{
	return (x);
}
int main()
{

	cout << "RK38" << endl;
	Algorithm<double> *RK = new RK38<double>(0.10, 2.0, -1.0);
	RK->Eqn.differential_equation = MyFunction;
	RK->Solve();
	
	cout << endl << "RK4" << endl;
	Algorithm<double> *ClassicRK = new RK4<double>(0.10, 2.0, -1.0);
	ClassicRK->Eqn.differential_equation = MyFunction;
	ClassicRK->Solve();

	cout << endl << "Forward Euler" << endl;
	Algorithm<double> *FEuler = new ForwardEuler<double>(0.10, 2.0, -1.0);
	FEuler->Eqn.differential_equation = MyFunction;
	FEuler->Solve();


	cout << endl << "RK45" << endl;
	Algorithm<double> *RKField = new RK45<double>(0.1, 2.0, -1.0);
	RKField->Eqn.differential_equation = MyFunction;
	RKField->Solve();

	
	cout << endl << "RKDP" << endl;
	Adaptive<double> *DormandPrince = new RKDP<double>(0.1, 2.0, -1.0);
	DormandPrince->Eqn.differential_equation = MyFunction;
	DormandPrince->Solve();
	
//	SimulSolve(*DormandPrince, *RKField);

	return 0;
}


