#include "stdafx.h"
#include "QODE.hpp"
double MyFunction(double x, double y);


int main()
{
	cout<< "RK38" << endl;
	Algorithm *RK = new RK38(0.10, 10, 1);
	RK->Eqn.differential_equation = MyFunction;
	RK->Solve();

	cout << endl << "RK4" << endl;
	Algorithm *ClassicRK = new RK4(0.10, 10, 1);
	ClassicRK->Eqn.differential_equation = MyFunction;
	ClassicRK->Solve();

	cout << endl << "Forward Euler" << endl;
	Algorithm *FEuler = new ForwardEuler(0.10, 10, 1);
	FEuler->Eqn.differential_equation = MyFunction;
	FEuler->Solve();


	cout << endl << "RK45" << endl;
	Algorithm *RKField = new RK45(1, 10, 1);
	RKField->Eqn.differential_equation = MyFunction;
	RKField->Solve();


	cout << endl << "RKDP" << endl;
	Algorithm *DormandPrince = new RKDP(1, 10, 1);
	DormandPrince->Eqn.differential_equation = MyFunction;
	DormandPrince->Solve();

	return 0;
}


double MyFunction(double x, double y)
{
	double result = x;
	return result;
}
