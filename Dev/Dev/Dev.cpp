#include "stdafx.h"
#include "../../QODES.hpp"
//double MyFunction(double x, double y);


int main()
{

	cout << "RK38" << endl;
	Algorithm<double> *RK = new RK38<double>(0.10, 10.0, 1.0);
	RK->Eqn.differential_equation = MyFunction;
	RK->Solve();
	
	cout << endl << "RK4" << endl;
	Algorithm<double> *ClassicRK = new RK4<double>(0.10, 10.0, 1.0);
	ClassicRK->Eqn.differential_equation = MyFunction;
	ClassicRK->Solve();

	cout << endl << "Forward Euler" << endl;
	Algorithm<double> *FEuler = new ForwardEuler<double>(0.10, 10, 1.0);
	FEuler->Eqn.differential_equation = MyFunction;
	FEuler->Solve();


	cout << endl << "RK45" << endl;
	Algorithm<double> *RKField = new RK45<double>(1.0, 10.0, 1.0);
	RKField->Eqn.differential_equation = MyFunction;
	RKField->Solve();


	cout << endl << "RKDP" << endl;
	Algorithm<double> *DormandPrince = new RKDP<double>(1.0, 10.0, 1.0);
	DormandPrince->Eqn.differential_equation = MyFunction;
	DormandPrince->Solve();
	
	return 0;
}


