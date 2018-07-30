#pragma once

#include <iostream>
#include <vector>
#include <math.h>

using std::vector;
using std::cout;
using std::endl;


// Abstract Base Class
class Algorithm
{
public:
	virtual float Solve() = 0; // Main Algorithm
	Algorithm()
	{
		step_size = 0;
		target_x = 0;
	};
	virtual void SetIC() = 0;
	void FillX();

	class Equation // Holds Equation-Based Data
	{
	public:
		int dimensions;
		float initial_condition; // IC @ x=0
		vector<double> x, y; // Containers for x and y
		double(*differential_equation)(double, double);
	}Eqn;

protected:
	float step_size; // Step Size Used
	float target_x; // x Value of Interest
	float initial_condition;
};
void Algorithm::FillX()
{
	for (double i = 0; i <= this->target_x; i += this->step_size)
	{
		this->Eqn.x.push_back(i);
	}
}

// Definitions for the RK 3/8 method
class RK38 : public Algorithm
{
public:
	RK38(float step, float final_x, float IC);
	float Solve();
	void SetIC();
private:
	double k1, k2, k3, k4;
};

RK38::RK38(float step, float final_x, float IC)
{
	this->step_size = step;
	this->target_x = final_x;
	this->initial_condition = IC;
};
void RK38::SetIC()
{
	this->Eqn.x.push_back(0);
	this->Eqn.y.push_back(initial_condition);
}
float RK38::Solve()
{
	SetIC();

	for (float i = 1; i<this->target_x; i+=step_size)
	{
		this->Eqn.x.push_back(this->Eqn.x.back() + step_size);
		k1 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(i - 1), this->Eqn.y.at(i - 1));
		k2 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(i - 1) + (this->step_size / 3), this->Eqn.y.at(i - 1) + (k1 / 3.0));
		k3 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(i - 1) + (2 * this->step_size) / 3.0, this->Eqn.y.at(i - 1) - (k1 / 3.0) + k2);
		k4 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(i - 1) + this->step_size, this->Eqn.y.at(i - 1) + k1 - k2 + k3);
		this->Eqn.y.push_back(this->Eqn.y.at(i - 1) + (1.0 / 8.0)*(k1 + 3 * k2 + 3 * k3 + k4));
		cout << "X: " << this->Eqn.x.back() << "    " << "Y: " << this->Eqn.y.back() << endl;
	}

	return 0;
}

// Definitions for the Class RK Method
class RK4 : public Algorithm
{
public:
	RK4(float, float, float);
	float Solve();
	void SetIC();
protected:
	double k1, k2, k3, k4;
};
RK4::RK4(float step, float final_x, float IC)
{
	this->step_size = step;
	this->target_x = final_x;
	this->initial_condition = IC;
}
void RK4::SetIC()
{
	this->Eqn.x.push_back(0);
	this->Eqn.y.push_back(initial_condition);
}
float RK4::Solve()
{
	SetIC();
	double result;
	for (int position = 1; Eqn.x.back() < target_x; ++position)
	{
		Eqn.x.push_back(Eqn.x.back() + step_size);

		k1 = step_size*Eqn.differential_equation(Eqn.x.at(position-1), 2); //Eqn.differential_equation(Eqn.x.at(position), Eqn.y.at(position-1));
		k2 = step_size*Eqn.differential_equation(Eqn.x.at(position) + 0.5*step_size, Eqn.y.at(position-1) + 0.5*k1);
		k3 = step_size*Eqn.differential_equation(Eqn.x.at(position) + 0.5*step_size, Eqn.y.at(position-1) + 0.5*k2);
		k4 = step_size*Eqn.differential_equation(Eqn.x.at(position) + step_size, Eqn.y.at(position-1) + k3);
		result = Eqn.y.back() + (1.0 / 6.0)*k1 + (1.0 / 3.0)*k2 + (1.0 / 3.0)*k3 + (1.0 / 6.0)*k4;
		Eqn.y.push_back(result);
		cout << "X: " << this->Eqn.x.back() << "    " << "Y: " << this->Eqn.y.back() << endl;
	}
	return 0;
}

// Definitions for the forward Euler Methods
class ForwardEuler : public Algorithm
{
public:
	ForwardEuler(float, float, float);
	float Solve();
	void SetIC();
private:
};
ForwardEuler::ForwardEuler(float step, float final_x, float IC)
{
	step_size = step;
	target_x = final_x;
	initial_condition = IC;
}
void ForwardEuler::SetIC()
{
	this->Eqn.x.push_back(0);
	this->Eqn.y.push_back(this->Eqn.initial_condition);
}
float ForwardEuler::Solve()
{
	cout << "Solving..." << endl;
	SetIC();
	while (this->Eqn.x.back() < target_x)
	{
		this->Eqn.x.push_back(this->Eqn.x.back() + step_size);
		this->Eqn.y.push_back(this->Eqn.y.back() + step_size*Eqn.differential_equation(this->Eqn.x.back(),this->Eqn.y.back()));
		cout << "X: " << this->Eqn.x.back() << "    " << "Y: " << this->Eqn.y.back() << endl;
	}
	return 0;
}

// Runge–Kutta–Fehlberg 
class RK45 : public Algorithm
{
public:
	RK45(float step, float final_x, float IC);
	float Solve();
	void SetIC();
private:
	void CalculateStepSize(double, double);
	double k1, k2, k3, k4, k5, k6;
	bool flag;
};
RK45::RK45(float step, float final_x, float IC)
{
	this->step_size = step;
	this->target_x = final_x;
	this->initial_condition = IC;
	this->Eqn.y.push_back(initial_condition);
	this->Eqn.x.push_back(0);
};
float RK45::Solve()
{
	double temp_y, temp_x = 0;
	double position = 1;
	for (double i = 1; this->Eqn.x.back()<= this->target_x+step_size; ++i)
	{
		k1 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1), this->Eqn.y.at(position - 1));
		k2 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + (this->step_size / 4), this->Eqn.y.at(position - 1) + (k1 / 4.0));
		k3 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + 3.0*(this->step_size / 8.0), this->Eqn.y.at(position - 1) + (3.0*k1 / 32.0) + k2*(9.0 / 32.0));
		k4 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + this->step_size, this->Eqn.y.at(position - 1) + (1932.0 / 2197.0)*k1 - (7200 / 2197)*k2 + (7296 / 2197)*k3);
		k5 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + this->step_size, this->Eqn.y.at(position - 1) + (439 / 216)*k1 - (8.0*k2) + (3680.0 / 513.0)*k3 - (845 / 4104)*k4);
		k6 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + 0.5*this->step_size, this->Eqn.y.at(position - 1) - (8.0 / 27)*k1 + 2.0*k2 - (3544 / 2565)*k3 + (1859 / 4104)*k4 - (11.0 / 40.0)*k5);
		temp_y = this->Eqn.y.at(position - 1) + (16.0 / 135.0)*k1 + (6656.0 / 12825.0)*k3 + (28561.0 / 56430.0)*k4 - (9.0 / 50.0)*k5 + (2.0 / 55.0)*k6;
		cout << "X: " << this->Eqn.x.back() << "    " << "Y: " << this->Eqn.y.back() << endl;

		if (flag == false)
		{
			temp_x += step_size;
			this->Eqn.x.push_back(temp_x);
			CalculateStepSize(temp_y, temp_x);
			this->Eqn.y.push_back(temp_y);
		}
		else
		{
			this->Eqn.y.pop_back();
			this->Eqn.y.push_back(temp_y);
			flag = false;
			position++;
		}
	}
	return 0;
}
void RK45::CalculateStepSize(double t_y, double t_x)
{
	double s;
	s = 0.84*pow(step_size / (2.0 * abs((t_x - t_y))), 0.25);
	if (s*this->Eqn.x.back() < 1)
		step_size = s*t_x;
	flag = 1;
}
void RK45::SetIC()
{
	cout << "Garbage" << endl;
}








class RKDP : public Algorithm
{
public:
	RKDP(float step, float final_x, float IC);
	float Solve();
	void SetIC();
private:
	void CalculateStepSize(double, double);
	double k1, k2, k3, k4, k5, k6, k7;
	bool flag;
};
RKDP::RKDP(float step, float final_x, float IC)
{
	this->step_size = step;
	this->target_x = final_x;
	this->initial_condition = IC;
	this->Eqn.y.push_back(initial_condition);
	this->Eqn.x.push_back(0);
};
float RKDP::Solve()
{
	double temp_y, temp_x = 0;
	double position = 1;
	for (double i = 1; this->Eqn.x.back()<= this->target_x+step_size; ++i)
	{

		k1 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1), this->Eqn.y.at(position - 1));
		k2 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + (this->step_size / 5.0), this->Eqn.y.at(position - 1) +k1/5.0);
		k3 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + 3.0*(this->step_size / 10.0), this->Eqn.y.at(position - 1) + (3.0*k1 / 40.0) + k2*(9.0 / 40.0));
		k4 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + this->step_size*4.0/5.0, this->Eqn.y.at(position - 1) + (44.0 / 45.0)*k1 - (56.0/15.0)*k2 + (32.0/9.0)*k3);
		k5 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + this->step_size*8.0/9.0, this->Eqn.y.at(position - 1) + (19372.0/6561.0)*k1 - (25360.0/2187.0)*k2 + (6448.0/6561.0)*k3 - (212.0/729.0)*k4);
		k6 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + this->step_size, this->Eqn.y.at(position - 1) + (9017.0/3168.0)*k1 - (355.0/33.0)*k2 - (46732.0/5247.0)*k3 + (49.0/176.0)*k4 - (5103.0/18656.0)*k5);
		k7 = this->step_size*this->Eqn.differential_equation(this->Eqn.x.at(position - 1) + this->step_size, this->Eqn.y.at(position - 1) + (35.0/384.0)*k1 + (500.0/1113.0)*k3 + (125.0/192.0)*k4 - (2187.0/6784.0)*k5+(11.0/84.0)*k6);
		temp_y = this->Eqn.y.at(position - 1) + (35.0/384.0)*k1 + (500.0/1113.0)*k3 + (125.0/192.0)*k4 - (2187.0/6784.0)*k5 + (11.0/84.0)*k6;
		cout << "X: " << this->Eqn.x.back() << "    " << "Y: " << this->Eqn.y.back() << endl;

		if (flag == false)
		{
			temp_x += step_size;
			this->Eqn.x.push_back(temp_x);
			CalculateStepSize(temp_y, temp_x);
			this->Eqn.y.push_back(temp_y);
		}
		else
		{
			this->Eqn.y.pop_back();
			this->Eqn.y.push_back(temp_y);
			flag = false;
			position++;
		}
	}
	return 0;
}
void RKDP::CalculateStepSize(double t_y, double t_x)
{
	double s, temp;
	temp = t_y + (5179.0 / 57600.0)*k1 + (7571.0 / 16695.0)*k3 + (393.0 / 640.0)*k4 - (92097.0 / 39200.0)*k5 + (1.0 / 40.0)*k7;
	s = pow(step_size / (2.0 * abs((t_x - t_y))),1.0/5.0);
	if ((s*this->Eqn.x.back() < 1.0) && (s*step_size > 0.001))
		step_size = s*step_size;

	flag = 1;
}
void RKDP::SetIC()
{
	cout << "Garbage check backburner" << endl;
}













// Backburner
// Make a new derived Algorithm derived class to provide an interface for adaptive methods
// -> virtual float CalculateStepSize()=0 
// Look at moving most of derived class constructors to Equation constructor
// Get rid of SetIC()
// ....And FillX()
// Template non adaptive methods
// PUT TARGET_X IN X????

/*
Sources

http://depa.fquim.unam.mx/amyd/archivero/DormandPrince_19856.pdf

http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html




*/
