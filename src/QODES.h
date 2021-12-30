#pragma once

#include <iostream>
#include <vector>
#include <math.h>

using std::vector;
using std::cout;
using std::endl;

// Abstract base class that all solution methods inherit from
template <class T>
class Algorithm
{
public:
	virtual T Solve() = 0; // Main Algorithm
	Algorithm() {
		m_step_size = 0;
		m_target_x = 0;
	};
	void FillX();

	class Equation // Represents an equation. It contains the function iteself, its initial conditions and all the values
	{
	public:
		int dimensions;
		T m_initial_condition; // IC @ x=0
		vector<T> x, y; // Containers for x and y
		T(*differential_equation)(T, T);
	} m_equation;

protected:
	T m_step_size; // Step Size Used
	T m_target_x; // x Value of Interest
	T m_initial_condition;
};

template <class T>
class RungeKutta : public Algorithm<T>
{
public:
	T k1, k2, k3, k4;
};

template<class T>
void Algorithm<T>::FillX()
{
	for (double i = 0; i <= m_target_x; i += m_step_size)
	{
		m_equation.x.push_back(i);
	}
	m_equation.x.push_back(m_target_x);
}

// Definitions for the RK 3/8 method
template <class T>
class RK38 : public RungeKutta<T>
{
public:
	RK38(T step, T final_x, T IC);
	T Solve();
};

template <class T>
RK38<T>::RK38(T step, T final_x, T IC)
{
	m_step_size = step;
	m_target_x = final_x;
	m_initial_condition = IC;
};

template <class T>
T RK38<T>::Solve()
{
	m_equation.x.push_back(0);
	m_equation.y.push_back(m_initial_condition);
	for (auto i = 1; m_equation.x.back() < m_target_x; ++i) {
		m_equation.x.push_back(m_equation.x.back() + m_step_size);
		k1 = m_step_size * m_equation.differential_equation(m_equation.x.at(i - 1), m_equation.y.at(i - 1));
		k2 = m_step_size * m_equation.differential_equation(m_equation.x.at(i - 1) + (m_step_size / 3.0), m_equation.y.at(i - 1) + (k1 / 3.0));
		k3 = m_step_size * m_equation.differential_equation(m_equation.x.at(i - 1) + (2.0 * m_step_size) / 3.0, m_equation.y.at(i - 1) - (k1 / 3.0) + k2);
		k4 = m_step_size * m_equation.differential_equation(m_equation.x.at(i - 1) + m_step_size, m_equation.y.at(i - 1) + k1 - k2 + k3);
		m_equation.y.push_back(m_equation.y.at(i - 1) + (1.0 / 8.0) * (k1 + 3 * k2 + 3 * k3 + k4));
		cout << "X: " << m_equation.x.back() << "    " << "Y: " << m_equation.y.back() << endl;
	}

	return 0;
}

// Definitions for the Class RK Method
template <class T>
class RK4 : public RungeKutta<T>
{
public:
	RK4(T a, T b, T c);
	T Solve();
};

template <class T>
RK4<T>::RK4(T step, T final_x, T IC)
{
	m_step_size = step;
	m_target_x = final_x;
	m_initial_condition = IC;
}

template <class T>
T RK4<T>::Solve()
{
	m_equation.x.push_back(0);
	m_equation.y.push_back(m_initial_condition);
	T result;
	for (auto position = 1; m_equation.x.back() < m_target_x; ++position) {
		m_equation.x.push_back(m_equation.x.back() + m_step_size);
		k1 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1), 2); //m_equation.differential_equation(m_equation.x.at(position), m_equation.y.at(position-1));
		k2 = m_step_size * m_equation.differential_equation(m_equation.x.at(position) + 0.5 * m_step_size, m_equation.y.at(position - 1) + 0.5 * k1);
		k3 = m_step_size * m_equation.differential_equation(m_equation.x.at(position) + 0.5 * m_step_size, m_equation.y.at(position - 1) + 0.5 * k2);
		k4 = m_step_size * m_equation.differential_equation(m_equation.x.at(position) + m_step_size, m_equation.y.at(position - 1) + k3);
		result = m_equation.y.back() + (1.0 / 6.0) * k1 + (1.0 / 3.0) * k2 + (1.0 / 3.0) * k3 + (1.0 / 6.0) * k4;
		m_equation.y.push_back(result);
		cout << "X: " << m_equation.x.back() << "    " << "Y: " << m_equation.y.back() << endl;
	}
	return 0;
}

// Definitions for the forward Euler Methods
template <class T>
class ForwardEuler : public RungeKutta<T>
{
public:
	ForwardEuler(T, T, T);
	T Solve();

};

template <class T>
ForwardEuler<T>::ForwardEuler(T step, T final_x, T IC)
{
	m_step_size = step;
	m_target_x = final_x;
	m_initial_condition = IC;
}

template <class T>
T ForwardEuler<T>::Solve()
{
	cout << "Solving..." << endl;
	m_equation.x.push_back(0);
	m_equation.y.push_back(m_equation.m_initial_condition);
	while (m_equation.x.back() < m_target_x) {
		m_equation.x.push_back(m_equation.x.back() + m_step_size);
		m_equation.y.push_back(m_equation.y.back() + m_step_size * m_equation.differential_equation(m_equation.x.back(), m_equation.y.back()));
		cout << "X: " << m_equation.x.back() << "    " << "Y: " << m_equation.y.back() << endl;
	}
	m_equation.x.push_back(m_target_x);
	m_equation.y.push_back(m_equation.y.back() + m_step_size * m_equation.differential_equation(m_equation.x.back(), m_equation.y.back()));
	return 0;
}


// Runge–Kutta–Fehlberg 
template <class T>
class RK45 : public RungeKutta<T>
{
public:
	RK45(T step, T final_x, T IC, T error_tolerance);
	T Solve();
	T ComputeDerivative(int);

protected:
	T k5, k6;
	bool m_stepsize_set;
	T m_error_tolerance;
};

template <class T>
RK45<T>::RK45(T step, T final_x, T IC, T error_tolerance)
{
	m_step_size = step;
	m_target_x = final_x;
	m_initial_condition = IC;
	m_error_tolerance = error_tolerance;
	m_equation.y.push_back(m_initial_condition);
	m_equation.x.push_back(0);
};


template <class T>
T RK45<T>::Solve()
{
	T new_y = 0;
	T new_x = 0;
	std::size_t position = 1;
	for (auto i = 1; m_equation.x.back() < m_target_x + m_step_size; ++i) {
		new_x += m_step_size;
		m_equation.x.push_back(new_x);
		k1 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1), m_equation.y.at(position - 1));
		k2 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + (m_step_size / 4), m_equation.y.at(position - 1) + (k1 / 4.0));
		k3 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + 3.0 * (m_step_size / 8.0), m_equation.y.at(position - 1) + (3.0 * k1 / 32.0) + k2 * (9.0 / 32.0));
		k4 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size, m_equation.y.at(position - 1) + (1932.0 / 2197.0) * k1 - (7200 / 2197) * k2 + (7296 / 2197) * k3);
		k5 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size, m_equation.y.at(position - 1) + (439 / 216) * k1 - (8.0 * k2) + (3680.0 / 513.0) * k3 - (845 / 4104) * k4);
		k6 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + 0.5 * m_step_size, m_equation.y.at(position - 1) - (8.0 / 27) * k1 + 2.0 * k2 - (3544 / 2565) * k3 + (1859 / 4104) * k4 - (11.0 / 40.0) * k5);
		cout << "X: " << m_equation.x.back() << "    " << "Y: " << m_equation.y.back() << endl;
		auto error = abs((1 / 360) * k1 - (128 / 4275) * k3 - (2197 / 75240) * k4 + (1 / 50) * k5 + (2 / 55) * k6);
		auto new_step_size = 0.9 * m_step_size * pow(m_error_tolerance / error, 1 / 5);
		if (error > m_error_tolerance) {
			m_step_size = new_step_size;
		}
		else {
			new_y = m_equation.y.at(position - 1) + (16.0 / 135.0) * k1 + (6656.0 / 12825.0) * k3 + (28561.0 / 56430.0) * k4 - (9.0 / 50.0) * k5 + (2.0 / 55.0) * k6;
			m_equation.y.push_back(new_y);
		}
		position++;
	}
	return 0;
}

// Runge Kutta Dormand Prince Definitions
template< class T>
class RKDP : public RungeKutta<T>
{
public:
	RKDP(T step, T final_x, T IC, T error_tolerance);
	T Solve();
protected:
	T k5, k6, k7;
	T m_error_tolerance;
};

template <class T>
RKDP<T>::RKDP(T step, T final_x, T IC, T error_tolerance)
{
	m_step_size = step;
	m_target_x = final_x;
	m_initial_condition = IC;
	m_error_tolerance = error_tolerance;
	k5 = 0;
	k6 = 0;
	k7 = 0;
	// Push back the initial conditions (x=0, y=IC)
	m_equation.y.push_back(m_initial_condition);
	m_equation.x.push_back(0);
};

template <class T>
T RKDP<T>::Solve()
{
	T new_y;
	T new_x= 0;
	std::size_t position = 1;
	for (auto i = 1; m_equation.x.back() <= m_target_x + m_step_size; ++i) {
		k1 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1), m_equation.y.at(position - 1));
		k2 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + (m_step_size / 5.0), m_equation.y.at(position - 1) + k1 / 5.0);
		k3 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + 3.0 * (m_step_size / 10.0), m_equation.y.at(position - 1) + (3.0 * k1 / 40.0) + k2 * (9.0 / 40.0));
		k4 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size * 4.0 / 5.0, m_equation.y.at(position - 1) + (44.0 / 45.0) * k1 - (56.0 / 15.0) * k2 + (32.0 / 9.0) * k3);
		k5 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size * 8.0 / 9.0, m_equation.y.at(position - 1) + (19372.0 / 6561.0) * k1 - (25360.0 / 2187.0) * k2 + (6448.0 / 6561.0) * k3 - (212.0 / 729.0) * k4);
		k6 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size, m_equation.y.at(position - 1) + (9017.0 / 3168.0) * k1 - (355.0 / 33.0) * k2 - (46732.0 / 5247.0) * k3 + (49.0 / 176.0) * k4 - (5103.0 / 18656.0) * k5);
		k7 = m_step_size * m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size, m_equation.y.at(position - 1) + (35.0 / 384.0) * k1 + (500.0 / 1113.0) * k3 + (125.0 / 192.0) * k4 - (2187.0 / 6784.0) * k5 + (11.0 / 84.0) * k6);
		new_y = m_equation.y.at(position - 1) + (35.0 / 384.0) * k1 + (500.0 / 1113.0) * k3 + (125.0 / 192.0) * k4 - (2187.0 / 6784.0) * k5 + (11.0 / 84.0) * k6;
		T s;
		// Order 5 calculation
		auto new_y_5 = new_y + (5179.0 / 57600.0) * k1 + (7571.0 / 16695.0) * k3 + (393.0 / 640.0) * k4 - (92097.0 / 39200.0) * k5 + (1.0 / 40.0) * k7;
		auto error = abs(new_y_5 - new_y);
		s = pow((m_error_tolerance * m_step_size) / (2.0 * abs(new_y_5 - new_y)), 1 / 5);
		m_step_size = m_step_size * s;
		new_x += m_step_size;
		m_equation.x.push_back(new_x);
		m_equation.y.push_back(new_y);
		position++;
		cout << "X: " << m_equation.x.back() << "    " << "Y: " << m_equation.y.back() << endl;
	}
	return 0;
}
