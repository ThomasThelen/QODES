#pragma once

#include <iostream>
#include <vector>
#include <math.h>

using std::vector;
using std::cout;
using std::endl;

// Abstract Base Class -> Is the interface for all solution methods ie all solution classes are derived from Algorithm
template <class T>
class Algorithm
{
public:
	virtual T Solve() = 0; // Main Algorithm
	Algorithm()
	{
		m_step_size= 0;
		m_target_x = 0;
	};
	virtual void SetIC() = 0;
	void FillX();

	class Equation // Represents an equation. It contains the function iteself, its initial conditions and all the values
	{
	public:
		int dimensions;
		T m_initial_condition; // IC @ x=0
		vector<T> x, y; // Containers for x and y
		T(*differential_equation)(T, T);
	}m_equation;

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

template <class T>
class Adaptive : public RungeKutta<T>
{
	virtual void CalculateStepSize(T, T)=0;
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
	void SetIC();
};

template <class T>
RK38<T>::RK38(T step, T final_x, T IC)
{
	m_step_size= step;
	m_target_x = final_x;
	m_initial_condition = IC;
};

template <class T>
void RK38<T>::SetIC()
{
	m_equation.x.push_back(0);
	m_equation.y.push_back(m_initial_condition);
}

template <class T>
T RK38<T>::Solve()
{
	SetIC();

	for (float i = 1; m_equation.x.back()<m_target_x; ++i)
	{
		m_equation.x.push_back(m_equation.x.back() + m_step_size);
		k1 = m_step_size*m_equation.differential_equation(m_equation.x.at(i - 1), m_equation.y.at(i - 1));
		k2 = m_step_size*m_equation.differential_equation(m_equation.x.at(i - 1) + (m_step_size/ 3.0), m_equation.y.at(i - 1) + (k1 / 3.0));
		k3 = m_step_size*m_equation.differential_equation(m_equation.x.at(i - 1) + (2.0 * m_step_size) / 3.0, m_equation.y.at(i - 1) - (k1 / 3.0) + k2);
		k4 = m_step_size*m_equation.differential_equation(m_equation.x.at(i - 1) + m_step_size, m_equation.y.at(i - 1) + k1 - k2 + k3);
		m_equation.y.push_back(m_equation.y.at(i - 1) + (1.0 / 8.0)*(k1 + 3 * k2 + 3 * k3 + k4));
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
	void SetIC();
};

template <class T>
RK4<T>::RK4(T step, T final_x, T IC)
{
	m_step_size= step;
	m_target_x = final_x;
	m_initial_condition = IC;
}

template <class T>
void RK4<T>::SetIC()
{
	m_equation.x.push_back(0);
	m_equation.y.push_back(m_initial_condition);
}

template <class T>
T RK4<T>::Solve()
{
	SetIC();
	T result;
	for (int position = 1; m_equation.x.back() < m_target_x; ++position)
	{
		m_equation.x.push_back(m_equation.x.back() + m_step_size);

		k1 = m_step_size*m_equation.differential_equation(m_equation.x.at(position-1), 2); //m_equation.differential_equation(m_equation.x.at(position), m_equation.y.at(position-1));
		k2 = m_step_size*m_equation.differential_equation(m_equation.x.at(position) + 0.5*m_step_size, m_equation.y.at(position-1) + 0.5*k1);
		k3 = m_step_size*m_equation.differential_equation(m_equation.x.at(position) + 0.5*m_step_size, m_equation.y.at(position-1) + 0.5*k2);
		k4 = m_step_size*m_equation.differential_equation(m_equation.x.at(position) + m_step_size, m_equation.y.at(position-1) + k3);
		result = m_equation.y.back() + (1.0 / 6.0)*k1 + (1.0 / 3.0)*k2 + (1.0 / 3.0)*k3 + (1.0 / 6.0)*k4;
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
	void SetIC();

};

template <class T>
ForwardEuler<T>::ForwardEuler(T step, T final_x, T IC)
{
	m_step_size= step;
	m_target_x = final_x;
	m_initial_condition = IC;
}

template <class T>
void ForwardEuler<T>::SetIC()
{
	m_equation.x.push_back(0);
	m_equation.y.push_back(m_equation.m_initial_condition);
}

template <class T>
T ForwardEuler<T>::Solve()
{
	cout << "Solving..." << endl;
	SetIC();
	while (m_equation.x.back() < m_target_x)
	{
		m_equation.x.push_back(m_equation.x.back() + m_step_size);
		m_equation.y.push_back(m_equation.y.back() + m_step_size*m_equation.differential_equation(m_equation.x.back(),m_equation.y.back()));
		cout << "X: " << m_equation.x.back() << "    " << "Y: " << m_equation.y.back() << endl;
	}
	m_equation.x.push_back(m_target_x);
	m_equation.y.push_back(m_equation.y.back() + m_step_size*m_equation.differential_equation(m_equation.x.back(), m_equation.y.back()));
	return 0;
}


// Runge–Kutta–Fehlberg 
template <class T>
class RK45 : public Adaptive<T>
{
public:
	RK45(T step, T final_x, T IC);
	T Solve();
	void ComputeDerivative(int, T&);
	void SetIC();

protected:
	void CalculateStepSize(T, T);
	T k5, k6;
	bool flag;
};

template <class T>
RK45<T>::RK45(T step, T final_x, T IC)
{
	m_step_size= step;
	m_target_x = final_x;
	m_initial_condition = IC;
	m_equation.y.push_back(m_initial_condition);
	m_equation.x.push_back(0);
};


template <class T>
void RK45<T>::ComputeDerivative(int position, T& temp_y)
{
	k1 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1), m_equation.y.at(position - 1));
	k2 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + (m_step_size/ 4), m_equation.y.at(position - 1) + (k1 / 4.0));
	k3 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + 3.0*(m_step_size/ 8.0), m_equation.y.at(position - 1) + (3.0*k1 / 32.0) + k2*(9.0 / 32.0));
	k4 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size, m_equation.y.at(position - 1) + (1932.0 / 2197.0)*k1 - (7200 / 2197)*k2 + (7296 / 2197)*k3);
	k5 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size, m_equation.y.at(position - 1) + (439 / 216)*k1 - (8.0*k2) + (3680.0 / 513.0)*k3 - (845 / 4104)*k4);
	k6 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + 0.5*m_step_size, m_equation.y.at(position - 1) - (8.0 / 27)*k1 + 2.0*k2 - (3544 / 2565)*k3 + (1859 / 4104)*k4 - (11.0 / 40.0)*k5);
	temp_y = m_equation.y.at(position - 1) + (16.0 / 135.0)*k1 + (6656.0 / 12825.0)*k3 + (28561.0 / 56430.0)*k4 - (9.0 / 50.0)*k5 + (2.0 / 55.0)*k6;


	cout << "X: " << m_equation.x.back() << "    " << "Y: " << m_equation.y.back() << endl;
}

template <class T>
T RK45<T>::Solve()
{
	T temp_y, temp_x = 0;
	T position = 1;
	for (double i = 1; m_equation.x.back()< m_target_x+m_step_size; ++i)
	{
		ComputeDerivative(position, temp_y);

		if (flag == false)
		{
			temp_x += m_step_size;
			m_equation.x.push_back(temp_x);
			CalculateStepSize(temp_y, temp_x);
			m_equation.y.push_back(temp_y);
		}
		else
		{
			m_equation.y.pop_back();
			m_equation.y.push_back(temp_y);
			flag = false;
			position++;
		}
	}
	m_equation.x.push_back(m_target_x);
	ComputeDerivative(position, temp_y);
	return 0;
}

template <class T>
void RK45<T>::CalculateStepSize(T t_y, T t_x)
{
	T s;
	s = 0.84*pow(m_step_size/ (2.0 * abs((t_x - t_y))), 0.25);
	if (s+m_equation.x.back() > m_target_x)
		m_step_size= s*0.001;
	flag = 1;
}

template<class T>
void RK45<T>::SetIC()
{
	cout << "Garbage" << endl;
}

// Runge Kutta Dormand Prince Definitions
template< class T>
class RKDP : public Adaptive<T>
{
public:
	RKDP(T step, T final_x, T IC);
	T Solve();
	void SetIC();
protected:
	void CalculateStepSize(T, T);
	T k5, k6, k7;
	bool flag;
};

template <class T>
RKDP<T>::RKDP(T step, T final_x, T IC)
{
	m_step_size= step;
	m_target_x = final_x;
	m_initial_condition = IC;
	m_equation.y.push_back(m_initial_condition);
	m_equation.x.push_back(0);
};

template <class T>
T RKDP<T>::Solve()
{
	T temp_y, temp_x = 0;
	T position = 1;
	for (double i = 1; m_equation.x.back()<= m_target_x+m_step_size; ++i)
	{

		k1 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1), m_equation.y.at(position - 1));
		k2 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + (m_step_size/ 5.0), m_equation.y.at(position - 1) +k1/5.0);
		k3 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + 3.0*(m_step_size/ 10.0), m_equation.y.at(position - 1) + (3.0*k1 / 40.0) + k2*(9.0 / 40.0));
		k4 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size*4.0/5.0, m_equation.y.at(position - 1) + (44.0 / 45.0)*k1 - (56.0/15.0)*k2 + (32.0/9.0)*k3);
		k5 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size*8.0/9.0, m_equation.y.at(position - 1) + (19372.0/6561.0)*k1 - (25360.0/2187.0)*k2 + (6448.0/6561.0)*k3 - (212.0/729.0)*k4);
		k6 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size, m_equation.y.at(position - 1) + (9017.0/3168.0)*k1 - (355.0/33.0)*k2 - (46732.0/5247.0)*k3 + (49.0/176.0)*k4 - (5103.0/18656.0)*k5);
		k7 = m_step_size*m_equation.differential_equation(m_equation.x.at(position - 1) + m_step_size, m_equation.y.at(position - 1) + (35.0/384.0)*k1 + (500.0/1113.0)*k3 + (125.0/192.0)*k4 - (2187.0/6784.0)*k5+(11.0/84.0)*k6);
		temp_y = m_equation.y.at(position - 1) + (35.0/384.0)*k1 + (500.0/1113.0)*k3 + (125.0/192.0)*k4 - (2187.0/6784.0)*k5 + (11.0/84.0)*k6;
		cout << "X: " << m_equation.x.back() << "    " << "Y: " << m_equation.y.back() << endl;

		if (flag == false)
		{
			temp_x += m_step_size;
			m_equation.x.push_back(temp_x);
			CalculateStepSize(temp_y, temp_x);
			m_equation.y.push_back(temp_y);
		}
		else
		{
			m_equation.y.pop_back();
			m_equation.y.push_back(temp_y);
			flag = false;
			position++;
		}
	}
	return 0;
}

template <class T>
void RKDP<T>::CalculateStepSize(T t_y, T t_x)
{
	T s, temp;
	temp = t_y + (5179.0 / 57600.0)*k1 + (7571.0 / 16695.0)*k3 + (393.0 / 640.0)*k4 - (92097.0 / 39200.0)*k5 + (1.0 / 40.0)*k7;
	s = pow(m_step_size/ (2.0 * abs((t_x - t_y))),1.0/5.0);
	if ((s*m_equation.x.back() < 1.0) && (s*m_step_size> 0.001))
	{
		m_step_size= s*m_step_size;
	}
	if (m_step_size+ m_equation.x.back() >= m_target_x)
	{
		m_step_size= m_step_size*.10;
	}
	flag = 1;
}

template <class T>
void RKDP<T>::SetIC()
{
	cout << "Garbage check backburner" << endl;
}
