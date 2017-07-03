#pragma once

class CPlaneForce
{
public:
	CPlaneForce(const int &Force_NO, double* Force_array);
	~CPlaneForce(void);

public:
	int	m_Force_NO;
	double m_Force_array[3];
};
