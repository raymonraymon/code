#pragma once

class CPlaneDisplacement
{
public:
	CPlaneDisplacement(const int &displacement_NO,double* displacement_array,bool* displacement_bool );
	~CPlaneDisplacement(void);
	void show(void);
public:
	int m_displacement_NO;
	double m_displacement_array[3];
	bool m_displacement_bool[3];
};
