#pragma once

class CSpaceForce
{
public:
	CSpaceForce(void);
	CSpaceForce(const int &ForceNo, double* array);
	~CSpaceForce(void);
public:
	int m_ForceNo;
	double m_array[6];
};
