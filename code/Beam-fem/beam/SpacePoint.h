#pragma once

class CSpacePoint
{
public:
	CSpacePoint(void);
	CSpacePoint(const int &No,const double &x,const double &y,const double &z);
	~CSpacePoint(void);
public:
	int    m_No;
	double m_x;
	double m_y;
	double m_z;
};
