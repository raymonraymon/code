#pragma once

class CPlanePoint
{
public:
	CPlanePoint(void);
	CPlanePoint(const int &pointNO,const double &x,const double &y);
	~CPlanePoint(void);
public:
	double m_x;		 //������
	double m_y;		 //������
	int m_pointNO;   //�ڵ���
};
