#pragma once

class CPlanePoint
{
public:
	CPlanePoint(void);
	CPlanePoint(const int &pointNO,const double &x,const double &y);
	~CPlanePoint(void);
public:
	double m_x;		 //横坐标
	double m_y;		 //纵坐标
	int m_pointNO;   //节点编号
};
