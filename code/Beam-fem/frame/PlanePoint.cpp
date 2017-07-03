//#include "StdAfx.h"
#include "PlanePoint.h"

CPlanePoint::CPlanePoint(void)
{
}
CPlanePoint::CPlanePoint(const int &pointNO,
						 const double &x,
						 const double &y
						 ):
m_pointNO(pointNO),
m_x(x),
m_y(y)
{
}

CPlanePoint::~CPlanePoint(void)
{
}
