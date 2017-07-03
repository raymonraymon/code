//#include "StdAfx.h"
#include "SpacePoint.h"

CSpacePoint::CSpacePoint(void)
{
}

CSpacePoint::CSpacePoint(const int &No,const double &x,const double &y,const double &z)
{
	m_No = No;
	m_x  = x;
	m_y  = y;
	m_z  = z;
}

CSpacePoint::~CSpacePoint(void)
{
}
