#include "StdAfx.h"
#include "PlaneForce.h"

#include "StdAfx.h"
CPlaneForce::CPlaneForce(const int &Force_NO, double* Force_array)
{
	m_Force_NO = Force_NO;
	for (int i=0;i<3; ++i)
	{
		m_Force_array[i] = Force_array[i];
	}
}

CPlaneForce::~CPlaneForce(void)
{
}
