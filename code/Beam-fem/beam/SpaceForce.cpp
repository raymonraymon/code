//#include "StdAfx.h"
#include "SpaceForce.h"

CSpaceForce::CSpaceForce(void)
{
}
//节点力类
CSpaceForce::CSpaceForce(const int &ForceNo, double* array)
{
	m_ForceNo = ForceNo;
	for (int i=0; i<6; ++i)
	{
		m_array[i] = array[i];
	}
}

CSpaceForce::~CSpaceForce(void)
{
}
