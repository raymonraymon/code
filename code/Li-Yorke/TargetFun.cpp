//#include "stdafx.h"
#include "TargetFun.h"


CTargetFun::CTargetFun()
{
}


CTargetFun::~CTargetFun()
{
}

/// <summary>
/// ͬ�׺�������ֵ.
/// </summary>
/// <param name="initpoint">��ֵ.</param>
void CTargetFun::Setinit(ublas::matrix<double>& initpoint)
{
	m_initpoint = initpoint;
}