#include "stdafx.h"
#include "TargetFun.h"


CTargetFun::CTargetFun()
{
}


CTargetFun::~CTargetFun()
{
}

/// <summary>
/// 同伦函数赋初值.
/// </summary>
/// <param name="initpoint">初值.</param>
void CTargetFun::Setinit(ublas::matrix<double>& initpoint)
{
	m_initpoint = initpoint;
}