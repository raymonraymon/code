/**
* @file 自定义分布荷载类的函数实现
*  
* @author 陈瑞山Chen ruishan
*/


//#include "StdAfx.h"
#include "Distribute_force.h"

/**
*@brief 保留空的构造函数以避免编译出错keep the empty constructor,avoiding debug error.
*
*/
CDistribute_force::CDistribute_force(void)
{
}

/**
*@brief 各自定义
*@param element_no 荷载所在的单元编号
*@param type g表示全局坐标系,l表示局部坐标系
*/
CDistribute_force::CDistribute_force(const int &element_no, distribute_type type,double* distribute)
{
	m_element_no  = element_no;
	m_type		  =	type;
	for (int i=0; i<6; ++i)
	{
		m_distribute[i] = distribute[i];
	}
}

CDistribute_force::~CDistribute_force(void)
{
}

