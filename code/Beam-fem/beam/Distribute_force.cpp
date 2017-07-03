/**
* @file �Զ���ֲ�������ĺ���ʵ��
*  
* @author ����ɽChen ruishan
*/


//#include "StdAfx.h"
#include "Distribute_force.h"

/**
*@brief �����յĹ��캯���Ա���������keep the empty constructor,avoiding debug error.
*
*/
CDistribute_force::CDistribute_force(void)
{
}

/**
*@brief ���Զ���
*@param element_no �������ڵĵ�Ԫ���
*@param type g��ʾȫ������ϵ,l��ʾ�ֲ�����ϵ
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

