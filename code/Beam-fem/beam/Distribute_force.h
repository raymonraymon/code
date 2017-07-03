/**
* @file �ռ�����Ԫ���Զ���ֲ�������Finite Element Method for Space beam!
* 
* @author ����ɽChen ruishan 
*/
#pragma once

/**
*@brief �����ֲ������࣬����һ����Ԫ�϶��������أ������û��Լ����ܼ���
*
*/

class CDistribute_force
{

public: 
	enum distribute_type {L,G}; ///<����ö������

public:
	CDistribute_force(void);
	CDistribute_force(const int &element_no, distribute_type type, double* distribute);
	~CDistribute_force(void);
	
public:
	int m_element_no;///<�������ڵ�Ԫ 
	distribute_type m_type;  ///<��������������ϵ 
	double m_distribute[6];  ///<����=3�������أ�3����غ���
};
