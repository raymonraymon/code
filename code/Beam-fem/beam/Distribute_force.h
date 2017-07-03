/**
* @file 空间梁单元的自定义分布荷载类Finite Element Method for Space beam!
* 
* @author 陈瑞山Chen ruishan 
*/
#pragma once

/**
*@brief 创建分布荷载类，方便一个单元上定义多项荷载，避免用户自己汇总计算
*
*/

class CDistribute_force
{

public: 
	enum distribute_type {L,G}; ///<定义枚举类型

public:
	CDistribute_force(void);
	CDistribute_force(const int &element_no, distribute_type type, double* distribute);
	~CDistribute_force(void);
	
public:
	int m_element_no;///<荷载所在单元 
	distribute_type m_type;  ///<荷载所在坐标体系 
	double m_distribute[6];  ///<荷载=3个力荷载，3个弯矩荷载
};
