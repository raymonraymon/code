#pragma once
#include "Displace.h"
#include "SpaceElement.h"
#include "SpaceForce.h"
#include "Distribute_force.h"
#include "SpacePoint.h"

#include <map>
#include <vector>
#include <iostream>

class CStructer
{
public:
	CStructer(void);
	CStructer(const std::map<int,CSpacePoint>	&m_points,
		const std::map<int,CSpaceElement>		&m_elements,
		const std::vector<CSpaceForce>			&v_forces,
		const std::vector<CDistribute_force>	&v_distribute_forces,
		const std::map<int,CDisplace>			&m_displaces);
	~CStructer(void);

	bool solve();
	void show_result();

private:
	bool								  local_vector_generate();
	boost::numeric::ublas::matrix<double> Assemble();
	boost::numeric::ublas::matrix<double> RHS_generate();
	boost::numeric::ublas::matrix<double> RHS_generate1();  //���Ӿ��Ⱥ��ص��Ҷ�����㺯��
	boost::numeric::ublas::matrix<double> Dirichlet_generate();

public:
	std::map<int,CSpacePoint>      m_m_points;
	std::map<int,CSpaceElement>	   m_m_elements;
	std::vector<CSpaceForce>       m_v_forces;
	std::vector<CDistribute_force> m_v_distribute_forces;
	std::map<int,CDisplace>		   m_m_displaces;
	std::map<int,CDisplace>		   m_m_result;

	std::map<int,std::vector<int>> m_point_local;//�ڵ㶨λ�������� 
	std::map<int,std::vector<int>> m_element_local;//��Ԫ��λ����

	/********************************************************/
	int							   m_Free_degree;   //���ɽڵ�ĸ���
	int							   m_Dirichlet_degree;   //ǿ��λ�ƽڵ����
	int							   m_total_degree; 

};
