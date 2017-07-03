#pragma once
#include "PlaneDisplacement.h"
#include "PlaneElement.h"
#include "PlanePoint.h"
#include "PlaneForce.h"
#include <map>
#include <vector>
#include <iostream>



class CPlaneStructer
{
public:
	CPlaneStructer(const std::map<int,CPlanePoint> &point_map,
		const std::map<int,CPlaneElement> &element_map,
		const std::vector<CPlaneForce> &force_vector,
		const std::map<int,CPlaneDisplacement> &displacement_map);
	~CPlaneStructer(void);

	bool solve();
	void show_result();

private:
	boost::numeric::ublas::matrix<double> Assemble();
	bool local_vector_generate();
	boost::numeric::ublas::matrix<double> RHS_generate();
	boost::numeric::ublas::matrix<double> Dirichlet_generate();

public:
	std::map<int,CPlanePoint> m_point_map;   //�ڵ����ݿ�
	std::map<int,CPlaneElement> m_element_map; //��Ԫ���ݿ�
	std::vector<CPlaneForce> m_force_vector;    //����
	std::map<int,CPlaneDisplacement> m_displacement_map;//ǿ��λ��

	std::map<int,CPlaneDisplacement> m_displacement_map_result; //����
	/*-----------------------------------------------------*/

	std::map<int,std::vector<int>> m_point_local;//�ڵ㶨λ�������� 
	std::map<int,std::vector<int>> m_element_local;//��Ԫ��λ����

	/********************************************************/
	int m_Free_degree;   //���ɽڵ�ĸ���
	int m_Dirichlet_degree;   //ǿ��λ�ƽڵ����
	int m_total_degree; 

};
