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
	std::map<int,CPlanePoint> m_point_map;   //节点数据库
	std::map<int,CPlaneElement> m_element_map; //单元数据库
	std::vector<CPlaneForce> m_force_vector;    //荷载
	std::map<int,CPlaneDisplacement> m_displacement_map;//强制位移

	std::map<int,CPlaneDisplacement> m_displacement_map_result; //最后解
	/*-----------------------------------------------------*/

	std::map<int,std::vector<int>> m_point_local;//节点定位向量集合 
	std::map<int,std::vector<int>> m_element_local;//单元定位向量

	/********************************************************/
	int m_Free_degree;   //自由节点的个数
	int m_Dirichlet_degree;   //强制位移节点个数
	int m_total_degree; 

};
