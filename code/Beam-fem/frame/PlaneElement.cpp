#include "StdAfx.h"
#include "PlaneElement.h"
#include <iostream>

using namespace std;
using namespace boost::numeric::ublas;

CPlaneElement::CPlaneElement(void)
{
}

CPlaneElement::	CPlaneElement(const int &elementNO,
							  const int &leftpoint,
							  const int &rightpoint,
							  const double &elasticity,
							  const double &area,
							  const double &inertia):
m_elementNO(elementNO),
m_leftpoint(leftpoint),
m_rightpoint(rightpoint),
m_elasticity(elasticity),
m_area(area),
m_inertia(inertia)

{

}

CPlaneElement::~CPlaneElement(void)
{
}

void CPlaneElement::setpoint(std::map<int,CPlanePoint> *pt)
{
	m_points = pt;
}
matrix<double> CPlaneElement::stiffmatrix()const
{
	matrix<double> tran_m (6, 6);
	tran_m = transfermatrix();

	symmetric_matrix<double> m (6, 6);  //如何通过点集的map指针调用计算坐标差。
	//通过找到左节点键对应的值 find

	std::map<int,CPlanePoint>::const_iterator pos1,pos2; //
	pos1=m_points->find(m_leftpoint);
	pos2=m_points->find(m_rightpoint);

	assert( pos1 != m_points->end() && pos2 != m_points->end() );

	double length =
		sqrt(
		pow(pos1->second.m_x - pos2->second.m_x,2)+
		pow(pos1->second.m_y - pos2->second.m_y,2)
		)
		;
	assert(length > 0);  

	m(0,0) = 1*m_elasticity*m_area/length;
	m(0,1) = 0;
	m(0,2) = 0;
	m(0,3) = -1*m_elasticity*m_area/length;
	m(0,4) = 0;
	m(0,5) = 0;

	m(1,1) = 12*m_elasticity*m_inertia/pow(length,3);
	m(1,2) = 6*m_elasticity*m_inertia/pow(length,2);
	m(1,3) = 0;
	m(1,4) = -12*m_elasticity*m_inertia/pow(length,3);
	m(1,5) = 6*m_elasticity*m_inertia/pow(length,2);

	m(2,2) = 4*m_elasticity*m_inertia/length;
	m(2,3) = 0;
	m(2,4) = -6*m_elasticity*m_inertia/pow(length,2);
	m(2,5) = 2*m_elasticity*m_inertia/length;

	m(3,3) = 1*m_elasticity*m_area/length;
	m(3,4) = 0;
	m(3,5) = 0;

	m(4,4) = 12*m_elasticity*m_inertia/pow(length,3);
	m(4,5) = -6*m_elasticity*m_inertia/pow(length,2);

	m(5,5) = 4*m_elasticity*m_inertia/length;
	matrix<double> temp (6, 6) ;
	temp = prod(m , tran_m );
	//cout<<"the element's stiffmatrix is:"<<endl;
	//cout << prod (trans (tran_m) , temp) << endl;
	return prod (trans (tran_m) , temp);
}


matrix<double> CPlaneElement::transfermatrix()const
{
	std::map<int,CPlanePoint>::const_iterator pos1,pos2; //
	pos1=m_points->find(m_leftpoint);
	pos2=m_points->find(m_rightpoint);
	
	matrix<double> m (6, 6);  
	m.clear();
	double length =
		sqrt(
		pow(pos1->second.m_x - pos2->second.m_x,2)+
		pow(pos1->second.m_y - pos2->second.m_y,2)
		);
	if (abs(length)>1e-10)
	{
	double s=(pos2->second.m_y-pos1->second.m_y)/length;
	double c=(pos2->second.m_x-pos1->second.m_x)/length;

	/*pos1->second.m_x;
	pos1->second.m_y;
	pos2->second.m_x;
	pos2->second.m_y;*/
	m(0,0) = c;
	m(0,1) = s;
	m(1,0) = -1*s;
	m(1,1) = c;
	m(2,2) = 1;
	m(3,3) = c;
	m(3,4) = s;
	m(4,3) = -1*s;
	m(4,4) = c;
	m(5,5) = 1;
	}
	else
	{
	    cout<<"有节点重合了,请重新输入"<<endl;
	}
	//cout<<"the element's transfermatrix is:"<<endl;
	//cout << m << endl;
	return m;
	
}
