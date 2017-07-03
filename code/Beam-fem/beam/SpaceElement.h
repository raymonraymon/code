#pragma once

#include "SpacePoint.h"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>      

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <math.h>
#include <map>
class CSpaceElement
{
public:
	CSpaceElement(void);
	CSpaceElement(const int		&No,
				  const int		&left,
				  const int		&right,
				  const double	&E,
				  const double	&G,
				  const double	&A,
				  const double	&Iy,
				  const double	&Iz,
				  const double	&J/*,
				  double* distribute,
				  double* global_distribute*/);
	~CSpaceElement(void);
	
	boost::numeric::ublas::matrix<double> get_stiffmatrix();   //得到局部坐标系下的单元刚度矩阵	
	void init(std::map<int,CSpacePoint> *pt);       //对单元赋值

	//boost::numeric::ublas::matrix<double> elementRHS();  //单元坐标系局部荷载生成的局部右端
	


	

public:
	int m_elementNO;   //单元编号
	int m_leftpoint;   //左节点编号
	int m_rightpoint;  //右节点编号

	double m_element_length ;   //为了算均匀荷载的积分
	boost::numeric::ublas::matrix<double> m_transfermatrix; //桥博转换矩阵

	boost::numeric::ublas::matrix<double> m_stiffmatrix; //局部坐标系下的单元刚度矩阵	

	boost::numeric::ublas::matrix<double> m_transfermatrix1;/*(12,12)*/
	//Kattan书上所用的转换矩阵

	
    //材料和几何属性
	double m_E;
	double m_G;
	double m_A;
	double m_Iy;
	double m_Iz;
	double m_J;

	//double m_distribute[3];//加载这个单元局部坐标系上的均匀荷载分三个方向
	//未来考虑扩展成+三个方向弯矩
	double m_distribute[6];
	//未来考虑扩展成整体坐标系下也可以。
	double m_global_distribute[6];

	std::map<int,CSpacePoint> *m_points;
	

};
