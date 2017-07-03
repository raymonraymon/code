#pragma once
#include "PlanePoint.h"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>                                
#include <math.h>
#include <map>

class CPlaneElement
{
public:
	CPlaneElement(void);
	CPlaneElement(const int &elementNO,
				  const int &leftpoint,
				  const int &rightpoint,
				  const double &elasticity,
				  const double &area,
		          const double &inertia);
	~CPlaneElement(void);
	boost::numeric::ublas::matrix<double> stiffmatrix()const;   //局部坐标系下的单元刚度矩阵	
	void setpoint(std::map<int,CPlanePoint> *pt);       //对单元赋值
private:

	boost::numeric::ublas::matrix<double> transfermatrix()const; //转换矩阵

public:
	int m_leftpoint;   //左节点编号
	int m_rightpoint;  //右节点编号
	int m_elementNO;   //单元编号
	
	//double m_theta;   //单元转角
	double m_elasticity;//弹性模量
	double m_area;      //截面面积
	double m_inertia;	//惯性矩

	std::map<int,CPlanePoint> *m_points;
};
