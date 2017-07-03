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
	boost::numeric::ublas::matrix<double> stiffmatrix()const;   //�ֲ�����ϵ�µĵ�Ԫ�նȾ���	
	void setpoint(std::map<int,CPlanePoint> *pt);       //�Ե�Ԫ��ֵ
private:

	boost::numeric::ublas::matrix<double> transfermatrix()const; //ת������

public:
	int m_leftpoint;   //��ڵ���
	int m_rightpoint;  //�ҽڵ���
	int m_elementNO;   //��Ԫ���
	
	//double m_theta;   //��Ԫת��
	double m_elasticity;//����ģ��
	double m_area;      //�������
	double m_inertia;	//���Ծ�

	std::map<int,CPlanePoint> *m_points;
};
