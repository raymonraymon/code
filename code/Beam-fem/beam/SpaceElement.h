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
	
	boost::numeric::ublas::matrix<double> get_stiffmatrix();   //�õ��ֲ�����ϵ�µĵ�Ԫ�նȾ���	
	void init(std::map<int,CSpacePoint> *pt);       //�Ե�Ԫ��ֵ

	//boost::numeric::ublas::matrix<double> elementRHS();  //��Ԫ����ϵ�ֲ��������ɵľֲ��Ҷ�
	


	

public:
	int m_elementNO;   //��Ԫ���
	int m_leftpoint;   //��ڵ���
	int m_rightpoint;  //�ҽڵ���

	double m_element_length ;   //Ϊ������Ⱥ��صĻ���
	boost::numeric::ublas::matrix<double> m_transfermatrix; //�Ų�ת������

	boost::numeric::ublas::matrix<double> m_stiffmatrix; //�ֲ�����ϵ�µĵ�Ԫ�նȾ���	

	boost::numeric::ublas::matrix<double> m_transfermatrix1;/*(12,12)*/
	//Kattan�������õ�ת������

	
    //���Ϻͼ�������
	double m_E;
	double m_G;
	double m_A;
	double m_Iy;
	double m_Iz;
	double m_J;

	//double m_distribute[3];//���������Ԫ�ֲ�����ϵ�ϵľ��Ⱥ��ط���������
	//δ��������չ��+�����������
	double m_distribute[6];
	//δ��������չ����������ϵ��Ҳ���ԡ�
	double m_global_distribute[6];

	std::map<int,CSpacePoint> *m_points;
	

};
