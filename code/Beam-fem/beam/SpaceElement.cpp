//#include "StdAfx.h"
#include "SpaceElement.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace boost::numeric::ublas;

//保留空的构造函数以便编译
CSpaceElement::CSpaceElement(void)
{
}
//实际赋初值
CSpaceElement::	CSpaceElement(const int		&No,
							  const int		&left,
							  const int		&right,
							  const double	&E,
							  const double	&G,
							  const double	&A,
							  const double	&Iy,
							  const double	&Iz,
							  const double	&J/*,
							  double* distribute,
							  double* global_distribute*/)
{
	m_elementNO = No;
	m_leftpoint = left;
	m_rightpoint= right;
	m_E         = E;
	m_G         = G;
	m_A         = A;
	m_Iy        = Iy;
	m_Iz        = Iz;
	m_J         = J;
	//for (int i=0; i<6; ++i)
	//{
	//	m_distribute[i] = distribute[i];
	//	m_global_distribute[i] = global_distribute[i];
	//}

}

CSpaceElement::~CSpaceElement(void)
{
}
//double CSpaceElement::element_length(void)const
//{
//	std::map<int,CSpacePoint>::const_iterator pos1,pos2; //
//	pos1 = m_points->find(m_leftpoint);
//	pos2 = m_points->find(m_rightpoint);
//
//	assert( pos1 != m_points->end() && pos2 != m_points->end() );
//
//	double len = sqrt(
//		pow(pos1->second.m_x - pos2->second.m_x,2) +
//		pow(pos1->second.m_y - pos2->second.m_y,2) +
//		pow(pos1->second.m_z - pos2->second.m_z,2)
//		);
//	return len;
//}

void CSpaceElement::init(std::map<int,CSpacePoint> *pt)
{
	m_points = pt;


	std::map<int,CSpacePoint>::const_iterator pos1,pos2; //
	pos1 = m_points->find(m_leftpoint);
	pos2 = m_points->find(m_rightpoint);
	assert( pos1 != m_points->end() && pos2 != m_points->end() );
	double x1 = pos1->second.m_x;
	double y1 = pos1->second.m_y;
	double z1 = pos1->second.m_z;
	double x2 = pos2->second.m_x;
	double y2 = pos2->second.m_y;
	double z2 = pos2->second.m_z;

	double len = sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) +	pow(z1 - z2,2));
	m_element_length = len;

	cout<<"length of the element is :"<<len<<endl;

	matrix<double> trans_m (12, 12);  
	trans_m.clear();

	matrix<double> L (3, 3);  
	L.clear();

	if (x1 == x2 && y1 == y2)
	{
		if (z2 > z1)
		{
			L(0,2) = 1;
			L(1,0) = 1;
			L(2,1) = 1;
		} 
		else
		{
			L(0,2) = -1;
			L(1,0) = 1;
			L(2,1) = -1;
		}
	} 
	else
	{
		double length =	m_element_length;//sqrt( pow(x2 - x1,2) + pow(y2 - y1,2) + pow(z2 - z1,2));
		assert(length > 0); 
		//文档中的e1
		double CXx = (x2 - x1)/length;
		double CYx = (y2 - y1)/length;
		double CZx = (z2 - z1)/length;
		//文档中的e3
		double D = sqrt(pow(CXx,2)+pow(CYx,2));
		double CXz = CYx/D;
		double CYz = -1*CXx/D;
		double CZz = 0;
		//文档中的e2
		double CXy = -1*CXx*CZx/D; /**之前这里算错了**/
		double CYy = -1*CYx*CZx/D;
		double CZy = D;


		L(0,0) = CXx;
		L(0,1) = CYx;
		L(0,2) = CZx;
		L(1,0) = CXy;
		L(1,1) = CYy;
		L(1,2) = CZy;
		L(2,0) = CXz;
		L(2,1) = CYz;
		L(2,2) = CZz;


	}

	cout<<"校验转换矩阵是否为单位阵："<<prod(trans(L),L)<<endl;


	for (int i=0;i<3;++i)		
	{
		for (int j=0;j<3;++j)
		{
			trans_m(i,j)     = L(i,j);
			trans_m(i+3,j+3) = L(i,j);
			trans_m(i+6,j+6) = L(i,j);
			trans_m(i+9,j+9) = L(i,j);
		}
	}
	m_transfermatrix = trans_m;

	ofstream outfile("result.csv", ios::out);

	if (!outfile)
	{
		cerr << "open error" << endl;
		exit(1);
	}
	cout<<setprecision(15);
	cout<<"the element's transfer matrix is:"<<endl;
	cout << subrange(trans_m,0,12,0,12) << endl;
	outfile << "the element's transfer matrix is:" << endl;
	outfile << subrange(trans_m, 0, 12, 0, 12) << endl;
	//cout<<setprecision(6);

	symmetric_matrix<double> m (12, 12);
	m.clear();

	double length = m_element_length;

	assert(length > 0); 
	double w1 =	m_E * m_A / length;
	double w2 = 12*m_E*m_Iz/pow(length,3);
	double w3 = 6*m_E*m_Iz/pow(length,2);
	double w4 = 4*m_E*m_Iz/length;
	double w5 = 2*m_E*m_Iz/length;
	double w6 = 12*m_E*m_Iy/pow(length,3);
	double w7 = 6*m_E*m_Iy/pow(length,2);
	double w8 = 4*m_E*m_Iy/length;
	double w9 = 2*m_E*m_Iy/length;
	double w10= m_G*m_J/length;

	m(0,0) = w1;
	m(0,6) =-1*w1;

	m(1,1) =w2;
	m(1,5) =w3;
	m(1,7) =-1*w2;
	m(1,11)=w3;

	m(2,2)  = w6;
	m(2,4)  = -1*w7;
	m(2,8)  = -1*w6;
	m(2,10) =-1*w7;

	m(3,3) = w10;
	m(3,9) = -1*w10;

	//m(4,2) = -1*w7;
	m(4,4) = w8;
	m(4,8) = w7;
	m(4,10) = w9;

	//m(5,1) = w3;
	m(5,5) = w4;
	m(5,7) =-1*w3;
	m(5,11) = w5;

	//m(6,0) = -1*w1;
	m(6,6) = w1;

	//m(7,1) =-1*w2;
	//m(7,5) =-1*w3;
	m(7,7) = w2;
	m(7,11) = -1*w3;

	//m(8,2) = -1*w6;
	//m(8,4) = w7;
	m(8,8) = w6;
	m(8,10) =w7;

	//m(9,3) = -1*w10;
	m(9,9) = w10;

	//m(10,2) = -1*w7;
	//m(10,4) = w9;
	//m(10,8) = w7;
	m(10,10) = w8;

	//m(11,1) = w3;
	//m(11,5) = w5;
	//m(11,7) = -1*w3;
	m(11,11) = w4;


	cout << "the element's local stiffmatrix is:" << endl;
	cout << subrange(m, 0, 12, 0, 12) << endl;
	outfile << "the element's local stiffmatrix is:" << endl;
	outfile << subrange(m, 0, 12, 0, 12) << endl;


	matrix<double> temp (12, 12) ;
	temp = prod(m , m_transfermatrix );


	m_stiffmatrix = prod (trans (m_transfermatrix) , temp);


	cout << "the element's stiffmatrix is:" << endl;
	//cout << subrange(m_stiffmatrix, 6, 12, 6, 12) << endl;

	outfile << "the element's stiffmatrix is:" << endl;
	//outfile << subrange(m_stiffmatrix, 6, 12, 6, 12) << endl;

	outfile.close();

}

matrix<double> CSpaceElement::get_stiffmatrix()
{
	return m_stiffmatrix;
}

//matrix<double> CSpaceElement::elementRHS()
//{
//	//matrix<double> tran_m (12, 12);
//	//tran_m = m_transfermatrix;
//
//	matrix<double> m (12, 1);
//	m.clear();
//
//	/************************************************************************/
//	/*                        把分布的右端统一变成矩阵来处理                */
//	/************************************************************************/
//
//	//把坐标转换矩阵乘上去即可.
//	matrix<double> k(6,1);k.clear();
//	matrix<double> k_local(6,1);k_local.clear();
//	matrix<double> k_global(6,1);k_global.clear();
//
//
//
//	for (int i=0;i<6;++i)
//	{
//		k_local(i,0)	= m_distribute[i];
//		k_global(i,0)	= m_global_distribute[i];
//	}
//	k = k_local + prod(subrange(m_transfermatrix,0,6,0,6),k_global);
//
//	double length = m_element_length;
//
//	assert( length > 0 ); 
//	//局部坐标系下的三个分布荷载	
//	matrix<double> m1 (12, 1);
//	m1.clear();
//	m1(0,0) = k(0,0)*length/2;
//	m1(6,0) = k(0,0)*length/2;
//
//	matrix<double> m2 (12, 1);
//	m2.clear();
//	m2(1,0)	= k(1,0)*length/2;
//	m2(5,0)	= k(1,0)*pow(length,2)/12;
//	m2(7,0)	= k(1,0)*length/2;
//	m2(11,0)= -1*k(1,0)*pow(length,2)/12;
//
//	matrix<double> m3 (12, 1);
//	m3.clear();
//	m3(2,0)	= k(2,0)*length/2;
//	m3(4,0)	= -1*k(2,0)*pow(length,2)/12;
//	m3(8,0)	= k(2,0)*length/2;
//	m3(10,0)= k(2,0)*pow(length,2)/12;
//	
//	//局部坐标系下的三个分布弯矩
//	matrix<double> m4 (12, 1);
//	m4.clear();
//	m4(3,0) = k(3,0)*length/2;
//	m4(9,0) = k(3,0)*length/2;
//
//	matrix<double> m5 (12, 1);
//	m5.clear();
//	m5(2,0)	= k(4,0);
//	m5(8,0) = -1*k(4,0);
//
//	matrix<double> m6 (12, 1);
//	m6.clear();
//	m6(1,0)	= -1*k(5,0);
//	m6(7,0) = k(5,0);
//
//	m = m1 + m2 + m3 + m4 + m5 + m6;
//
//	//for (int i=0;i<12;++i)
//	//{
//	//	m(i,0) = m1(i,0)+m2(i,0)+m3(i,0);
//	//}
//
//	//cout<<"the NO" <<m_elementNO<<"'s element's RHS is:"<<endl;
//	//cout << m << endl;
//	//cout << prod (trans (m_transfermatrix) , m) << endl;
//	return prod (trans (m_transfermatrix) , m);
//}

//matrix<double> CSpaceElement::transfermatrix1()const
//{
//	std::map<int,CSpacePoint>::const_iterator pos1,pos2; //
//	pos1=m_points->find(m_leftpoint);
//	pos2=m_points->find(m_rightpoint);
//
//	matrix<double> m (12, 12);  
//	m.clear();
//	
//	matrix<double> L (3, 3);  
//	L.clear();
//	double x1 = pos1->second.m_x;
//	double y1 = pos1->second.m_y;
//	double z1 = pos1->second.m_z;
//	double x2 = pos2->second.m_x;
//	double y2 = pos2->second.m_y;
//	double z2 = pos2->second.m_z;
//
//
//	if (abs(x1-x2)<1E-7 && abs(y1-y2)<1E-7)
//	{
//		if (z2 > z1)
//		{
//			L(0,2) = 1;
//			L(1,1) = 1;
//			L(2,0) = -1;
//		} 
//		else
//		{
//			L(0,2) = -1;
//			L(1,1) = 1;
//			L(2,0) = 1;
//		}
//	} 
//	else
//	{
//		assert( pos1 != m_points->end() && pos2 != m_points->end() );
//
//		double length =	sqrt( pow(x2 - x1,2) + pow(y2 - y1,2) + pow(z2 - z1,2));
//		assert(length > 0); 
//		double CXx = (x2 - x1)/length;
//		double CYx = (y2 - y1)/length;
//		double CZx = (z2 - z1)/length;
//		double D = sqrt(pow(CXx,2)+pow(CYx,2));
//		double CXy =-1 * CYx/D;
//		double CYy = CXx/D;
//		double CZy = 0;
//		double CXz = -1 * CXx*CZx/D;
//		double CYz = -1 * CYx*CZx/D;
//		double CZz = D;
//		L(0,0) = CXx;
//		L(0,1) = CYx;
//		L(0,2) = CZx;
//		L(1,0) = CXy;
//		L(1,1) = CYy;
//		L(1,2) = CZy;
//		L(2,0) = CXz;
//		L(2,1) = CYz;
//		L(2,2) = CZz;
//	}
//	for (int i=0;i<3;++i)		
//	{
//		for (int j=0;j<3;++j)
//		{
//			m(i,j)     = L(i,j);
//			m(i+3,j+3) = L(i,j);
//			m(i+6,j+6) = L(i,j);
//			m(i+9,j+9) = L(i,j);
//		}
//	}
//	return m;
//}

//matrix<double> CSpaceElement::transfermatrix()const  //桥博的转换关系写法
//{
//	std::map<int,CSpacePoint>::const_iterator pos1,pos2; //
//	pos1=m_points->find(m_leftpoint);
//	pos2=m_points->find(m_rightpoint);
//	assert( pos1 != m_points->end() && pos2 != m_points->end() );
//	matrix<double> m (12, 12);  
//	m.clear();
//
//	matrix<double> L (3, 3);  
//	L.clear();
//	double x1 = pos1->second.m_x;
//	double y1 = pos1->second.m_y;
//	double z1 = pos1->second.m_z;
//	double x2 = pos2->second.m_x;
//	double y2 = pos2->second.m_y;
//	double z2 = pos2->second.m_z;
//
//
//	if (x1 == x2 && y1 == y2)
//	{
//		if (z2 > z1)
//		{
//			L(0,2) = 1;
//			L(1,0) = 1;
//			L(2,1) = 1;
//		} 
//		else
//		{
//			L(0,2) = -1;
//			L(1,0) = 1;
//			L(2,1) = -1;
//		}
//	} 
//	else
//	{
//		double length =	element_length();//sqrt( pow(x2 - x1,2) + pow(y2 - y1,2) + pow(z2 - z1,2));
//		assert(length > 0); 
//		//文档中的e1
//		double CXx = (x2 - x1)/length;
//		double CYx = (y2 - y1)/length;
//		double CZx = (z2 - z1)/length;
//		//文档中的e3
//		double D = sqrt(pow(CXx,2)+pow(CYx,2));
//		double CXz = CYx/D;
//		double CYz = -1*CXx/D;
//		double CZz = 0;
//		//文档中的e2
//		double CXy = CYx*CZx/D;
//		double CYy = -1 * CYx*CZx/D;
//		double CZy = D;
//
//
//		L(0,0) = CXx;
//		L(0,1) = CYx;
//		L(0,2) = CZx;
//		L(1,0) = CXy;
//		L(1,1) = CYy;
//		L(1,2) = CZy;
//		L(2,0) = CXz;
//		L(2,1) = CYz;
//		L(2,2) = CZz;
//	}
//	for (int i=0;i<3;++i)		
//	{
//		for (int j=0;j<3;++j)
//		{
//			m(i,j)     = L(i,j);
//			m(i+3,j+3) = L(i,j);
//			m(i+6,j+6) = L(i,j);
//			m(i+9,j+9) = L(i,j);
//		}
//	}
//	return m;
//}