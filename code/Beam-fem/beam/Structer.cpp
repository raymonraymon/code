// ***********************************************************************
// Assembly         : SpaceFrame
// Author           : chenruishan
// Created          : 09-28-2014
//
// Last Modified By : chenruishan
// Last Modified On : 09-28-2014
// ***********************************************************************
// <copyright file="Structer.cpp" company="">
//     Copyright (c) . All rights reserved.
// </copyright>
// <summary>结构类
//
//
//
//
//
//
//
//
//
//
//
//
//</summary>
// ***********************************************************************
//#include "StdAfx.h"
#include "Structer.h"
//#include "SFrame\linalg.h"
using namespace std;
using namespace boost::numeric::ublas;


/// <summary>
/// Initializes a new instance of the <see cref="CStructer"/> class.
/// </summary>
CStructer::CStructer(void)
{
}

/// <summary>
/// Initializes a new instance of the <see cref="CStructer"/> class.
/// </summary>
/// <param name="m_points">The m_points.节点集合</param>
/// <param name="m_elements">The m_elements.单元集合</param>
/// <param name="v_forces">The v_forces.节点荷载向量</param>
/// <param name="v_distribute_forces">The v_distribute_forces.单元荷载向量</param>
/// <param name="m_displaces">The m_displaces.强制位移</param>
CStructer::CStructer(const std::map<int,CSpacePoint>		&m_points,
					 const std::map<int,CSpaceElement>		&m_elements,
					 const std::vector<CSpaceForce>			&v_forces,
					 const std::vector<CDistribute_force>	&v_distribute_forces,
					 const std::map<int,CDisplace>			&m_displaces)
{
	m_m_points				= m_points;
	m_m_elements			= m_elements;
	m_v_forces				= v_forces;
	m_v_distribute_forces	= v_distribute_forces;
	m_m_displaces			= m_displaces;
}

/// <summary>
/// Finalizes an instance of the <see cref="CStructer"/> class.
/// </summary>
CStructer::~CStructer(void)
{
}

/// <summary>
/// Local_vector_generates this instance.
///节点定位向量生成
///
///并利用节点定位向量生成单元定位向量
///
/// </summary>
/// <returns>bool.</returns>
bool CStructer::local_vector_generate()
{
	m_point_local.clear();///<保证可以重入
	m_element_local.clear();

	m_total_degree = 6 * (int)m_m_points.size();
	m_Dirichlet_degree = 0; 

	for (std::map<int,CDisplace>::iterator pos1 = m_m_displaces.begin(); pos1 != m_m_displaces.end(); ++pos1)
	{
		for (int i=0;i<6;++i)
		{
			if(pos1->second.m_arraybool[i])
			{
				++m_Dirichlet_degree;
			}
		}
	}
	cout<<"m_Dirichlet_degree:"<<m_Dirichlet_degree<<endl;
	/*m_Dirichlet_degree = 3 * (int)m_displacement_map.size();*/
	m_Free_degree = m_total_degree - m_Dirichlet_degree;


	int i=0,j=0;
	for (std::map<int,CSpacePoint>::iterator pos = m_m_points.begin();pos != m_m_points.end();++pos)
	{
		/**为了做位移是否是true的判断**/
		std::map<int,CDisplace>::iterator pos_bool ;
		pos_bool = m_m_displaces.find(pos->first);
		/*************/
		std::vector<int>  Position1(6);
		if (pos_bool == m_m_displaces.end())
		{
			for (int k=0;k<6;++k)
			{
				Position1[k]=j+1;
				++j;
			}
		} 
		else
		{
			for (int k=0;k<6;++k)
			{
				if (pos_bool->second.m_arraybool[k])
				{
					Position1[k]=-1*(i+1);
					++i;
				} 
				else
				{
					Position1[k]=j+1;
					++j;
				}

			}
		}
		m_point_local.insert(make_pair(pos->first,Position1));
	}
	/*节点定位向量结束*/


	for (std::map<int,CSpaceElement>::iterator pos = m_m_elements.begin();pos!= m_m_elements.end();++pos)
	{
		std::map<int,std::vector<int>>::iterator pos_point =m_point_local.find(pos->second.m_leftpoint);
		std::vector<int> positon_left;
		if (pos_point != m_point_local.end())
		{
			positon_left = pos_point->second;
		} 
		else
		{
			cout<<"input is wrong!"<<endl;
		}

		std::map<int,std::vector<int>>::iterator pos_point1 =m_point_local.find(pos->second.m_rightpoint);
		std::vector<int> positon_right;
		if (pos_point1 != m_point_local.end())
		{
			positon_right = pos_point1->second;
		} 
		else
		{
			cout<<"input is wrong!"<<endl;
		}


		std::vector<int> position_elem;

		for (int i=0;i < (int)positon_left.size(); ++i)
		{
			position_elem.push_back(positon_left[i]);
		}
		for (int i=0;i < (int)positon_right.size(); ++i)
		{
			position_elem.push_back(positon_right[i]);
		}				
		m_element_local.insert(make_pair(pos->first,position_elem));
	}

	return true;
}


/// <summary>
/// Assembles this instance.
/// </summary>
/// <returns>boost.numeric.ublas.matrix{double}.</returns>
matrix<double> CStructer::Assemble()
{
	for (std::map<int,CSpaceElement>::iterator pos = m_m_elements.begin();pos != m_m_elements.end();++pos)
	{
		pos->second.init(&m_m_points);
	}

	symmetric_matrix<double> m (m_total_degree, m_total_degree);   //对称矩阵
	m.clear();//保证重入
	int I,J;
	for (std::map<int,CSpaceElement>::iterator pos = m_m_elements.begin();pos != m_m_elements.end();++pos)
	{
		symmetric_matrix<double> mm = pos->second.get_stiffmatrix();
		for(unsigned int i=0;i < mm.size1();++i)
			for(unsigned int j=i;j < mm.size2();++j)
			{

				//m(?,?) += pos->second.stiffmatrix(i,j);
				if (m_element_local[pos->first][i] > 0)
				{
					I = m_element_local[pos->first][i]-1;				
				}
				else
				{
					I = -1*m_element_local[pos->first][i]-1+m_Free_degree;			
				}

				if (m_element_local[pos->first][j] > 0)
				{
					J = m_element_local[pos->first][j]-1;
				}
				else
				{
					J = -1*m_element_local[pos->first][j]-1+m_Free_degree;			
				}


				m(I,J) += mm(i,j);
			}
	}

	//cout<<"the structer's stiffmatrix is:"<<endl;
	//cout << m << endl;
	return m;
}

//matrix<double> CStructer::RHS_generate()
//{
//	//确定右边向量，求解方程的时候再截断成自由度那么长
//	matrix<double> RHS (m_total_degree,1);
//	RHS.clear();
//	int m = 0;
//	for (std::vector<CSpaceForce>::const_iterator pos = m_v_forces.begin(); pos != m_v_forces.end(); ++pos)
//	{
//
//		for (int i=0;i<6;++i)
//		{
//			if (m_point_local[pos->m_ForceNo][i]>0)
//			{
//				m = m_point_local[pos->m_ForceNo][i] - 1;
//			}
//			else
//			{
//				m = -1 * m_point_local[pos->m_ForceNo][i] - 1 + m_Free_degree;		
//			}
//
//			RHS(m,0) += pos->m_array[i];   //逐个向量相加
//
//		}
//	}
//	//集中荷载完成
//	//分布荷载完成
//	int I;
//	for (std::map<int,CSpaceElement>::iterator pos = m_m_elements.begin();pos != m_m_elements.end();++pos)
//	{
//		matrix<double> mm = pos->second.elementRHS();
//		for(unsigned int i=0;i < mm.size1();++i)
//			{
//				if (m_element_local[pos->first][i] > 0)
//				{
//					I = m_element_local[pos->first][i]-1;				
//				}
//				else
//				{
//					I = -1*m_element_local[pos->first][i]-1 + m_Free_degree;			
//				}
//				RHS(I,0) += mm(i,0);
//			}
//	}
//
//	//右边向量完成
//	cout<<"the Right hand side is:"<<endl;
//	cout << RHS << endl;
//	return RHS;
//}

/// <summary>
/// Rhes the s_generate1.
/// </summary>
/// <returns>boost.numeric.ublas.matrix{double}.</returns>
matrix<double> CStructer::RHS_generate1()
{
	//确定右边向量，求解方程的时候再截断成自由度那么长
	matrix<double> RHS (m_total_degree,1);
	RHS.clear();
	int m = 0;
	for (std::vector<CSpaceForce>::const_iterator pos = m_v_forces.begin(); pos != m_v_forces.end(); ++pos)
	{

		for (int i=0;i<6;++i)
		{
			if (m_point_local[pos->m_ForceNo][i]>0)
			{
				m = m_point_local[pos->m_ForceNo][i] - 1;
			}
			else
			{
				m = -1 * m_point_local[pos->m_ForceNo][i] - 1 + m_Free_degree;		
			}

			RHS(m,0) += pos->m_array[i];   //逐个向量相加

		}
	}

	cout<<"the Right  hand side for point force is:"<<endl;
	cout << RHS << endl;
	//集中荷载的右边向量完成
	/************************************************************************/
	/*				分布荷载  涉及到坐标转换矩阵吗？                        */
	/************************************************************************/
	int n1 = 0,n2 = 0;

	for (std::vector<CDistribute_force>::const_iterator distribute_pos 
		= m_v_distribute_forces.begin(); distribute_pos != m_v_distribute_forces.end();
		++distribute_pos)
	{

		/**为了做位移是否是true的判断**/
		std::map<int,CSpaceElement>::iterator pos_bool ;
		pos_bool = m_m_elements.find(distribute_pos->m_element_no);
		/*************/
		
		assert(pos_bool != m_m_elements.end());
		//CSpaceElement E1 = m_m_elements[distribute_pos->m_element_no];//危险做法

		matrix<double>	transm(12,12);
		transm = pos_bool->second.m_transfermatrix;

		matrix<double> dis(6,1);
	
		for (int i=0;i<6;++i)
		{
			dis(i,0) = distribute_pos->m_distribute[i];
		}
			/*dis(0,0) = distribute_pos->m_distribute[0];
			dis(1,0) = distribute_pos->m_distribute[1];
			dis(2,0) = distribute_pos->m_distribute[2];
			dis(3,0) = distribute_pos->m_distribute[3];
			dis(4,0) = distribute_pos->m_distribute[4];
			dis(5,0) = distribute_pos->m_distribute[5];*/
		if (distribute_pos->m_type == CDistribute_force::G)
			{
				matrix<double> dis1(6,1);
				dis1 = prod(/*trans*/(subrange(transm,0,6,0,6)),dis);  /**整体到局部**/
				dis  =dis1;
			}
	double	len;
	len = pos_bool->second.m_element_length;

	//分布荷载的插值会影响到单元的两个端点

	//分布荷载变成单元集中力，可以像单元集中力一样使用

	//double array_left[6]  = {q_x*len/2, q_y*len/2, q_z*len/2,
	//	-1*q_x*pow(len,2)/12, q_z*pow(len,2)/12, q_y*pow(len,2)/12};

	//double array_right[6] = {q_x*len/2, q_y*len/2, q_z*len/2,
	//	q_x*pow(len,2)/12, -1*q_z*pow(len,2)/12,-1*q_y*pow(len,2)/12};
	
	//这样就把y方向上的荷载等效成了y方向的位移和z方向的转角所受力

	double array_left[6]  = {dis(0,0)*len/2,dis(1,0)*len/2-dis(5,0), dis(2,0)*len/2+dis(4,0),
		dis(3,0)*len/2, -1*dis(2,0)*pow(len,2)/12, dis(1,0)*pow(len,2)/12};

	double array_right[6] = {dis(0,0)*len/2, dis(1,0)*len/2+dis(5,0), dis(2,0)*len/2-dis(4,0),
		dis(3,0)*len/2, dis(2,0)*pow(len,2)/12, -1*dis(1,0)*pow(len,2)/12};

	matrix<double>  matrix_left(6,1);
	matrix<double>  matrix_right(6,1);
	matrix<double>  matrix_left1(6,1);
	matrix<double>  matrix_right1(6,1);

	for (int i=0;i<6;++i)
	{
		matrix_left1(i,0) = array_left[i];
		matrix_right1(i,0) = array_right[i];
	}

	matrix_left=prod(trans(subrange(transm,0,6,0,6)),matrix_left1);/**局部到整体**/
	matrix_right=prod(trans(subrange(transm,0,6,0,6)),matrix_right1); 

	
	
	for (int i=0; i<6; ++i)
	{
		if (m_point_local[pos_bool->second.m_leftpoint][i]>0)
		{
			n1 = m_point_local[pos_bool->second.m_leftpoint][i] - 1;
		}
		else
		{
			n1 = -1 * m_point_local[pos_bool->second.m_leftpoint][i] - 1 + m_Free_degree;		
		}

		RHS(n1,0) += matrix_left(i,0);   //逐个向量相加

	//}

	//for (int i=0;i<6;++i)
	//{
		if (m_point_local[pos_bool->second.m_rightpoint][i]>0)
		{
			n2 = m_point_local[pos_bool->second.m_rightpoint][i] - 1;
		}
		else
		{
			n2 = -1 * m_point_local[pos_bool->second.m_rightpoint][i] - 1 + m_Free_degree;		
		}

		RHS(n2,0) += matrix_right(i,0);   //逐个向量相加

	}


	}
	

	cout<<"the Right hand side of distribute force and point force is:"<<endl;
	cout << RHS << endl;
	return RHS;
}

/// <summary>
/// Dirichlet_generates this instance.
/// </summary>
/// <returns>boost.numeric.ublas.matrix{double}.</returns>
matrix<double> CStructer::Dirichlet_generate()
{
	//引入位移边界条件
	matrix<double>  dirichlet (m_Dirichlet_degree,1);
	dirichlet.clear();
	int d = 0;
	for (std::map<int,CDisplace>::const_iterator pos = m_m_displaces.begin(); pos != m_m_displaces.end(); ++pos)
	{
		for(int i=0;i<6;++i)
		{
			if (pos->second.m_arraybool[i])
			{
				d=-1 * m_point_local[pos->second.m_No][i]-1;//+m_FreeNodes;
				dirichlet(d,0)=pos->second.m_array[i];
			}
		}
	}
	cout<<"the dirichlet boundary is:"<<endl;
	cout << dirichlet << endl;
	//位移边界条件结束
	return dirichlet;
}

/// <summary>
/// Solves this instance.
/// </summary>
/// <returns>bool.</returns>
bool CStructer::solve()
{
	if (local_vector_generate() == false)
	{
		return false;
	}
	symmetric_matrix<double> K = Assemble();
	matrix<double> Kfree = subrange(K,0,m_Free_degree,0,m_Free_degree);
	//matrix<double> RHS = RHS_generate();
	
	matrix<double> dirichlet = Dirichlet_generate();  

	matrix<double> RHS = RHS_generate1();  //增加一个局部坐标系下均匀荷载
	
	//cout<<"Kfree:"<<Kfree<<endl;
	
	//matrix<double> answer(6,1);
	//answer(0,0) = -0.00000000001722;
	//answer(1,0) = 0.00000000001859;
	//answer(2,0) = -0.000000007097;
	//answer(3,0) = -0.0000000007195;
	//answer(4,0) = -0.000000001746;
	//answer(5,0) = -0.000000000002509;


	//cout<<"Kfree*answer:"<<prod(Kfree,answer)<<endl;


	matrix<double> Kab =subrange(K,0,m_Free_degree,m_Free_degree,m_total_degree);

	matrix<double> Unknown=trans(subrange(RHS,0,m_Free_degree,0,1) - prod (Kab, dirichlet));//(m_displacement_set.size(),0);

	/************************************************************************/
	/*                     调用第三方Fortran库						        */
	/************************************************************************/
	
	//if (linalg::Solve(Kfree,Unknown) == 0)
	//{
	//	return false;
	//}
	//Unknown = trans(Unknown);

	/************************************************************************/
	/*     利用ublas库LU分解(当其所有顺序主子式都不为0时)来求解		        */
	/************************************************************************/

    permutation_matrix<double> P(m_Free_degree);
    boost::numeric::ublas::vector<double> x(m_Free_degree);
    boost::numeric::ublas::vector<double> v(m_Free_degree);

    for (int j = 0; j < m_Free_degree; ++j)
    {
        v(j) = Unknown(0, j);
    }

    //try
    {
        lu_factorize(Kfree, P);
        x = v;
        lu_substitute(Kfree, P, x);
        cout << "解向量: " << x << endl;
        for (int j = 0; j < m_Free_degree; ++j)
        {
            Unknown(0, j) = x(j);
        }
        Unknown = trans(Unknown);
    }
//catch (CMemoryException* e)
//{
//	std::cout << "内存操作失败" << std::endl; 
//}
//catch (CFileException* e)
//{
//	std::cout << "文件操作失败." << std::endl; 
//}
//catch (CException* e)
//{
//	std::cout << "奇异矩阵(无解或无穷解)." << std::endl; 
//}






/*************************************************/
	for (std::map<int,std::vector<int>>::iterator pos = m_point_local.begin(); pos != m_point_local.end(); ++pos)
	{
		//double displace_bool[3] = {true,true,true};
		//if (pos->second[0]>0)
		//{
		//	double displace1[3] = {Unknown(pos->second[0]-1,0),Unknown(pos->second[1]-1,0),Unknown(pos->second[2]-1,0)};
		//	
		//	CPlaneDisplacement Dis1(pos->first,displace1,displace_bool);
		//	m_displacement_map_result.insert(make_pair(pos->first,Dis1));
		//	
		//} 
		//else
		//{
		//	double displace1[3] = {dirichlet(-1*pos->second[0]-1,0),dirichlet(-1*pos->second[1]-1,0),dirichlet(-1*pos->second[2]-1,0)};
		//	CPlaneDisplacement Dis1(pos->first,displace1,displace_bool);
		//	m_displacement_map_result.insert(make_pair(pos->first,Dis1));
		//	
		//}
		/************************************************************************/
		/* 以上是假设强制位移全部为true的情况，接下来是针对有false的情形        */
		/************************************************************************/
		bool displace_bool[6] = {true,true,true,true,true,true};
		double displace1[6];

		for (int k=0;k<6;++k)
		{
			if (pos->second[k]>0)
			{
				displace1[k]=Unknown(pos->second[k]-1,0);
			}
			else
			{
				displace1[k]=dirichlet(-1*pos->second[k]-1,0);
			}
		}
		CDisplace Dis1(pos->first,displace1,displace_bool);
		m_m_result.insert(make_pair(pos->first,Dis1));
	}

	return true;
}
/************************************************************************/
/* 面向对象的编程的重点之一就是代码的重用                               */
/************************************************************************/
/// <summary>
/// Show_results this instance.
/// </summary>
void CStructer::show_result()
{
	for (std::map<int,CDisplace>::iterator pos = m_m_result.begin(); pos != m_m_result.end(); ++pos)
	{
		pos->second.show();

	}
}	