#include "StdAfx.h"
#include "PlaneStructer.h"
#include "SFrame\linalg.h"
using namespace std;
using namespace boost::numeric::ublas;

CPlaneStructer::CPlaneStructer(const std::map<int,CPlanePoint> &point_map,
							   const std::map<int,CPlaneElement> &element_map,
							   const std::vector<CPlaneForce> &force_vector,
							   const std::map<int,CPlaneDisplacement> &displacement_map):
m_point_map(point_map),
m_element_map(element_map),
m_force_vector(force_vector),
m_displacement_map(displacement_map)
{
}

CPlaneStructer::~CPlaneStructer(void)
{
}

bool CPlaneStructer::local_vector_generate()
{
	m_point_local.clear();//保证可以重入
	m_element_local.clear();

	m_total_degree = 3 * (int)m_point_map.size();
	m_Dirichlet_degree = 0; 

	for (std::map<int,CPlaneDisplacement>::iterator pos1 = m_displacement_map.begin(); pos1 != m_displacement_map.end(); ++pos1)
	{
		for (int i=0;i<3;++i)
		{
			if(pos1->second.m_displacement_bool[i])
			{
				++m_Dirichlet_degree;
			}
		}
	}
	cout<<"m_Dirichlet_degree:"<<m_Dirichlet_degree<<endl;
	/*m_Dirichlet_degree = 3 * (int)m_displacement_map.size();*/
	m_Free_degree = m_total_degree - m_Dirichlet_degree;


	int i=0,j=0;
	for (std::map<int,CPlanePoint>::iterator pos = m_point_map.begin();pos != m_point_map.end();++pos)
	{
		/**为了做位移是否是true的判断**/
		std::map<int,CPlaneDisplacement>::iterator pos_bool ;
		pos_bool = m_displacement_map.find(pos->first);
		/*************/
		std::vector<int>  Position1(3);
		if (m_displacement_map.find(pos->first) == m_displacement_map.end())
		{
			for (int k=0;k<3;++k)
			{
				Position1[k]=j+1;
				++j;
			}
			//Position1[0]=j+1;
			//Position1[1]=j+2;
			//Position1[2]=j+3;
			//++(++(++j));
		} 
		else
		{
			for (int k=0;k<3;++k)
			{
				//Position1[k]=-1*(i+1);
				//++i;
				if (pos_bool->second.m_displacement_bool[k])
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


	for (std::map<int,CPlaneElement>::iterator pos = m_element_map.begin();pos!= m_element_map.end();++pos)
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

matrix<double> CPlaneStructer::Assemble()
{
	for (std::map<int,CPlaneElement>::iterator pos = m_element_map.begin();pos != m_element_map.end();++pos)
	{
		pos->second.setpoint(&m_point_map);
	}

	//声明自由度n*n一样大的矩阵
	symmetric_matrix<double> m (m_total_degree, m_total_degree);   //对称矩阵
	m.clear();
	int I,J;
	for (std::map<int,CPlaneElement>::const_iterator pos = m_element_map.begin();pos != m_element_map.end();++pos)
	{
		symmetric_matrix<double> mm = pos->second.stiffmatrix();
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

	cout<<"the structer's stiffmatrix is:"<<endl;
	cout << m << endl;
	return m;
}

matrix<double> CPlaneStructer::RHS_generate()
{
	//确定右边向量，求解方程的时候再截短成自由度那么长
	matrix<double> RHS (m_total_degree,1);
	RHS.clear();
	int m = 0;
	for (std::vector<CPlaneForce>::const_iterator pos = m_force_vector.begin(); pos != m_force_vector.end(); ++pos)
	{

		for (int i=0;i<3;++i)
		{
			if (m_point_local[pos->m_Force_NO][i]>0)
			{
				m = m_point_local[pos->m_Force_NO][i] - 1;
			}
			else
			{
				m = -1 * m_point_local[pos->m_Force_NO][i] - 1 + m_Free_degree;		
			}

			RHS(m,0) += pos->m_Force_array[i];   //逐个向量相加

		}
	}
	//右边向量完成
	cout<<"the Right hand side is:"<<endl;
	cout << RHS << endl;
	return RHS;
}

matrix<double> CPlaneStructer::Dirichlet_generate()
{
	//引入位移边界条件
	matrix<double>  dirichlet (m_Dirichlet_degree,1);
	dirichlet.clear();
	int d = 0;
	for (std::map<int,CPlaneDisplacement>::const_iterator pos = m_displacement_map.begin(); pos != m_displacement_map.end(); ++pos)
	{
		for(int i=0;i<3;++i)
		{
			if (pos->second.m_displacement_bool[i])
			{
				d=-1 * m_point_local[pos->second.m_displacement_NO][i]-1;//+m_FreeNodes;
				dirichlet(d,0)=pos->second.m_displacement_array[i];
			}
		}
	}
	cout<<"the dirichlet boundary is:"<<endl;
	cout << dirichlet << endl;
	//位移边界条件结束
	return dirichlet;
}

bool CPlaneStructer::solve()
{
	if (local_vector_generate() == false)
	{
		return false;
	}
	symmetric_matrix<double> K = Assemble();
	matrix<double> RHS = RHS_generate();
	matrix<double> dirichlet = Dirichlet_generate();  
	matrix<double> Kfree = subrange(K,0,m_Free_degree,0,m_Free_degree);

	matrix<double> Kab =subrange(K,0,m_Free_degree,m_Free_degree,m_total_degree);

	matrix<double> Unknown=trans(subrange(RHS,0,m_Free_degree,0,1) - prod (Kab, dirichlet));//(m_displacement_set.size(),0);

	if (linalg::Solve(Kfree,Unknown) == 0)
	{
		return false;
	}
	Unknown = trans(Unknown);
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
		bool displace_bool[3] = {true,true,true};
		double displace1[3];

		for (int k=0;k<3;++k)
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
		CPlaneDisplacement Dis1(pos->first,displace1,displace_bool);
		m_displacement_map_result.insert(make_pair(pos->first,Dis1));
	}

	return true;
}
/************************************************************************/
/* 面向对象的编程的重点之一就是代码的重用                               */
/************************************************************************/
void CPlaneStructer::show_result()
{
	for (std::map<int,CPlaneDisplacement>::iterator pos = m_displacement_map_result.begin(); pos != m_displacement_map_result.end(); ++pos)
	{
		pos->second.show();

	}
}