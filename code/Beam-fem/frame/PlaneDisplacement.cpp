#include "StdAfx.h"
#include "PlaneDisplacement.h"
using namespace std;

CPlaneDisplacement::CPlaneDisplacement(const int &displacement_NO,double* displacement_array,bool* displacement_bool)
{
	m_displacement_NO = displacement_NO;
	//for (int i=0;i<3; ++i)
	//{
	//	m_displacement_array[i] = displacement_array[i];
	//	m_displacement_bool[i]=displacement_bool[i];
	//}
	memcpy(m_displacement_array,displacement_array,sizeof(double)*3);
	memcpy(m_displacement_bool,displacement_bool,sizeof(bool)*3);
	 
}

CPlaneDisplacement::~CPlaneDisplacement(void)
{
}

void CPlaneDisplacement::show(void)
{
	cout<<"the No of the node:"<<m_displacement_NO<<endl;
	for (int i=0;i<3;++i)
	{
		if(m_displacement_bool[i])
		{
			cout<<"the NO" << i+1<<" displacement :"<<m_displacement_array[i]<<endl;
		}
	}
}