//#include "StdAfx.h"
#include "Displace.h"

#include <iomanip>
#include <iostream>
using namespace std;

/// <summary>
/// Initializes a new instance of the <see cref="CDisplace"/> class.
/// </summary>
CDisplace::CDisplace(void)
{
}

/// <summary>
/// Initializes a new instance of the <see cref="CDisplace"/> class.
/// </summary>
/// <param name="No">The no.</param>
/// <param name="array">The array.</param>
/// <param name="arraybool">The arraybool.</param>
CDisplace::CDisplace(const int &No, double* array, bool* arraybool)
{
	m_No  = No;
	for (int i=0; i<6; ++i)
	{
		m_array[i]	   = array[i];
		m_arraybool[i] = arraybool[i];
	}
}

/// <summary>
/// Finalizes an instance of the <see cref="CDisplace"/> class.
/// </summary>
CDisplace::~CDisplace(void)
{
}

/// <summary>
/// Shows this instance.
/// </summary>
void CDisplace::show(void)
{	
	cout<<setprecision(11);//为了跟sap比较

	cout<<"the No of the node:"<<m_No<<endl;
	for (int i=0;i<3;++i)
	{
		if(m_arraybool[i])
		{
			cout<<"第" << i+1<<" 个位移 :"<<m_array[i]<<endl;
		}
	}

	for (int i=3;i<6;++i)
	{
		if(m_arraybool[i])
		{
			cout<<"第" << i-2<<" 个转角 :"<<m_array[i]<<endl;
		}
	}
	cout<<setprecision(6);
}
