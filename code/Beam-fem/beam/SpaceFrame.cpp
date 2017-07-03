// ***********************************************************************
// Assembly         : SpaceFrame
// Author           : chenruishan
// Created          : 09-28-2014
//
// Last Modified By : chenruishan
// Last Modified On : 09-29-2014
// ***********************************************************************
// <copyright file="SpaceFrame.cpp" company="">
//     Copyright (c) . All rights reserved.
// </copyright>
// <summary>
//这是空间梁计算的入口程序
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
//
//
//
//
//</summary>
// ***********************************************************************

//#include "stdafx.h"
#include "SpaceFrame.h"
#include "Structer.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// The one and only application object

/// <summary>
/// The application
/// </summary>
CWinApp theApp;

using namespace std;

void example1();

void example2();
void example3();
void example4();
/// <summary>
/// Mains the specified argc.
/// </summary>
/// <param name="argc">The argc.</param>
/// <param name="argv">The argv.</param>
/// <param name="envp">The envp.</param>
/// <returns>int.</returns>
int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int nRetCode = 0;

	// initialize MFC and print and error on failure
	if (!AfxWinInit(::GetModuleHandle(NULL), NULL, ::GetCommandLine(), 0))
	{
		// TODO: change error code to suit your needs
		_tprintf(_T("Fatal Error: MFC initialization failed\n"));
		nRetCode = 1;
	}
	else
	{
		example4();

	}

	return nRetCode;
}


void example1()
{
		cout<<"kattan book测试"<<endl;
		/***************************************************************************/
		/*如何把txt文件的内容直接读入到文件中,这样用户就不用打开文件就能运行程序了 */
		/***************************************************************************/
		// TODO: code your application's behavior here.
		CSpacePoint P1(1,0,0,0);
		CSpacePoint P2(2,3,0,1);
		CSpacePoint P3(3,0,0,-3);
		CSpacePoint P4(4,0,-4,0);

		std::map<int,CSpacePoint> point_map;

		point_map.insert(make_pair(1,P1));
		point_map.insert(make_pair(2,P2));
		point_map.insert(make_pair(3,P3));
		point_map.insert(make_pair(4,P4));

		/************************************************************************/
		/*                想办法 从文件自动输入                                 */
		/************************************************************************/
		/*******************************************************/
		CSpaceElement E1(1,1,2,21E10,84E9,2E-2,2E-4,1E-4,5E-5);
		CSpaceElement E2(2,1,3,21E10,84E9,2E-2,2E-4,1E-4,5E-5);
		CSpaceElement E3(3,1,4,21E10,84E9,2E-2,2E-4,1E-4,5E-5);

		std::map<int,CSpaceElement> element_map;

		element_map.insert(make_pair(1,E1));
		element_map.insert(make_pair(2,E2));
		element_map.insert(make_pair(3,E3));



		double f1[6] = {0,0,0,0,0,0};
		CSpaceForce F1(1,f1);
		std::vector<CSpaceForce> force_vector;
		force_vector.push_back(F1);
		/************************************************************************/
		/*                               增加了均匀荷载                         */
		/************************************************************************/
		double distribute_f1[6] = {0,0,0,0,0,10};
		CDistribute_force dis_F1(1,CDistribute_force::G,distribute_f1);

		//double distribute_f2[6] = {0,0,0,0,0,0};
		//CDistribute_force dis_F2(1,G,distribute_f2);

		std::vector<CDistribute_force> distribute_forces;
		distribute_forces.push_back(dis_F1);
		//distribute_forces.push_back(dis_F2);


		double displace2[6] = {0,0,0,0,0,0};
		bool displace_bool2[6] = {true,true,true,true,true,true};
		CDisplace Dis2(2,displace2,displace_bool2);

		double displace3[6] = {0,0,0,0,0,0};
		bool displace_bool3[6] = {true,true,true,true,true,true};
		CDisplace Dis3(3,displace3,displace_bool3);

		double displace4[6] = {0,0,0,0,0,0};
		bool displace_bool4[6] = {true,true,true,true,true,true};
		CDisplace Dis4(4,displace4,displace_bool4);

		std::map<int,CDisplace> displacement_map;
		displacement_map.insert(make_pair(2,Dis2));
		displacement_map.insert(make_pair(3,Dis3));
		displacement_map.insert(make_pair(4,Dis4));
		

		CStructer str1(point_map,element_map,force_vector,distribute_forces,displacement_map);

		if(str1.solve())
		{
			str1.show_result();
		}
		else
		{
			cout<<"It does not work!"<<endl;
		}

}

void example2()
{
	cout<<"一般化测试"<<endl;
	/***************************************************************************/
	/*如何把txt文件的内容直接读入到文件中,这样用户就不用打开文件就能运行程序了 */
	/***************************************************************************/
	// TODO: code your application's behavior here.
	CSpacePoint P1(1,0,1,2);
	CSpacePoint P2(2,3,5,14);
	CSpacePoint P3(3,6,8,10);


	std::map<int,CSpacePoint> point_map;

	point_map.insert(make_pair(1,P1));
	point_map.insert(make_pair(2,P2));
	point_map.insert(make_pair(3,P3));

	/*******************************************************/
	CSpaceElement E1(1,1,2,21E10,84E9,2E-2,2E-4,1E-4,5E-5);
	CSpaceElement E2(2,2,3,21E10,84E9,2E-2,2E-4,1E-4,5E-5);

	std::map<int,CSpaceElement> element_map;

	element_map.insert(make_pair(1,E1));
	element_map.insert(make_pair(2,E2));


	double f1[6] = {1,2,3,-1,-2,-3};
	CSpaceForce F1(3,f1);
	std::vector<CSpaceForce> force_vector;
	force_vector.push_back(F1);

	double distribute_f1[6] = {-1,-2,-3,1,2,3};
	CDistribute_force dis_F1(2,CDistribute_force::L,distribute_f1);

	std::vector<CDistribute_force> distribute_forces;
	distribute_forces.push_back(dis_F1);

	double displace1[6] = {0,0,0,0,0,0};
	bool displace_bool1[6] = {true,true,true,true,true,true};
	CDisplace Dis1(1,displace1,displace_bool1);

	double displace2[6] = {1,2,3,0,0,0};
	bool displace_bool2[6] = {true,true,true,false,false,false};
	CDisplace Dis2(2,displace2,displace_bool2);

	std::map<int,CDisplace> displacement_map;
	displacement_map.insert(make_pair(1,Dis1));
	displacement_map.insert(make_pair(2,Dis2));

	CStructer str1(point_map,element_map,force_vector,distribute_forces,displacement_map);

	if(str1.solve())
	{
		str1.show_result();
	}
	else
	{
		cout<<"It does not work!"<<endl;
	}

}

void example3()
{
	cout << "翘曲测试" << endl;
	/***************************************************************************/
	/*如何把txt文件的内容直接读入到文件中,这样用户就不用打开文件就能运行程序了 */
	/***************************************************************************/
	// TODO: code your application's behavior here.
	CSpacePoint P1(1, 0, 0, 0);
	CSpacePoint P2(2, 50, 30, 40);



	std::map<int, CSpacePoint> point_map;

	point_map.insert(make_pair(1, P1));
	point_map.insert(make_pair(2, P2));


	/*******************************************************/
	CSpaceElement E1(1, 1, 2, 3.25E10, 1.3E10, 2.24, 1.4379, 4.4459, 3.3418);

	std::map<int, CSpaceElement> element_map;

	element_map.insert(make_pair(1, E1));



	double f1[6] = { 0, 0, 0, 1000000000, 0, 0 };
	CSpaceForce F1(2, f1);
	std::vector<CSpaceForce> force_vector;
	force_vector.push_back(F1);

	double distribute_f1[6] = { -1, -2, -3, 1, 2, 3 };
	CDistribute_force dis_F1(2, CDistribute_force::L, distribute_f1);

	std::vector<CDistribute_force> distribute_forces;
	//distribute_forces.push_back(dis_F1);

	double displace1[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	bool displace_bool1[6] = { true, true, true, true, true, true };
	CDisplace Dis1(1, displace1, displace_bool1);

	double displace2[6] = { 0, 0, 0, 0, 0, 0 };
	bool displace_bool2[6] = { false, false, false, false, false, false };
	CDisplace Dis2(2, displace2, displace_bool2);

	std::map<int, CDisplace> displacement_map;
	displacement_map.insert(make_pair(1, Dis1));
	displacement_map.insert(make_pair(2, Dis2));

	CStructer str1(point_map, element_map, force_vector, distribute_forces, displacement_map);

	if (str1.solve())
	{
		str1.show_result();
	}
	else
	{
		cout << "It does not work!" << endl;
	}

}

void example4()
{
	cout << "两单元翘曲测试" << endl;
	/***************************************************************************/
	/*如何把txt文件的内容直接读入到文件中,这样用户就不用打开文件就能运行程序了 */
	/***************************************************************************/
	// TODO: code your application's behavior here.
	CSpacePoint P1(1, 0, 0, 0);
	CSpacePoint P2(2, 50, 30, 40);
	CSpacePoint P3(3, 80, 40, 10);


	std::map<int, CSpacePoint> point_map;

	point_map.insert(make_pair(1, P1));
	point_map.insert(make_pair(2, P2));
	point_map.insert(make_pair(3, P3));

	/*******************************************************/
	CSpaceElement E1(1, 1, 2, 3.25E10, 1.3E10, 2.24, 1.4379, 4.4459, 3.3418);
	CSpaceElement E2(2, 2, 3, 3.25E10, 1.3E10, 2.24, 1.4379, 4.4459, 3.3418);

	std::map<int, CSpaceElement> element_map;

	element_map.insert(make_pair(1, E1));
	element_map.insert(make_pair(2, E2));


	double f2[6] = { 3440, 1790, 6530, 1000, 8540, 4750 };
	CSpaceForce F2(2, f2);
	double f3[6] = { 2220, 7081, 4850, 1300, 5410, 4710 };
	CSpaceForce F3(3, f3);
	std::vector<CSpaceForce> force_vector;
	force_vector.push_back(F2);
	force_vector.push_back(F3);

	double distribute_f1[6] = { -1, -2, -3, 1, 2, 3 };
	CDistribute_force dis_F1(2, CDistribute_force::L, distribute_f1);

	std::vector<CDistribute_force> distribute_forces;
	//distribute_forces.push_back(dis_F1);

	double displace1[6] = { 0, 0, 0, 0, 0, 0 };
	bool displace_bool1[6] = { true, true, true, true, true, true };
	CDisplace Dis1(1, displace1, displace_bool1);

	double displace2[6] = { 1, 2, 3, 0, 0, 0 };
	bool displace_bool2[6] = { false, false, false, false, false, false };
	CDisplace Dis2(2, displace2, displace_bool2);

	std::map<int, CDisplace> displacement_map;
	displacement_map.insert(make_pair(1, Dis1));
	displacement_map.insert(make_pair(2, Dis2));

	CStructer str1(point_map, element_map, force_vector, distribute_forces, displacement_map);

	if (str1.solve())
	{
		str1.show_result();
	}
	else
	{
		cout << "It does not work!" << endl;
	}

}