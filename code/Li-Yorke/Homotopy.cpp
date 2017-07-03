// Homotopy.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Homotopy.h"
//#include "EulerNewton.h"
#include "LiYorke.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif



// The one and only application object

CWinApp theApp;

using namespace std;
using namespace boost::numeric::ublas;

bool LYexample1();
bool LYexample2();
bool LYexample3(const double & x);
//bool ENexample1();
//bool ENexample2();
//bool ENexample3();


int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int nRetCode = 0;

	HMODULE hModule = ::GetModuleHandle(NULL);

	if (hModule != NULL)
	{
		// initialize MFC and print and error on failure
		if (!AfxWinInit(hModule, NULL, ::GetCommandLine(), 0))
		{
			// TODO: change error code to suit your needs
			_tprintf(_T("Fatal Error: MFC initialization failed\n"));
			nRetCode = 1;
		}
		else
		{


			// TODO: code your application's behavior here.
			
			{
				if (LYexample1() == false)
				{
					cout << "失败!" << endl;
				}
				else
				{
					cout << "成功!" << endl;
				}
			}

		}
	}
	else
	{
		// TODO: change error code to suit your needs
		_tprintf(_T("Fatal Error: GetModuleHandle failed\n"));
		nRetCode = 1;
	}

	return nRetCode;
}


bool LYexample1()
{
	//非线性数值分析 黄象鼎  例4.9.1
	BrouwerFixedPointFun fun1;

	matrix<double> initpoint1(3,1);

	double eps = 0.0E-3;
	initpoint1(0, 0) = 0.54232000 + eps;
	initpoint1(1, 0) = -0.321000 - eps;

	//initpoint1(0, 0) = -0.707027000;
	//initpoint1(1, 0) = -0.688818000;
	initpoint1(2, 0) = 0.0;

	fun1.Setinit(initpoint1);
	
	std::map<int, double> delta;
	std::map<int, double> angle;
	std::map<int, matrix<double>> pointsOnRoad;

	CLiYorke BrouwerFixedPoint;

	if (BrouwerFixedPoint.Solve(
		&fun1,
		initpoint1,
		1E-5,//牛顿收敛的准则
		0.01,//初始步长
		10,
		delta,
		angle,
		pointsOnRoad) == false)
	{
		cout << "求解失败！" << endl;
		return false;
	}
	return true;
}

bool LYexample2()
{
	//同伦方法引论 王则柯 例一
	//Problem1 fun1;
	minusProblem1 fun1;

	matrix<double> initpoint1(3, 1);
	initpoint1(0, 0) = 0.50;
	initpoint1(1, 0) = 1.0;
	initpoint1(2, 0) = 0.0;

	fun1.Setinit(initpoint1);

	std::map<int, double> delta;
	std::map<int, double> angle;
	std::map<int, matrix<double>> pointsOnRoad;

	CLiYorke BrouwerFixedPoint;

	if (BrouwerFixedPoint.Solve(
		&fun1,
		initpoint1,
		1E-3,//牛顿收敛的准则
		0.01,//初始步长
		20,
		delta,
		angle,
		pointsOnRoad) == false)
	{
		cout << "求解失败！" << endl;
		return false;
	}
	return true;
}

bool LYexample3(const double & x)
{

	matrix<double> initpoint1(2, 1);
	initpoint1(0, 0) = x;
	initpoint1(1, 0) = 0.0;

	Snap fun1;

	fun1.Setinit(initpoint1);

	std::map<int, double> delta;
	std::map<int, double> angle;
	std::map<int, matrix<double>> pointsOnRoad;
	double r = 200.0;

	double initstepsize = min(0.1, 0.1*r);

	CLiYorke snapback;

	if (snapback.Solve(
		&fun1,
		initpoint1,
		1E-3,//牛顿收敛的准则
		initstepsize,//初始步长
		r,
		delta,
		angle,
		pointsOnRoad) == false)
	{
		cout << "求解失败！" << endl;

		//cout << "初值为" << initpoint1(0, 0) << endl;
		
		return false;
	}

	//cout << "初值为" << initpoint1(0, 0) << endl;
	return true;
}


//bool ENexample1()
//{
//	BrouwerFixedPointFun fun1;
//
//	matrix<double> initpoint1(3, 1);
//
//	double eps = 0.0E-3;
//	initpoint1(0, 0) = 0.54232000 + eps;
//	initpoint1(1, 0) = -0.321000 - eps;
//
//	//initpoint1(0, 0) = -0.707027000;
//	//initpoint1(1, 0) = -0.688818000;
//	initpoint1(2, 0) = 0.0;
//
//	fun1.Setinit(initpoint1);
//
//	std::map<int, double> delta;
//	std::map<int, double> angle;
//	std::map<int, matrix<double>> pointsOnRoad;
//
//	CEulerNewton BrouwerFixedPoint;
//
//	if (BrouwerFixedPoint.Solve(
//		&fun1,
//		initpoint1,
//		1E-5,//牛顿收敛的准则
//		0.01,//初始步长
//		10,
//		delta,
//		angle,
//		pointsOnRoad) == false)
//	{
//		cout << "求解失败！" << endl;
//		return false;
//	}
//	return true;
//}
//
//bool ENexample2()
//{
//	//同伦方法引论 王则柯 例一
//	//Problem1 fun1;
//	minusProblem1 fun1;
//
//	matrix<double> initpoint1(3, 1);
//	initpoint1(0, 0) = 0.5;
//	initpoint1(1, 0) = 1.0;
//
//	initpoint1(2, 0) = 0.0;
//	fun1.Setinit(initpoint1);
//
//	std::map<int, double> delta;
//	std::map<int, double> angle;
//	std::map<int, matrix<double>> pointsOnRoad;
//
//	CEulerNewton BrouwerFixedPoint;
//
//	if (BrouwerFixedPoint.Solve(
//		&fun1,
//		initpoint1,
//		1E-2,//牛顿收敛的准则
//		0.5,//初始步长
//		20,
//		delta,
//		angle,
//		pointsOnRoad) == false)
//	{
//		cout << "求解失败！" << endl;
//		return false;
//	}
//	return true;
//}
//
//bool ENexample3()
//{
//
//	Snap fun1;
//
//	matrix<double> initpoint1(2, 1);
//	initpoint1(0, 0) = 2.0;
//	initpoint1(1, 0) = 0.0;
//	fun1.Setinit(initpoint1);
//	std::map<int, double> delta;
//	std::map<int, double> angle;
//	std::map<int, matrix<double>> pointsOnRoad;
//	double r = 100.0;
//
//	double initstepsize = min(0.01, 0.1*r);
//
//	CEulerNewton snapback;
//
//	if (snapback.Solve(
//		&fun1,
//		initpoint1,
//		1E-7,//牛顿收敛的准则
//		initstepsize,//初始步长
//		r,
//		delta,
//		angle,
//		pointsOnRoad) == false)
//	{
//		cout << "求解失败！" << endl;
//
//		cout << "初值为" << initpoint1(0, 0) << endl;
//
//		return false;
//	}
//
//	cout << "初值为" << initpoint1(0, 0) << endl;
//	return true;
//}