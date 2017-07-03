// Homotopy.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include "Homotopy.h"
//#include "EulerNewton.h"
#include "LiYorke.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif



// The one and only application object



using namespace std;
using namespace boost::numeric::ublas;

bool LYexample1();
bool LYexample2();
bool LYexample3(const double & x);
//bool ENexample1();
//bool ENexample2();
//bool ENexample3();


int main()
{
    int nRetCode = 0;

    if (LYexample3(0.5) == false)
    {
        cout << "ʧ��!" << endl;
    }
    else
    {
        cout << "�ɹ�!" << endl;
    }


    return nRetCode;
}


bool LYexample1()
{
	//��������ֵ���� ����  ��4.9.1
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
		1E-5,//ţ��������׼��
		0.01,//��ʼ����
		10,
		delta,
		angle,
		pointsOnRoad) == false)
	{
		cout << "���ʧ�ܣ�" << endl;
		return false;
	}
	return true;
}

bool LYexample2()
{
	//ͬ�׷������� ����� ��һ
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
		1E-3,//ţ��������׼��
		0.01,//��ʼ����
		20,
		delta,
		angle,
		pointsOnRoad) == false)
	{
		cout << "���ʧ�ܣ�" << endl;
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
		1E-3,//ţ��������׼��
		initstepsize,//��ʼ����
		r,
		delta,
		angle,
		pointsOnRoad) == false)
	{
		cout << "���ʧ�ܣ�" << endl;

		//cout << "��ֵΪ" << initpoint1(0, 0) << endl;
		
		return false;
	}

	//cout << "��ֵΪ" << initpoint1(0, 0) << endl;
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
//		1E-5,//ţ��������׼��
//		0.01,//��ʼ����
//		10,
//		delta,
//		angle,
//		pointsOnRoad) == false)
//	{
//		cout << "���ʧ�ܣ�" << endl;
//		return false;
//	}
//	return true;
//}
//
//bool ENexample2()
//{
//	//ͬ�׷������� ����� ��һ
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
//		1E-2,//ţ��������׼��
//		0.5,//��ʼ����
//		20,
//		delta,
//		angle,
//		pointsOnRoad) == false)
//	{
//		cout << "���ʧ�ܣ�" << endl;
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
//		1E-7,//ţ��������׼��
//		initstepsize,//��ʼ����
//		r,
//		delta,
//		angle,
//		pointsOnRoad) == false)
//	{
//		cout << "���ʧ�ܣ�" << endl;
//
//		cout << "��ֵΪ" << initpoint1(0, 0) << endl;
//
//		return false;
//	}
//
//	cout << "��ֵΪ" << initpoint1(0, 0) << endl;
//	return true;
//}