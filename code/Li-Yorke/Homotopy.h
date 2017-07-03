#pragma once
#include <math.h>
#include <vector>
#include "resource.h"
#include "TargetFun.h"

namespace ublas = boost::numeric::ublas;


/// <summary>
/// CTargetFun的派生类 BrouwerFixedPointFun.
/// </summary>
class BrouwerFixedPointFun:public CTargetFun
{
public:
	BrouwerFixedPointFun();
	~BrouwerFixedPointFun();
	virtual ublas::matrix<double> Homomap(ublas::matrix<double>& point) const;
	virtual ublas::matrix<double> HomomapJacobi(ublas::matrix<double>& point) const;
	virtual ublas::matrix<double> HomomapT(ublas::matrix<double>& point) const;
	virtual double AdditionalEqn(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const;
	virtual ublas::matrix<double> AdditionalEqnPartial(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const;
};

BrouwerFixedPointFun::BrouwerFixedPointFun()
{
}

BrouwerFixedPointFun::~BrouwerFixedPointFun()
{
}
ublas::matrix<double> BrouwerFixedPointFun::Homomap(ublas::matrix<double>& point) const
{


	ublas::matrix<double> a(2, 1);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double y0 = m_initpoint(1, 0);
	double x = point(0, 0);
	double y = point(1, 0);
	double t = point(2, 0);

	a(0, 0) = x - ((1 - t)*x0 + t*(sin(100 * x) - cos(100 * y)));
	a(1, 0) = y - ((1 - t)*y0 + t*(sin(100 * y) - cos(100 * x)));

	//a(0, 0) = x - sin(100 * x) + cos(100 * y) + (t - 1)*(x0 - sin(100 * x0) + cos(100 * y0));
	//a(1, 0) = y + cos(100 * x) - sin(100 * y) + (t - 1)*(y0 + cos(100 * x0) - sin(100 * y0));

	return a;
}

ublas::matrix<double> BrouwerFixedPointFun::HomomapJacobi(ublas::matrix<double>& point) const
{
	ublas::matrix<double> a(2, 2);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double y0 = m_initpoint(1, 0);
	double x = point(0, 0);
	double y = point(1, 0);
	double t = point(2, 0);
	a(0, 0) = 1 - 100 * t*cos(100 * x);
	a(1, 0) = -100 * t*sin(100 * x);
	a(0, 1) = -100 * t*sin(100 * y);
	a(1, 1) = 1 - 100 * t*cos(100 * y);

	/*a(0, 0) = 1 - 100*cos(100*x);
	a(1, 0) = -100*sin(100*x);
	a(0, 1) = -100*sin(100*y);
	a(1, 1) = 1 - 100*cos(100*y);*/
	return a;
}

ublas::matrix<double> BrouwerFixedPointFun::HomomapT(ublas::matrix<double>& point) const
{

	ublas::matrix<double> a(2, 1);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double y0 = m_initpoint(1, 0);
	double x = point(0, 0);
	double y = point(1, 0);
	double t = point(2, 0);

	a(0, 0) = x0 + cos(100 * y) - sin(100 * x);
	a(1, 0) = y0 + cos(100 * x) - sin(100 * y);

	//a(0, 0) = x0 - sin(100 * x0) + cos(100 * y0);
	//a(1, 0) = y0 + cos(100 * x0) - sin(100 * y0);

	return a;
}

double BrouwerFixedPointFun::AdditionalEqn(
	ublas::matrix<double>& point,
	ublas::matrix<double>& prepoint,
	ublas::matrix<double>& tanvect,
	ublas::matrix<double>& verticalpoint,
	double& delta,
	double& sigma) const
{
	assert(point.size1() == tanvect.size1() && tanvect.size1() == verticalpoint.size1());
	assert(point.size2() == 1 && tanvect.size2() == 1 && verticalpoint.size2() == 1);
	ublas::matrix<double> tempvec1(tanvect);
	ublas::matrix<double> tempvec2(point - verticalpoint);

	ublas::matrix_column<ublas::matrix<double> > vec1column(tempvec1, 0);
	ublas::matrix_column<ublas::matrix<double> > vec2column(tempvec2, 0);
	return inner_prod(vec1column, vec2column);
}


ublas::matrix<double> BrouwerFixedPointFun::AdditionalEqnPartial(
	ublas::matrix<double>& point,
	ublas::matrix<double>& prepoint,
	ublas::matrix<double>& tanvect,
	ublas::matrix<double>& verticalpoint,
	double& delta,
	double& sigma) const
{
	ublas::matrix<double> a(tanvect);
	return a;
}

/// <summary>
/// CTargetFun的派生类 Problem1.
/// </summary>
class Problem1 :public CTargetFun
{
public:
	Problem1();
	~Problem1();
	virtual ublas::matrix<double> Homomap(ublas::matrix<double>& point) const;
	virtual ublas::matrix<double> HomomapJacobi(ublas::matrix<double>& point) const;
	virtual ublas::matrix<double> HomomapT(ublas::matrix<double>& point) const;
	virtual double AdditionalEqn(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const;
	virtual ublas::matrix<double> AdditionalEqnPartial(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const;

};

Problem1::Problem1()
{
}

Problem1::~Problem1()
{
}
ublas::matrix<double> Problem1::Homomap(ublas::matrix<double>& point) const
{
	ublas::matrix<double> a(2, 1);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double y0 = m_initpoint(1, 0);
	double x = point(0, 0);
	double y = point(1, 0);
	double t = point(2, 0);

	a(0, 0) = (t)*(x*x + x*y - 3) + (1 - t)*(x - x0);
	a(1, 0) = (t)*(2 * sin(x) + y - 7) + (1 - t)*(y - y0);
	return a;
}

ublas::matrix<double> Problem1::HomomapJacobi(ublas::matrix<double>& point) const
{
	ublas::matrix<double> a(2, 2);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double y0 = m_initpoint(1, 0);
	double x = point(0, 0);
	double y = point(1, 0);
	double t = point(2, 0);
	a(0, 0) = (t)*(2 * x + y) + 1-t;
	a(1, 0) = (t) * 2 * cos(x);
	a(0, 1) = (t)*x;
	a(1, 1) = 1;
	return a;
}

ublas::matrix<double> Problem1::HomomapT(ublas::matrix<double>& point) const
{

	ublas::matrix<double> a(2, 1);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double y0 = m_initpoint(1, 0);
	double x = point(0, 0);
	double y = point(1, 0);
	double t = point(2, 0);

	a(0, 0) = 1 * (x*x + x*y - 3) -1*(x - x0);
	a(1, 0) = 1 * (2 * sin(x) + y - 7) - 1 * (y - y0);
	return a;
}

double Problem1::AdditionalEqn(
	ublas::matrix<double>& point,
	ublas::matrix<double>& prepoint,
	ublas::matrix<double>& tanvect,
	ublas::matrix<double>& verticalpoint,
	double& delta,
	double& sigma) const  //超圆柱面切割
{
	assert(point.size1() == tanvect.size1() && tanvect.size1() == verticalpoint.size1());
	assert(point.size2() == 1 && tanvect.size2() == 1 && verticalpoint.size2() == 1);
	
	ublas::matrix<double> radisofCylinder(point - prepoint);
	int n = (int)radisofCylinder.size1();
	ublas::matrix_range<ublas::matrix<double> > sub(radisofCylinder, ublas::range(0, n), ublas::range(0, 1));
	
	sub *= sqrt(sigma);

	return pow(ublas::norm_frobenius(sub), 2.0) - delta*delta;
}

ublas::matrix<double> Problem1::AdditionalEqnPartial(
	ublas::matrix<double>& point,
	ublas::matrix<double>& prepoint,
	ublas::matrix<double>& tanvect,
	ublas::matrix<double>& verticalpoint,
	double& delta,
	double& sigma) const //超圆柱面切割
{
	ublas::matrix<double> radisofCylinder(point - prepoint);
	radisofCylinder *= 2.0;
	int n = (int)radisofCylinder.size1();
	ublas::matrix_range<ublas::matrix<double> > sub(radisofCylinder, ublas::range(0, n), ublas::range(0, 1));
	sub *= sigma;

	radisofCylinder(n - 1, 0) = 0.0;

	return radisofCylinder;
}

/// <summary>
/// CTargetFun的派生类 minusProblem1.
/// </summary>
class minusProblem1 :public CTargetFun
{
public:
	minusProblem1();
	~minusProblem1();
	virtual ublas::matrix<double> Homomap(ublas::matrix<double>& point) const;
	virtual ublas::matrix<double> HomomapJacobi(ublas::matrix<double>& point) const;
	virtual ublas::matrix<double> HomomapT(ublas::matrix<double>& point) const;
	virtual double AdditionalEqn(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const;
	virtual ublas::matrix<double> AdditionalEqnPartial(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const;

};

minusProblem1::minusProblem1()
{
}

minusProblem1::~minusProblem1()
{
}

ublas::matrix<double> minusProblem1::Homomap(ublas::matrix<double>& point) const
{
	ublas::matrix<double> a(2, 1);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double y0 = m_initpoint(1, 0);
	double x = point(0, 0);
	double y = point(1, 0);
	double t = point(2, 0);

	a(0, 0) = (-1*t)*(x*x + x*y - 3) + (1 - t)*(x - x0);
	a(1, 0) = (t)*(2 * sin(x) + y - 7) + (1 - t)*(y - y0);
	return a;
}

ublas::matrix<double> minusProblem1::HomomapJacobi(ublas::matrix<double>& point) const
{
	ublas::matrix<double> a(2, 2);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double y0 = m_initpoint(1, 0);
	double x = point(0, 0);
	double y = point(1, 0);
	double t = point(2, 0);
	a(0, 0) = (-1*t)*(2 * x + y) + 1 - t;
	a(0, 1) = (-1*t)*x;

	a(1, 0) = (t)* 2 * cos(x);	
	a(1, 1) = 1;
	return a;
}

ublas::matrix<double> minusProblem1::HomomapT(ublas::matrix<double>& point) const
{

	ublas::matrix<double> a(2, 1);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double y0 = m_initpoint(1, 0);
	double x = point(0, 0);
	double y = point(1, 0);
	double t = point(2, 0);

	a(0, 0) = -1 * (x*x + x*y - 3) - 1 * (x - x0);
	a(1, 0) = 1 * (2 * sin(x) + y - 7) - 1 * (y - y0);
	return a;
}

double minusProblem1::AdditionalEqn(
	ublas::matrix<double>& point,
	ublas::matrix<double>& prepoint,
	ublas::matrix<double>& tanvect,
	ublas::matrix<double>& verticalpoint,
	double& delta,
	double& sigma) const //超球面切割
{

	assert(point.size1() == tanvect.size1() && tanvect.size1() == verticalpoint.size1());
	assert(point.size2() == 1 && tanvect.size2() == 1 && verticalpoint.size2() == 1);
	
	ublas::matrix<double> radisofsphere(point - prepoint);
	int n = (int)radisofsphere.size1();
	ublas::matrix_range<ublas::matrix<double> > sub(radisofsphere, ublas::range(0, n), ublas::range(0, 1));
	sub *= sqrt(sigma);

	return pow(ublas::norm_frobenius(radisofsphere), 2.0) - delta*delta;
}

ublas::matrix<double> minusProblem1::AdditionalEqnPartial(
	ublas::matrix<double>& point,
	ublas::matrix<double>& prepoint,
	ublas::matrix<double>& tanvect,
	ublas::matrix<double>& verticalpoint,
	double& delta,
	double& sigma) const //超球面切割
{

	ublas::matrix<double> radisofsphere(point - prepoint);
	radisofsphere *= 2.0;
	int n = (int)radisofsphere.size1();
	ublas::matrix_range<ublas::matrix<double> > sub(radisofsphere, ublas::range(0, n), ublas::range(0, 1));
	sub *= sigma;
	return radisofsphere;
}
/// <summary>
/// CTargetFun的派生类 Snap.
/// </summary>
class Snap :public CTargetFun
{
public:
	Snap();
	~Snap();
	virtual ublas::matrix<double> Homomap(ublas::matrix<double>& point) const;
	virtual ublas::matrix<double> HomomapJacobi(ublas::matrix<double>& point) const;
	virtual ublas::matrix<double> HomomapT(ublas::matrix<double>& point) const;
	virtual double AdditionalEqn(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const;
	virtual ublas::matrix<double> AdditionalEqnPartial(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const;
public:
	double para = 1.0;

};

Snap::Snap()
{
}

Snap::~Snap()
{
}


ublas::matrix<double> Snap::Homomap(ublas::matrix<double>& point) const
{

	ublas::matrix<double> a(1, 1);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double x = point(0, 0);
	double t = point(1, 0);
	double phi = para * (1 + t);

	a(0, 0) = (2 - 4 * x + pow(x, 2)) - (1 - t)*(2 - 4 * x0 + pow(x0, 2));

	//a(0, 0) = (2 - 4 * x + pow(x, 1)) - (1 - t)*(2 - 4 * x0 + pow(x0, 1));

	a(0, 0) = 8 - 11 * pow(x, 2) + 6 * pow(x, 4) - pow(x, 6) + (t - 1)*(8 - 11 * pow(x0, 2) + 6 * pow(x0, 4) - pow(x0, 6));
	
	a(0, 0) = 8 - (11 * pow(x, 2)) / 100. + (3 * pow(x, 4)) / 5000. - pow(x, 6) / 1.e6 + (t - 1)*(8 - (11 * pow(x0, 2)) / 100. + (3 * pow(x0, 4)) / 5000. - pow(x0, 6) / 1.e6);

	a(0, 0) = x * sin(x) - 10 + (t - 1)*(x0*sin(x0) - 10);

	a(0, 0) = 0.3*x + sin(x) + (t - 1)*(0.3*x0 + sin(x0));//牛顿同伦

	a(0, 0) = (1 - t)*(x - x0) + t*(0.3*x + sin(x));//不动点同伦

	//a(0, 0) = 0.3*x + sin(x) - (1 - t)*(0.3*x0 + sin(x0) - phi*(x - x0));//自适应同伦
	return a;
}


ublas::matrix<double> Snap::HomomapJacobi(ublas::matrix<double>& point) const
{

	ublas::matrix<double> a(1, 1);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double x = point(0, 0);
	double t = point(1, 0);
	double phi = para * (1 + t);

	a(0, 0) = 2 * (-2 + x);
	//a(0, 0) = -3.0;
	a(0, 0) = -22 * x + 24 * pow(x, 3) - 6 * pow(x, 5);
	a(0, 0) = (-11 * x) / 50. + (3 * pow(x, 3)) / 1250. - (3 * pow(x, 5)) / 500000.;
	a(0, 0) = x*cos(x) + sin(x);

	a(0, 0) = 0.3 + cos(x);//牛顿同伦

	a(0, 0) = 1 - t + 0.3*t + t*cos(x);//不动点同伦

	//a(0, 0) = 0.3 + cos(x) + (1 - t)* phi;//自适应同伦
	return a;
}


ublas::matrix<double> Snap::HomomapT(ublas::matrix<double>& point) const
{
	ublas::matrix<double> a(1, 1);
	a.clear();
	double x0 = m_initpoint(0, 0);
	double x = point(0, 0);
	double t = point(1, 0);
	double phi = para * (1 + t);

	a(0, 0) = 2 - 4 * x0 + pow(x0, 2);
	//a(0, 0) = 2 - 4 * x0 + pow(x0, 1);
	a(0, 0) = 8 - 11 * pow(x0, 2) + 6 * pow(x0, 4) - pow(x0, 6);
	a(0, 0) = 8 - (11 * pow(x0, 2)) / 100. + (3 * pow(x0, 4)) / 5000. - pow(x0, 6) / 1.e6;
	a(0, 0) = x0*sin(x0) - 10;

	a(0, 0) = 0.3*x0 + sin(x0);//牛顿同伦

	a(0, 0) = -1 * (x - x0) + 0.3*x + sin(x);//不动点同伦

	//a(0, 0) = 0.3*x0 + sin(x0) - phi*(x - x0);//自适应同伦
 	return a;
}


double Snap::AdditionalEqn(
	ublas::matrix<double>& point,
	ublas::matrix<double>& prepoint,
	ublas::matrix<double>& tanvect,
	ublas::matrix<double>& verticalpoint,
	double& delta,
	double& sigma) const
{
	//assert(point.size1() == tanvect.size1() && tanvect.size1() == verticalpoint.size1());
	//assert(point.size2() == 1 && tanvect.size2() == 1 && verticalpoint.size2() == 1);
	//ublas::matrix<double> tempvec1(tanvect);
	//ublas::matrix<double> tempvec2(point - verticalpoint);

	//ublas::matrix_column<ublas::matrix<double> > vec1column(tempvec1, 0);
	//ublas::matrix_column<ublas::matrix<double> > vec2column(tempvec2, 0);
	//return inner_prod(vec1column, vec2column);

	//assert(point.size1() == tanvect.size1() && tanvect.size1() == verticalpoint.size1());
	//assert(point.size2() == 1 && tanvect.size2() == 1 && verticalpoint.size2() == 1);

	//ublas::matrix<double> radisofsphere(point - prepoint);
	//int n = (int)radisofsphere.size1();
	//ublas::matrix_range<ublas::matrix<double> > sub(radisofsphere, ublas::range(0, n), ublas::range(0, 1));
	//sub *= sqrt(sigma);

	//return pow(ublas::norm_frobenius(radisofsphere), 2.0) - delta*delta;


	assert(point.size1() == tanvect.size1() && tanvect.size1() == verticalpoint.size1());
	assert(point.size2() == 1 && tanvect.size2() == 1 && verticalpoint.size2() == 1);

	ublas::matrix<double> radisofCylinder(point - prepoint);
	int n = (int)radisofCylinder.size1();
	ublas::matrix_range<ublas::matrix<double> > sub(radisofCylinder, ublas::range(0, n), ublas::range(0, 1));

	sub *= sqrt(sigma);

	return pow(ublas::norm_frobenius(sub), 2.0) - delta*delta;
}

ublas::matrix<double> Snap::AdditionalEqnPartial(
	ublas::matrix<double>& point,
	ublas::matrix<double>& prepoint,
	ublas::matrix<double>& tanvect,
	ublas::matrix<double>& verticalpoint,
	double& delta,
	double& sigma) const
{
	//ublas::matrix<double> a(tanvect);
	//return a;

	//ublas::matrix<double> radisofsphere(point - prepoint);
	//radisofsphere *= 2.0;
	//int n = (int)radisofsphere.size1();
	//ublas::matrix_range<ublas::matrix<double> > sub(radisofsphere, ublas::range(0, n), ublas::range(0, 1));
	//sub *= sigma;
	//return radisofsphere;

	ublas::matrix<double> radisofCylinder(point - prepoint);
	radisofCylinder *= 2.0;
	int n = (int)radisofCylinder.size1();
	ublas::matrix_range<ublas::matrix<double> > sub(radisofCylinder, ublas::range(0, n), ublas::range(0, 1));
	sub *= sigma;

	radisofCylinder(n - 1, 0) = 0.0;

	return radisofCylinder;
}