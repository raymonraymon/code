#pragma once
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

namespace ublas = boost::numeric::ublas;

/// <summary>
/// 同伦函数基类(含校正时用到的约束),实际使用时需要派生.
/// </summary>
class CTargetFun
{
public:
	CTargetFun();
	~CTargetFun();
	void Setinit(ublas::matrix<double>& initpoint);
	/// <summary>
	/// 某点处的同伦函数值.
	/// </summary>
	/// <param name="point">输入点.</param>
	/// <returns>返回函数值<double>.</returns>
	virtual ublas::matrix<double> Homomap(ublas::matrix<double>& point) const = 0;
	/// <summary>
	/// 某点处的同伦关于自变量x导函数值.
	/// </summary>
	/// <param name="point">输入点.</param>
	/// <returns>返回导函数值<double>.</returns>
	virtual ublas::matrix<double> HomomapJacobi(ublas::matrix<double>& point) const = 0;
	/// <summary>
	/// 某点处的同伦关于t导函数值.
	/// </summary>
	/// <param name="point">输入点.</param>
	/// <returns>返回导函数值<double>.</returns>
	virtual ublas::matrix<double> HomomapT(ublas::matrix<double>& point) const = 0;
	/// <summary>
	/// 补充方程.参数尽可能多的罗列上以便应对不同形式的切割方式.
	/// </summary>
	/// <param name="point">当前点.</param>
	/// <param name="prepoint">前一个点.</param>
	/// <param name="tanvect">当前点的切向量.</param>
	/// <param name="verticalpoint">垂点.</param>
	/// <param name="delta">超圆柱面切割时对应底面圆的半径.</param>
	/// <param name="sigma">超球面切割时对应半径.</param>
	/// <returns>double.</returns>
	virtual double AdditionalEqn(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const = 0; 
	/// <summary>
	/// 补充方程关于x,t的导数.参数尽可能多的罗列上以便应对不同形式的切割方式.
	/// </summary>
	/// <param name="point">当前点.</param>
	/// <param name="prepoint">前一个点.</param>
	/// <param name="tanvect">当前点的切向量.</param>
	/// <param name="verticalpoint">垂点.</param>
	/// <param name="delta">超圆柱面切割时对应底面圆的半径.</param>
	/// <param name="sigma">超球面切割时对应半径.</param>
	/// <returns>double.</returns>
	virtual ublas::matrix<double> AdditionalEqnPartial(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const = 0;
protected:
	ublas::matrix<double> m_initpoint;
};

