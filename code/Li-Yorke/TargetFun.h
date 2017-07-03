#pragma once
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

namespace ublas = boost::numeric::ublas;

/// <summary>
/// ͬ�׺�������(��У��ʱ�õ���Լ��),ʵ��ʹ��ʱ��Ҫ����.
/// </summary>
class CTargetFun
{
public:
	CTargetFun();
	~CTargetFun();
	void Setinit(ublas::matrix<double>& initpoint);
	/// <summary>
	/// ĳ�㴦��ͬ�׺���ֵ.
	/// </summary>
	/// <param name="point">�����.</param>
	/// <returns>���غ���ֵ<double>.</returns>
	virtual ublas::matrix<double> Homomap(ublas::matrix<double>& point) const = 0;
	/// <summary>
	/// ĳ�㴦��ͬ�׹����Ա���x������ֵ.
	/// </summary>
	/// <param name="point">�����.</param>
	/// <returns>���ص�����ֵ<double>.</returns>
	virtual ublas::matrix<double> HomomapJacobi(ublas::matrix<double>& point) const = 0;
	/// <summary>
	/// ĳ�㴦��ͬ�׹���t������ֵ.
	/// </summary>
	/// <param name="point">�����.</param>
	/// <returns>���ص�����ֵ<double>.</returns>
	virtual ublas::matrix<double> HomomapT(ublas::matrix<double>& point) const = 0;
	/// <summary>
	/// ���䷽��.���������ܶ���������Ա�Ӧ�Բ�ͬ��ʽ���иʽ.
	/// </summary>
	/// <param name="point">��ǰ��.</param>
	/// <param name="prepoint">ǰһ����.</param>
	/// <param name="tanvect">��ǰ���������.</param>
	/// <param name="verticalpoint">����.</param>
	/// <param name="delta">��Բ�����и�ʱ��Ӧ����Բ�İ뾶.</param>
	/// <param name="sigma">�������и�ʱ��Ӧ�뾶.</param>
	/// <returns>double.</returns>
	virtual double AdditionalEqn(
		ublas::matrix<double>& point,
		ublas::matrix<double>& prepoint,
		ublas::matrix<double>& tanvect,
		ublas::matrix<double>& verticalpoint,
		double& delta,
		double& sigma) const = 0; 
	/// <summary>
	/// ���䷽�̹���x,t�ĵ���.���������ܶ���������Ա�Ӧ�Բ�ͬ��ʽ���иʽ.
	/// </summary>
	/// <param name="point">��ǰ��.</param>
	/// <param name="prepoint">ǰһ����.</param>
	/// <param name="tanvect">��ǰ���������.</param>
	/// <param name="verticalpoint">����.</param>
	/// <param name="delta">��Բ�����и�ʱ��Ӧ����Բ�İ뾶.</param>
	/// <param name="sigma">�������и�ʱ��Ӧ�뾶.</param>
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

