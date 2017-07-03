// ***********************************************************************
// Assembly         : Homotopy
// Author           : chenruishan
// Created          : 04-23-2015
//
// Last Modified By : chenruishan
// Last Modified On : 04-29-2015
// ***********************************************************************
// <copyright file="LiYorke.h" company="china">
//     Copyright (c) china. All rights reserved.
// </copyright>
// <summary>Li-Yorke算法实现</summary>
//***********************************************************************
/*! \page Li-Yorke Li-Yorke算法实现
*Li-Yorke算法实现\n
*有一个主要的对外接口：<see cref="CLiYorke::Solve()"/>求解。\n
*\section cited 参考文献
*[1] 黄象鼎,曾钟钢,马亚南, 非线性数值分析[O], 武汉: 武汉大学出版社, 2004.\n
*[2] 张丽琴,王家映, 地球物理资料非线性反演方法讲座(七) 同伦反演法[J], 工程地球物理学报, 2008, 5:509∼515.\n
*[3] Arie de Niet, Step-size control and corrector methods in numerical continuation of ocean circulation and fill-reducing orderings in multilevel ILU methods [Master’s Thesis] Department of Mathematics University of Groningen, (Aug., 2002)\n
*/
//***********************************************************************
#pragma once

#include <vector>
#include <map>
#include <math.h>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include "TargetFun.h"
#include "SFrame\linalg.h"

namespace ublas = boost::numeric::ublas;

/// <summary>
/// Li-Yorke算法求解同伦延拓问题.
/// </summary>
class CLiYorke
{
private:
	/// <summary>
	/// 指向目标函数抽象类的指针
	/// </summary>
	const CTargetFun *m_pt; 
	/// <summary>
	/// 迭代初始值
	/// </summary>
	ublas::matrix<double> m_initPoint;
	/// <summary>
	/// 牛顿迭代的误差	
	/// </summary>
	double m_epsilon;
	/// <summary>
	/// 初始步长
	/// </summary>
	double m_initStepsize;
	/// <summary>
	/// 无限性半径
	/// </summary>
	double m_infRadius;
	/// <summary>
	/// 拟弧长参数	
	/// </summary>
	double m_sigma; 
	/// <summary>
	/// 上一步的切向量
	/// </summary>
	ublas::matrix<double> m_lastTanVec;
	/// <summary>
	/// 机器精度 
	/// </summary>
	double m_precision;
	/// <summary>
	/// 雅可比行列式的符号
	/// </summary>
	int m_sgndexHx;
	/// <summary>
	/// 最优迭代次数
	/// </summary>
	int m_optIter;

	//std::map<int, double> m_delta;//步长
	//std::map<int, double> m_angle;//步长
	//std::map<int, std::vector<double>> m_points;//道路上的点

public:
	CLiYorke();
	~CLiYorke();

private:
	//主干子函数
	//1计算切向量

	bool TangentVec1(const ublas::matrix<double>& point,
		ublas::matrix<double>& TangentVector);//迭代求解切向量,张丽琴

	//黄象鼎书中计算切向量的方法
	//计算切向量时需要用到的5个子函数
	//1.0计算符号
	int Tansign(const double& val);
	//1.1计算jacobi矩阵在初值点point处的行列式
	double	TanDet(const ublas::matrix<double>& Mat);
	//1.2计算主元列序号置换的符号的 相反数
	void Tanswap(int& a, int& b);
	int	TanBubbleCount(std::vector<int> &noSortedVector);
	//1.3计算Householder中间结果的主元所在列序号
	int	TanMaxColumn(const ublas::matrix<double>& matrixAfterTansform);
	//1.4Householder变换
	bool TanHouseholder(ublas::matrix<double>& A,
		const int& num,//第num次变换
		const int& k,
		ublas::matrix<double>& Q);
	bool TangentVec(const ublas::matrix<double>& point,
		ublas::matrix<double>& TangentVector);//QR分解精确求解切向量,黄象鼎改良版

	bool FirstTangentVec(const ublas::matrix<double>& point,
		ublas::matrix<double>& TangentVector);//deNiet硕士论文

	//2预估
	bool EulerPredict(const ublas::matrix<double>& currentPoint,
		const ublas::matrix<double>& tangentVector,
		double& delta, //最后一步要传回计算出来的步长，所以不能const         
		ublas::matrix<double>& PredictRes);
	//3.1子函数 牛顿迭代中遇到的矩阵求解
	bool NewtonItergen(
		const ublas::matrix<double>& movingPoint,
		const ublas::matrix<double>& previousPoint,
		const ublas::matrix<double>& predictPoint,
		double& stepsize,
		ublas::matrix<double>& jacobi,
		ublas::matrix<double>& rhs);


	//3牛顿迭代
	bool NewtonIter(
		const ublas::matrix<double>& previousPoint,
		const ublas::matrix<double>& predictPoint,
		double& stepsize,
		int &k,//迭代次数
		ublas::matrix<double>& ConvergentPoint);

	//4计算前后切向量的夹角
	bool Angle(const ublas::matrix<double>& vec1,
		const ublas::matrix<double>& vec2,
		double& ang);

	bool StepSolve(const ublas::matrix<double>& StartPoint,
		double& stepsize,
		ublas::matrix<double>& endPoint,
		double& stepangle);
	//对外接口 
	//用函数指针读入用户输入函数。
	//用纯虚函数基类来派生读入。
public:
	bool Solve(
		const CTargetFun *pt,
		const ublas::matrix<double>& initPoint,
		const double& epsilon,
		const double& initDelta,
		const double& infRadius,
		std::map<int, double>& stepsizemap,
		std::map<int, double>& anglemap,
		std::map<int, ublas::matrix<double>>& pointsmap);
};

