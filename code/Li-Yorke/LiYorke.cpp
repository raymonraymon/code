//#include "stdafx.h"
#include "LiYorke.h"
#include <iomanip>
#include <iosfwd>

#define PI 3.141592653589793

using namespace std;
using namespace boost::numeric::ublas;

/// <summary>
/// 构造 <see cref="CLiYorke"/> 类.
/// </summary>
CLiYorke::CLiYorke()
{
}


/// <summary>
/// 析构 <see cref="CLiYorke"/> 类.
/// </summary>
CLiYorke::~CLiYorke()
{
}

/// <summary>
/// LI-YORKE算法求解同伦延拓.
/// </summary>
/// <param name="pt">同伦函数的指针.</param>
/// <param name="initPoint">起点.</param>
/// <param name="epsilon">误差限.</param>
/// <param name="initDelta">初始步长.</param>
/// <param name="infRadius">无限性半径.</param>
/// <param name="stepsizemap">求解路径上的步长.</param>
/// <param name="anglemap">求解路径上的夹角.</param>
/// <param name="pointsmap">求解路径上的点.</param>
/// <returns>bool.</returns>
bool CLiYorke::Solve(
	const CTargetFun *pt,
	const matrix<double>& initPoint,
	const double& epsilon,
	const double& initStepsize,
	const double& infRadius,
	map<int, double>& stepsizemap,
	map<int, double>& anglemap,
	map<int, matrix<double>>& pointsmap)
{
	cout << setprecision(20) << endl;

	m_pt = pt;
	
	m_initPoint = initPoint;
	matrix<double> initPt(m_initPoint);

	matrix<double> Matjacobi0 = m_pt->HomomapJacobi(initPt);
	m_sgndexHx = Tansign(TanDet(Matjacobi0));

	m_epsilon = epsilon;
	m_initStepsize = initStepsize;
	m_infRadius = infRadius;

	m_sigma = 1.0 / 1.0;//m_sigma =1/dim
	m_precision = 1E-7;
	m_optIter = 5;

	/************************************************************************/
	/*                           初始切向量赋初值【0，1】                     */
	/************************************************************************/
	int dim = (int)initPt.size1();
	matrix<double> lastTanVec(dim, 1);
	lastTanVec.clear();
	for (int i = 0; i < dim - 1; ++i)
	{
		lastTanVec(i, 0) = 0.0;
	}
	lastTanVec(dim - 1, 0) = 1.0; //张丽琴推荐用【0，1】代替初始切向量，但因为只需要计算一次，用deNiet的方法可以保证精确。

	if (FirstTangentVec(initPt, lastTanVec) == false)
	{
		return false;
	}
	//cout << "deNiet第一个切向量:" << lastTanVec << endl;
	
	m_lastTanVec = lastTanVec;

	double maxAngle = 18.0;
	double minAngle = 6.0;

	//matrix<double> Matjacobi0 = m_pt->HomomapJacobi(initPt);

	int n = (int)m_initPoint.size1() - 1;

	matrix<double> endPt(initPt); //赋给初值
	double ang;

	double stepsize = m_initStepsize;//初始化步长
	int k = 0;

	pointsmap.insert(make_pair(k, endPt));


	cout << "temp道路上的点：" << endl;

	while (abs(endPt(n, 0) - 1.0) > m_precision || 
		norm_frobenius(m_pt->Homomap(endPt)) > m_precision)
	{
		if (k > 1050)
		{
			//cout << "预估校正超过最大迭代次数！" << endl;
			return false;
		}
		if (StepSolve(initPt, stepsize, endPt, ang) == false)
		{
			return false;
		}

		if (ang >= 90.0) //大于90度,切向量计算有误
		{
			cout << "角度大于90度！" << endl;
			return false;
		}

		if (ang >= maxAngle) //大于18度，上步结果不加入输出map,折回预估
		{
			stepsize /= 2.0;
		}
		else  //上步结果加入输出map，
		{
			++k;
			initPt = endPt;
			stepsizemap.insert(make_pair(k, stepsize));
			anglemap.insert(make_pair(k, ang));
			pointsmap.insert(make_pair(k, endPt));

            cout << setw(2) << k << ",";
			for (int i = 0; i < n + 1; ++i)
			{
                cout << setw(20) << endPt(i, 0) << ",";
			}
            cout << endl;

			if (ang < minAngle)//介于6-18度，步长不变继续前进
			{
				stepsize *= 2.0;
			}
		}
		//cout << "endPt[t]:" << endPt(n, 0) << endl;
	}
	//cout << "迭代次数:" << k << endl;

	double error = norm_frobenius(m_pt->Homomap(initPt));
	//cout << "误差:" << error << endl;
	assert(error < m_precision);


    cout << setprecision(ios::scientific) << setprecision(20);
    cout << "道路上的点：" << endl;
	for (std::map<int, matrix<double>>::iterator pos = pointsmap.begin();
		pos != pointsmap.end();
		++pos)
	{
        cout << setw(2) << pos->first << ",";
		for (int i = 0; i < n + 1; ++i)
		{
            cout << setw(20) << pos->second(i, 0) << ",";
		}
        cout << endl;
	}
	return true;
}

/// <summary>
/// 判断正负号和零.
/// </summary>
/// <param name="val">输入值.</param>
/// <returns>int.</returns>
int CLiYorke::Tansign(const double& val)
{
	if (abs(val) <= m_precision)
	{
		return 0;
	}
	else if (val < -1 * m_precision)
	{
		return -1;
	}
	else
	{
		return 1;
	}
}

/// <summary>
/// 递归计算行列式.
/// </summary>
/// <param name="Mat">待求矩阵.</param>
/// <returns>行列式值.</returns>
double	CLiYorke::TanDet(const matrix<double>& Mat)
{

	assert(Mat.size1() == Mat.size2());
	int m = Mat.size1();
	int len = m - 1;/*余子式的阶*/
	double s = 0;

	/*按照定义，初始化一个子行列式数组的空间*/
	ublas::matrix<double> p(Mat.size1() - 1, Mat.size1() - 1);

	/*阶为1，按照定义计算*/
	if (1 == m)
		return Tansign(Mat(0, 0));

	for (int k = 0; k < m; ++k)
	{
		for (int i = 0; i < len; ++i)
			for (int j = 0; j < len; ++j)
			{
				if (i < k)
					p(i, j) = Mat(i, j + 1);/*初始化子行列式的值*/
				if (i >= k)
					p(i, j) = Mat(i + 1, j + 1);
			}
		s += (double)pow(-1.0, k) * Mat(k, 0) * TanDet(p);/*递归计算*/
	}
	//cout << "行列式：" << s << endl;
	return s;
}

/// <summary>
/// 交换两个int值.
/// </summary>
/// <param name="a">待交换的输入1.</param>
/// <param name="b">待交换的输入2.</param>
void CLiYorke::Tanswap(int& a, int& b)
{
	int temp;
	temp = a;
	a = b;
	b = temp;
}
/// <summary>
/// 计算交换顺序的次数.
/// </summary>
/// <param name="noSortedVector">输入排序前的数组(对应列号)\f$ sgn(a_{i,k_{i}}) \f$.</param>
/// <returns>int.</returns>
int	CLiYorke::TanBubbleCount(std::vector<int> &noSortedVector)
{
	int length = noSortedVector.size();
	int k = 0;
	for (int i = 0; i < length; i++)
	{
		for (int j = length - 1; j > i; j--)
		{
			if (noSortedVector[j] < noSortedVector[j - 1])
			{
				Tanswap(noSortedVector[j], noSortedVector[j - 1]);
				k++;
			}

		}
	}
	return k;
}

/// <summary>
/// 计算2范数最大的子列的列号.
/// </summary>
/// <param name="matrixAfterTansform">输入矩阵.</param>
/// <returns>int.</returns>
int	CLiYorke::TanMaxColumn(const matrix<double>& matrixAfterTansform)
{
	int m = matrixAfterTansform.size1();
	int n = matrixAfterTansform.size2();

	std::vector<double> a(n);
	a.clear();

	for (int i = 0; i < n; ++i)
	{
		a.push_back(0.0);
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			a[i] += pow(abs(matrixAfterTansform(j, i)), 2);
		}
	}

	int k = 0;
	for (int i = 0; i < n; ++i)
	{
		if (a[i]>a[k])
		{
			k = i; //如有两列同为最大范数，返回后一个列
		}
	}
	return k;
}
/// <summary>
/// householder变换.
/// </summary>
/// <param name="A">输入矩阵.</param>
/// <param name="num">第num次变换.</param>
/// <param name="k">A的列数减去2范数最大的子列的列号.</param>
/// <param name="Q">输出矩阵.</param>
/// <returns>bool.</returns>
bool CLiYorke::TanHouseholder(matrix<double>& A,
	const int& num,//第num次变换
	const int& k,
	matrix<double>& Q)
{
	int n = (int)A.size1();
	int m = (int)A.size2();

	matrix<double> xminusy(n, 1);

	for (int i = 0; i < num - 1; i++)
	{
		xminusy(i, 0) = 0;
	}
	for (int i = num - 1; i < n; i++)
	{
		xminusy(i, 0) = A(i, m - k);
	}

	//matrix_range<matrix<double> > xminusy1(xminusy, range(num, n), range(0, 1));

	xminusy(num - 1, 0) += Tansign(xminusy(num - 1, 0)) *  norm_frobenius(xminusy);

	double d = norm_frobenius(xminusy);
	if (abs(d) < m_precision)
	{
		return false;
	}

	for (int i = 0; i < n; i++)
	{
		xminusy(i, 0) = xminusy(i, 0) / d;
	}
	matrix<double> H(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			H(i, j) = -2 * xminusy(i, 0) * xminusy(j, 0);
		}
		H(i, i)++;
	} //H_k=I-2vvT 

	matrix<double> temp1(n, m);
	temp1 = prod(H, A);
	A = temp1;

	matrix<double> temp2(n, n);
	temp2 = prod(Q, H);
	Q = temp2;
	return true;
}

/// <summary>
/// QR分解精确求解切向量,黄象鼎改良版.
/// </summary>
/// <param name="point">输入点坐标.</param>
/// <param name="TangentVector">输入点的切向量.</param>
/// <returns>bool.</returns>
bool CLiYorke::TangentVec(const matrix<double>& point,
	matrix<double>& TangentVector)
{
	matrix<double> pt(point);

	matrix<double> Matjacobi = m_pt->HomomapJacobi(pt);
	matrix<double> MatHt = m_pt->HomomapT(pt);

	int n = (int)Matjacobi.size1();
	matrix<double> MatJacobiT(n, n + 1);

	matrix_range<matrix<double> > mvr1(MatJacobiT, range(0, n), range(0, n));
	mvr1 = Matjacobi;

	//std::cout << "mvr1:" << mvr1 << std::endl;
	matrix_range<matrix<double> > mvr2(MatJacobiT, range(0, n), range(n, n + 1));
	mvr2 = MatHt;

	//std::cout << "mvr2:" << mvr2 << std::endl;
	assert(m_sigma > 0.0 && m_sigma <= 1.0);

	matrix<double> QRF(MatJacobiT);

	//std::cout << "QRF:" << QRF << std::endl;
	int m = n + 1;

	matrix<double> Q(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Q(i, j) = 0;
		}
		Q(i, i) = 1;
	}

	//for (int i = n; i >= 2; i--)
	//{
	//	if (householder(AA, i, QQ) == false)
	//	{
	//		return false;
	//	}		
	//}

	matrix<double> BB1;
	std::vector<int> maxCol;//记录主元排列
	for (int i = 0; i < n - 1; ++i)
	{
		BB1 = QRF;
		for (int k = 0; k < i; ++k)
		{
			for (int j = 0; j < m; ++j)
			{
				BB1(k, j) = 0;
			}
		}
		maxCol.push_back(TanMaxColumn(BB1));

		if (TanHouseholder(QRF, i + 1, m - maxCol[i], Q) == false)
		{
			return false;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	//cout << "Q:" << Q << endl;
	//cout << "Q_norm:" << pow(norm_frobenius(prod(Q, trans(Q))),2) << endl;
	//cout << "R:" << QRF << endl;
	//cout << "MatJacobiT:" << MatJacobiT << endl;
	//cout << "error:" << norm_frobenius(prod(Q, QRF) - MatJacobiT) << endl;

	assert(norm_frobenius(prod(Q, QRF) - MatJacobiT) < m_precision);

	matrix<double> CC(QRF);
	for (int i = 0; i < n - 1; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			CC(i, j) = 0; //还要加入一个循环
		}
	}

	maxCol.push_back(TanMaxColumn(CC));

	//剩下的那列的序号怎么得到？
	int lastColumn;
	for (int i = 0; i < m; ++i)
	{
		std::vector<int>::iterator result = find(maxCol.begin(), maxCol.end(), i); //查找3
		if (result == maxCol.end())//没找到
		{
			lastColumn = i;
			break;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	//解增广矩阵得到初始切向量
	matrix<double> extendedMatrix(n + 1, m);
	extendedMatrix.clear();
	matrix_range<matrix<double> > mrcopyA(extendedMatrix, range(0, n), range(0, m));
	mrcopyA = QRF;

	for (int j = 0; j < m; ++j)
	{
		extendedMatrix(n, j) = 0.0;
	}
	extendedMatrix(n, lastColumn) = 1;

	matrix<double> theRightHand(n + 1, 1);
	theRightHand.clear();
	for (int i = 0; i < n; ++i)
	{
		theRightHand(i, 0) = 0.0;
	}
	//theRightHand(n, 0) = ?;//两个符号函数的乘积,
	int k = 1;
	for (int i = 0; i < n; ++i)
	{
		//cout << "QRF("<<i<<","<< maxCol[i]<<"):" << QRF(i, maxCol[i]) << endl;

		k *= Tansign(QRF(i, maxCol[i]));	//找到每列零上面的那个元
	}

	//matrix<double> Matjacobi0 = m_pt->HomomapJacobi(pt);
	//m_sgndexHx = Tansign(TanDet(Matjacobi0));

	int powtimes = TanBubbleCount(maxCol);
	double sgntime = pow(-1.0, powtimes + 1);

	theRightHand(n, 0) = 1;// sgntime * k *  m_sgndexHx;
	/*如果用与前一个切向量的夹角来选择下一个切向量的方向，这一步就可以省略*/

	matrix<double> unKnown = trans(theRightHand);


	//if (linalg::Solve(extendedMatrix, unKnown) == 0)
	//{
	//	return false;
	//}
    //unKnown = trans(unKnown);


    int m_Free_degree = n + 1;

    permutation_matrix<double> P(m_Free_degree);
    boost::numeric::ublas::vector<double> x(m_Free_degree);
    boost::numeric::ublas::vector<double> v(m_Free_degree);

    for (int j = 0; j < m_Free_degree; ++j)
    {
        v(j) = unKnown(0, j);
    }

    //try
    {
        lu_factorize(extendedMatrix, P);
        x = v;
        lu_substitute(extendedMatrix, P, x);
        cout << "解向量: " << x << endl;
        for (int j = 0; j < m_Free_degree; ++j)
        {
            unKnown(0, j) = x(j);
        }
        unKnown = trans(unKnown);
    }

	/************************************************************************/
	/* 记下上一个切向量                                                       */
	/************************************************************************/

	

	matrix_column<matrix<double> > vec1column(unKnown, 0);
	matrix_column<matrix<double> > vec2column(m_lastTanVec, 0);
	double anglecos = inner_prod(vec1column, vec2column);

	//////////////////////////////////////////////////////////////////////////
	if (anglecos < 0.0) //保证切向量中的t永远为正。黄象鼎的算法实现有误需要查
	{
		unKnown *= -1.0;
	}
	//////////////////////////////////////////////////////////////////////////

	//引入拟弧长法参数m_sigma
	matrix<double> unKnownSigma;
	unKnownSigma = sqrt(m_sigma)*unKnown;
	unKnownSigma(lastColumn, 0) = unKnown(lastColumn, 0);
	double normUnknown = norm_frobenius(unKnownSigma);
	if (normUnknown < m_precision)
	{
		return false;
	}
	TangentVector = unKnown / normUnknown;
	//cout << "切向量：" << TangentVector << endl;
	return true;

}

/// <summary>
/// 迭代求解切向量，不是精确的切向量，但相比黄象鼎书中的方法，更简单，精度也够.\n
/// 来自张丽琴的文章.
/// </summary>
/// <param name="point">输入点坐标.</param>
/// <param name="TangentVector">输入点的切向量.</param>
/// <returns>bool.</returns>
bool CLiYorke::TangentVec1(const matrix<double>& point,
	matrix<double>& TangentVector)//
{
	matrix<double> pt(point);

	matrix<double> Matjacobi = m_pt->HomomapJacobi(pt);
	matrix<double> MatHt = m_pt->HomomapT(pt);

	int n = (int)Matjacobi.size1();
	matrix<double> MatJacobiT(n + 1, n + 1);

	matrix_range<matrix<double> > mvr1(MatJacobiT, range(0, n), range(0, n));
	mvr1 = Matjacobi;

	//std::cout << "mvr1:" << mvr1 << std::endl;
	matrix_range<matrix<double> > mvr2(MatJacobiT, range(0, n), range(n, n + 1));
	mvr2 = MatHt;

	matrix_range<matrix<double> > mvr3(MatJacobiT, range(n, n + 1), range(0, n + 1));
	mvr3 = trans(m_lastTanVec);


	matrix<double> theRightHand(n + 1, 1);
	theRightHand.clear();
	for (int i = 0; i < n; ++i)
	{
		theRightHand(i, 0) = 0.0;
	}
	theRightHand(n, 0) = 1.0;

	matrix<double> unKnown = trans(theRightHand);
	//if (linalg::Solve(MatJacobiT, unKnown) == 0)
	//{
	//	return false;
	//}
    //unKnown = trans(unKnown);

    int m_Free_degree = n + 1;

    permutation_matrix<double> P(m_Free_degree);
    boost::numeric::ublas::vector<double> x(m_Free_degree);
    boost::numeric::ublas::vector<double> v(m_Free_degree);

    for (int j = 0; j < m_Free_degree; ++j)
    {
        v(j) = unKnown(0, j);
    }

    //try
    {
        lu_factorize(MatJacobiT, P);
        x = v;
        lu_substitute(MatJacobiT, P, x);
        cout << "解向量: " << x << endl;
        for (int j = 0; j < m_Free_degree; ++j)
        {
            unKnown(0, j) = x(j);
        }
        unKnown = trans(unKnown);
    }



	matrix<double> unKnownSigma;
	unKnownSigma = sqrt(m_sigma)*unKnown;
	unKnownSigma(n, 0) = unKnown(n, 0);
	double normUnknown = norm_frobenius(unKnownSigma);
	if (normUnknown < m_precision)
	{
		return false;
	}
	TangentVector = unKnown / normUnknown;

	//cout << "切向量：" << TangentVector << endl;
	return true;
}
/// <summary>
/// deNiet的硕士论文中提到的第一个切向量的计算方法.
/// </summary>
/// <param name="point">输入点坐标.</param>
/// <param name="TangentVector">第一个切向量.</param>
/// <returns>bool.</returns>
bool CLiYorke::FirstTangentVec(const matrix<double>& point,
	matrix<double>& TangentVector)
{
	matrix<double> pt(point);
	matrix<double> Matjacobi = m_pt->HomomapJacobi(pt);
	matrix<double> MatHt = m_pt->HomomapT(pt);
	matrix<double> unKnown = trans(-1 * MatHt);
	//cout << "Matjacobi：" << Matjacobi << endl;
	//if (linalg::Solve(Matjacobi, unKnown) == 0)
	//{
	//	return false;
	//}
	//unKnown = trans(unKnown);

    int m_Free_degree = Matjacobi.size1();

    permutation_matrix<double> P(m_Free_degree);
    boost::numeric::ublas::vector<double> x(m_Free_degree);
    boost::numeric::ublas::vector<double> v(m_Free_degree);

    for (int j = 0; j < m_Free_degree; ++j)
    {
        v(j) = unKnown(0, j);
    }

    //try
    {
        lu_factorize(Matjacobi, P);
        x = v;
        lu_substitute(Matjacobi, P, x);
        cout << "解向量: " << x << endl;
        for (int j = 0; j < m_Free_degree; ++j)
        {
            unKnown(0, j) = x(j);
        }
        unKnown = trans(unKnown);
    }



	double t = 1 / sqrt(1 + m_sigma*pow(norm_frobenius(unKnown), 2));
	unKnown = t * unKnown;
	TangentVector = point;
	int n = (int)TangentVector.size1();
	matrix_range<matrix<double> > mvr1(TangentVector, range(0, n - 1), range(0, 1));
	mvr1 = unKnown;
	TangentVector(n - 1, 0) = t;

	return true;
}

/// <summary>
/// 欧拉法预估.
/// </summary>
/// <param name="currentPoint">出发点.</param>
/// <param name="tangentVector">出发点的切向量.</param>
/// <param name="delta">步长.</param>
/// <param name="PredictRes">预估点.</param>
/// <returns>bool.</returns>
bool CLiYorke::EulerPredict(const matrix<double>& curPoint,
	const matrix<double>& tanVec,
	double& delta, //最后一步要传回计算出来的步长，所以不能const
	matrix<double>& PredictRes)
{
	assert(curPoint.size1() == tanVec.size1() && curPoint.size2() == tanVec.size2());
	//assert(delta > m_precision);

	PredictRes = curPoint + delta * tanVec;

	//应对最后一步
	int n = (int)curPoint.size1();

	//cout <<"中间t值："<< PredictRes(n - 1, 0) << endl;

	//if (PredictRes(n - 1, 0) < 0.0)
	//{
	//	return false;
	//}

	if (PredictRes(n - 1, 0) > 1.0)
	{
		//cout << "到达t=1平面！" << endl;
		if (abs(tanVec(n - 1, 0)) < m_precision)
		{
			return false;
		} 
		else
		{
			delta = (1.000000000 - curPoint(n - 1, 0)) / tanVec(n - 1, 0);
			PredictRes = curPoint + delta * tanVec;
		}
		
	}

	//cout << "当前点：" << curPoint << endl;
	//cout << "切向量：" << tanVec << endl;
	//cout << "步长：" << delta << endl;
	//cout << "欧拉预测点：" << PredictRes << endl;

	return true;
}

/// <summary>
/// 牛顿迭代校正的矩阵生成.参考黄象鼎书中4.28式
/// </summary>
/// <param name="movingPoint">.</param>
/// <param name="previousPoint">上一步的求解点，满足同伦方程.</param>
/// <param name="predictPoint">预估点.</param>
/// <param name="stepsize">步长.</param>
/// <param name="jacobi">雅可比矩阵.</param>
/// <param name="rhs">右端向量.</param>
/// <returns>bool.</returns>
bool CLiYorke::NewtonItergen(
	const matrix<double>& movingPoint,
	const matrix<double>& previousPoint,
	const matrix<double>& predictPoint,
	double& stepsize,
	matrix<double>& jacobi,
	matrix<double>& rhs)
{
	assert(movingPoint.size1() == previousPoint.size1() && movingPoint.size2() == previousPoint.size2());
	assert(previousPoint.size1() == predictPoint.size1() && previousPoint.size2() == predictPoint.size2());
	//补充校验参数的维数

	matrix<double> previousPt(previousPoint);
	matrix<double> predictPt(predictPoint);
	matrix<double> movingPt(movingPoint);

	/*生成jabobi矩阵*/
	matrix<double> jacobi00 = m_pt->HomomapJacobi(movingPt);
	matrix<double> jacobi01 = m_pt->HomomapT(movingPt);

	int n = (int)jacobi00.size1();

	if (predictPt(n, 0) < 0.0 || predictPt(n, 0) > 1.0)
	{
		cout << "预测点的t值必须大于零小于等于1！"<< endl;
		return false;
	}

	matrix<double> Tanvector;
	if (TangentVec(previousPt, Tanvector) == false)
	{
		return false;
	}
	matrix<double> jacobi1(1, n + 1);

	double delta = stepsize;//辅助函数中用到的参数.

	double abstemp = abs(predictPoint(n, 0) - 1.0);

	if (abstemp < m_precision)
	{
		for (int i = 0; i < n; ++i)
		{
			jacobi1(0, i) = 0;
		}
		jacobi1(0, n) = 1;
	}
	else
	{
		jacobi1 = trans(m_pt->AdditionalEqnPartial(movingPt, previousPt, Tanvector, predictPt, delta, m_sigma));
	}

	matrix<double> jacobitemp(n + 1, n + 1);//Target1

	matrix_range<matrix<double> > mvr00(jacobitemp, range(0, n), range(0, n));
	mvr00 = jacobi00;

	matrix_range<matrix<double> > mvr01(jacobitemp, range(0, n), range(n, n + 1));
	mvr01 = jacobi01;

	matrix_range<matrix<double> > mvr1(jacobitemp, range(n, n + 1), range(0, n + 1));
	mvr1 = jacobi1;
	/*生成jabobi矩阵end*/

	/*生成右端矩阵*/
	matrix<double> rhs0 = m_pt->Homomap(movingPt);

	double rhs1;
	if (abstemp < m_precision)
	{
		rhs1 = movingPt(n, 0) - predictPt(n, 0);
	}
	else
	{
		rhs1 = m_pt->AdditionalEqn(movingPt, previousPt, Tanvector, predictPt, delta, m_sigma);
	}
	matrix<double> rhstemp(n + 1, 1);//Target2
	matrix_range<matrix<double> > rhsr01(rhstemp, range(0, n), range(0, 1));
	rhsr01 = rhs0;
	rhstemp(n, 0) = rhs1;
	/*生成右端矩阵end*/

	jacobi.assign_temporary(jacobitemp);
	rhs.assign_temporary(rhstemp);
	return true;
}


/// <summary>
/// 牛顿迭代.
/// </summary>
/// <param name="previousPoint">上一步的求解点，满足同伦方程.</param>
/// <param name="predictPoint">预估点.</param>
/// <param name="stepsize">步长.</param>
/// <param name="k">迭代次数.</param>
/// <param name="ConvergentPoint">迭代收敛点.</param>
/// <returns>bool.</returns>
bool CLiYorke::NewtonIter(
	const matrix<double>& previousPoint,
	const matrix<double>& predictPoint,
	double& stepsize,
	int& k,
	matrix<double>& ConvergentPoint)
{
	k = 0;
	matrix<double> Mjacobi;
	matrix<double> Mrhs;
	matrix<double> previousPt(previousPoint);
	matrix<double> predictPt(predictPoint);
	if (NewtonItergen(predictPt, previousPt, predictPt, stepsize, Mjacobi, Mrhs) == false)
	{
		return false;
	}
	matrix<double> unKnown = trans(Mrhs);
	//if (linalg::Solve(Mjacobi, unKnown) == 0)
	//{
	//	return false;
	//}

    int m_Free_degree = Mjacobi.size1();

    permutation_matrix<double> P(m_Free_degree);
    boost::numeric::ublas::vector<double> x(m_Free_degree);
    boost::numeric::ublas::vector<double> v(m_Free_degree);

    for (int j = 0; j < m_Free_degree; ++j)
    {
        v(j) = unKnown(0, j);
    }

    //try
    {
        lu_factorize(Mjacobi, P);
        x = v;
        lu_substitute(Mjacobi, P, x);
        cout << "解向量: " << x << endl;
        for (int j = 0; j < m_Free_degree; ++j)
        {
            unKnown(0, j) = x(j);
        }
        //unKnown = trans(unKnown);
    }


	matrix<double> interResult1 = predictPt - trans(unKnown);
	matrix<double> interResult2 = predictPt;

	std::vector<double> sigma;
	sigma.clear();
	sigma.push_back(norm_frobenius(interResult1 - interResult2) / norm_frobenius(interResult2));

	/*计算好第一步*/
	while (sigma[k] > m_epsilon)
	{
		k++;
		//cout << k << endl;
		//cout << norm_frobenius(interResult1 - interResult2) << endl;
		if (k > m_optIter) //5可以定义为成员变量，最优迭代次数
		{
			cout << "超过最大迭代次数,步长减半！" << endl;
			return false;
		}
		interResult2 = interResult1;
		if (NewtonItergen(interResult1, previousPt, predictPt, stepsize, Mjacobi, Mrhs) == false)
		{
			return false;
		}
		unKnown = trans(Mrhs);
		//if (linalg::Solve(Mjacobi, unKnown) == 0)
		//{
		//	return false;
		//}

        int m_Free_degree = Mjacobi.size1();

        permutation_matrix<double> P(m_Free_degree);
        boost::numeric::ublas::vector<double> x(m_Free_degree);
        boost::numeric::ublas::vector<double> v(m_Free_degree);

        for (int j = 0; j < m_Free_degree; ++j)
        {
            v(j) = unKnown(0, j);
        }

        //try
    {
        lu_factorize(Mjacobi, P);
        x = v;
        lu_substitute(Mjacobi, P, x);
        cout << "解向量: " << x << endl;
        for (int j = 0; j < m_Free_degree; ++j)
        {
            unKnown(0, j) = x(j);
        }
        //unKnown = trans(unKnown);
    }

		interResult1 = interResult2 - trans(unKnown);
		sigma.push_back(norm_frobenius(interResult1 - interResult2) / norm_frobenius(interResult2));

		if (sigma[k] > sigma[k - 1] / 3)
		{
			//cout << "牛顿迭代提前判断失败！" << endl;
			return false;
		}
		if (sigma[k] * sigma[k] / (sigma[k - 1] - sigma[k]) < m_epsilon)
		{
			break;
		}
	}

	//cout << "牛顿迭代的收敛点:" << interResult1 << endl;


	assert(interResult1.size1() == m_initPoint.size1()
		&& interResult1.size2() == m_initPoint.size2());
	if (norm_frobenius(interResult1 - m_initPoint) > m_infRadius)
	{
		cout << "超过无限性半径！" << endl;
		return false;
	}

	if (interResult1((int)interResult1.size1() - 1, 0) < 0.0)
	{
		cout << "牛顿迭代超出t=0平面，请修改输入条件" << endl;
		return false;
	}
	if (interResult1((int)interResult1.size1() - 1, 0) > 1.0)
	{
		cout << "牛顿迭代超出t=1平面，请修改输入条件" << endl;
		return false;
	}

	ConvergentPoint = interResult1;
	return true;
}

/// <summary>
/// 计算两个切向量的夹角.
/// </summary>
/// <param name="vec1">切向量1.</param>
/// <param name="vec2">切向量2.</param>
/// <param name="ang">夹角.</param>
/// <returns>bool.</returns>
bool CLiYorke::Angle(const matrix<double>& vec1,
	const matrix<double>& vec2,
	double& ang)
{
	//保证输入两个非零列向量
	assert(vec1.size2() == 1 && vec2.size2() == 1);
	assert(vec1.size1() == vec2.size1());
	double norm1 = norm_frobenius(vec1);
	double norm2 = norm_frobenius(vec2);
	assert(norm1 > m_precision && norm2 > m_precision);

	//int size = vec1.size1();
	//double innerproduct = 0.0;
	//for (int i = 0; i < size;++i)
	//{
	//	innerproduct += vec1(i, 0)*vec2(i, 0);
	//}

	matrix<double> tempvec1(vec1);
	matrix<double> tempvec2(vec2);

	matrix_column<matrix<double> > vec1column(tempvec1, 0);
	matrix_column<matrix<double> > vec2column(tempvec2, 0);
	double innerproduct = inner_prod(vec1column, vec2column);

	double a = innerproduct / (norm1 * norm2);
	//cout << "tempvec1:" << tempvec1 << endl;
	//cout << "tempvec2:" << tempvec2 << endl;
	//cout << "a:" << a << endl;
	//cout << "(180 / PI) * acos(a):" << (180 / PI) * acos(a) << endl;

	if (abs(a) > 1.0)
	{
		cout << "内积达到数值精度"<< endl;
		return false;
	}
	ang = (180 / PI) * acos(a);
	return true;
}

/// <summary>
/// 单步求解.
/// </summary>
/// <param name="StartPoint">起始点.</param>
/// <param name="stepsize">步长.</param>
/// <param name="endPoint">终止点.</param>
/// <param name="stepangle">夹角.</param>
/// <returns>bool.</returns>
bool CLiYorke::StepSolve(const matrix<double>& StartPoint,
	double& stepsize,
	matrix<double>& endPoint,
	double& stepangle)
{
	if (abs(stepsize) < m_precision)
	{
		cout << "步长达到机器精度！" << endl;
		return false;
	}

	matrix<double> initPt(StartPoint);
	matrix<double> tanvec1;
	if (TangentVec(initPt, tanvec1) == false)
	{
		return false;
	}
	//cout << "出发点:" << initPt << endl;
	matrix<double> pdictpt;
	matrix<double> tanvec2;
	//cout << "切向量:" << tanvec1 << endl;


	//cout << "步长:" << stepsize << endl;
	if (EulerPredict(initPt, tanvec1, stepsize, pdictpt) == false)
	{
		return false;
	}
	//cout << "预测点:" << pdictpt << endl;
	/*如果牛顿迭代不收敛，则缩小步长转回预估步*/

	int num = 0;
	int k;

	while (NewtonIter(initPt, pdictpt, stepsize, k, endPoint) == false)
		//如果小于最优迭代步就收敛，需要加大步长。如何实现？
	{
		//t从0.9(到达t==1之前)经过牛顿迭代到1之外，
		num++;
		if (num > 20)
		{
			return false;
		}
		stepsize /= 2.0;
		if (EulerPredict(initPt, tanvec1, stepsize, pdictpt) == false)
		{
			return false;
		}
		//cout << "重新的预测点:" << pdictpt << endl;
	}
	//k是第一次或者反复迭代中最后一次成功迭代中的次数。
	double p = (double)m_optIter / (k + 1);
	//cout << "小于最优迭代次数的倍数:" << p << endl;

	//if (p > 2.0) 
	//{
	//	stepsize *= 1.5;//sqrt(2.0);//如此判断导致假收敛
	//}

	//cout << "收敛点:" << endPoint << endl;

	if (TangentVec(endPoint, tanvec2) == false)
	{
		return false;
	}
	m_lastTanVec = tanvec2; //重置
	//cout << "tanvec2:" << tanvec2 << endl;
	if (Angle(tanvec1, tanvec2, stepangle) == false)
	{
		return false;
	}
	//cout << "stepangle:" << stepangle << endl;
	//cout << "stepsize:" << stepsize << endl;
	return true;
}

