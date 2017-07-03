#include "stdafx.h"
#include "LiYorke.h"

#define PI 3.141592653589793

using namespace std;
using namespace boost::numeric::ublas;

/// <summary>
/// ���� <see cref="CLiYorke"/> ��.
/// </summary>
CLiYorke::CLiYorke()
{
}


/// <summary>
/// ���� <see cref="CLiYorke"/> ��.
/// </summary>
CLiYorke::~CLiYorke()
{
}

/// <summary>
/// LI-YORKE�㷨���ͬ������.
/// </summary>
/// <param name="pt">ͬ�׺�����ָ��.</param>
/// <param name="initPoint">���.</param>
/// <param name="epsilon">�����.</param>
/// <param name="initDelta">��ʼ����.</param>
/// <param name="infRadius">�����԰뾶.</param>
/// <param name="stepsizemap">���·���ϵĲ���.</param>
/// <param name="anglemap">���·���ϵļн�.</param>
/// <param name="pointsmap">���·���ϵĵ�.</param>
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
	/*                           ��ʼ����������ֵ��0��1��                     */
	/************************************************************************/
	int dim = (int)initPt.size1();
	matrix<double> lastTanVec(dim, 1);
	lastTanVec.clear();
	for (int i = 0; i < dim - 1; ++i)
	{
		lastTanVec(i, 0) = 0.0;
	}
	lastTanVec(dim - 1, 0) = 1.0; //�������Ƽ��á�0��1�������ʼ������������Ϊֻ��Ҫ����һ�Σ���deNiet�ķ������Ա�֤��ȷ��

	if (FirstTangentVec(initPt, lastTanVec) == false)
	{
		return false;
	}
	//cout << "deNiet��һ��������:" << lastTanVec << endl;
	
	m_lastTanVec = lastTanVec;

	double maxAngle = 18.0;
	double minAngle = 6.0;

	//matrix<double> Matjacobi0 = m_pt->HomomapJacobi(initPt);

	int n = (int)m_initPoint.size1() - 1;

	matrix<double> endPt(initPt); //������ֵ
	double ang;

	double stepsize = m_initStepsize;//��ʼ������
	int k = 0;

	pointsmap.insert(make_pair(k, endPt));

	ofstream outfile1("LiYorke_temp.csv", ios::out);
	if (!outfile1)
	{
		cerr << "open error" << endl;
		//exit(1);
		return false;
	}
	outfile1 << "temp��·�ϵĵ㣺" << endl;

	while (abs(endPt(n, 0) - 1.0) > m_precision || 
		norm_frobenius(m_pt->Homomap(endPt)) > m_precision)
	{
		if (k > 1050)
		{
			//cout << "Ԥ��У������������������" << endl;
			return false;
		}
		if (StepSolve(initPt, stepsize, endPt, ang) == false)
		{
			return false;
		}

		if (ang >= 90.0) //����90��,��������������
		{
			cout << "�Ƕȴ���90�ȣ�" << endl;
			return false;
		}

		if (ang >= maxAngle) //����18�ȣ��ϲ�������������map,�ۻ�Ԥ��
		{
			stepsize /= 2.0;
		}
		else  //�ϲ�����������map��
		{
			++k;
			initPt = endPt;
			stepsizemap.insert(make_pair(k, stepsize));
			anglemap.insert(make_pair(k, ang));
			pointsmap.insert(make_pair(k, endPt));

			outfile1 << setw(2) << k << ",";
			for (int i = 0; i < n + 1; ++i)
			{
			outfile1 << setw(20) << endPt(i, 0) << ",";
			}
			outfile1 << endl;

			if (ang < minAngle)//����6-18�ȣ������������ǰ��
			{
				stepsize *= 2.0;
			}
		}
		//cout << "endPt[t]:" << endPt(n, 0) << endl;
	}
	//cout << "��������:" << k << endl;

	double error = norm_frobenius(m_pt->Homomap(initPt));
	//cout << "���:" << error << endl;
	assert(error < m_precision);

	ofstream outfile("LiYorke.csv", ios::out);
	if (!outfile)
	{
		cerr << "open error" << endl;
		//exit(1);
		return false;
	}
	outfile << setprecision(ios::scientific) << setprecision(20);
	outfile << "��·�ϵĵ㣺" << endl;
	for (std::map<int, matrix<double>>::iterator pos = pointsmap.begin();
		pos != pointsmap.end();
		++pos)
	{
		outfile << setw(2) << pos->first << ",";
		for (int i = 0; i < n + 1; ++i)
		{
			outfile << setw(20) << pos->second(i, 0) << ",";
		}
		outfile << endl;
	}
	return true;
}

/// <summary>
/// �ж������ź���.
/// </summary>
/// <param name="val">����ֵ.</param>
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
/// �ݹ��������ʽ.
/// </summary>
/// <param name="Mat">�������.</param>
/// <returns>����ʽֵ.</returns>
double	CLiYorke::TanDet(const matrix<double>& Mat)
{

	assert(Mat.size1() == Mat.size2());
	int m = Mat.size1();
	int len = m - 1;/*����ʽ�Ľ�*/
	double s = 0;

	/*���ն��壬��ʼ��һ��������ʽ����Ŀռ�*/
	ublas::matrix<double> p(Mat.size1() - 1, Mat.size1() - 1);

	/*��Ϊ1�����ն������*/
	if (1 == m)
		return Tansign(Mat(0, 0));

	for (int k = 0; k < m; ++k)
	{
		for (int i = 0; i < len; ++i)
			for (int j = 0; j < len; ++j)
			{
				if (i < k)
					p(i, j) = Mat(i, j + 1);/*��ʼ��������ʽ��ֵ*/
				if (i >= k)
					p(i, j) = Mat(i + 1, j + 1);
			}
		s += (double)pow(-1.0, k) * Mat(k, 0) * TanDet(p);/*�ݹ����*/
	}
	//cout << "����ʽ��" << s << endl;
	return s;
}

/// <summary>
/// ��������intֵ.
/// </summary>
/// <param name="a">������������1.</param>
/// <param name="b">������������2.</param>
void CLiYorke::Tanswap(int& a, int& b)
{
	int temp;
	temp = a;
	a = b;
	b = temp;
}
/// <summary>
/// ���㽻��˳��Ĵ���.
/// </summary>
/// <param name="noSortedVector">��������ǰ������(��Ӧ�к�)\f$ sgn(a_{i,k_{i}}) \f$.</param>
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
/// ����2�����������е��к�.
/// </summary>
/// <param name="matrixAfterTansform">�������.</param>
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
			k = i; //��������ͬΪ����������غ�һ����
		}
	}
	return k;
}
/// <summary>
/// householder�任.
/// </summary>
/// <param name="A">�������.</param>
/// <param name="num">��num�α任.</param>
/// <param name="k">A��������ȥ2�����������е��к�.</param>
/// <param name="Q">�������.</param>
/// <returns>bool.</returns>
bool CLiYorke::TanHouseholder(matrix<double>& A,
	const int& num,//��num�α任
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
/// QR�ֽ⾫ȷ���������,���󶦸�����.
/// </summary>
/// <param name="point">���������.</param>
/// <param name="TangentVector">������������.</param>
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
	std::vector<int> maxCol;//��¼��Ԫ����
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
			CC(i, j) = 0; //��Ҫ����һ��ѭ��
		}
	}

	maxCol.push_back(TanMaxColumn(CC));

	//ʣ�µ����е������ô�õ���
	int lastColumn;
	for (int i = 0; i < m; ++i)
	{
		std::vector<int>::iterator result = find(maxCol.begin(), maxCol.end(), i); //����3
		if (result == maxCol.end())//û�ҵ�
		{
			lastColumn = i;
			break;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	//���������õ���ʼ������
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
	//theRightHand(n, 0) = ?;//�������ź����ĳ˻�,
	int k = 1;
	for (int i = 0; i < n; ++i)
	{
		//cout << "QRF("<<i<<","<< maxCol[i]<<"):" << QRF(i, maxCol[i]) << endl;

		k *= Tansign(QRF(i, maxCol[i]));	//�ҵ�ÿ����������Ǹ�Ԫ
	}

	//matrix<double> Matjacobi0 = m_pt->HomomapJacobi(pt);
	//m_sgndexHx = Tansign(TanDet(Matjacobi0));

	int powtimes = TanBubbleCount(maxCol);
	double sgntime = pow(-1.0, powtimes + 1);

	theRightHand(n, 0) = 1;// sgntime * k *  m_sgndexHx;
	/*�������ǰһ���������ļн���ѡ����һ���������ķ�����һ���Ϳ���ʡ��*/

	matrix<double> unKnown = trans(theRightHand);
	if (linalg::Solve(extendedMatrix, unKnown) == 0)
	{
		return false;
	}
	/************************************************************************/
	/* ������һ��������                                                       */
	/************************************************************************/

	unKnown = trans(unKnown);

	matrix_column<matrix<double> > vec1column(unKnown, 0);
	matrix_column<matrix<double> > vec2column(m_lastTanVec, 0);
	double anglecos = inner_prod(vec1column, vec2column);

	//////////////////////////////////////////////////////////////////////////
	if (anglecos < 0.0) //��֤�������е�t��ԶΪ�������󶦵��㷨ʵ��������Ҫ��
	{
		unKnown *= -1.0;
	}
	//////////////////////////////////////////////////////////////////////////

	//�����⻡��������m_sigma
	matrix<double> unKnownSigma;
	unKnownSigma = sqrt(m_sigma)*unKnown;
	unKnownSigma(lastColumn, 0) = unKnown(lastColumn, 0);
	double normUnknown = norm_frobenius(unKnownSigma);
	if (normUnknown < m_precision)
	{
		return false;
	}
	TangentVector = unKnown / normUnknown;
	//cout << "��������" << TangentVector << endl;
	return true;

}

/// <summary>
/// ������������������Ǿ�ȷ��������������Ȼ������еķ��������򵥣�����Ҳ��.\n
/// ���������ٵ�����.
/// </summary>
/// <param name="point">���������.</param>
/// <param name="TangentVector">������������.</param>
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
	if (linalg::Solve(MatJacobiT, unKnown) == 0)
	{
		return false;
	}

	unKnown = trans(unKnown);
	matrix<double> unKnownSigma;
	unKnownSigma = sqrt(m_sigma)*unKnown;
	unKnownSigma(n, 0) = unKnown(n, 0);
	double normUnknown = norm_frobenius(unKnownSigma);
	if (normUnknown < m_precision)
	{
		return false;
	}
	TangentVector = unKnown / normUnknown;

	//cout << "��������" << TangentVector << endl;
	return true;
}
/// <summary>
/// deNiet��˶ʿ�������ᵽ�ĵ�һ���������ļ��㷽��.
/// </summary>
/// <param name="point">���������.</param>
/// <param name="TangentVector">��һ��������.</param>
/// <returns>bool.</returns>
bool CLiYorke::FirstTangentVec(const matrix<double>& point,
	matrix<double>& TangentVector)
{
	matrix<double> pt(point);
	matrix<double> Matjacobi = m_pt->HomomapJacobi(pt);
	matrix<double> MatHt = m_pt->HomomapT(pt);
	matrix<double> unKnown = trans(-1 * MatHt);
	//cout << "Matjacobi��" << Matjacobi << endl;
	if (linalg::Solve(Matjacobi, unKnown) == 0)
	{
		return false;
	}
	unKnown = trans(unKnown);
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
/// ŷ����Ԥ��.
/// </summary>
/// <param name="currentPoint">������.</param>
/// <param name="tangentVector">�������������.</param>
/// <param name="delta">����.</param>
/// <param name="PredictRes">Ԥ����.</param>
/// <returns>bool.</returns>
bool CLiYorke::EulerPredict(const matrix<double>& curPoint,
	const matrix<double>& tanVec,
	double& delta, //���һ��Ҫ���ؼ�������Ĳ��������Բ���const
	matrix<double>& PredictRes)
{
	assert(curPoint.size1() == tanVec.size1() && curPoint.size2() == tanVec.size2());
	//assert(delta > m_precision);

	PredictRes = curPoint + delta * tanVec;

	//Ӧ�����һ��
	int n = (int)curPoint.size1();

	//cout <<"�м�tֵ��"<< PredictRes(n - 1, 0) << endl;

	//if (PredictRes(n - 1, 0) < 0.0)
	//{
	//	return false;
	//}

	if (PredictRes(n - 1, 0) > 1.0)
	{
		//cout << "����t=1ƽ�棡" << endl;
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

	//cout << "��ǰ�㣺" << curPoint << endl;
	//cout << "��������" << tanVec << endl;
	//cout << "������" << delta << endl;
	//cout << "ŷ��Ԥ��㣺" << PredictRes << endl;

	return true;
}

/// <summary>
/// ţ�ٵ���У���ľ�������.�ο���������4.28ʽ
/// </summary>
/// <param name="movingPoint">.</param>
/// <param name="previousPoint">��һ�������㣬����ͬ�׷���.</param>
/// <param name="predictPoint">Ԥ����.</param>
/// <param name="stepsize">����.</param>
/// <param name="jacobi">�ſɱȾ���.</param>
/// <param name="rhs">�Ҷ�����.</param>
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
	//����У�������ά��

	matrix<double> previousPt(previousPoint);
	matrix<double> predictPt(predictPoint);
	matrix<double> movingPt(movingPoint);

	/*����jabobi����*/
	matrix<double> jacobi00 = m_pt->HomomapJacobi(movingPt);
	matrix<double> jacobi01 = m_pt->HomomapT(movingPt);

	int n = (int)jacobi00.size1();

	if (predictPt(n, 0) < 0.0 || predictPt(n, 0) > 1.0)
	{
		cout << "Ԥ����tֵ���������С�ڵ���1��"<< endl;
		return false;
	}

	matrix<double> Tanvector;
	if (TangentVec(previousPt, Tanvector) == false)
	{
		return false;
	}
	matrix<double> jacobi1(1, n + 1);

	double delta = stepsize;//�����������õ��Ĳ���.

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
	/*����jabobi����end*/

	/*�����Ҷ˾���*/
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
	/*�����Ҷ˾���end*/

	jacobi.assign_temporary(jacobitemp);
	rhs.assign_temporary(rhstemp);
	return true;
}


/// <summary>
/// ţ�ٵ���.
/// </summary>
/// <param name="previousPoint">��һ�������㣬����ͬ�׷���.</param>
/// <param name="predictPoint">Ԥ����.</param>
/// <param name="stepsize">����.</param>
/// <param name="k">��������.</param>
/// <param name="ConvergentPoint">����������.</param>
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
	if (linalg::Solve(Mjacobi, unKnown) == 0)
	{
		return false;
	}
	matrix<double> interResult1 = predictPt - trans(unKnown);
	matrix<double> interResult2 = predictPt;

	std::vector<double> sigma;
	sigma.clear();
	sigma.push_back(norm_frobenius(interResult1 - interResult2) / norm_frobenius(interResult2));

	/*����õ�һ��*/
	while (sigma[k] > m_epsilon)
	{
		k++;
		//cout << k << endl;
		//cout << norm_frobenius(interResult1 - interResult2) << endl;
		if (k > m_optIter) //5���Զ���Ϊ��Ա���������ŵ�������
		{
			cout << "��������������,�������룡" << endl;
			return false;
		}
		interResult2 = interResult1;
		if (NewtonItergen(interResult1, previousPt, predictPt, stepsize, Mjacobi, Mrhs) == false)
		{
			return false;
		}
		unKnown = trans(Mrhs);
		if (linalg::Solve(Mjacobi, unKnown) == 0)
		{
			return false;
		}
		interResult1 = interResult2 - trans(unKnown);
		sigma.push_back(norm_frobenius(interResult1 - interResult2) / norm_frobenius(interResult2));

		if (sigma[k] > sigma[k - 1] / 3)
		{
			//cout << "ţ�ٵ�����ǰ�ж�ʧ�ܣ�" << endl;
			return false;
		}
		if (sigma[k] * sigma[k] / (sigma[k - 1] - sigma[k]) < m_epsilon)
		{
			break;
		}
	}

	//cout << "ţ�ٵ�����������:" << interResult1 << endl;


	assert(interResult1.size1() == m_initPoint.size1()
		&& interResult1.size2() == m_initPoint.size2());
	if (norm_frobenius(interResult1 - m_initPoint) > m_infRadius)
	{
		cout << "���������԰뾶��" << endl;
		return false;
	}

	if (interResult1((int)interResult1.size1() - 1, 0) < 0.0)
	{
		cout << "ţ�ٵ�������t=0ƽ�棬���޸���������" << endl;
		return false;
	}
	if (interResult1((int)interResult1.size1() - 1, 0) > 1.0)
	{
		cout << "ţ�ٵ�������t=1ƽ�棬���޸���������" << endl;
		return false;
	}

	ConvergentPoint = interResult1;
	return true;
}

/// <summary>
/// ���������������ļн�.
/// </summary>
/// <param name="vec1">������1.</param>
/// <param name="vec2">������2.</param>
/// <param name="ang">�н�.</param>
/// <returns>bool.</returns>
bool CLiYorke::Angle(const matrix<double>& vec1,
	const matrix<double>& vec2,
	double& ang)
{
	//��֤������������������
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
		cout << "�ڻ��ﵽ��ֵ����"<< endl;
		return false;
	}
	ang = (180 / PI) * acos(a);
	return true;
}

/// <summary>
/// �������.
/// </summary>
/// <param name="StartPoint">��ʼ��.</param>
/// <param name="stepsize">����.</param>
/// <param name="endPoint">��ֹ��.</param>
/// <param name="stepangle">�н�.</param>
/// <returns>bool.</returns>
bool CLiYorke::StepSolve(const matrix<double>& StartPoint,
	double& stepsize,
	matrix<double>& endPoint,
	double& stepangle)
{
	if (abs(stepsize) < m_precision)
	{
		cout << "�����ﵽ�������ȣ�" << endl;
		return false;
	}

	matrix<double> initPt(StartPoint);
	matrix<double> tanvec1;
	if (TangentVec(initPt, tanvec1) == false)
	{
		return false;
	}
	//cout << "������:" << initPt << endl;
	matrix<double> pdictpt;
	matrix<double> tanvec2;
	//cout << "������:" << tanvec1 << endl;


	//cout << "����:" << stepsize << endl;
	if (EulerPredict(initPt, tanvec1, stepsize, pdictpt) == false)
	{
		return false;
	}
	//cout << "Ԥ���:" << pdictpt << endl;
	/*���ţ�ٵ���������������С����ת��Ԥ����*/

	int num = 0;
	int k;

	while (NewtonIter(initPt, pdictpt, stepsize, k, endPoint) == false)
		//���С�����ŵ���������������Ҫ�Ӵ󲽳������ʵ�֣�
	{
		//t��0.9(����t==1֮ǰ)����ţ�ٵ�����1֮�⣬
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
		//cout << "���µ�Ԥ���:" << pdictpt << endl;
	}
	//k�ǵ�һ�λ��߷������������һ�γɹ������еĴ�����
	double p = (double)m_optIter / (k + 1);
	//cout << "С�����ŵ��������ı���:" << p << endl;

	//if (p > 2.0) 
	//{
	//	stepsize *= 1.5;//sqrt(2.0);//����жϵ��¼�����
	//}

	//cout << "������:" << endPoint << endl;

	if (TangentVec(endPoint, tanvec2) == false)
	{
		return false;
	}
	m_lastTanVec = tanvec2; //����
	//cout << "tanvec2:" << tanvec2 << endl;
	if (Angle(tanvec1, tanvec2, stepangle) == false)
	{
		return false;
	}
	//cout << "stepangle:" << stepangle << endl;
	//cout << "stepsize:" << stepsize << endl;
	return true;
}

