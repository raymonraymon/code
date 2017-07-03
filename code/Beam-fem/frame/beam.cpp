/************************************************************************/
/*                 ���ƽ��ռ�Ԫ������Ԫ����                           */
/* author:����ɽ   date:2014-09-10 ver:1.0                              */
/* �ö�λ����(�����Ƶ����ɶȴ�1��ʼ�������Ƶ����ɶȴ�-1��ʼ)�����      */
/*                ��Ŀ�����ˡ�MATLAB����Ԫ������Ӧ�á�һ���е�8.1       */
/************************************************************************/

// beam.cpp : Defines the entry point for the console application.
//


#include "beam.h"
#include "PlaneStructer.h"

#define PI 3.14159265358979323846264338328



#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// The one and only application object



using namespace std;

int main()
{
	int nRetCode = 0;

		// TODO: code your application's behavior here.
		CPlanePoint P1(1,0.0,0.0);
		CPlanePoint	P2(2,0.0,3.0);
		CPlanePoint	P3(3,4.0,3.0);
		CPlanePoint	P4(4,4.0,0.0);
		std::map<int,CPlanePoint> point_map;
		point_map.insert(make_pair(1,P1));
		point_map.insert(make_pair(2,P2));
		point_map.insert(make_pair(3,P3));
		point_map.insert(make_pair(4,P4));
		/*�ڵ����ݿ����*/
		CPlaneElement E1(1,1,2,21e7,2e-2,5e-5);//��Ԫ��ţ����ҽڵ��ţ��������ϲ���
		CPlaneElement E2(2,2,3,21e7,2e-2,5e-5);
		CPlaneElement E3(3,3,4,21e7,2e-2,5e-5);
		std::map<int,CPlaneElement> element_map;
		element_map.insert(make_pair(1,E1));
		element_map.insert(make_pair(2,E2));
		element_map.insert(make_pair(3,E3));
		/*��Ԫ���ݿ����*/

		double f2[3] = {-20,0,0};
		double f3[3] = {0,0,12};

		CPlaneForce F2(2,f2);
		CPlaneForce F3(3,f3);
		std::vector<CPlaneForce> force_vector;
		force_vector.push_back(F2);
		force_vector.push_back(F3);

		/*�������*/

		double displace1[3] = {0,0,0};
		bool displace_bool1[3] = {true,true,true};
		CPlaneDisplacement Dis1(1,displace1,displace_bool1);

		double displace3[3] = {4,5,6};
		bool displace_bool3[3] = {true,true,true};
		CPlaneDisplacement Dis3(3,displace3,displace_bool3);

		double displace4[3] = {0,0,0};
		bool displace_bool4[3] = {true,true,true};
		CPlaneDisplacement Dis4(4,displace4,displace_bool4);


		std::map<int,CPlaneDisplacement> displacement_map;
		displacement_map.insert(make_pair(1,Dis1));
		//displacement_map.insert(make_pair(3,Dis3));
		displacement_map.insert(make_pair(4,Dis4));
		/*�߽��������*/

		CPlaneStructer str1(point_map,element_map,force_vector,displacement_map);

		if(str1.solve())
		{
			str1.show_result();
		}
		else
		{
			cout<<"It does not work!"<<endl;
		}


	return nRetCode;
}
