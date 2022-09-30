#include "stdafx.h"
#include "Myheader.h"
#include "Matrix.h"
#include "Lambda.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <assert.h>

using namespace std;

//声明ReadSPP_KIN函数，读取spp文件（针对本人自定义的spp文件格式，可修改）
bool ReadSPP_KIN(string filename, observe_spp& sppfile);

//声明OnCycleDetect函数，周跳探测
void OnCycleDetect(Obs_epoch& tempEpoch, Obs_epoch& preEpoch);

//声明EpochCommonPro函数，对基站和流动站卫星观测数据进行共视处理
void EpochCommonPro(Obs_epoch tempData1, Obs_epoch tempData2, Obs_epoch& useData);

//声明SelectRef函数，选择参考卫星
void SelectRef(Obs_epoch tempData, int& RefNo);

//声明InitState函数，初始化状态方程
void InitState(CMatrix& QX);

//声明GetGRCQValue函数，生成转换矩阵
bool GetGRCQValue(int PRN_X[], int PRN_old[], int num_sat, int num_sat_old, int Ref, int Ref_old, CMatrix& TT_Matrix);

//声明GetDDR_Peo函数，依据误差传播定律获取噪声矩阵（伪距）
void GetDDR_Peo(Obs_epoch GPSData, vector<int> useprn, int GPSRef, CMatrix& R);

//声明GetDDR_Leo函数，依据误差传播定律获取噪声矩阵（载波）
void GetDDR_Leo(Obs_epoch GPSData, vector<int> useprn, int GPSRef, CMatrix& R);

//声明Error_Correction函数，误差改正
double Error_Correction(Sat& tempSat1, Sat& RefSat1, Sat& tempSat2, Sat& RefSat2);

//重载ProGRCTrans函数，拼接矩阵
void ProGRCTrans(CMatrix& TT_G, CMatrix& TT_C, CMatrix& TT_B, CMatrix& TT_E, CMatrix& TT_ALL)
{
	int Row_G = TT_G.Row;
	int Col_G = TT_G.Col;
	int Row_C = TT_C.Row;
	int Col_C = TT_C.Col;
	int Row_B = TT_B.Row;
	int Col_B = TT_B.Col;
	int Row_E = TT_E.Row;
	int Col_E = TT_E.Col;

	int Tot_Row = Row_G + Row_C + Row_B + Row_E;
	int Tot_Col = Col_G + Col_C + Col_B + Col_E;
	TT_ALL.SetSize(Tot_Row, Tot_Col);

	int i = 0;
	int j = 0;

	//GPS部分赋值给新的转换矩阵
	for (i = 0; i < Row_G; i++)
		for (j = 0; j < Col_G; j++)
			TT_ALL[i][j] = TT_G[i][j];

	//BDS-2部分赋值给新的转换矩阵
	for (i = 0; i < Row_C; i++)
		for (j = 0; j < Col_C; j++)
			TT_ALL[Row_G + i][Col_G + j] = TT_C[i][j];

	//BDS-3部分赋值给新的转换矩阵
	for (i = 0; i < Row_B; i++)
		for (j = 0; j < Col_B; j++)
			TT_ALL[Row_G + Row_C + i][Col_G + Col_C + j] = TT_B[i][j];

	//Galileo部分赋值给新的转换矩阵
	for (i = 0; i < Row_E; i++)
		for (j = 0; j < Col_E; j++)
			TT_ALL[Row_G + Row_C + Row_B + i][Col_G + Col_C + Col_B + j] = TT_E[i][j];
}

//重载CheckNewSats函数，检测新卫星
void CheckNewSats(int* SatNum, int* oldSatNum, int SatPRN[4][40], int oldSatPRN[4][40], int NewSatPRN[4][40], int* NewSatNum)
{
	for (int i = 0; i < 4; i++)
	{
		int newsatsnum = 0;
		for (int j = 0; j < SatNum[i]; j++)
		{
			bool bFind = false;
			for (int k = 0; k < oldSatNum[i]; k++)
			{
				if (SatPRN[i][j] == oldSatPRN[i][k])
				{
					bFind = true;
					break;
				}
			}

			if (!bFind)
			{
				NewSatPRN[i][newsatsnum] = SatPRN[i][j];
				newsatsnum++;
			}
		}
		NewSatNum[i] = newsatsnum;
	}
}

//重载UpdataState函数，更新协方差
void UpdataState(CMatrix& QX, CMatrix& X, int* SatNum, int* oldSatNum, int SatPRN[4][40], int oldSatPRN[4][40], int* RefSat, int* OldRef)
{
	CMatrix T[4], ALL_T, Matrix_t;
	int sys = 0, oldsys = 0, ALL_NUM, oldALL_NUM;

	if (SatNum[0] >= 2)
	{
		sys++;
	}
	if (SatNum[1] >= 2)
	{
		sys++;
	}
	if (SatNum[2] >= 2)
	{
		sys++;
	}
	if (SatNum[3] >= 2)
	{
		sys++;
	}

	if (oldSatNum[0] >= 2)
	{
		oldsys++;
	}
	if (oldSatNum[1] >= 2)
	{
		oldsys++;
	}
	if (oldSatNum[2] >= 2)
	{
		oldsys++;
	}
	if (oldSatNum[3] >= 2)
	{
		oldsys++;
	}

	ALL_NUM = SatNum[0] + SatNum[1] + SatNum[2] + SatNum[3] - sys;
	oldALL_NUM = oldSatNum[0] + oldSatNum[1] + oldSatNum[2] + oldSatNum[3] - oldsys;

	//L1或者L1变化矩阵
	for (int i = 0; i < 4; i++)
	{
		if (SatNum[i] > 1)
		{
			if (oldSatNum[i] > 0)
			{
				T[i].SetSize(SatNum[i] - 1, oldSatNum[i] - 1);//当前历元GPS卫星数及上一历元GPS卫星数
				GetGRCQValue(SatPRN[i], oldSatPRN[i], SatNum[i], oldSatNum[i], RefSat[i], OldRef[i], T[i]);//生成转换矩阵
			}
			else
			{
				T[i].SetSize(SatNum[i] - 1, 0);
			}
		}
		else
		{
			T[i].SetSize(0, 0);
			if (oldSatNum[i] != 0)
			{
				T[i].SetSize(0, oldSatNum[i] - 1);
			}
		}
	}

	ALL_T.SetSize(ALL_NUM, oldALL_NUM);
	ProGRCTrans(T[0], T[1], T[2], T[3], ALL_T);//拼接矩阵

	Matrix_t.SetSize(ALL_NUM + 3, oldALL_NUM + 3);
	for (int j = 0; j < 3; j++)
	{
		Matrix_t[j][j] = 1.0;
	}

	for (int ii = 0; ii < ALL_NUM; ii++)
	{
		for (int jj = 0; jj < oldALL_NUM; jj++)
		{
			Matrix_t[ii + 3][jj + 3] = ALL_T[ii][jj];
		}
	}

	QX = Matrix_t * QX * Matrix_t.T();
	X = Matrix_t * X;
}

//声明KalmanFilter，卡尔曼滤波
void KalmanFilter(CMatrix& B, CMatrix& L, CMatrix& R, CMatrix& Qw, CMatrix& Q, CMatrix& X);

//声明Resamb_LAMBDA函数，Lambda算法固定模糊度
void Resamb_LAMBDA(CMatrix& Q, CMatrix& X, CMatrix& aXb, double& ratio, CMatrix& fixX);

//函数功能：通过读取spp文件，进行多系统RTK解算
bool MSMF_RTK(string filename1, string filename2)
{
	//默认第一个站为基准站，第二个为流动站
	observe_spp sppfile1;
	observe_spp sppfile2;
	ReadSPP_KIN(filename1, sppfile1);
	ReadSPP_KIN(filename2, sppfile2);

	//投影矩阵：东北天坐标系（此处简化直接采用头文件坐标）
	CMatrix TT(3, 3);
	TT[0][0] = -sin(sppfile2.B) * cos(sppfile2.L);
	TT[0][1] = -sin(sppfile2.B) * sin(sppfile2.L);
	TT[0][2] = cos(sppfile2.B);

	TT[1][0] = -sin(sppfile2.L);
	TT[1][1] = cos(sppfile2.L);
	TT[1][2] = 0;

	TT[2][0] = cos(sppfile2.B) * cos(sppfile2.L);
	TT[2][1] = cos(sppfile2.B) * sin(sppfile2.L);
	TT[2][2] = sin(sppfile2.B);

	bool bInit = false;//初始化

	//输出差分定位结果用
	string name = filename2.substr(0, 3) + "_MSMF_RTK.txt";
	const char* output_filename = name.c_str();
	FILE* fp3;
	fopen_s(&fp3, output_filename, "w");

	int OldSatPRN[4][40] = { 0 };
	int OldRefPRN[4] = { -1,-1,-1,-1 };
	int OldSatNum[4] = { 0 };

	CMatrix Q, X;

	//根据流动站匹配基准站对应历元数据
	int posk = 0;   //用以记录基站数据的匹配
	//流动站逐历元解算
	for (int i = 0; i < sppfile2.liyuan_num; i++)
	{
		Obs_epoch tempData;

		//流动站当前时刻
		Obs_epoch tempData2 = sppfile2.epoch[i];

		if (i >= 1)
		{
			Obs_epoch preData2 = sppfile2.epoch[i - 1];
			OnCycleDetect(tempData2, preData2);     //判断此此历元与上一个历元中的卫星载波测量值是否有周跳,如果有周跳就将此历元中的卫星记为不可用
		}

		//基准站数据相应时刻必须在前
		int j = 0;
		for (j = posk; j < sppfile1.liyuan_num; j++)
		{
			if (sppfile1.epoch[j].GPSTIME <= tempData2.GPSTIME)
				continue;
			else
				break;
		}

		//回溯一个历元：posk即为要找的最近的历元且满足posk的时刻小于tempData2的时间
		if (j > 0) posk = --j;
		else
		{
			posk = j;
			continue;
		}

		//判断时间差：两个是否相隔超过2s
		if (fabs(sppfile1.epoch[j].GPSTIME - tempData2.GPSTIME) >= 2)
		{
			printf("未匹配到基站数据\n");
			continue;
		}

		//如果基准站数据滞后
		if (sppfile1.epoch[j].GPSTIME > tempData2.GPSTIME)
		{
			printf("基站数据匹配错误\n");
			exit(1);
		}

		//获得基站参与解算的历元数据
		Obs_epoch tempData1 = sppfile1.epoch[posk];

		//对基站进行周跳探测
		if (posk >= 1)
		{
			Obs_epoch preData1 = sppfile1.epoch[posk - 1];
			OnCycleDetect(tempData1, preData1);
		}

		//差分定位共视处理
		EpochCommonPro(tempData1, tempData2, tempData);

		//分离GPS、BDS-2、BDS-3、Galileo卫星，便于后续单系统、多系统处理
		Obs_epoch GPSdata;
		Obs_epoch BDS2data;
		Obs_epoch BDS3data;
		Obs_epoch Galileodata;

		int SatPRN[4][40] = { 0 };		//当前历元卫星号 依次存储
		int SatNum[4] = { 0 };	//当前历元各系统卫星数目
		int RefPRN[4] = { 0 };  //当前历元各系统参考卫星号

		int NewSatNum[4] = { 0 };
		int NewSatPRN[4][40] = { 0 };

		int k = 0;
		int m = 0;
		int gpsnum, bds2num, bds3num, galileonum;
		gpsnum = bds2num = bds3num = galileonum = 0;
		for (m = 0; m < tempData.sat_num; m++)
		{
			if (tempData.sat1[m].sattype == "G")
			{
				GPSdata.sat1.push_back(tempData.sat1[m]);
				GPSdata.sat2.push_back(tempData.sat2[m]);
				SatPRN[0][gpsnum] = tempData.sat1[m].numofsat;
				gpsnum++;
			}

			if (tempData.sat1[m].sattype == "C")
			{
				BDS2data.sat1.push_back(tempData.sat1[m]);
				BDS2data.sat2.push_back(tempData.sat2[m]);
				SatPRN[1][bds2num] = tempData.sat1[m].numofsat;
				bds2num++;
			}

			if (tempData.sat1[m].sattype == "B")
			{
				BDS3data.sat1.push_back(tempData.sat1[m]);
				BDS3data.sat2.push_back(tempData.sat2[m]);
				SatPRN[2][bds3num] = tempData.sat1[m].numofsat;
				bds3num++;
			}

			if (tempData.sat1[m].sattype == "E")
			{
				Galileodata.sat1.push_back(tempData.sat1[m]);
				Galileodata.sat2.push_back(tempData.sat2[m]);
				SatPRN[3][galileonum] = tempData.sat1[m].numofsat;
				galileonum++;
			}
		}

		//各类型卫星数
		SatNum[0] = gpsnum;
		SatNum[1] = bds2num;
		SatNum[2] = bds3num;
		SatNum[3] = galileonum;

		GPSdata.sat_num = GPSdata.sat1.size();
		BDS2data.sat_num = BDS2data.sat1.size();
		BDS3data.sat_num = BDS3data.sat1.size();
		Galileodata.sat_num = Galileodata.sat1.size();
		tempData.GPSTIME = tempData2.GPSTIME;

		//卫星数至少要4颗
		if ((GPSdata.sat_num + BDS2data.sat_num + BDS3data.sat_num + Galileodata.sat_num) < 4)
		{
			printf("卫星数少于4颗");
			bInit = false;
			continue;
		}

		//各系统至少要2颗卫星,区分参考星和普通星
		if (GPSdata.sat_num < 2)
		{
			printf("GPS卫星数少于2颗\n");
			bInit = false;
			continue;
		}
		if (BDS2data.sat_num < 2)
		{
			printf("BDS2卫星数少于2颗\n");
			bInit = false;
			continue;
		}
		if (BDS3data.sat_num < 2)
		{
			printf("BDS3卫星数少于2颗\n");
			bInit = false;
			continue;
		}
		if (Galileodata.sat_num < 2)
		{
			printf("Galileo卫星数少于2颗\n");
			bInit = false;
			continue;
		}

		//采用站间星间二次差分模型，选择参考卫星
		int GPSref = -1;
		int BDS2ref = -2;
		int BDS3ref = -3;
		int Galref = -4;

		GPSref = OldRefPRN[0];
		BDS2ref = OldRefPRN[1];
		BDS3ref = OldRefPRN[2];
		Galref = OldRefPRN[3];

		//选择最大高度角卫星作为参考卫星，不同类型卫星分开看
		SelectRef(GPSdata, GPSref);
		SelectRef(BDS2data, BDS2ref);
		SelectRef(BDS3data, BDS3ref);
		SelectRef(Galileodata, Galref);

		RefPRN[0] = GPSref;  //GPS参考卫星序号
		RefPRN[1] = BDS2ref; //BDS-2参考卫星序号
		RefPRN[2] = BDS3ref; //BDS-3参考卫星序号
		RefPRN[3] = Galref;	 //Galileo参考卫星序号

		//当前历元GPS、BDS-2、BDS-3、Galileo卫星序列
		int PRN_X_G[32];
		int PRN_X_B2[20];
		int PRN_X_B3[40];
		int PRN_X_E[40];

		//数列初始化
		memset(PRN_X_G, 0, sizeof(PRN_X_G));
		memset(PRN_X_B2, 0, sizeof(PRN_X_B2));
		memset(PRN_X_B3, 0, sizeof(PRN_X_B3));
		memset(PRN_X_E, 0, sizeof(PRN_X_E));

		//储存历元中的卫星序列
		for (k = 0; k < GPSdata.sat_num; k++)
			PRN_X_G[k] = GPSdata.sat1[k].numofsat;
		for (k = 0; k < BDS2data.sat_num; k++)
			PRN_X_B2[k] = BDS2data.sat1[k].numofsat;
		for (k = 0; k < BDS3data.sat_num; k++)
			PRN_X_B3[k] = BDS3data.sat1[k].numofsat;
		for (k = 0; k < Galileodata.sat_num; k++)
			PRN_X_E[k] = Galileodata.sat1[k].numofsat;

		//各卫星个数
		int num_X_G = GPSdata.sat_num;
		int num_X_B2 = BDS2data.sat_num;
		int num_X_B3 = BDS3data.sat_num;
		int num_X_E = Galileodata.sat_num;

		//各系统双差方程个数
		int numgps_x = num_X_G - 1;
		int numbds2_x = num_X_B2 - 1;
		int numbds3_x = num_X_B3 - 1;
		int numgal_x = num_X_E - 1;

		//多系统融合双差方程个数 P1 + L1 载波观测方程加伪距观测方程
		int num_x_p = (numgps_x + numbds2_x + numbds3_x + numgal_x) * 2;

		//待估参数个数：ΔX，ΔY，ΔZ，每颗卫星的模糊度
		int num_x = numgps_x + numbds2_x + numbds3_x + numgal_x + 3;

		//模糊度个数
		int num_n = numgps_x + numbds2_x + numbds3_x + numgal_x;

		//定位解算相关的矩阵
		//1-多系统
		CMatrix B(num_x_p, num_x);                //系数矩阵  
		CMatrix L(num_x_p, 1);                    //观测值矩阵
		CMatrix R(num_x_p, num_x_p);              //误差矩阵（逆阵即为权矩阵）
		CMatrix Qw(num_x, num_x);                 //动态噪声阵

		//2-单GPS
		CMatrix RG_P(numgps_x, numgps_x);         //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RG_L(numgps_x, numgps_x);         //误差矩阵（逆阵即为权矩阵） 载波

		//3-单BDS2
		CMatrix RC2_P(numbds2_x, numbds2_x);      //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RC2_L(numbds2_x, numbds2_x);      //误差矩阵（逆阵即为权矩阵） 载波

		//4-单BDS3
		CMatrix RC3_P(numbds3_x, numbds3_x);      //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RC3_L(numbds3_x, numbds3_x);      //误差矩阵（逆阵即为权矩阵） 载波

		//5-单Galileo
		CMatrix RE_P(numgal_x, numgal_x);         //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RE_L(numgal_x, numgal_x);         //误差矩阵（逆阵即为权矩阵） 载波

		//位置的动态噪声阵，如果没有周跳，模糊度不变，所以模糊度的噪声为0，只有位置噪声
		for (int idx = 0; idx < 3; idx++) {
			Qw[idx][idx] = 30.0 * 30.0;
		}

		//获取GPS和BDS的误差矩阵
		if (!bInit)
		{
			Q.SetSize(num_x, num_x);
			X.SetSize(num_x, 1);
			InitState(Q);
		}
		else
		{
			//检测新卫星
			CheckNewSats(SatNum, OldSatNum, SatPRN, OldSatPRN, NewSatPRN, NewSatNum);

			//X.MyTRACE();

			UpdataState(Q, X, SatNum, OldSatNum, SatPRN, OldSatPRN, RefPRN, OldRefPRN);
		}

		//各类型卫星序号存储容器
		vector<int> GPS_PRN;
		vector<int> BDS2_PRN;
		vector<int> BDS3_PRN;
		vector<int> Gal_PRN;

		for (k = 0; k < GPSdata.sat_num; k++)
			GPS_PRN.push_back(PRN_X_G[k]);
		for (k = 0; k < BDS2data.sat_num; k++)
			BDS2_PRN.push_back(PRN_X_B2[k]);
		for (k = 0; k < BDS3data.sat_num; k++)
			BDS3_PRN.push_back(PRN_X_B3[k]);
		for (k = 0; k < Galileodata.sat_num; k++)
			Gal_PRN.push_back(PRN_X_E[k]);

		//伪距
		GetDDR_Peo(GPSdata, GPS_PRN, GPSref, RG_P);
		GetDDR_Peo(BDS2data, BDS2_PRN, BDS2ref, RC2_P);
		GetDDR_Peo(BDS3data, BDS3_PRN, BDS3ref, RC3_P);
		GetDDR_Peo(Galileodata, Gal_PRN, Galref, RE_P);

		//载波
		GetDDR_Leo(GPSdata, GPS_PRN, GPSref, RG_L);
		GetDDR_Leo(BDS2data, BDS2_PRN, BDS2ref, RC2_L);
		GetDDR_Leo(BDS3data, BDS3_PRN, BDS3ref, RC3_L);
		GetDDR_Leo(Galileodata, Gal_PRN, Galref, RE_L);

		//融合系统的R矩阵，分块对角矩阵
		//伪距
		for (int k1 = 0; k1 < numgps_x; k1++)
			for (int k2 = 0; k2 < numgps_x; k2++)
			{
				R[k1][k2] = RG_P[k1][k2];
			}

		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R[k1 + numgps_x][k2 + numgps_x] = RC2_P[k1][k2];
			}

		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R[k1 + numgps_x + numbds2_x][k2 + numgps_x + numbds2_x] = RC3_P[k1][k2];
			}

		for (int k1 = 0; k1 < numgal_x; k1++)
			for (int k2 = 0; k2 < numgal_x; k2++)
			{
				R[k1 + numgps_x + numbds2_x + numbds3_x][k2 + numgps_x + numbds2_x + numbds3_x] = RE_P[k1][k2];
			}
		//载波
		for (int k1 = 0; k1 < numgps_x; k1++)
			for (int k2 = 0; k2 < numgps_x; k2++)
			{
				R[k1 + numgps_x + numbds2_x + numbds3_x + numgal_x][k2 + numgps_x + numbds2_x + numbds3_x + numgal_x] = RG_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R[k1 + numgps_x * 2 + numbds2_x + numbds3_x + numgal_x][k2 + numgps_x * 2 + numbds2_x + numbds3_x + numgal_x] = RC2_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R[k1 + numgps_x * 2 + numbds2_x * 2 + numbds3_x + numgal_x][k2 + numgps_x * 2 + numbds2_x * 2 + numbds3_x + numgal_x] = RC3_L[k1][k2];
			}
		for (int k1 = 0; k1 < numgal_x; k1++)
			for (int k2 = 0; k2 < numgal_x; k2++)
			{
				R[k1 + numgps_x * 2 + numbds2_x * 2 + numbds3_x * 2 + numgal_x][k2 + numgps_x * 2 + numbds2_x * 2 + numbds3_x * 2 + numgal_x] = RE_L[k1][k2];
			}

		//初始坐标信息：基站采用spp文件头中坐标（固定），流动站采用每个历元的概率坐标（考虑实际应用中其可能为动态情况）
		//本组数据中：流动站为静态，因此其文件头中坐标也读取下来便于比较
		//基站坐标：头文件坐标
		double POSITION_X1 = sppfile1.APP_X;
		double POSITION_Y1 = sppfile1.APP_Y;
		double POSITION_Z1 = sppfile1.APP_Z;

		//流动站坐标：是每个历元单点定位解出的坐标
		double POSITION_X2 = tempData2.posX;
		double POSITION_Y2 = tempData2.posY;
		double POSITION_Z2 = tempData2.posZ;

		//流动站比较值（可认为真值，头文件坐标）
		double P_X2 = sppfile2.APP_X;
		double P_Y2 = sppfile2.APP_Y;
		double P_Z2 = sppfile2.APP_Z;

		//找到参考卫星下标
		//初始化
		int ref_gps = -1;
		int ref_bds2 = -1;
		int ref_bds3 = -1;
		int ref_gal = -1;

		for (k = 0; k < GPSdata.sat_num; k++)
		{
			if (GPSdata.sat1[k].numofsat == GPSref)
			{
				ref_gps = k;
				break;
			}
		}

		for (k = 0; k < BDS2data.sat_num; k++)
		{
			if (BDS2data.sat1[k].numofsat == BDS2ref)
			{
				ref_bds2 = k;
				break;
			}
		}

		for (k = 0; k < BDS3data.sat_num; k++)
		{
			if (BDS3data.sat1[k].numofsat == BDS3ref)
			{
				ref_bds3 = k;
				break;
			}
		}

		for (k = 0; k < Galileodata.sat_num; k++)
		{
			if (Galileodata.sat1[k].numofsat == Galref)
			{
				ref_gal = k;
				break;
			}
		}

		int tp = 0;  //记录单系统方程数
		int tw = 0;  //记录双系统方程数
		double DDerror = 0;  //初始化误差改正项
		//GPS
		for (k = 0; k < GPSdata.sat_num; k++)
		{
			bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = GPSdata.sat1[k].numofsat;
			//参考卫星跳过
			if (PRN == GPSref)
			{
				continue;
			}

			if (ref_gps < 0)
			{
				assert(0);
			}

			Sat tempSat1 = GPSdata.sat1[k];
			Sat RefSat1 = GPSdata.sat1[ref_gps];
			Sat tempSat2 = GPSdata.sat2[k];
			Sat RefSat2 = GPSdata.sat2[ref_gps];

			double length_O1, length_R1, length_O2, length_R2;

			//计算卫星到测站距离（4个距离）
			length_O1 = sqrt(pow((tempSat1.POS_X - POSITION_X1), 2) + pow((tempSat1.POS_Y - POSITION_Y1), 2) + pow((tempSat1.POS_Z - POSITION_Z1), 2));
			length_R1 = sqrt(pow((RefSat1.POS_X - POSITION_X1), 2) + pow((RefSat1.POS_Y - POSITION_Y1), 2) + pow((RefSat1.POS_Z - POSITION_Z1), 2));

			length_O2 = sqrt(pow((tempSat2.POS_X - POSITION_X2), 2) + pow((tempSat2.POS_Y - POSITION_Y2), 2) + pow((tempSat2.POS_Z - POSITION_Z2), 2));
			length_R2 = sqrt(pow((RefSat2.POS_X - POSITION_X2), 2) + pow((RefSat2.POS_Y - POSITION_Y2), 2) + pow((RefSat2.POS_Z - POSITION_Z2), 2));

			//三个坐标系数
			double Cof1 = (tempSat2.POS_X - POSITION_X2) / length_O2 - (RefSat2.POS_X - POSITION_X2) / length_R2;
			double Cof2 = (tempSat2.POS_Y - POSITION_Y2) / length_O2 - (RefSat2.POS_Y - POSITION_Y2) / length_R2;
			double Cof3 = (tempSat2.POS_Z - POSITION_Z2) / length_O2 - (RefSat2.POS_Z - POSITION_Z2) / length_R2;

			B[tw][0] = Cof1;
			B[tw][1] = Cof2;
			B[tw][2] = Cof3;         //伪距在前

			B[tw + num_n][0] = Cof1;
			B[tw + num_n][1] = Cof2;
			B[tw + num_n][2] = Cof3;  //载波在后
			B[tw + num_n][3 + tw] = lambda_L1;

			//误差改正
			DDerror = Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2);

			//双差形式的几何项和误差
			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;
			double deltp1 = (tempSat2.data[0] - RefSat2.data[0]) - (tempSat1.data[0] - RefSat1.data[0]);  //伪距
			double deltL1 = ((tempSat2.data[2] - RefSat2.data[2]) - (tempSat1.data[2] - RefSat1.data[2])) * lambda_L1; //载波

			L[tw][0] = deltp1 - DDgeo;         //伪距在前
			L[tw + num_n][0] = deltL1 - DDgeo; //载波在后

			for (int ii = 0; ii < SatNum[0]; ii++)
			{
				if (OldRefPRN[0] == SatPRN[0][ii])
				{
					Init = true;
					break;
				}
				else
					Init = false;
			}

			//若旧参考星被剔除，需初始化模糊度
			if (!Init)
			{
				X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1;
				Q[3 + tw][3 + tw] = 100.0 * 100.0;
			}

			if (!bInit)
			{
				X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1;       //初始化模糊度
			}

			for (int ii = 0; ii < NewSatNum[0]; ii++)
			{
				if (NewSatPRN[0][ii] == PRN)
				{
					X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1;   //初始化模糊度
					Q[3 + tw][3 + tw] = 100.0 * 100.0;
					break;
				}
			}
			tw++;
		}

		//BDS-2
		for (k = 0; k < BDS2data.sat_num; k++)
		{
			bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = BDS2data.sat1[k].numofsat;
			//参考卫星跳过
			if (PRN == BDS2ref)
			{
				continue;
			}

			if (ref_bds2 < 0)
			{
				assert(0);
			}

			Sat tempSat1 = BDS2data.sat1[k];
			Sat RefSat1 = BDS2data.sat1[ref_bds2];
			Sat tempSat2 = BDS2data.sat2[k];
			Sat RefSat2 = BDS2data.sat2[ref_bds2];

			double length_O1, length_R1, length_O2, length_R2;

			//计算卫星到测站距离（4个距离）
			length_O1 = sqrt(pow((tempSat1.POS_X - POSITION_X1), 2) + pow((tempSat1.POS_Y - POSITION_Y1), 2) + pow((tempSat1.POS_Z - POSITION_Z1), 2));
			length_R1 = sqrt(pow((RefSat1.POS_X - POSITION_X1), 2) + pow((RefSat1.POS_Y - POSITION_Y1), 2) + pow((RefSat1.POS_Z - POSITION_Z1), 2));

			length_O2 = sqrt(pow((tempSat2.POS_X - POSITION_X2), 2) + pow((tempSat2.POS_Y - POSITION_Y2), 2) + pow((tempSat2.POS_Z - POSITION_Z2), 2));
			length_R2 = sqrt(pow((RefSat2.POS_X - POSITION_X2), 2) + pow((RefSat2.POS_Y - POSITION_Y2), 2) + pow((RefSat2.POS_Z - POSITION_Z2), 2));

			//三个坐标系数
			double Cof1 = (tempSat2.POS_X - POSITION_X2) / length_O2 - (RefSat2.POS_X - POSITION_X2) / length_R2;
			double Cof2 = (tempSat2.POS_Y - POSITION_Y2) / length_O2 - (RefSat2.POS_Y - POSITION_Y2) / length_R2;
			double Cof3 = (tempSat2.POS_Z - POSITION_Z2) / length_O2 - (RefSat2.POS_Z - POSITION_Z2) / length_R2;

			//系数矩阵
			B[tw][0] = Cof1;
			B[tw][1] = Cof2;
			B[tw][2] = Cof3;
			B[tw + num_n][0] = Cof1;
			B[tw + num_n][1] = Cof2;
			B[tw + num_n][2] = Cof3;
			B[tw + num_n][3 + tw] = lambda_L1_C;

			//误差改正
			DDerror = Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2);

			//双差形式的几何项和误差
			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;
			double deltp1 = (tempSat2.data[0] - RefSat2.data[0]) - (tempSat1.data[0] - RefSat1.data[0]);  //伪距
			double deltL1 = ((tempSat2.data[2] - RefSat2.data[2]) - (tempSat1.data[2] - RefSat1.data[2])) * lambda_L1_C;  //载波
			L[tw][0] = deltp1 - DDgeo;         //伪距在前
			L[tw + num_n][0] = deltL1 - DDgeo; //载波在后

			for (int ii = 0; ii < SatNum[1]; ii++)
			{
				if (OldRefPRN[1] == SatPRN[1][ii])
				{
					Init = true;
					break;
				}
				else
					Init = false;
			}

			//若旧参考星被剔除，需初始化模糊度
			if (!Init)
			{
				X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1_C;
				Q[3 + tw][3 + tw] = 100.0 * 100.0;
			}

			if (!bInit)
			{
				X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1_C;       //初始化模糊度
			}

			for (int ii = 0; ii < NewSatNum[1]; ii++)
			{
				if (NewSatPRN[1][ii] == PRN)
				{
					X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1_C;   //初始化模糊度
					Q[3 + tw][3 + tw] = 100.0 * 100.0;
					break;
				}
			}
			tw++;
		}

		//BDS-3
		if (BDS3data.sat_num > 1)
		{
			for (k = 0; k < BDS3data.sat_num; k++)
			{
				bool Init = true;                      //记录是否有旧参考星被剔除
				int PRN = BDS3data.sat1[k].numofsat;
				//参考卫星跳过
				if (PRN == BDS3ref)
				{
					continue;
				}

				if (ref_bds3 < 0)
				{
					assert(0);
				}

				Sat tempSat1 = BDS3data.sat1[k];
				Sat RefSat1 = BDS3data.sat1[ref_bds3];
				Sat tempSat2 = BDS3data.sat2[k];
				Sat RefSat2 = BDS3data.sat2[ref_bds3];

				double length_O1, length_R1, length_O2, length_R2;

				//计算卫星到测站距离（4个距离）
				length_O1 = sqrt(pow((tempSat1.POS_X - POSITION_X1), 2) + pow((tempSat1.POS_Y - POSITION_Y1), 2) + pow((tempSat1.POS_Z - POSITION_Z1), 2));
				length_R1 = sqrt(pow((RefSat1.POS_X - POSITION_X1), 2) + pow((RefSat1.POS_Y - POSITION_Y1), 2) + pow((RefSat1.POS_Z - POSITION_Z1), 2));

				length_O2 = sqrt(pow((tempSat2.POS_X - POSITION_X2), 2) + pow((tempSat2.POS_Y - POSITION_Y2), 2) + pow((tempSat2.POS_Z - POSITION_Z2), 2));
				length_R2 = sqrt(pow((RefSat2.POS_X - POSITION_X2), 2) + pow((RefSat2.POS_Y - POSITION_Y2), 2) + pow((RefSat2.POS_Z - POSITION_Z2), 2));

				//三个坐标系数
				double Cof1 = (tempSat2.POS_X - POSITION_X2) / length_O2 - (RefSat2.POS_X - POSITION_X2) / length_R2;
				double Cof2 = (tempSat2.POS_Y - POSITION_Y2) / length_O2 - (RefSat2.POS_Y - POSITION_Y2) / length_R2;
				double Cof3 = (tempSat2.POS_Z - POSITION_Z2) / length_O2 - (RefSat2.POS_Z - POSITION_Z2) / length_R2;

				//系数矩阵
				B[tw][0] = Cof1;
				B[tw][1] = Cof2;
				B[tw][2] = Cof3;
				B[tw + num_n][0] = Cof1;
				B[tw + num_n][1] = Cof2;
				B[tw + num_n][2] = Cof3;
				B[tw + num_n][3 + tw] = lambda_L1_B;

				//误差改正
				DDerror = Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2);

				//双差形式的几何项和误差
				double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;
				double deltp1 = (tempSat2.data[0] - RefSat2.data[0]) - (tempSat1.data[0] - RefSat1.data[0]);  //伪距
				double deltL1 = ((tempSat2.data[2] - RefSat2.data[2]) - (tempSat1.data[2] - RefSat1.data[2])) * lambda_L1_B;  //载波
				L[tw][0] = deltp1 - DDgeo;         //伪距在前
				L[tw + num_n][0] = deltL1 - DDgeo; //载波在后

				for (int ii = 0; ii < SatNum[2]; ii++)
				{
					if (OldRefPRN[2] == SatPRN[2][ii])
					{
						Init = true;
						break;
					}
					else
						Init = false;
				}

				//若旧参考星被剔除，需初始化模糊度
				if (!Init)
				{
					X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1_B;
					Q[3 + tw][3 + tw] = 100.0 * 100.0;
				}

				if (!bInit)
				{
					X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1_B;       //初始化模糊度
				}

				for (int ii = 0; ii < NewSatNum[2]; ii++)
				{
					if (NewSatPRN[2][ii] == PRN)
					{
						X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1_B;   //初始化模糊度
						Q[3 + tw][3 + tw] = 100.0 * 100.0;
						break;
					}
				}
				tw++;
			}
		}

		//Galileo
		for (k = 0; k < Galileodata.sat_num; k++)
		{
			bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = Galileodata.sat1[k].numofsat;
			//参考卫星跳过
			if (PRN == Galref)
			{
				continue;
			}

			if (ref_gal < 0)
			{
				assert(0);
			}

			Sat tempSat1 = Galileodata.sat1[k];
			Sat RefSat1 = Galileodata.sat1[ref_gal];
			Sat tempSat2 = Galileodata.sat2[k];
			Sat RefSat2 = Galileodata.sat2[ref_gal];

			double length_O1, length_R1, length_O2, length_R2;

			//计算卫星到测站距离（4个距离）
			length_O1 = sqrt(pow((tempSat1.POS_X - POSITION_X1), 2) + pow((tempSat1.POS_Y - POSITION_Y1), 2) + pow((tempSat1.POS_Z - POSITION_Z1), 2));
			length_R1 = sqrt(pow((RefSat1.POS_X - POSITION_X1), 2) + pow((RefSat1.POS_Y - POSITION_Y1), 2) + pow((RefSat1.POS_Z - POSITION_Z1), 2));

			length_O2 = sqrt(pow((tempSat2.POS_X - POSITION_X2), 2) + pow((tempSat2.POS_Y - POSITION_Y2), 2) + pow((tempSat2.POS_Z - POSITION_Z2), 2));
			length_R2 = sqrt(pow((RefSat2.POS_X - POSITION_X2), 2) + pow((RefSat2.POS_Y - POSITION_Y2), 2) + pow((RefSat2.POS_Z - POSITION_Z2), 2));

			//三个坐标系数
			double Cof1 = (tempSat2.POS_X - POSITION_X2) / length_O2 - (RefSat2.POS_X - POSITION_X2) / length_R2;
			double Cof2 = (tempSat2.POS_Y - POSITION_Y2) / length_O2 - (RefSat2.POS_Y - POSITION_Y2) / length_R2;
			double Cof3 = (tempSat2.POS_Z - POSITION_Z2) / length_O2 - (RefSat2.POS_Z - POSITION_Z2) / length_R2;

			//系数矩阵
			B[tw][0] = Cof1;
			B[tw][1] = Cof2;
			B[tw][2] = Cof3;
			B[tw + num_n][0] = Cof1;
			B[tw + num_n][1] = Cof2;
			B[tw + num_n][2] = Cof3;
			B[tw + num_n][3 + tw] = lambda_L1_E;

			//误差改正
			DDerror = Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2);

			//双差形式的几何项和误差
			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;
			double deltp1 = (tempSat2.data[0] - RefSat2.data[0]) - (tempSat1.data[0] - RefSat1.data[0]);  //伪距
			double deltL1 = ((tempSat2.data[2] - RefSat2.data[2]) - (tempSat1.data[2] - RefSat1.data[2])) * lambda_L1_E;  //载波
			L[tw][0] = deltp1 - DDgeo;         //伪距在前
			L[tw + num_n][0] = deltL1 - DDgeo; //载波在后

			for (int ii = 0; ii < SatNum[3]; ii++)
			{
				if (OldRefPRN[3] == SatPRN[3][ii])
				{
					Init = true;
					break;
				}
				else
					Init = false;
			}

			//若旧参考星被剔除，需初始化模糊度
			if (!Init)
			{
				X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1_E;
				Q[3 + tw][3 + tw] = 100.0 * 100.0;
			}

			if (!bInit)
			{
				X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1_E;       //初始化模糊度
			}

			for (int ii = 0; ii < NewSatNum[3]; ii++)
			{
				if (NewSatPRN[3][ii] == PRN)
				{
					X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1_E;   //初始化模糊度
					Q[3 + tw][3 + tw] = 100.0 * 100.0;
					break;
				}
			}
			tw++;
		}

		if ((GPSdata.sat_num + BDS2data.sat_num + BDS3data.sat_num + Galileodata.sat_num) >= 4)
		{
			//权矩阵PG、法方程矩阵NG(逆矩阵即为参数方差)、参数求解(相对概率坐标的改正数)

			//B系数矩阵  L观测值矩阵  R误差矩阵  Qw动态噪声矩阵  Q协方差矩阵  X初始化模糊度

			//kalman滤波
			KalmanFilter(B, L, R, Qw, Q, X);

			double ratio = 1.0;
			CMatrix aXb(3, 1);
			CMatrix fix(num_n, 1);

			//lambda算法固定模糊度
			Resamb_LAMBDA(Q, X, aXb, ratio, fix);

			//恢复出估计的坐标值
			CMatrix X_POS(3, 1);
			X_POS[0][0] = POSITION_X2 - aXb[0][0];
			X_POS[1][0] = POSITION_Y2 - aXb[1][0];
			X_POS[2][0] = POSITION_Z2 - aXb[2][0];

			//真实坐标写成矩阵
			CMatrix XT(3, 1);
			XT[0][0] = P_X2;
			XT[1][0] = P_Y2;
			XT[2][0] = P_Z2;

			//与准确值作差得到XYZ坐标系下的误差
			CMatrix ErrorX = X_POS - XT;

			//将误差投影到NEU东北天坐标系
			CMatrix ErrorX_NEU = TT * ErrorX;

			fprintf(fp3, "EpNum: %2d, %2d, %2d, %10.6f, %2d, %2d, %2d, %2d, %15.8g,%15.8g, %15.8g, %15.8g\n", i + 1, tempData2.hour, tempData2.minute, tempData2.second,
				GPSdata.sat_num, BDS2data.sat_num, BDS3data.sat_num, Galileodata.sat_num, ratio, ErrorX_NEU[0][0], ErrorX_NEU[1][0], ErrorX_NEU[2][0]);
		}

		for (int idx = 0; idx < 4; idx++)
		{
			OldRefPRN[idx] = RefPRN[idx];
			for (int ii = 0; ii < OldSatNum[idx]; ii++)
			{
				OldSatPRN[idx][ii] = 0;
			}
			OldSatNum[idx] = SatNum[idx];
			for (int ii = 0; ii < SatNum[idx]; ii++)
			{
				OldSatPRN[idx][ii] = SatPRN[idx][ii];
			}
			for (int ii = 0; ii < SatNum[idx]; ii++)
			{
				SatPRN[idx][ii] = 0;
			}
			SatNum[idx] = 0;
			for (int ii = 0; ii < NewSatNum[idx]; ii++)
			{
				NewSatPRN[idx][ii] = 0;
			}
			NewSatNum[idx] = 0;
		}
		bInit = true;  //初始化成功
	}
	fclose(fp3);
	return true;
}