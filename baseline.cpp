#include "stdafx.h"
#include "Myheader.h"
#include "Matrix.h"
#include "Lambda.h"
#include <ctime>
#include <cmath>
#include <set>
#include <vector>
#include <fstream>
#include <string>
#include <Windows.h>
#include <assert.h>
#include <map>

using namespace std;

//声明getGPSweek函数功能，得到GPS周
double getGPSweek(double& gpstime);

//声明EpochCommonPro函数，对基站和流动站卫星观测数据进行共视处理
void EpochCommonPro(Obs_epoch tempData1, Obs_epoch tempData2, Obs_epoch& useData);

//声明SelectRef函数，选择参考卫星
void SelectRef(Obs_epoch tempData, int& RefNo);

//声明ReadSPP_KIN函数，读取spp文件（针对本人自定义的spp文件格式，可修改）
bool ReadSPP_KIN(string filename, observe_spp& sppfile);

//声明GetDDR_Leo函数，依据误差传播定律获取噪声矩阵（载波）
void GetDDR_Leo(Obs_epoch GPSData, vector<int> useprn, int GPSRef, CMatrix& R);

//声明GetDDR_Peo函数，依据误差传播定律获取噪声矩阵（伪距）
void GetDDR_Peo(Obs_epoch GPSData, vector<int> useprn, int GPSRef, CMatrix& R);

//函数功能：初始化状态方程,中长基线用
void InitState_long(CMatrix& QX)
{
	int num_x = QX.Row;
	for (int idx = 0; idx < num_x; idx++)
	{
		if (idx <= 2)
		{
			QX[idx][idx] = 30.0 * 30.0;    //初始化位置
		}
		else
		{
			QX[idx][idx] = 5.0 * 5.0;      //初始化模糊度
		}
	}
}

//函数功能：初始化状态方程,三频用
void InitState_Tri(CMatrix& QX)
{
	int num_x = QX.Row;
	int num_s = (num_x - 3) / 3;
	for (int idx = 0; idx < num_x; idx++)
	{
		if (idx <= 2)
		{
			QX[idx][idx] = 30.0 * 30.0;    //初始化位置
		}
		else if (idx <= 2 + num_s)
		{
			QX[idx][idx] = 5.0 * 5.0;      //初始化超宽巷模糊度
		}
		else if(idx<=2+num_s*2)
		{
			QX[idx][idx] = 30.0 * 30.0;    //初始化宽巷模糊度
		}
		else if(idx<=2+num_s*3)
		{
			QX[idx][idx] = 100.0 * 100.0;    //初始化基础模糊度
		}
	}
}

//声明InitState函数，初始化Q矩阵
void InitState(CMatrix& QX);

//重载CheckNewSats函数，检测新卫星，BDS2和BDS3
void CheckNewSats(int* SatNum, int* oldSatNum, int SatPRN[2][30], int oldSatPRN[2][30], int NewSatPRN[2][30], int* NewSatNum)
{
	for (int i = 0; i < 2; i++)
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

//声明GetGRCQValue函数，生成转换矩阵
bool GetGRCQValue(int PRN_X[], int PRN_old[], int num_sat, int num_sat_old, int Ref, int Ref_old, CMatrix& TT_Matrix);

//声明ProGRCTrans函数，拼接矩阵
void ProGRCTrans(CMatrix& TT_G, CMatrix& TT_C, CMatrix& TT_ALL);

//重载UpdataState函数，更新协方差，BDS2和BDS3
void UpdataState(CMatrix& QX, CMatrix& X, int* SatNum, int* oldSatNum, int SatPRN[2][30], int oldSatPRN[2][30], int* RefSat, int* OldRef)
{
	CMatrix T[2], ALL_T, Matrix_t;
	int sys = 0, oldsys = 0, ALL_NUM, oldALL_NUM;
	if (SatNum[0] >= 2)
	{
		sys++;
	}
	if (SatNum[1] >= 2)
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
	ALL_NUM = SatNum[0] + SatNum[1] - sys;
	oldALL_NUM = oldSatNum[0] + oldSatNum[1] - oldsys;
	//L1或者L1变化矩阵
	for (int i = 0; i < 2; i++)
	{
		if (SatNum[i] > 1)
		{
			if (oldSatNum[i] > 0)
			{
				T[i].SetSize(SatNum[i] - 1, oldSatNum[i] - 1);//当前历元GPS卫星数及上一历元GPS卫星数
				GetGRCQValue(SatPRN[i], oldSatPRN[i], SatNum[i], oldSatNum[i], RefSat[i], OldRef[i], T[i]);
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
	ProGRCTrans(T[0], T[1], ALL_T);

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

//重载UpdataState函数，更新协方差，BDS2和BDS3，三频合算使用
void UpTridataState(CMatrix& QX, CMatrix& X, int* SatNum, int* oldSatNum, int SatPRN[2][30], int oldSatPRN[2][30], int* RefSat, int* OldRef)
{
	CMatrix T[2], Single_T, double_T, ALL_T, Matrix_t;
	int sys = 0, oldsys = 0, ALL_NUM, oldALL_NUM;
	if (SatNum[0] >= 2)
	{
		sys++;
	}
	if (SatNum[1] >= 2)
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
	ALL_NUM = SatNum[0] + SatNum[1] - sys;
	oldALL_NUM = oldSatNum[0] + oldSatNum[1] - oldsys;
	//L1或者L1变化矩阵
	for (int i = 0; i < 2; i++)
	{
		if (SatNum[i] > 1)
		{
			if (oldSatNum[i] > 0)
			{
				T[i].SetSize(SatNum[i] - 1, oldSatNum[i] - 1);//当前历元GPS卫星数及上一历元GPS卫星数
				GetGRCQValue(SatPRN[i], oldSatPRN[i], SatNum[i], oldSatNum[i], RefSat[i], OldRef[i], T[i]);
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

	Single_T.SetSize(ALL_NUM, oldALL_NUM);
	double_T.SetSize(2 * ALL_NUM, 2 * oldALL_NUM);
	ALL_T.SetSize(3 * ALL_NUM, 3 * oldALL_NUM);
	ProGRCTrans(T[0], T[1], Single_T);
	ProGRCTrans(Single_T, Single_T, double_T);
	ProGRCTrans(Single_T, double_T, ALL_T);
	
	Matrix_t.SetSize(3 * ALL_NUM + 3, 3 * oldALL_NUM + 3);
	for (int j = 0; j < 3; j++)
	{
		Matrix_t[j][j] = 1.0;
	}

	for (int ii = 0; ii < 3 * ALL_NUM; ii++)
	{
		for (int jj = 0; jj < 3 * oldALL_NUM; jj++)
		{
			Matrix_t[ii + 3][jj + 3] = ALL_T[ii][jj];
		}
	}

	QX = Matrix_t * QX * Matrix_t.T();
	X = Matrix_t * X;
}

//声明Error_Correction函数，误差改正
double Error_Correction(Sat& tempSat1, Sat& RefSat1, Sat& tempSat2, Sat& RefSat2);

//声明KalmanFilter，卡尔曼滤波
void KalmanFilter(CMatrix& B, CMatrix& L, CMatrix& R, CMatrix& Qw, CMatrix& Q, CMatrix& X);

//声明Resamb_LAMBDA函数，Lambda算法固定模糊度
void Resamb_LAMBDA(CMatrix& Q, CMatrix& X, CMatrix& aXb, double& ratio, CMatrix& fixX);

//系数对应的频率为L1,L2,L5
double GPS(Sat& tempSat, double a, double b, double c)
{
	const double GPS_TriP = (a * f1G * tempSat.data[0] + b * f2G * tempSat.data[1] + c * f5G * tempSat.data[4]) / (a * f1G + b * f2G + c * f5G);
	return GPS_TriP;
}
//系数对应的频率为B1I,B2I,B3I
double BDS2(Sat& tempSat, double a, double b, double c)
{
	const double BDS2_TriP = (a * f1C * tempSat.data[0] + b * f2C * tempSat.data[1] + c * f3C * tempSat.data[4]) / (a * f1C + b * f2C + c * f3C);
	return BDS2_TriP;
}
//系数对应的频率为B1I,B2a,B3I,需要看ssp文件各频率观测值顺序
double BDS3(Sat& tempSat, double a, double b, double c)
{ 
	const double BDS3_TriP = (a * f1B * tempSat.data[0] + b * f2B * tempSat.data[1] + c * f3B * tempSat.data[4]) / (a * f1B + b * f2B + c * f3B);
	return BDS3_TriP;
}
//函数功能：三频伪距线性组合，操作数0为GPS，1为BDS2，2为BDS2，单位为米
double (*oper_TriP[])(Sat& tempSat, double a, double b, double c) = { GPS,BDS2,BDS3 };

//系数对应的频率为L1,L2,L5
double GPS(Sat& tempSat, int i, int j, int k)
{
	const double GPS_TriPhi = i * tempSat.data[2] + j * tempSat.data[3] + k * tempSat.data[5];
	return GPS_TriPhi;
}

//系数对应的频率为B1I,B2I,B3I
double BDS2(Sat& tempSat, int i, int j, int k)
{
	const double BDS2_TriPhi= i * tempSat.data[2] + j * tempSat.data[3] + k * tempSat.data[5];
	return BDS2_TriPhi;
}

//系数对应的频率为B1I,B2a,B3I,需要看ssp文件各频率观测值顺序
double BDS3(Sat& tempSat, int i, int j, int k)
{
	const double BDS3_TriPhi = i * tempSat.data[2] + j * tempSat.data[3] + k * tempSat.data[5];
	return BDS3_TriPhi;
}

//函数功能：三频载波线性组合，操作数0为GPS，1为BDS2，单位为周
double (*oper_TriPhi[])(Sat& tempSat, int i, int j, int k) = { GPS,BDS2,BDS3 };

//声明OnCycleDetect函数，周跳探测
void OnCycleDetect(Obs_epoch& tempEpoch, Obs_epoch& preEpoch);

//函数功能：三频周跳探测
void TriOnCycleDetect(Obs_epoch& tempEpoch, Obs_epoch& preEpoch)
{
	double slip1 = 0;
	double slip2 = 0;
	double slip3 = 0;
	double slip4 = 0;
	double slip5 = 0;
	for (int i = 0; i < tempEpoch.sat_num; i++) {
		Sat* sat = &tempEpoch.sat[i];
		if (sat->judge_use)
			continue;
		if (sat->sattype != "C" && sat->sattype != "B")
			continue;
		for (int j = 0; j < preEpoch.sat_num; j++) {
			if ((tempEpoch.sat[i].sattype == preEpoch.sat[j].sattype) && (tempEpoch.sat[i].numofsat == preEpoch.sat[j].numofsat)) {
				Sat* presat = &preEpoch.sat[j];
				if (presat->judge_use)
					continue;
				if (presat->sattype != "C" && presat->sattype != "B")
					continue;
				//BDS-2
				if (presat->sattype == "C") {
					slip1 = sat->data[2] - presat->data[2] - (FREQ1_BDS / FREQ2_BDS) * (sat->data[3] - presat->data[3]);//电离层残差
					slip2 = (-(sat->data[2] - sat->data[3]) + 1 / lambda_L_wide_C * lambda_L_narrow_C * (sat->data[0] / lambda_L1_C + sat->data[1] / lambda_L2_C))
						- (-(presat->data[2] - presat->data[3]) + 1 / lambda_L_wide_C * lambda_L_narrow_C * (presat->data[0] / lambda_L1_C + presat->data[1] / lambda_L2_C)); //MW

					slip3 = sat->data[2] - presat->data[2] - (FREQ1_BDS / FREQ3_BDS) * (sat->data[5] - presat->data[5]);//电离层残差
					slip4 = (-(sat->data[2] - sat->data[5]) + 1 / lambda_L13_wide_C * lambda_L13_narrow_C * (sat->data[0] / lambda_L1_C + sat->data[4] / lambda_L3_C))
						- (-(presat->data[2] - presat->data[5]) + 1 / lambda_L13_wide_C * lambda_L13_narrow_C * (presat->data[0] / lambda_L1_C + presat->data[4] / lambda_L3_C)); //MW

					slip5 = sat->data[3] - presat->data[3] - (FREQ2_BDS / FREQ3_BDS) * (sat->data[5] - presat->data[5]);//电离层残差
				}
				//BDS-3
				else if (presat->sattype == "B") {
					slip1 = sat->data[2] - presat->data[2] - (f1B / f2B) * (sat->data[3] - presat->data[3]);//电离层残差
					slip2 = (-(sat->data[2] - sat->data[3]) + 1 / lambda_L_wide_B * lambda_L_narrow_B * (sat->data[0] / lambda_L1_B + sat->data[1] / lambda_L2_B))
						- (-(presat->data[2] - presat->data[3]) + 1 / lambda_L_wide_B * lambda_L_narrow_B * (presat->data[0] / lambda_L1_B + presat->data[1] / lambda_L2_B));//MW

					slip3 = sat->data[2] - presat->data[2] - (f1B / f3B) * (sat->data[5] - presat->data[5]);//电离层残差
					slip4 = (-(sat->data[2] - sat->data[5]) + 1 / lambda_L13_wide_B * lambda_L13_narrow_B * (sat->data[0] / lambda_L1_B + sat->data[4] / lambda_L3_B))
						- (-(presat->data[2] - presat->data[5]) + 1 / lambda_L13_wide_B * lambda_L13_narrow_B * (presat->data[0] / lambda_L1_B + presat->data[4] / lambda_L3_B));//MW

					slip5 = sat->data[3] - presat->data[3] - (f2B/ f3B) * (sat->data[5] - presat->data[5]);//电离层残差
				}
				if (fabs(slip1) > 0.05 || fabs(slip2) > 10 || fabs(slip3) > 0.05 || fabs(slip4) > 10 || fabs(slip5) > 0.05) {
					sat->judge_use = 1;
				}
				break;
			}
		}
	}
}

//函数功能，部分模糊度固定，按超宽巷，宽巷和窄巷（基础模糊度）顺序固定
void ParAR_LAMBDA(CMatrix& X, CMatrix& Q, CMatrix& aXb, int& satnum, double ratio[3], CMatrix& fixX)
{
	CMatrix Qaa, Qbb, Qab, Xa, Xb, fix;
	Qbb.SetSize(3, 3);
	Qaa.SetSize(satnum, satnum);
	Qab.SetSize(3, satnum);
	Xb.SetSize(3, 1);
	Xa.SetSize(satnum, 1);

	//取出分块矩阵,超宽巷和基线位置解
	for (int ii = 0; ii < 3; ii++)
		for (int jj = 0; jj < 3; jj++)
		{
			Qbb[ii][jj] = Q[ii][jj];            //位置解方差
		}
	for (int ii = 0; ii < satnum; ii++)
		for (int jj = 0; jj < satnum; jj++)
		{
			Qaa[ii][jj] = Q[3 + ii][3 + jj];    //超宽巷模糊度方差
		}
	for (int ii = 0; ii < 3; ii++)
		for (int jj = 0; jj < satnum; jj++)
		{
			Qab[ii][jj] = Q[ii][jj + 3];        //位置解和超宽巷模糊度的协方差
		}
	for (int ii = 0; ii < 3; ii++)
		Xb[ii][0] = X[ii][0];          //基线位置解
	for (int ii = 0; ii < satnum; ii++)
		Xa[ii][0] = X[ii + 3][0];      //超宽巷模糊度

	//LAMBDA解算
	//超宽巷固定
	CLambda Elambda;
	Elambda.afloat.SetSize(satnum, 1);
	Elambda.Qahat.SetSize(satnum, satnum);
	Elambda.afloat = Xa;
	Elambda.Qahat = Qaa;
	//lambda算法固定超宽巷模糊度
	Elambda.lambda2(Elambda.afloat, Elambda.Qahat);
	//超宽巷固定解
	CMatrix Efix;
	CMatrix EaXb;
	EaXb.SetSize(3, 1);
	Efix.SetSize(satnum, 1);
	for (int ii = 0; ii < satnum; ii++)
		Efix[ii][0] = Elambda.afixed[ii][0];
	//回代,得到超宽巷基线解
	EaXb = Xb - Qab * Qaa.InvertGaussJordan() * (Xa - Efix);
	Qbb = Qbb - Qab * Qaa.InvertGaussJordan() * Qab.T();

	/*test*/
	// Efix.MyTRACE();
	// cout << endl;
	// Xa.MyTRACE();
	// cout << endl;
	// cout << Elambda.ratio << endl;
	/*test*/

	//检验超宽巷是否固定
	double Eratio;
	Eratio = Elambda.ratio;
	if (Eratio < 3.0)  //未固定
	{
		ratio[0] = Eratio;
		aXb = EaXb;
	}
	else               //已固定
	{
		CMatrix EQbb, EQab, EXb;
		EQbb.SetSize(2 * satnum, 2 * satnum);
		EQab.SetSize(satnum, 2 * satnum);
		EXb.SetSize(2 * satnum, 1);
		//宽巷固定
		CLambda Wlambda;
		Wlambda.afloat.SetSize(satnum, 1);
		Wlambda.Qahat.SetSize(satnum, satnum);
		//取出分块矩阵,用固定的超宽巷模糊度约束宽巷和窄巷模糊度
		for (int ii = 0; ii < 2 * satnum; ii++)
			for (int jj = 0; jj < 2 * satnum; jj++)
				EQbb[ii][jj] = Q[ii + 3 + satnum][jj + 3 + satnum];       //宽巷和窄巷模糊度方差
		for (int ii = 0; ii < satnum; ii++)
			for (int jj = 0; jj < 2 * satnum; jj++)
				EQab[ii][jj] = Q[ii + 3][jj + 3 + satnum];                //超宽巷和其他模糊度的协方差
		for (int ii = 0; ii < 2 * satnum; ii++)
			EXb[ii][0] = X[ii + 3 + satnum][0];                           //宽巷和窄巷模糊度
		//回代约束宽巷和窄巷模糊度
		EXb = EXb - EQab.T() * Qaa.InvertGaussJordan() * (Xa - Efix);
		EQbb = EQbb - EQab.T() * Qaa.InvertGaussJordan() * EQab;
		//取出宽巷模糊度及其方差
		for (int ii = 0; ii < satnum; ii++)
			for (int jj = 0; jj < satnum; jj++)
				Wlambda.Qahat[ii][jj] = EQbb[ii][jj];
		for (int ii = 0; ii < satnum; ii++)
			Wlambda.afloat[ii][0] = EXb[ii][0];
		CMatrix wXa = Wlambda.afloat;     //记录宽巷模糊度的浮点解
		CMatrix wQaa = Wlambda.Qahat;     //记录宽巷模糊度的方差-协方差阵
		//lambda算法固定宽巷模糊度
		Wlambda.lambda2(Wlambda.afloat, Wlambda.Qahat);
		//宽巷固定解
		CMatrix WaXb, Wfix;
		WaXb.SetSize(3, 1);
		Wfix.SetSize(satnum, 1);
		for (int ii = 0; ii < satnum; ii++)
			Wfix[ii][0] = Wlambda.afixed[ii][0];
		//取出基线位置解和宽巷模糊度的方差、协方差
		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < satnum; jj++)
				Qab[ii][jj] = Q[ii][jj + 3 + satnum];
		//回代,得到宽巷基线解
		WaXb = Xb - Qab * wQaa.InvertGaussJordan() * (wXa - Wfix);
		Qbb = Qbb - Qab * wQaa.InvertGaussJordan() * Qab.T();
		//检验宽巷是否固定
		double Wratio = Wlambda.ratio;

		/*test*/
		/*Wfix.MyTRACE();
		cout << endl;
		wXa.MyTRACE();
		cout << endl;
		cout << Wlambda.ratio << endl;*/
		/*test*/

		if (Wratio < 3.0) //未固定
		{
			ratio[1] = Wratio;
			aXb = WaXb;
		}
		else              //宽巷已固定
		{
			CMatrix WQbb, WQab, WXb;
			WQbb.SetSize(satnum, satnum);
			WQab.SetSize(satnum, satnum);
			WXb.SetSize(satnum, 1);
			//窄巷（基础模糊度）固定
			CLambda Nlambda;
			Nlambda.afloat.SetSize(satnum, 1);
			Nlambda.Qahat.SetSize(satnum, satnum);
			//取出分块矩阵,用固定的宽巷模糊度约束窄巷模糊度
			for (int ii = 0; ii < satnum; ii++)
				for (int jj = 0; jj < satnum; jj++)
					WQbb[ii][jj] = EQbb[ii + satnum][jj + satnum];
			for (int ii = 0; ii < satnum; ii++)
				for (int jj = 0; jj < satnum; jj++)
					WQab[ii][jj] = EQbb[ii][jj + satnum];
			for (int ii = 0; ii < satnum; ii++)
				WXb[ii][0] = EXb[ii + satnum][0];
			//回代约束窄巷模糊度
			WXb = WXb - WQab.T() * wQaa.InvertGaussJordan() * (wXa - Wfix);   //窄巷模糊度的浮点解
			WQbb = WQbb - WQab.T() * wQaa.InvertGaussJordan() * WQab;     //窄巷模糊度的方差-协方差阵
			Nlambda.afloat = WXb;
			Nlambda.Qahat = WQbb;
			//lambda算法固定窄巷模糊度
			Nlambda.lambda2(Nlambda.afloat, Nlambda.Qahat);
			//窄巷固定解
			CMatrix NaXb, Nfix;
			NaXb.SetSize(3, 1);
			Nfix.SetSize(satnum, 1);
			for (int ii = 0; ii < satnum; ii++)
				Nfix[ii][0] = Nlambda.afixed[ii][0];
			//取出基线位置解和窄巷模糊度的方差、协方差
			for (int ii = 0; ii < 3; ii++)
				for (int jj = 0; jj < satnum; jj++)
					Qab[ii][jj] = Q[ii][jj + 3 + 2 * satnum];
			//回代,得到窄巷基线解
			NaXb = Xb - Qab * WQbb.InvertGaussJordan() * (WXb - Nfix);
			Qbb = Qbb - Qab * WQbb.InvertGaussJordan() * Qab.T();
			ratio[2] = Nlambda.ratio;
			aXb = NaXb;

			/*test*/
			/*Nfix.MyTRACE();
			cout << endl;
			WXb.MyTRACE();
			cout << endl;
			cout << Nlambda.ratio << endl;*/
			/*test*/
		}
	}
}

//函数功能，TCAR法固定模糊度,采用几何相关模型,固定站
void TCAR_fix(string& filename1, string& filename2)
{
	//默认第一个站为基准站，第二个为流动站
	observe_spp sppfile1;
	observe_spp sppfile2;
	ReadSPP_KIN(filename1, sppfile1);
	ReadSPP_KIN(filename2, sppfile2);
	bool bInit = false; //初始化
	bool Init = true;   //周跳后初始化

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

	//输出差分定位结果用
	string name = filename2.substr(0, 3) + "_TriRTK.txt";
	const char* output_filename = name.c_str();
	FILE* fp;
	fopen_s(&fp, output_filename, "w");

	int OldSatPRN[2][30] = {{0}};
	int OldRefPRN[2] = { -1,-1 };
	int OldSatNum[2] = { 0 };

	CMatrix Q_EWL1, Q_EWL2, Q_WL, X_EWL1, X_EWL2, X_WL;

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
			//TriOnCycleDetect(tempData2, preData2);     //判断此此历元与上一个历元中的卫星载波测量值是否有周跳,如果有周跳就将此历元中的卫星记为不可用
			OnCycleDetect(tempData2, preData2);
		}

		//基准站数据相应时刻必须在前
		int j = 0;
		for (j = posk; j < sppfile1.liyuan_num; j++)
		{
			if (sppfile1.epoch[j].GPSTIME <= tempData2.GPSTIME)
				// ReSharper disable once CppRedundantControlFlowJump
				continue;
			break;
		}

		//回溯一个历元：posk即为要找的最近的历元且满足posk的时刻小于tempData2的时间
		if (j > 0) posk = --j;
		else
		{
			posk = j;
			continue;
		}

		//判断时间差：两个是否相隔超过历元间隔时间
		if (fabs(sppfile1.epoch[j].GPSTIME - tempData2.GPSTIME) >= (sppfile1.INTERVAL + 1))
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
			//TriOnCycleDetect(tempData1, preData1);
			OnCycleDetect(tempData1, preData1);
		}

		//差分定位共视处理
		EpochCommonPro(tempData1, tempData2, tempData);

		//分离BDS-2、BDS-3卫星，便于后续单系统、多系统处理
		Obs_epoch bds2_data;
		Obs_epoch bds3_data;

		int SatPRN[2][30] = {{0}};		//当前历元卫星号 依次存储
		int SatNum[2] = { 0 };	//当前历元各系统卫星数目
		int RefPRN[2] = { 0 };  //当前历元各系统参考卫星号

		int NewSatNum[2] = { 0 };
		int NewSatPRN[2][30] = {{0}};

		int k = 0;
		int m = 0;
		int bds2num, bds3num;
		bds2num = bds3num = 0;
		for (m = 0; m < tempData.sat_num; m++)
		{
			if (tempData.sat1[m].sattype == "C")
			{
				bds2_data.sat1.push_back(tempData.sat1[m]);
				bds2_data.sat2.push_back(tempData.sat2[m]);
				SatPRN[0][bds2num] = tempData.sat1[m].numofsat;
				bds2num++;
			}

			if (tempData.sat1[m].sattype == "B")
			{
				bds3_data.sat1.push_back(tempData.sat1[m]);
				bds3_data.sat2.push_back(tempData.sat2[m]);
				SatPRN[1][bds3num] = tempData.sat1[m].numofsat;
				bds3num++;
			}
		}

		//各类型卫星数
		SatNum[0] = bds2num;
		SatNum[1] = bds3num;

		bds2_data.sat_num = bds2_data.sat1.size();
		bds3_data.sat_num = bds3_data.sat1.size();
		tempData.GPSTIME = tempData2.GPSTIME;

		//记录新卫星
		if (bInit == true)
		{
			for (int idx = 0; idx < 2; idx++)
			{
				set<int>s;
				for (int ii = 0; ii < OldSatNum[idx]; ii++)
					s.insert(OldSatPRN[idx][ii]);
				for (int ii = 0; ii < SatNum[idx]; ii++)
				{
					if (s.find(SatPRN[idx][ii]) == s.end())
					{
						NewSatNum[idx]++;
						for (int iii = 0; iii < NewSatNum[idx]; iii++)
							NewSatPRN[idx][iii] = SatPRN[idx][ii];
					}
				}
			}
		}

		//卫星数至少要4颗
		if (bds2_data.sat_num + bds3_data.sat_num < 4)
		{
			printf("卫星数少于4颗\n");
			bInit = false;
			continue;
		}

		//各系统至少要2颗卫星,区分参考星和普通星
		if (bds2_data.sat_num < 2)
		{
			printf("BDS2卫星数少于2颗\n");
			Init = false;
			continue;
		}
		if (bds3_data.sat_num < 2)
		{
			printf("BDS3卫星数少于2颗\n");
			Init = false;
			continue;
		}

		//采用站间星间二次差分模型，选择参考卫星
		int BDS2ref = -1;
		int BDS3ref = -2;

		BDS2ref = OldRefPRN[0];
		BDS3ref = OldRefPRN[1];

		//选择最大高度角卫星作为参考卫星，不同类型卫星分开看
		SelectRef(bds2_data, BDS2ref);
		SelectRef(bds3_data, BDS3ref);

		RefPRN[0] = BDS2ref; //BDS-2参考卫星序号
		RefPRN[1] = BDS3ref; //BDS-3参考卫星序号

		//当前历元BDS-2、BDS-3卫星序列
		int PRN_X_B2[20];
		int PRN_X_B3[30];

		//数列初始化
		memset(PRN_X_B2, 0, sizeof(PRN_X_B2));
		memset(PRN_X_B3, 0, sizeof(PRN_X_B3));

		//储存历元中的卫星序列
		for (k = 0; k < bds2_data.sat_num; k++)
			PRN_X_B2[k] = bds2_data.sat1[k].numofsat;
		for (k = 0; k < bds3_data.sat_num; k++)
			PRN_X_B3[k] = bds3_data.sat1[k].numofsat;

		//各卫星个数
		int num_X_B2 = bds2_data.sat_num;
		int num_X_B3 = bds3_data.sat_num;

		//各系统双差方程个数
		int numbds2_x = num_X_B2 - 1;
		int numbds3_x = num_X_B3 - 1;

		//多系统融合双差方程个数（载波观测方程+伪距观测方程)，松组合
		int num_x_p = (numbds2_x + numbds3_x) * 2;

		//待估参数个数：ΔX，ΔY，ΔZ，每颗卫星的模糊度
		int num_x = numbds2_x + numbds3_x + 3;

		//模糊度个数
		int num_n = numbds2_x + numbds3_x;

		//定位解算相关的矩阵
		//多系统
		CMatrix B_EWL1(num_x_p, num_x);           //系数矩阵,EWL1
		CMatrix B_EWL2(num_x_p, num_x);           //系数矩阵,EWL2
		CMatrix B_WL(num_n, num_x);               //系数矩阵（只有载波观测值）,WL

		CMatrix L_EWL1(num_x_p, 1);               //观测值矩阵,EWL1
		CMatrix L_EWL2(num_x_p, 1);               //观测值矩阵,EWL2
		CMatrix L_WL(num_n, 1);                   //观测值矩阵（只有载波观测值）,WL

		CMatrix R_EWL1(num_x_p, num_x_p);         //误差矩阵（逆阵即为权矩阵）
		CMatrix R_EWL2(num_x_p, num_x_p);         //误差矩阵（逆阵即为权矩阵）
		CMatrix R_WL(num_n, num_n);               //误差矩阵（逆阵即为权矩阵）

		CMatrix Qw_EWL1(num_x, num_x);            //动态噪声阵,EWL1
		CMatrix Qw_EWL2(num_x, num_x);            //动态噪声阵,EWL2
		CMatrix Qw_WL(num_x, num_x);              //动态噪声阵,WL
		
		//单BDS2
		CMatrix RC2_P(numbds2_x, numbds2_x);      //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RC2_L(numbds2_x, numbds2_x);      //误差矩阵（逆阵即为权矩阵） 载波

		//单BDS3
		CMatrix RC3_P(numbds3_x, numbds3_x);      //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RC3_L(numbds3_x, numbds3_x);      //误差矩阵（逆阵即为权矩阵） 载波

		//位置的动态噪声阵，如果没有周跳，模糊度不变，所以模糊度的噪声为0，只有位置噪声
		for (int idx = 0; idx < 3; idx++) {
			Qw_EWL1[idx][idx] = 100.0 * 100.0;
			Qw_EWL2[idx][idx] = 100.0 * 100.0;
			Qw_WL[idx][idx] = 30.0 * 30.0;
		}

		//更新BDS2和BDS2的协方差矩阵Q和结果矩阵X
		if (!bInit)
		{
			Q_EWL1.SetSize(num_x, num_x);
			Q_EWL2.SetSize(num_x, num_x);
			Q_WL.SetSize(num_x, num_x);
			X_EWL1.SetSize(num_x, 1);
			X_EWL2.SetSize(num_x, 1);
			X_WL.SetSize(num_x, 1);
			InitState_long(Q_EWL1);
			InitState_long(Q_EWL2);
			InitState_long(Q_WL);
		}
		else
		{
			//检测新卫星
			CheckNewSats(SatNum, OldSatNum, SatPRN, OldSatPRN, NewSatPRN, NewSatNum);

			//更新协方差和X矩阵
			UpdataState(Q_EWL1, X_EWL1, SatNum, OldSatNum, SatPRN, OldSatPRN, RefPRN, OldRefPRN);
			UpdataState(Q_EWL2, X_EWL2, SatNum, OldSatNum, SatPRN, OldSatPRN, RefPRN, OldRefPRN);
			UpdataState(Q_WL, X_WL, SatNum, OldSatNum, SatPRN, OldSatPRN, RefPRN, OldRefPRN);
		}

		//各类型卫星序号存储容器
		vector<int> BDS2_PRN;
		vector<int> BDS3_PRN;

		for (k = 0; k < bds2_data.sat_num; k++)
			BDS2_PRN.push_back(PRN_X_B2[k]);
		for (k = 0; k < bds3_data.sat_num; k++)
			BDS3_PRN.push_back(PRN_X_B3[k]);

		//伪距
		GetDDR_Peo(bds2_data, BDS2_PRN, BDS2ref, RC2_P);
		GetDDR_Peo(bds3_data, BDS3_PRN, BDS3ref, RC3_P);
		//载波
		GetDDR_Leo(bds2_data, BDS2_PRN, BDS2ref, RC2_L);
		GetDDR_Leo(bds3_data, BDS3_PRN, BDS3ref, RC3_L);

		//融合系统的R矩阵，分块对角矩阵
		//EWL1
		//伪距
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R_EWL1[k1][k2] = RC2_P[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R_EWL1[k1 + numbds2_x][k2 + numbds2_x] = RC3_P[k1][k2];
			}
		//载波
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R_EWL1[k1 + numbds2_x + numbds3_x][k2 + numbds2_x + numbds3_x] = RC2_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R_EWL1[k1 + numbds2_x * 2 + numbds3_x][k2 + numbds2_x * 2 + numbds3_x] = RC3_L[k1][k2];
			}

		//EWL2
		//伪距
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R_EWL2[k1][k2] = RC2_P[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R_EWL2[k1 + numbds2_x][k2 + numbds2_x] = RC3_P[k1][k2];
			}
		//载波
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R_EWL2[k1 + numbds2_x + numbds3_x][k2 + numbds2_x + numbds3_x] = RC2_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R_EWL2[k1 + numbds2_x * 2 + numbds3_x][k2 + numbds2_x * 2 + numbds3_x] = RC3_L[k1][k2];
			}

		//WL
		//宽巷
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R_WL[k1][k2] = RC2_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R_WL[k1 + numbds2_x][k2 + numbds2_x] = RC3_L[k1][k2];
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
		int ref_bds2 = -1;
		int ref_bds3 = -1;

		for (k = 0; k < bds2_data.sat_num; k++)
		{
			if (bds2_data.sat1[k].numofsat == BDS2ref)
			{
				ref_bds2 = k;
				break;
			}
		}

		for (k = 0; k < bds3_data.sat_num; k++)
		{
			if (bds3_data.sat1[k].numofsat == BDS3ref)
			{
				ref_bds3 = k;
				break;
			}
		}
		
		int tp = 0;  //记录单系统方程数
		int tw = 0;  //记录多系统方程数
		double DDerror = 0;  //初始化误差改正项
		//BDS-2
		for (k = 0; k < bds2_data.sat_num; k++)
		{
			//bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = bds2_data.sat1[k].numofsat;
			//参考卫星跳过
			if (PRN == BDS2ref)
			{
				continue;
			}

			if (ref_bds2 < 0)
			{
				assert(0);
			}

			Sat tempSat1 = bds2_data.sat1[k];
			Sat RefSat1 = bds2_data.sat1[ref_bds2];
			Sat tempSat2 = bds2_data.sat2[k];
			Sat RefSat2 = bds2_data.sat2[ref_bds2];

			//误差改正
			DDerror = Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2);

			//计算卫星到测站距离（4个距离）
			double length_O1, length_R1, length_O2, length_R2;
			length_O1 = sqrt(pow((tempSat1.POS_X - POSITION_X1), 2) + pow((tempSat1.POS_Y - POSITION_Y1), 2) + pow((tempSat1.POS_Z - POSITION_Z1), 2));
			length_R1 = sqrt(pow((RefSat1.POS_X - POSITION_X1), 2) + pow((RefSat1.POS_Y - POSITION_Y1), 2) + pow((RefSat1.POS_Z - POSITION_Z1), 2));

			length_O2 = sqrt(pow((tempSat2.POS_X - POSITION_X2), 2) + pow((tempSat2.POS_Y - POSITION_Y2), 2) + pow((tempSat2.POS_Z - POSITION_Z2), 2));
			length_R2 = sqrt(pow((RefSat2.POS_X - POSITION_X2), 2) + pow((RefSat2.POS_Y - POSITION_Y2), 2) + pow((RefSat2.POS_Z - POSITION_Z2), 2));

			//三个坐标系数
			double Cof1 = (tempSat2.POS_X - POSITION_X2) / length_O2 - (RefSat2.POS_X - POSITION_X2) / length_R2;
			double Cof2 = (tempSat2.POS_Y - POSITION_Y2) / length_O2 - (RefSat2.POS_Y - POSITION_Y2) / length_R2;
			double Cof3 = (tempSat2.POS_Z - POSITION_Z2) / length_O2 - (RefSat2.POS_Z - POSITION_Z2) / length_R2;

			//系数矩阵,EWL1
			B_EWL1[tw][0] = Cof1;
			B_EWL1[tw][1] = Cof2;
			B_EWL1[tw][2] = Cof3;
			B_EWL1[tw + num_n][0] = Cof1;
			B_EWL1[tw + num_n][1] = Cof2;
			B_EWL1[tw + num_n][2] = Cof3;
			B_EWL1[tw + num_n][3 + tw] = lambda_longwideC1;

			//系数矩阵,EWL2
			B_EWL2[tw][0] = Cof1;
			B_EWL2[tw][1] = Cof2;
			B_EWL2[tw][2] = Cof3;
			B_EWL2[tw + num_n][0] = Cof1;
			B_EWL2[tw + num_n][1] = Cof2;
			B_EWL2[tw + num_n][2] = Cof3;
			B_EWL2[tw + num_n][3 + tw] = lambda_longwideC2;

			//系数矩阵,WL
			B_WL[tw][0] = Cof1;
			B_WL[tw][1] = Cof2;
			B_WL[tw][2] = Cof3;
			B_WL[tw][3 + tw] = lambda_wideC;

			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;     //站星距
			
			//L矩阵,EWL1
			double deltP_EWL1 = oper_TriP[1](tempSat2, 0, 1, 1) - oper_TriP[1](RefSat2, 0, 1, 1) - (oper_TriP[1](tempSat1, 0, 1, 1) - oper_TriP[1](RefSat1, 0, 1, 1));
			double deltL_EWL1 = oper_TriPhi[1](tempSat2, 0, 1, -1) - oper_TriPhi[1](RefSat2, 0, 1, -1) - (oper_TriPhi[1](tempSat1, 0, 1, -1) - oper_TriPhi[1](RefSat1, 0, 1, -1));
			L_EWL1[tw][0] = deltP_EWL1 - DDgeo;
			L_EWL1[tw + num_n][0] = deltL_EWL1 * lambda_longwideC1 - DDgeo;

			//L矩阵,EWL1
			double deltP_EWL2 = oper_TriP[1](tempSat2, 1, 0, 0) - oper_TriP[1](RefSat2, 1, 0, 0) - (oper_TriP[1](tempSat1, 1, 0, 0) - oper_TriP[1](RefSat1, 1, 0, 0));
			double deltL_EWL2 = oper_TriPhi[1](tempSat2, 1, -5, 4) - oper_TriPhi[1](RefSat2, 1, -5, 4) - (oper_TriPhi[1](tempSat1, 1, -5, 4) - oper_TriPhi[1](RefSat1, 1, -5, 4));
			L_EWL2[tw][0] = deltP_EWL2 - DDgeo;
			L_EWL2[tw + num_n][0] = deltL_EWL2 * lambda_longwideC2 - DDgeo;

			//L矩阵,WL
			double deltL_WL= oper_TriPhi[1](tempSat2, 1, 0, -1) - oper_TriPhi[1](RefSat2, 1, 0, -1) - (oper_TriPhi[1](tempSat1, 1, 0, -1) - oper_TriPhi[1](RefSat1, 1, 0, -1));
			L_WL[tw][0] = deltL_WL * lambda_wideC - DDgeo;

			////判断旧参考星是否被剔除
			//for (int ii = 0; ii < SatNum[0]; ii++)
			//{
			//	if (OldRefPRN[0] == SatPRN[0][ii])
			//	{
			//		Init = true;
			//		break;
			//	}
			//	else
			//		Init = false;
			//}

			//若旧参考星被剔除，需初始化模糊度
			//if (!Init)
			//{
			//	X[3 + tw][0] = N_WC[tp][0];       //初始化模糊度
			//	Q[3 + tw][3 + tw] = 100.0 * 100.0;
			//}
			
			//if (!bInit)
			//{
			//	//初始化模糊度
			//	X_EWL1[3 + tw][0] = (L_EWL1[tw + num_n][0] - L_EWL1[tw][0]) / lambda_longwideC1;       
			//	X_EWL2[3 + tw][0] = (L_EWL2[tw + num_n][0] - L_EWL2[tw][0]) / lambda_longwideC2;
			//	X_WL[3 + tw][0] = 5 * X_EWL1[3 + tw][0] + X_EWL2[3 + tw][0];
			//}
			//else
			//{
			//	for (int ii = 0; ii < NewSatNum[0]; ii++)
			//	{
			//		if (NewSatPRN[0][ii] == PRN)
			//		{
			//			//初始化模糊度
			//			X_EWL1[3 + tw][0] = (L_EWL1[tw + num_n][0] - L_EWL1[tw][0]) / lambda_longwideC1;       
			//			X_EWL2[3 + tw][0] = (L_EWL2[tw + num_n][0] - L_EWL2[tw][0]) / lambda_longwideC2;
			//			X_WL[3 + tw][0] = 5 * X_EWL1[3 + tw][0] + X_EWL2[3 + tw][0];
			//			Q_EWL1[3 + tw][3 + tw] = 100.0 * 100.0;
			//			Q_EWL2[3 + tw][3 + tw] = 100.0 * 100.0;
			//			Q_WL[3 + tw][3 + tw] = 100.0 * 100.0;
			//			break;
			//		}
			//	}
			//}
			tw++;
		}

		tp = 0;           //单系统方程数初始化
		DDerror = 0;
		//BDS-3
		for (k = 0; k < bds3_data.sat_num; k++)
		{
			//bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = bds3_data.sat1[k].numofsat;
			//参考卫星跳过
			if (PRN == BDS3ref)
			{
				continue;
			}

			if (ref_bds3 < 0)
			{
				assert(0);
			}

			Sat tempSat1 = bds3_data.sat1[k];
			Sat RefSat1 = bds3_data.sat1[ref_bds3];
			Sat tempSat2 = bds3_data.sat2[k];
			Sat RefSat2 = bds3_data.sat2[ref_bds3];

			//误差改正
			DDerror = Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2);

			//计算卫星到测站距离（4个距离）
			double length_O1, length_R1, length_O2, length_R2;
			length_O1 = sqrt(pow((tempSat1.POS_X - POSITION_X1), 2) + pow((tempSat1.POS_Y - POSITION_Y1), 2) + pow((tempSat1.POS_Z - POSITION_Z1), 2));
			length_R1 = sqrt(pow((RefSat1.POS_X - POSITION_X1), 2) + pow((RefSat1.POS_Y - POSITION_Y1), 2) + pow((RefSat1.POS_Z - POSITION_Z1), 2));

			length_O2 = sqrt(pow((tempSat2.POS_X - POSITION_X2), 2) + pow((tempSat2.POS_Y - POSITION_Y2), 2) + pow((tempSat2.POS_Z - POSITION_Z2), 2));
			length_R2 = sqrt(pow((RefSat2.POS_X - POSITION_X2), 2) + pow((RefSat2.POS_Y - POSITION_Y2), 2) + pow((RefSat2.POS_Z - POSITION_Z2), 2));

			//三个坐标系数
			double Cof1 = (tempSat2.POS_X - POSITION_X2) / length_O2 - (RefSat2.POS_X - POSITION_X2) / length_R2;
			double Cof2 = (tempSat2.POS_Y - POSITION_Y2) / length_O2 - (RefSat2.POS_Y - POSITION_Y2) / length_R2;
			double Cof3 = (tempSat2.POS_Z - POSITION_Z2) / length_O2 - (RefSat2.POS_Z - POSITION_Z2) / length_R2;

			//系数矩阵,EWL1
			B_EWL1[tw][0] = Cof1;
			B_EWL1[tw][1] = Cof2;
			B_EWL1[tw][2] = Cof3;
			B_EWL1[tw + num_n][0] = Cof1;
			B_EWL1[tw + num_n][1] = Cof2;
			B_EWL1[tw + num_n][2] = Cof3;
			B_EWL1[tw + num_n][3 + tw] = lambda_longwideB1;

			//系数矩阵,EWL2
			B_EWL2[tw][0] = Cof1;
			B_EWL2[tw][1] = Cof2;
			B_EWL2[tw][2] = Cof3;
			B_EWL2[tw + num_n][0] = Cof1;
			B_EWL2[tw + num_n][1] = Cof2;
			B_EWL2[tw + num_n][2] = Cof3;
			B_EWL2[tw + num_n][3 + tw] = lambda_longwideB2;

			//系数矩阵,WL
			B_WL[tw][0] = Cof1;
			B_WL[tw][1] = Cof2;
			B_WL[tw][2] = Cof3;
			B_WL[tw][3 + tw] = lambda_wideB;

			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;        //站星距

			//L矩阵,EWL1
			double deltP_EWL1 = oper_TriP[2](tempSat2, 0, 1, 1) - oper_TriP[2](RefSat2, 0, 1, 1) - (oper_TriP[2](tempSat1, 0, 1, 1) - oper_TriP[2](RefSat1, 0, 1, 1));
			double deltL_EWL1 = oper_TriPhi[2](tempSat2, 0, -1, 1) - oper_TriPhi[2](RefSat2, 0, -1, 1) - (oper_TriPhi[2](tempSat1, 0, -1, 1) - oper_TriPhi[2](RefSat1, 0, -1, 1));
			L_EWL1[tw][0] = deltP_EWL1 - DDgeo;
			L_EWL1[tw + num_n][0] = deltL_EWL1 * lambda_longwideB1 - DDgeo;

			//L矩阵,EWL1
			double deltP_EWL2 = oper_TriP[2](tempSat2, 1, 0, 0) - oper_TriP[2](RefSat2, 1, 0, 0) - (oper_TriP[2](tempSat1, 1, 0, 0) - oper_TriP[2](RefSat1, 1, 0, 0));
			double deltL_EWL2 = oper_TriPhi[2](tempSat2, 1, 3, -4) - oper_TriPhi[2](RefSat2, 1, 3, -4) - (oper_TriPhi[2](tempSat1, 1, 3, -4) - oper_TriPhi[2](RefSat1, 1, 3, -4));
			L_EWL2[tw][0] = deltP_EWL2 - DDgeo;
			L_EWL2[tw + num_n][0] = deltL_EWL2 * lambda_longwideB2 - DDgeo;

			//L矩阵,WL
			double deltL_WL = oper_TriPhi[2](tempSat2, 1, 0, -1) - oper_TriPhi[2](RefSat2, 1, 0, -1) - (oper_TriPhi[2](tempSat1, 1, 0, -1) - oper_TriPhi[2](RefSat1, 1, 0, -1));
			L_WL[tw][0] = deltL_WL * lambda_wideB - DDgeo;

			////判断旧参考星是否被剔除
			//for (int ii = 0; ii < SatNum[1]; ii++)
			//{
			//	if (OldRefPRN[1] == SatPRN[1][ii])
			//	{
			//		Init = true;
			//		break;
			//	}
			//	else
			//		Init = false;
			//}

			//若旧参考星被剔除，需初始化模糊度
			//if (!Init)
			//{
			//	X[3 + tw][0] = N_WB[tp][0];       //初始化模糊度
			//	Q[3 + tw][3 + tw] = 10.0 * 10.0;
			//}

			//if (!bInit)
			//{
			//	//初始化模糊度
			//	X_EWL1[3 + tw][0] = (L_EWL1[tw + num_n][0] - L_EWL1[tw][0]) / lambda_longwideB1;       
			//	X_EWL2[3 + tw][0] = (L_EWL2[tw + num_n][0] - L_EWL2[tw][0]) / lambda_longwideB2;
			//	X_WL[3 + tw][0] = 3 * X_EWL1[3 + tw][0] + X_EWL2[3 + tw][0];
			//}
			//else
			//{
			//	for (int ii = 0; ii < NewSatNum[1]; ii++)
			//	{
			//		if (NewSatPRN[1][ii] == PRN)
			//		{
			//			//初始化模糊度
			//			X_EWL1[3 + tw][0] = (L_EWL1[tw + num_n][0] - L_EWL1[tw][0]) / lambda_longwideB1;       
			//			X_EWL2[3 + tw][0] = (L_EWL2[tw + num_n][0] - L_EWL2[tw][0]) / lambda_longwideB2;
			//			X_WL[3 + tw][0] = 3 * X_EWL1[3 + tw][0] + X_EWL2[3 + tw][0];
			//			Q_EWL1[3 + tw][3 + tw] = 100.0 * 100.0;
			//			Q_EWL2[3 + tw][3 + tw] = 100.0 * 100.0;
			//			Q_WL[3 + tw][3 + tw] = 100.0 * 100.0;
			//			break;
			//		}
			//	}
			//}
			tw++;
		}

		tp = 0;
		if ((bds2_data.sat_num + bds3_data.sat_num) >= 5)
		{
			//权矩阵PG、法方程矩阵NG(逆矩阵即为参数方差)、参数求解(相对概率坐标的改正数)

			//B系数矩阵  L观测值矩阵  R误差矩阵  Qw动态噪声矩阵  Q协方差矩阵  X初始化模糊度

			//kalman滤波
			KalmanFilter(B_EWL1, L_EWL1, R_EWL1, Qw_EWL1, Q_EWL1, X_EWL1);
			KalmanFilter(B_EWL2, L_EWL2, R_EWL2, Qw_EWL2, Q_EWL2, X_EWL2);
			KalmanFilter(B_WL, L_WL, R_WL, Qw_WL, Q_WL, X_WL);

			/*B_EWL1.MyTRACE();
			cout << endl;*/

			double ratio_EWL1 = 1.0;
			double ratio_EWL2 = 1.0;
			double ratio_WL = 1.0;
			CMatrix aXb_EWL1(3, 1);
			CMatrix aXb_EWL2(3, 1);
			CMatrix aXb_WL(3, 1);
			CMatrix fix_EWL1(num_n, 1);
			CMatrix fix_EWL2(num_n, 1);
			CMatrix fix_WL(num_n, 1);

			//lambda算法固定模糊度
			Resamb_LAMBDA(Q_EWL1, X_EWL1, aXb_EWL1, ratio_EWL1, fix_EWL1);
			Resamb_LAMBDA(Q_EWL2, X_EWL2, aXb_EWL2, ratio_EWL2, fix_EWL2);
			Resamb_LAMBDA(Q_WL, X_WL, aXb_WL, ratio_WL, fix_WL);
			X_WL.MyTRACE();
			cout << endl;

			tp = 3;
			for (int i1 = 0; i1 < bds2_data.sat_num; tp++, i1++)
			{
				if (abs(X_WL[tp][0] - 5 * X_EWL1[tp][0] - X_EWL2[tp][0]) > 0.5 && ratio_WL < 2 && ratio_EWL1>100 && ratio_EWL2 > 100)
				{
					X_WL[tp][0] = 5 * X_EWL1[tp][0] + X_EWL2[tp][0];
				}
			}

			for (int i2 = 0; i2 < bds3_data.sat_num; tp++, i2++)
			{
				if (abs(X_WL[tp][0] - 3 * X_EWL1[tp][0] - X_EWL2[tp][0]) > 0.5 && ratio_WL < 2 && ratio_EWL1>100 && ratio_EWL2 > 100)
				{
					X_WL[tp][0] = 3 * X_EWL1[tp][0] + X_EWL2[tp][0];
				}
			}

			//恢复出估计的坐标值
			CMatrix X_POS(3, 1);
			X_POS[0][0] = POSITION_X2 - aXb_WL[0][0];
			X_POS[1][0] = POSITION_Y2 - aXb_WL[1][0];
			X_POS[2][0] = POSITION_Z2 - aXb_WL[2][0];

			//真实坐标写成矩阵
			CMatrix XT(3, 1);
			XT[0][0] = P_X2;
			XT[1][0] = P_Y2;
			XT[2][0] = P_Z2;

			//与准确值作差得到XYZ坐标系下的误差
			CMatrix ErrorX = X_POS - XT;

			//将误差投影到NEU东北天坐标系
			CMatrix ErrorX_NEU = TT * ErrorX;

			fprintf(fp, "EpNum: %2d, %2d, %2d, %10.6f, %2d, %2d,%15.8g,%15.8g, %15.8g, %15.8g\n", j + 1, tempData2.hour, tempData2.minute, tempData2.second, bds2_data.sat_num, bds3_data.sat_num, ratio_WL,
				ErrorX_NEU[0][0], ErrorX_NEU[1][0], ErrorX_NEU[2][0]);

			for (int idx = 0; idx < 2; idx++)
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
			}
			bInit = true;  //初始化成功
			Init = true;
		}
	}
}

//函数功能，TCAR法固定模糊度,采用几何相关模型,流动站
void TCAR_move(string& filename1, string& filename2)
{
	//默认第一个站为基准站，第二个为流动站
	observe_spp sppfile1;
	observe_spp sppfile2;
	ReadSPP_KIN(filename1, sppfile1);
	ReadSPP_KIN(filename2, sppfile2);
	bool bInit = false; //初始化
	bool Init1 = true;  //单系统初始化
	bool Init2 = true;

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

	//输出差分定位结果用
	string name = filename2.substr(0, 3) + "_TriRTK.txt";
	string error = filename2.substr(0, 3) + "_TriDebug.txt";
	const char* output_filename = name.c_str();
	const char* out_file = error.c_str();
	FILE* fp;
	FILE* fp1;
	fopen_s(&fp, output_filename, "w");
	fopen_s(&fp1, out_file, "w");

	int OldSatPRN[2][30] = { {0} };
	int OldRefPRN[2] = { -1,-1 };
	int OldSatNum[2] = { 0 };

	CMatrix Q_EWL1, Q_EWL2, Q_WL, X_EWL1, X_EWL2, X_WL;
	map<string, double> preN_WL, preN_EWL1, preN_EWL2;

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
			TriOnCycleDetect(tempData2, preData2);     //判断此此历元与上一个历元中的卫星载波测量值是否有周跳,如果有周跳就将此历元中的卫星记为不可用
			//OnCycleDetect(tempData2, preData2);
		}

		//基准站数据相应时刻必须在前
		int j = 0;
		for (j = posk; j < sppfile1.liyuan_num; j++)
		{
			if (sppfile1.epoch[j].GPSTIME <= tempData2.GPSTIME)
				// ReSharper disable once CppRedundantControlFlowJump
				continue;
			break;
		}

		//回溯一个历元：posk即为要找的最近的历元且满足posk的时刻小于tempData2的时间
		if (j > 0) posk = --j;
		else
		{
			posk = j;
			continue;
		}

		//判断时间差：两个是否相隔超过历元间隔时间
		if (fabs(sppfile1.epoch[j].GPSTIME - tempData2.GPSTIME) >= sppfile1.INTERVAL)
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
			TriOnCycleDetect(tempData1, preData1);
			//OnCycleDetect(tempData1, preData1);
		}

		//差分定位共视处理
		EpochCommonPro(tempData1, tempData2, tempData);

		//分离BDS-2、BDS-3卫星，便于后续单系统、多系统处理
		Obs_epoch bds2_data;
		Obs_epoch bds3_data;

		int SatPRN[2][30] = { {0} };		//当前历元卫星号 依次存储
		int SatNum[2] = { 0 };	//当前历元各系统卫星数目
		int RefPRN[2] = { 0 };  //当前历元各系统参考卫星号

		int NewSatNum[2] = { 0 };
		int NewSatPRN[2][30] = { {0} };

		int k = 0;
		int m = 0;
		int bds2num, bds3num;
		bds2num = bds3num = 0;
		for (m = 0; m < tempData.sat_num; m++)
		{
			if (tempData.sat1[m].sattype == "C")
			{
				bds2_data.sat1.push_back(tempData.sat1[m]);
				bds2_data.sat2.push_back(tempData.sat2[m]);
				SatPRN[0][bds2num] = tempData.sat1[m].numofsat;
				bds2num++;
			}

			if (tempData.sat1[m].sattype == "B")
			{
				bds3_data.sat1.push_back(tempData.sat1[m]);
				bds3_data.sat2.push_back(tempData.sat2[m]);
				SatPRN[1][bds3num] = tempData.sat1[m].numofsat;
				bds3num++;
			}
		}

		//各类型卫星数
		SatNum[0] = bds2num;
		SatNum[1] = bds3num;

		bds2_data.sat_num = bds2_data.sat1.size();
		bds3_data.sat_num = bds3_data.sat1.size();
		tempData.GPSTIME = tempData2.GPSTIME;

		//记录新卫星
		if (bInit == true)
		{
			for (int idx = 0; idx < 2; idx++)
			{
				set<int>s;
				for (int ii = 0; ii < OldSatNum[idx]; ii++)
					s.insert(OldSatPRN[idx][ii]);
				for (int ii = 0; ii < SatNum[idx]; ii++)
				{
					if (s.find(SatPRN[idx][ii]) == s.end())
					{
						NewSatNum[idx]++;
						for (int iii = 0; iii < NewSatNum[idx]; iii++)
							NewSatPRN[idx][iii] = SatPRN[idx][ii];
					}
				}
			}
		}

		//卫星数至少要4颗
		if (bds2_data.sat_num + bds3_data.sat_num < 4)
		{
			printf("卫星数少于4颗\n");
			bInit = false;
			continue;
		}

		//各系统至少要2颗卫星,区分参考星和普通星
		if (bds2_data.sat_num < 2)
		{
			printf("BDS2卫星数少于2颗\n");
			Init1 = false;
			continue;
		}
		if (bds3_data.sat_num < 2)
		{
			printf("BDS3卫星数少于2颗\n");
			Init2 = false;
			continue;
		}

		//采用站间星间二次差分模型，选择参考卫星
		int BDS2ref = -1;
		int BDS3ref = -2;

		BDS2ref = OldRefPRN[0];
		BDS3ref = OldRefPRN[1];

		//选择最大高度角卫星作为参考卫星，不同类型卫星分开看
		SelectRef(bds2_data, BDS2ref);
		SelectRef(bds3_data, BDS3ref);

		RefPRN[0] = BDS2ref; //BDS-2参考卫星序号
		RefPRN[1] = BDS3ref; //BDS-3参考卫星序号

		//当前历元BDS-2、BDS-3卫星序列
		int PRN_X_B2[20];
		int PRN_X_B3[30];

		//数列初始化
		memset(PRN_X_B2, 0, sizeof(PRN_X_B2));
		memset(PRN_X_B3, 0, sizeof(PRN_X_B3));

		//储存历元中的卫星序列
		for (k = 0; k < bds2_data.sat_num; k++)
			PRN_X_B2[k] = bds2_data.sat1[k].numofsat;
		for (k = 0; k < bds3_data.sat_num; k++)
			PRN_X_B3[k] = bds3_data.sat1[k].numofsat;

		//各卫星个数
		int num_X_B2 = bds2_data.sat_num;
		int num_X_B3 = bds3_data.sat_num;

		//各系统双差方程个数
		int numbds2_x = num_X_B2 - 1;
		int numbds3_x = num_X_B3 - 1;

		//多系统融合双差方程个数（载波观测方程+伪距观测方程)，松组合
		int num_x_p = (numbds2_x + numbds3_x) * 2;

		//待估参数个数：ΔX，ΔY，ΔZ，每颗卫星的模糊度
		int num_x = numbds2_x + numbds3_x + 3;

		//模糊度个数
		int num_n = numbds2_x + numbds3_x;

		//定位解算相关的矩阵
		//多系统
		CMatrix B_EWL1(num_x_p, num_x);           //系数矩阵,EWL1
		CMatrix B_EWL2(num_x_p, num_x);           //系数矩阵,EWL2
		CMatrix B_WL(num_n, num_x);               //系数矩阵（只有载波观测值）,WL

		CMatrix L_EWL1(num_x_p, 1);               //观测值矩阵,EWL1
		CMatrix L_EWL2(num_x_p, 1);               //观测值矩阵,EWL2
		CMatrix L_WL(num_n, 1);                   //观测值矩阵（只有载波观测值）,WL

		CMatrix R_EWL1(num_x_p, num_x_p);         //误差矩阵（逆阵即为权矩阵）
		CMatrix R_EWL2(num_x_p, num_x_p);         //误差矩阵（逆阵即为权矩阵）
		CMatrix R_WL(num_n, num_n);               //误差矩阵（逆阵即为权矩阵）

		CMatrix Qw_EWL1(num_x, num_x);            //动态噪声阵,EWL1
		CMatrix Qw_EWL2(num_x, num_x);            //动态噪声阵,EWL2
		CMatrix Qw_WL(num_x, num_x);              //动态噪声阵,WL

		//单BDS2
		CMatrix RC2_P(numbds2_x, numbds2_x);      //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RC2_L(numbds2_x, numbds2_x);      //误差矩阵（逆阵即为权矩阵） 载波

		//单BDS3
		CMatrix RC3_P(numbds3_x, numbds3_x);      //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RC3_L(numbds3_x, numbds3_x);      //误差矩阵（逆阵即为权矩阵） 载波

		//位置的动态噪声阵，如果没有周跳，模糊度不变，所以模糊度的噪声为0，只有位置噪声
		for (int idx = 0; idx < 3; idx++) {
			Qw_EWL1[idx][idx] = 50.0 * 50.0;
			Qw_EWL2[idx][idx] = 50.0 * 50.0;
			Qw_WL[idx][idx] = 10.0 * 10.0;
		}

		//更新BDS2和BDS2的协方差矩阵Q和结果矩阵X
		if (!bInit)
		{
			Q_EWL1.SetSize(num_x, num_x);
			Q_EWL2.SetSize(num_x, num_x);
			Q_WL.SetSize(num_x, num_x);
			X_EWL1.SetSize(num_x, 1);
			X_EWL2.SetSize(num_x, 1);
			X_WL.SetSize(num_x, 1);
			InitState_long(Q_EWL1);
			InitState_long(Q_EWL2);
			InitState(Q_WL);
		}
		else
		{
			//检测新卫星
			CheckNewSats(SatNum, OldSatNum, SatPRN, OldSatPRN, NewSatPRN, NewSatNum);

			//更新协方差和X矩阵
			UpdataState(Q_EWL1, X_EWL1, SatNum, OldSatNum, SatPRN, OldSatPRN, RefPRN, OldRefPRN);
			UpdataState(Q_EWL2, X_EWL2, SatNum, OldSatNum, SatPRN, OldSatPRN, RefPRN, OldRefPRN);
			UpdataState(Q_WL, X_WL, SatNum, OldSatNum, SatPRN, OldSatPRN, RefPRN, OldRefPRN);
		}

		//各类型卫星序号存储容器
		vector<int> BDS2_PRN;
		vector<int> BDS3_PRN;

		for (k = 0; k < bds2_data.sat_num; k++)
			BDS2_PRN.push_back(PRN_X_B2[k]);
		for (k = 0; k < bds3_data.sat_num; k++)
			BDS3_PRN.push_back(PRN_X_B3[k]);

		//伪距
		GetDDR_Peo(bds2_data, BDS2_PRN, BDS2ref, RC2_P);
		GetDDR_Peo(bds3_data, BDS3_PRN, BDS3ref, RC3_P);
		//载波
		GetDDR_Leo(bds2_data, BDS2_PRN, BDS2ref, RC2_L);
		GetDDR_Leo(bds3_data, BDS3_PRN, BDS3ref, RC3_L);

		//融合系统的R矩阵，分块对角矩阵
		//EWL1
		//伪距
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R_EWL1[k1][k2] = 2 * RC2_P[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R_EWL1[k1 + numbds2_x][k2 + numbds2_x] = 2 * RC3_P[k1][k2];
			}
		//载波
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R_EWL1[k1 + numbds2_x + numbds3_x][k2 + numbds2_x + numbds3_x] = 2 * RC2_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R_EWL1[k1 + numbds2_x * 2 + numbds3_x][k2 + numbds2_x * 2 + numbds3_x] = 2 * RC3_L[k1][k2];
			}

		//EWL2
		//伪距
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R_EWL2[k1][k2] = 3 * RC2_P[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R_EWL2[k1 + numbds2_x][k2 + numbds2_x] = RC3_P[k1][k2];
			}
		//载波
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R_EWL2[k1 + numbds2_x + numbds3_x][k2 + numbds2_x + numbds3_x] = 42 * RC2_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R_EWL2[k1 + numbds2_x * 2 + numbds3_x][k2 + numbds2_x * 2 + numbds3_x] = 26 * RC3_L[k1][k2];
			}

		//WL
		//宽巷
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R_WL[k1][k2] = 2 * RC2_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R_WL[k1 + numbds2_x][k2 + numbds2_x] = 2 * RC3_L[k1][k2];
			}

		//初始坐标信息：基站采用spp文件头中坐标（固定），流动站采用每个历元的概率坐标（考虑实际应用中其可能为动态情况）
		//本组数据中：流动站为动态，因此需要更新其坐标
		//基站坐标：头文件坐标
		double POSITION_X1 = sppfile1.APP_X;
		double POSITION_Y1 = sppfile1.APP_Y;
		double POSITION_Z1 = sppfile1.APP_Z;
		/*double POSITION_X1 = tempData1.posX;
		double POSITION_Y1 = tempData1.posY;
		double POSITION_Z1 = tempData1.posZ;*/

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
		int ref_bds2 = -1;
		int ref_bds3 = -1;

		for (k = 0; k < bds2_data.sat_num; k++)
		{
			if (bds2_data.sat1[k].numofsat == BDS2ref)
			{
				ref_bds2 = k;
				break;
			}
		}

		for (k = 0; k < bds3_data.sat_num; k++)
		{
			if (bds3_data.sat1[k].numofsat == BDS3ref)
			{
				ref_bds3 = k;
				break;
			}
		}

		int tp = 0;  //记录单系统方程数
		int tw = 0;  //记录多系统方程数
		double DDerror = 0;  //初始化误差改正项
		//BDS-2
		for (k = 0; k < bds2_data.sat_num; k++)
		{
			//bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = bds2_data.sat1[k].numofsat;
			string nameofsat = "C" + to_string(PRN) + "-" + to_string(BDS2ref);
			//参考卫星跳过
			if (PRN == BDS2ref)
			{
				continue;
			}

			if (ref_bds2 < 0)
			{
				assert(0);
			}

			Sat tempSat1 = bds2_data.sat1[k];
			Sat RefSat1 = bds2_data.sat1[ref_bds2];
			Sat tempSat2 = bds2_data.sat2[k];
			Sat RefSat2 = bds2_data.sat2[ref_bds2];

			//误差改正
			DDerror = Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2);

			//计算卫星到测站距离（4个距离）
			double length_O1, length_R1, length_O2, length_R2;
			length_O1 = sqrt(pow((tempSat1.POS_X - POSITION_X1), 2) + pow((tempSat1.POS_Y - POSITION_Y1), 2) + pow((tempSat1.POS_Z - POSITION_Z1), 2));
			length_R1 = sqrt(pow((RefSat1.POS_X - POSITION_X1), 2) + pow((RefSat1.POS_Y - POSITION_Y1), 2) + pow((RefSat1.POS_Z - POSITION_Z1), 2));

			length_O2 = sqrt(pow((tempSat2.POS_X - POSITION_X2), 2) + pow((tempSat2.POS_Y - POSITION_Y2), 2) + pow((tempSat2.POS_Z - POSITION_Z2), 2));
			length_R2 = sqrt(pow((RefSat2.POS_X - POSITION_X2), 2) + pow((RefSat2.POS_Y - POSITION_Y2), 2) + pow((RefSat2.POS_Z - POSITION_Z2), 2));

			//三个坐标系数
			double Cof1 = (tempSat2.POS_X - POSITION_X2) / length_O2 - (RefSat2.POS_X - POSITION_X2) / length_R2;
			double Cof2 = (tempSat2.POS_Y - POSITION_Y2) / length_O2 - (RefSat2.POS_Y - POSITION_Y2) / length_R2;
			double Cof3 = (tempSat2.POS_Z - POSITION_Z2) / length_O2 - (RefSat2.POS_Z - POSITION_Z2) / length_R2;

			//系数矩阵,EWL1
			B_EWL1[tw][0] = Cof1;
			B_EWL1[tw][1] = Cof2;
			B_EWL1[tw][2] = Cof3;
			B_EWL1[tw + num_n][0] = Cof1;
			B_EWL1[tw + num_n][1] = Cof2;
			B_EWL1[tw + num_n][2] = Cof3;
			B_EWL1[tw + num_n][3 + tw] = lambda_longwideC1;

			//系数矩阵,EWL2
			B_EWL2[tw][0] = Cof1;
			B_EWL2[tw][1] = Cof2;
			B_EWL2[tw][2] = Cof3;
			B_EWL2[tw + num_n][0] = Cof1;
			B_EWL2[tw + num_n][1] = Cof2;
			B_EWL2[tw + num_n][2] = Cof3;
			B_EWL2[tw + num_n][3 + tw] = lambda_longwideC2;

			//系数矩阵,WL
			B_WL[tw][0] = Cof1;
			B_WL[tw][1] = Cof2;
			B_WL[tw][2] = Cof3;
			B_WL[tw][3 + tw] = lambda_wideC;

			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;     //站星距
			
			//L矩阵,EWL1
			double deltP_EWL1 = oper_TriP[1](tempSat2, 0, 1, 1) - oper_TriP[1](RefSat2, 0, 1, 1) - (oper_TriP[1](tempSat1, 0, 1, 1) - oper_TriP[1](RefSat1, 0, 1, 1));
			double deltL_EWL1 = oper_TriPhi[1](tempSat2, 0, -1, 1) - oper_TriPhi[1](RefSat2, 0, -1, 1) - (oper_TriPhi[1](tempSat1, 0, -1, 1) - oper_TriPhi[1](RefSat1, 0, -1, 1));
			L_EWL1[tw][0] = deltP_EWL1 - DDgeo;
			L_EWL1[tw + num_n][0] = deltL_EWL1 * lambda_longwideC1 - DDgeo;

			//L矩阵,EWL2
			double deltP_EWL2 = oper_TriP[1](tempSat2, 1, 1, 1) - oper_TriP[1](RefSat2, 1, 1, 1) - (oper_TriP[1](tempSat1, 1, 1, 1) - oper_TriP[1](RefSat1, 1, 1, 1));
			double deltL_EWL2 = oper_TriPhi[1](tempSat2, 1, 4, -5) - oper_TriPhi[1](RefSat2, 1, 4, -5) - (oper_TriPhi[1](tempSat1, 1, 4, -5) - oper_TriPhi[1](RefSat1, 1, 4, -5));
			L_EWL2[tw][0] = deltP_EWL2 - DDgeo;
			L_EWL2[tw + num_n][0] = deltL_EWL2 * lambda_longwideC2 - DDgeo;

			//L矩阵,WL
			double deltL_WL = oper_TriPhi[1](tempSat2, 1, -1, 0) - oper_TriPhi[1](RefSat2, 1, -1, 0) - (oper_TriPhi[1](tempSat1, 1, -1, 0) - oper_TriPhi[1](RefSat1, 1, -1, 0));
			L_WL[tw][0] = deltL_WL * lambda_wideC - DDgeo;
			
			//判断旧参考星是否被剔除
			/*for (int ii = 0; ii < SatNum[0]; ii++)
			{
				if (OldRefPRN[0] == SatPRN[0][ii])
				{
					Init = true;
					break;
				}
				else
					Init = false;
			}*/

			/*if (i >= 1)
			{
				if (fabs(X_WL[tw + 3][0] - preN_WL[nameofsat]) > 0.2)
				{
					double var = (X_WL[3 + tw][0] + preN_WL[nameofsat]) / 2;
					X_WL[3 + tw][0] = var;
				}
			}*/

			//初始化模糊度
			if (fabs(X_EWL1[3 + tw][0] - 0.0) < 1.0 || !Init1)
			{
				X_EWL1[3 + tw][0] = (L_EWL1[tw + num_n][0] - L_EWL1[tw][0]) / lambda_longwideC1;
				X_EWL2[3 + tw][0] = (L_EWL2[tw + num_n][0] - L_EWL2[tw][0]) / lambda_longwideC2;
				X_WL[3 + tw][0] = 5 * X_EWL1[3 + tw][0] + X_EWL2[3 + tw][0];
				Q_EWL1[3 + tw][3 + tw] = 10.0 * 10.0;
				Q_EWL2[3 + tw][3 + tw] = 10.0 * 10.0;
				Q_WL[3 + tw][3 + tw] = 30.0 * 30.0;
			}

			if (!bInit)
			{
				//初始化模糊度
				X_EWL1[3 + tw][0] = (L_EWL1[tw + num_n][0] - L_EWL1[tw][0]) / lambda_longwideC1;
				X_EWL2[3 + tw][0] = (L_EWL2[tw + num_n][0] - L_EWL2[tw][0]) / lambda_longwideC2;
				X_WL[3 + tw][0] = 5 * X_EWL1[3 + tw][0] + X_EWL2[3 + tw][0];
			}
			else
			{
				for (int ii = 0; ii < NewSatNum[0]; ii++)
				{
					if (NewSatPRN[0][ii] == PRN)
					{
						//初始化模糊度
						X_EWL1[3 + tw][0] = (L_EWL1[tw + num_n][0] - L_EWL1[tw][0]) / lambda_longwideC1;
						X_EWL2[3 + tw][0] = (L_EWL2[tw + num_n][0] - L_EWL2[tw][0]) / lambda_longwideC2;
						X_WL[3 + tw][0] = 5 * X_EWL1[3 + tw][0] + X_EWL2[3 + tw][0];
						Q_EWL1[3 + tw][3 + tw] = 5.0 * 5.0;
						Q_EWL2[3 + tw][3 + tw] = 5.0 * 5.0;
						Q_WL[3 + tw][3 + tw] = 30.0 * 30.0;
						break;
					}
				}
			}
			tw++;
		}

		tp = 0;           //单系统方程数初始化
		DDerror = 0;
		//BDS-3
		for (k = 0; k < bds3_data.sat_num; k++)
		{
			//bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = bds3_data.sat1[k].numofsat;
			string nameofsat = "B" + to_string(PRN) + "-" + to_string(BDS3ref);
			//参考卫星跳过
			if (PRN == BDS3ref)
			{
				continue;
			}

			if (ref_bds3 < 0)
			{
				assert(0);
			}

			Sat tempSat1 = bds3_data.sat1[k];
			Sat RefSat1 = bds3_data.sat1[ref_bds3];
			Sat tempSat2 = bds3_data.sat2[k];
			Sat RefSat2 = bds3_data.sat2[ref_bds3];

			//误差改正
			DDerror = Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2);

			//计算卫星到测站距离（4个距离）
			double length_O1, length_R1, length_O2, length_R2;
			length_O1 = sqrt(pow((tempSat1.POS_X - POSITION_X1), 2) + pow((tempSat1.POS_Y - POSITION_Y1), 2) + pow((tempSat1.POS_Z - POSITION_Z1), 2));
			length_R1 = sqrt(pow((RefSat1.POS_X - POSITION_X1), 2) + pow((RefSat1.POS_Y - POSITION_Y1), 2) + pow((RefSat1.POS_Z - POSITION_Z1), 2));

			length_O2 = sqrt(pow((tempSat2.POS_X - POSITION_X2), 2) + pow((tempSat2.POS_Y - POSITION_Y2), 2) + pow((tempSat2.POS_Z - POSITION_Z2), 2));
			length_R2 = sqrt(pow((RefSat2.POS_X - POSITION_X2), 2) + pow((RefSat2.POS_Y - POSITION_Y2), 2) + pow((RefSat2.POS_Z - POSITION_Z2), 2));

			//三个坐标系数
			double Cof1 = (tempSat2.POS_X - POSITION_X2) / length_O2 - (RefSat2.POS_X - POSITION_X2) / length_R2;
			double Cof2 = (tempSat2.POS_Y - POSITION_Y2) / length_O2 - (RefSat2.POS_Y - POSITION_Y2) / length_R2;
			double Cof3 = (tempSat2.POS_Z - POSITION_Z2) / length_O2 - (RefSat2.POS_Z - POSITION_Z2) / length_R2;

			//系数矩阵,EWL1
			B_EWL1[tw][0] = Cof1;
			B_EWL1[tw][1] = Cof2;
			B_EWL1[tw][2] = Cof3;
			B_EWL1[tw + num_n][0] = Cof1;
			B_EWL1[tw + num_n][1] = Cof2;
			B_EWL1[tw + num_n][2] = Cof3;
			B_EWL1[tw + num_n][3 + tw] = lambda_longwideB1;

			//系数矩阵,EWL2
			B_EWL2[tw][0] = Cof1;
			B_EWL2[tw][1] = Cof2;
			B_EWL2[tw][2] = Cof3;
			B_EWL2[tw + num_n][0] = Cof1;
			B_EWL2[tw + num_n][1] = Cof2;
			B_EWL2[tw + num_n][2] = Cof3;
			B_EWL2[tw + num_n][3 + tw] = lambda_longwideB2;

			//系数矩阵,WL
			B_WL[tw][0] = Cof1;
			B_WL[tw][1] = Cof2;
			B_WL[tw][2] = Cof3;
			B_WL[tw][3 + tw] = lambda_wideB;

			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;        //站星距

			//L矩阵,EWL1
			double deltP_EWL1 = oper_TriP[2](tempSat2, 0, 1, 1) - oper_TriP[2](RefSat2, 0, 1, 1) - (oper_TriP[2](tempSat1, 0, 1, 1) - oper_TriP[2](RefSat1, 0, 1, 1));
			double deltL_EWL1 = oper_TriPhi[2](tempSat2, 0, -1, 1) - oper_TriPhi[2](RefSat2, 0, -1, 1) - (oper_TriPhi[2](tempSat1, 0, -1, 1) - oper_TriPhi[2](RefSat1, 0, -1, 1));
			L_EWL1[tw][0] = deltP_EWL1 - DDgeo;
			L_EWL1[tw + num_n][0] = deltL_EWL1 * lambda_longwideB1 - DDgeo;

			//L矩阵,EWL2
			double deltP_EWL2 = oper_TriP[2](tempSat2, 1, 0, 0) - oper_TriP[2](RefSat2, 1, 0, 0) - (oper_TriP[2](tempSat1, 1, 0, 0) - oper_TriP[2](RefSat1, 1, 0, 0));
			double deltL_EWL2 = oper_TriPhi[2](tempSat2, 1, 3, -4) - oper_TriPhi[2](RefSat2, 1, 3, -4) - (oper_TriPhi[2](tempSat1, 1, 3, -4) - oper_TriPhi[2](RefSat1, 1, 3, -4));
			L_EWL2[tw][0] = deltP_EWL2 - DDgeo;
			L_EWL2[tw + num_n][0] = deltL_EWL2 * lambda_longwideB2 - DDgeo;

			//L矩阵,WL
			double deltL_WL = oper_TriPhi[2](tempSat2, 1, 0, -1) - oper_TriPhi[2](RefSat2, 1, 0, -1) - (oper_TriPhi[2](tempSat1, 1, 0, -1) - oper_TriPhi[2](RefSat1, 1, 0, -1));
			L_WL[tw][0] = deltL_WL * lambda_wideB - DDgeo;

			//判断旧参考星是否被剔除
			/*for (int ii = 0; ii < SatNum[1]; ii++)
			{
				if (OldRefPRN[1] == SatPRN[1][ii])
				{
					Init1 = true;
					break;
				}
				else
					Init1 = false;
			}*/

			/*if (i >= 1)
			{
				if (fabs(X_WL[tw + 3][0] - preN_WL[nameofsat]) > 0.2)
				{
					double var = (X_WL[3 + tw][0] + preN_WL[nameofsat]) / 2;
					X_WL[3 + tw][0] = var;
				}
			}*/

			//初始化模糊度
			if (fabs(X_EWL1[3 + tw][0] - 0.0) < 1.0 || !Init2)
			{
				X_EWL1[3 + tw][0] = (L_EWL1[tw + num_n][0] - L_EWL1[tw][0]) / lambda_longwideB1;
				X_EWL2[3 + tw][0] = (L_EWL2[tw + num_n][0] - L_EWL2[tw][0]) / lambda_longwideB2;
				X_WL[3 + tw][0] = 3 * X_EWL1[3 + tw][0] + X_EWL2[3 + tw][0];
				Q_EWL1[3 + tw][3 + tw] = 10.0 * 10.0;
				Q_EWL2[3 + tw][3 + tw] = 10.0 * 10.0;
				Q_WL[3 + tw][3 + tw] = 30.0 * 30.0;
			}
			

			if (!bInit)
			{
				//初始化模糊度
				X_EWL1[3 + tw][0] = (L_EWL1[tw + num_n][0] - L_EWL1[tw][0]) / lambda_longwideB1;
				X_EWL2[3 + tw][0] = (L_EWL2[tw + num_n][0] - L_EWL2[tw][0]) / lambda_longwideB2;
				X_WL[3 + tw][0] = 3 * X_EWL1[3 + tw][0] + X_EWL2[3 + tw][0];
			}
			else
			{
				for (int ii = 0; ii < NewSatNum[1]; ii++)
				{
					if (NewSatPRN[1][ii] == PRN)
					{
						//初始化模糊度
						X_EWL1[3 + tw][0] = (L_EWL1[tw + num_n][0] - L_EWL1[tw][0]) / lambda_longwideB1;
						X_EWL2[3 + tw][0] = (L_EWL2[tw + num_n][0] - L_EWL2[tw][0]) / lambda_longwideB2;
						X_WL[3 + tw][0] = 3 * X_EWL1[3 + tw][0] + X_EWL2[3 + tw][0];
						Q_EWL1[3 + tw][3 + tw] = 5.0 * 5.0;
						Q_EWL2[3 + tw][3 + tw] = 5.0 * 5.0;
						Q_WL[3 + tw][3 + tw] = 30.0 * 30.0;
						break;
					}
				}
			}
			tw++;
		}


		if (bds2_data.sat_num + bds3_data.sat_num >= 4)
		{
			//B系数矩阵  L观测值矩阵  R误差矩阵  Qw动态噪声矩阵  Q协方差矩阵  X位置参数和模糊度

			double ratio_EWL1 = 1.0;
			double ratio_EWL2 = 1.0;
			double ratio_WL = 1.0;
			CMatrix aXb_EWL1(3, 1);
			CMatrix aXb_EWL2(3, 1);
			CMatrix aXb_WL(3, 1);
			CMatrix fix_EWL1(num_n, 1);
			CMatrix fix_EWL2(num_n, 1);
			CMatrix fix_WL(num_n, 1);

			//kalman滤波
			KalmanFilter(B_EWL1, L_EWL1, R_EWL1, Qw_EWL1, Q_EWL1, X_EWL1);
			KalmanFilter(B_EWL2, L_EWL2, R_EWL2, Qw_EWL2, Q_EWL2, X_EWL2);
			KalmanFilter(B_WL, L_WL, R_WL, Qw_WL, Q_WL, X_WL);

			//lambda算法固定模糊度
			Resamb_LAMBDA(Q_EWL1, X_EWL1, aXb_EWL1, ratio_EWL1, fix_EWL1);
			Resamb_LAMBDA(Q_EWL2, X_EWL2, aXb_EWL2, ratio_EWL2, fix_EWL2);
			Resamb_LAMBDA(Q_WL, X_WL, aXb_WL, ratio_WL, fix_WL);

			/*for (int var1 = 0; var1 < bds2num + bds3num - 2; var1++)
			{
				if (fabs(X_EWL1[var1 + 3][0] - fix_EWL1[var1][0]) > 0.2)
				{
					X_EWL1[var1 + 3][0] = fix_EWL1[var1][0];
				}
			}

			for (int var2 = 0; var2 < bds2num + bds3num - 2; var2++)
			{
				if (fabs(X_EWL2[var2 + 3][0] - fix_EWL2[var2][0]) > 0.2)
				{
					X_EWL2[var2 + 3][0] = fix_EWL2[var2][0];
				}
			}*/

			int r = X_EWL1.Row;
			CMatrix tran1(r, r);
			CMatrix tran2(r, r);
			for (int i1 = 0; i1 < bds2num - 1; i1++)
				tran1[i1 + 3][i1 + 3] = 5;
			for (int i2 = 0; i2 < bds3num - 1; i2++)
				tran1[i2 + 2 + bds2num][i2 + 2 + bds2num] = 3;
			for (int i3 = 0; i3 < bds2num + bds3num - 2; i3++)
				tran2[i3 + 3][i3 + 3] = 1;
			double x = X_WL[0][0]; double y = X_WL[1][0]; double z = X_WL[2][0];

			if (ratio_WL < 2 /*&& ratio_EWL1 > 2 && ratio_EWL2 > 2*/)
			{
				/*for (int ii = 0; ii < bds2num - 1; ii++)
					X_WL[3 + ii][0] = 5 * fix_EWL1[ii][0] + fix_EWL2[ii][0];
				for (int ij = 0; ij < bds3num - 1; ij++)
					X_WL[2 + bds2num + ij][0] = 3 * fix_EWL1[ij + bds2num - 1][0] + fix_EWL2[ij + bds2num - 1][0];*/
				X_WL = tran1 * X_EWL1 + tran2 * X_EWL2;
				X_WL[0][0] = x;
				X_WL[1][0] = y;
				X_WL[2][0] = z;
				Q_WL = tran1 * Q_EWL1 * tran1.T() + tran2 * Q_EWL2 * tran2.T();
				KalmanFilter(B_WL, L_WL, R_WL, Qw_WL, Q_WL, X_WL);
				Resamb_LAMBDA(Q_WL, X_WL, aXb_WL, ratio_WL, fix_WL);
			}

			/*X_EWL1.MyTRACE();
			cout << ratio_EWL1 << endl;*/
			/*Q_EWL1.MyTRACE();
			cout << endl;
			fix_EWL1.MyTRACE();
			cout << endl;*/
			/*Q_EWL2.MyTRACE();
			cout << endl;*/
			X_EWL2.MyTRACE();
			cout << ratio_EWL2 << endl;
			fix_EWL2.MyTRACE();
			cout << endl;
			//L_EWL2.MyTRACE();
			//cout << endl;
			/*X_WL.MyTRACE();
			cout << ratio_WL << endl;*/

			//恢复出估计的坐标值
			CMatrix X_POS(3, 1);
			if (ratio_WL > 2)
			{
				X_POS[0][0] = POSITION_X2 - aXb_WL[0][0];
				X_POS[1][0] = POSITION_Y2 - aXb_WL[1][0];
				X_POS[2][0] = POSITION_Z2 - aXb_WL[2][0];
			}
			else
			{
				X_POS[0][0] = POSITION_X2 - X_WL[0][0];
				X_POS[1][0] = POSITION_Y2 - X_WL[1][0];
				X_POS[2][0] = POSITION_Z2 - X_WL[2][0];
			}

			//真实坐标写成矩阵
			CMatrix XT(3, 1);
			XT[0][0] = P_X2;
			XT[1][0] = P_Y2;
			XT[2][0] = P_Z2;

			//与准确值作差得到XYZ坐标系下的误差
			//CMatrix ErrorX = X_POS - XT;

			//将误差投影到NEU东北天坐标系
			//CMatrix ErrorX_NEU = TT * ErrorX;

			if (i == 0)
				fprintf(fp, "%% GPST                        x-ecef(m)       y-ecef(m)       z-ecef(m)   Q  ns  ratio\n");
			int Q = 2;
			if (ratio_WL > 2)
				Q = 1;
			//fprintf(fp, "EpNum: %2d, %2d, %2d, %10.6f, %2d, %2d,%15.8g,%15.8g, %15.8g, %15.8g\n", j + 1, tempData2.hour, tempData2.minute, tempData2.second, bds2_data.sat_num, bds3_data.sat_num, ratio_WL,
				//ErrorX_NEU[0][0], ErrorX_NEU[1][0], ErrorX_NEU[2][0]);
			fprintf(fp, "%4d/%2d/%2d %6.2f  %10.4f  %10.4f  %10.4f  %3d %3d  %4.3f\n", tempData2.year, tempData2.month, tempData2.day, tempData2.GPSTIME,
				X_POS[0][0], X_POS[1][0], X_POS[2][0], Q, bds2num+bds3num, ratio_WL);

			for (int idx = 0; idx < 2; idx++)
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
			}

			//记录上历元的模糊度
			int a = 0;
			for (int i1 = 0; i1 < bds2_data.sat_num; i1++)
			{
				int PRN = bds2_data.sat2[i1].numofsat;
				//参考星跳过
				if (PRN == BDS2ref)
					continue;
				string nameofsat = "C" + to_string(PRN) + "-" + to_string(BDS2ref);
				preN_WL[nameofsat] = X_WL[a + 3][0];
				preN_EWL1[nameofsat] = X_EWL1[a + 3][0];
				preN_EWL2[nameofsat] = X_EWL2[a + 3][0];
				a++;
				fprintf(fp1, "%2d %4s  %10.4f  %10.4f  %10.4f\n", i + 1, nameofsat, preN_EWL1[nameofsat], preN_EWL2[nameofsat], preN_WL[nameofsat]);
			}

			for (int i2 = 0; i2 < bds3_data.sat_num; i2++)
			{
				int PRN = bds3_data.sat2[i2].numofsat;
				//参考星跳过
				if (PRN == BDS3ref)
					continue;
				string nameofsat = "B" + to_string(PRN) + "-" + to_string(BDS3ref);
				preN_WL[nameofsat] = X_WL[a + 3][0];
				preN_EWL1[nameofsat] = X_EWL1[a + 3][0];
				preN_EWL2[nameofsat] = X_EWL2[a + 3][0];
				a++;
				fprintf(fp1, "%2d %4s  %10.4f  %10.4f  %10.4f\n", i + 1, nameofsat, preN_EWL1[nameofsat], preN_EWL2[nameofsat], preN_WL[nameofsat]);
			}

			bInit = true;  //初始化成功
			Init1 = true;
			Init2 = true;
		}

	}
	fclose(fp);
	fclose(fp1);
}

//函数功能，TCAR法固定模糊度,采用几何相关模型,合起来计算超宽巷，宽巷和窄巷（基础模糊度）
void TCAR_test(string& filename1, string& filename2)
{
	//默认第一个站为基准站，第二个为流动站
	observe_spp sppfile1;
	observe_spp sppfile2;
	ReadSPP_KIN(filename1, sppfile1);
	ReadSPP_KIN(filename2, sppfile2);
	bool bInit = false; //初始化

	bool Initgps = true;  //单系统初始化
	bool Initbds2 = true;
	bool Initbds3 = true;
	bool Initgal = true;

	//输出差分定位结果用
	string name = filename2.substr(0, 3) + "_TriRTKtest.txt";
	string error = filename2.substr(0, 3) + "_TriDebugtest.txt";
	const char* output_filename = name.c_str();
	const char* out_file = error.c_str();
	FILE* fp;
	FILE* fp1;
	fopen_s(&fp, output_filename, "w");
	fopen_s(&fp1, out_file, "w");

	int OldSatPRN[2][30] = { {0} };
	int OldRefPRN[2] = { -1,-1 };
	int OldSatNum[2] = { 0 };

	CMatrix Q, X;
	map<string, double> preN_EWL, preN_WL, preN;

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
			TriOnCycleDetect(tempData2, preData2);     //判断此此历元与上一个历元中的卫星载波测量值是否有周跳,如果有周跳就将此历元中的卫星记为不可用
			//OnCycleDetect(tempData2, preData2);
		}

		//基准站数据相应时刻必须在前
		int j = 0;
		for (j = posk; j < sppfile1.liyuan_num; j++)
		{
			if (sppfile1.epoch[j].GPSTIME <= tempData2.GPSTIME)
				// ReSharper disable once CppRedundantControlFlowJump
				continue;
			break;
		}

		//回溯一个历元：posk即为要找的最近的历元且满足posk的时刻小于tempData2的时间
		if (j > 0) posk = --j;
		else
		{
			posk = j;
			continue;
		}

		//判断时间差：两个是否相隔超过历元间隔时间
		if (fabs(sppfile1.epoch[j].GPSTIME - tempData2.GPSTIME) >= sppfile1.INTERVAL)
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
			TriOnCycleDetect(tempData1, preData1);
			//OnCycleDetect(tempData1, preData1);
		}

		//差分定位共视处理
		EpochCommonPro(tempData1, tempData2, tempData);

		//分离BDS-2、BDS-3卫星，便于后续单系统、多系统处理
		Obs_epoch bds2_data;
		Obs_epoch bds3_data;

		int SatPRN[2][30] = { {0} };		//当前历元卫星号 依次存储
		int SatNum[2] = { 0 };	//当前历元各系统卫星数目
		int RefPRN[2] = { 0 };  //当前历元各系统参考卫星号

		int NewSatNum[2] = { 0 };
		int NewSatPRN[2][30] = { {0} };

		int k = 0;
		int m = 0;
		int bds2num, bds3num;
		bds2num = bds3num = 0;
		for (m = 0; m < tempData.sat_num; m++)
		{
			if (tempData.sat1[m].sattype == "C")
			{
				bds2_data.sat1.push_back(tempData.sat1[m]);
				bds2_data.sat2.push_back(tempData.sat2[m]);
				SatPRN[0][bds2num] = tempData.sat1[m].numofsat;
				bds2num++;
			}

			if (tempData.sat1[m].sattype == "B")
			{
				bds3_data.sat1.push_back(tempData.sat1[m]);
				bds3_data.sat2.push_back(tempData.sat2[m]);
				SatPRN[1][bds3num] = tempData.sat1[m].numofsat;
				bds3num++;
			}
		}

		//各类型卫星数
		SatNum[0] = bds2num;
		SatNum[1] = bds3num;

		bds2_data.sat_num = bds2_data.sat1.size();
		bds3_data.sat_num = bds3_data.sat1.size();
		tempData.GPSTIME = tempData2.GPSTIME;

		//记录新卫星
		if (bInit == true)
		{
			for (int idx = 0; idx < 2; idx++)
			{
				set<int>s;
				for (int ii = 0; ii < OldSatNum[idx]; ii++)
					s.insert(OldSatPRN[idx][ii]);
				for (int ii = 0; ii < SatNum[idx]; ii++)
				{
					if (s.find(SatPRN[idx][ii]) == s.end())
					{
						NewSatNum[idx]++;
						for (int iii = 0; iii < NewSatNum[idx]; iii++)
							NewSatPRN[idx][iii] = SatPRN[idx][ii];
					}
				}
			}
		}

		//卫星数至少要5颗
		if (bds2_data.sat_num + bds3_data.sat_num < 5)
		{
			printf("卫星数少于5颗\n");
			bInit = false;
			continue;
		}

		//各系统至少要2颗卫星,区分参考星和普通星
		if (bds2_data.sat_num < 2)
		{
			printf("BDS2卫星数少于2颗\n");
			Initbds2 = false;
			continue;
		}
		if (bds3_data.sat_num < 2)
		{
			printf("BDS3卫星数少于2颗\n");
			Initbds3 = false;
			continue;
		}

		//采用站间星间二次差分模型，选择参考卫星
		int BDS2ref = -1;
		int BDS3ref = -2;

		BDS2ref = OldRefPRN[0];
		BDS3ref = OldRefPRN[1];

		//选择最大高度角卫星作为参考卫星，不同类型卫星分开看
		SelectRef(bds2_data, BDS2ref);
		SelectRef(bds3_data, BDS3ref);

		RefPRN[0] = BDS2ref; //BDS-2参考卫星序号
		RefPRN[1] = BDS3ref; //BDS-3参考卫星序号

		//当前历元BDS-2、BDS-3卫星序列
		int PRN_X_B2[20];
		int PRN_X_B3[30];

		//数列初始化
		memset(PRN_X_B2, 0, sizeof(PRN_X_B2));
		memset(PRN_X_B3, 0, sizeof(PRN_X_B3));

		//储存历元中的卫星序列
		for (k = 0; k < bds2_data.sat_num; k++)
			PRN_X_B2[k] = bds2_data.sat1[k].numofsat;
		for (k = 0; k < bds3_data.sat_num; k++)
			PRN_X_B3[k] = bds3_data.sat1[k].numofsat;

		//各卫星个数
		int num_X_B2 = bds2_data.sat_num;
		int num_X_B3 = bds3_data.sat_num;

		//各系统双差方程个数
		int numbds2_x = num_X_B2 - 1;
		int numbds3_x = num_X_B3 - 1;

		//多系统融合双差方程个数（载波观测方程+伪距观测方程)，松组合
		int num_x_p = (numbds2_x + numbds3_x) * 6;

		//待估参数个数：ΔX，ΔY，ΔZ，每颗卫星的模糊度
		int num_x = 3 * numbds2_x + 3 * numbds3_x + 3;

		//模糊度个数
		int num_n = 3 * numbds2_x + 3 * numbds3_x;

		//各组合模糊度个数
		int num_s = numbds2_x + numbds3_x;

		//定位解算相关的矩阵
		//多系统
		CMatrix B(num_x_p, num_x);             //系数矩阵
		CMatrix L(num_x_p, 1);                 //观测值矩阵
		CMatrix R(num_x_p, num_x_p);           //误差矩阵（逆阵即为权矩阵）
		CMatrix Qw(num_x, num_x);              //动态噪声阵

		//单BDS2
		CMatrix RC2_P(numbds2_x, numbds2_x);      //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RC2_L(numbds2_x, numbds2_x);      //误差矩阵（逆阵即为权矩阵） 载波

		//单BDS3
		CMatrix RC3_P(numbds3_x, numbds3_x);      //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RC3_L(numbds3_x, numbds3_x);      //误差矩阵（逆阵即为权矩阵） 载波

		//位置的动态噪声阵，如果没有周跳，模糊度不变，所以模糊度的噪声为0，只有位置噪声
		for (int idx = 0; idx < 3; idx++) {
			Qw[idx][idx] = 30.0 * 30.0;
		}

		//更新BDS2和BDS3的协方差矩阵Q和结果矩阵X
		//初始化
		if (!bInit)
		{
			Q.SetSize(num_x, num_x);
			X.SetSize(num_x, 1);
			InitState_Tri(Q);
		}
		else
		{
			//检测新卫星
			CheckNewSats(SatNum, OldSatNum, SatPRN, OldSatPRN, NewSatPRN, NewSatNum);

			//更新协方差和X矩阵
			UpTridataState(Q, X, SatNum, OldSatNum, SatPRN, OldSatPRN, RefPRN, OldRefPRN);
		}

		//各类型卫星序号存储容器
		vector<int> BDS2_PRN;
		vector<int> BDS3_PRN;

		for (k = 0; k < bds2_data.sat_num; k++)
			BDS2_PRN.push_back(PRN_X_B2[k]);
		for (k = 0; k < bds3_data.sat_num; k++)
			BDS3_PRN.push_back(PRN_X_B3[k]);

		//伪距
		GetDDR_Peo(bds2_data, BDS2_PRN, BDS2ref, RC2_P);
		GetDDR_Peo(bds3_data, BDS3_PRN, BDS3ref, RC3_P);
		//载波
		GetDDR_Leo(bds2_data, BDS2_PRN, BDS2ref, RC2_L);
		GetDDR_Leo(bds3_data, BDS3_PRN, BDS3ref, RC3_L);

		//融合系统的R矩阵，分块对角矩阵
		//伪距P1
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R[k1][k2] = RC2_P[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R[k1 + numbds2_x][k2 + numbds2_x] = RC3_P[k1][k2];
			}
		//伪距P2
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R[k1 + num_s][k2 + num_s] = RC2_P[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R[k1 + numbds2_x + num_s][k2 + numbds2_x + num_s] = RC3_P[k1][k2];
			}
		//伪距P3
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R[k1 + 2 * num_s][k2 + 2 * num_s] = RC2_P[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R[k1 + numbds2_x + 2 * num_s][k2 + numbds2_x + 2 * num_s] = RC3_P[k1][k2];
			}
		//载波EWL,BDS2/BDS3(0,-1,1)
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R[k1 + 3 * num_s][k2 + 3 * num_s] = 2 * RC2_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R[k1 + 3 * num_s + numbds2_x][k2 + 3 * num_s + numbds2_x] = 2 * RC3_L[k1][k2];
			}
		//载波WL,BDS2/BDS3(1,-1,0)
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R[k1 + 4 * num_s][k2 + 4 * num_s] = 2 * RC2_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R[k1 + 4 * num_s + numbds2_x][k2 + 4 * num_s + numbds2_x] = 2 * RC3_L[k1][k2];
			}
		//载波基础模糊度,BDS2/BDS3(1,0,0)
		for (int k1 = 0; k1 < numbds2_x; k1++)
			for (int k2 = 0; k2 < numbds2_x; k2++)
			{
				R[k1 + 5 * num_s][k2 + 5 * num_s] = RC2_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds3_x; k1++)
			for (int k2 = 0; k2 < numbds3_x; k2++)
			{
				R[k1 + 5 * num_s + numbds2_x][k2 + 5 * num_s + numbds2_x] = RC3_L[k1][k2];
			}

		//初始坐标信息：基站采用spp文件头中坐标（固定），流动站采用每个历元的概率坐标（考虑实际应用中其可能为动态情况）
		//本组数据中：流动站为动态，因此需要更新其坐标
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
		int ref_bds2 = -1;
		int ref_bds3 = -1;

		for (k = 0; k < bds2_data.sat_num; k++)
		{
			if (bds2_data.sat1[k].numofsat == BDS2ref)
			{
				ref_bds2 = k;
				break;
			}
		}

		for (k = 0; k < bds3_data.sat_num; k++)
		{
			if (bds3_data.sat1[k].numofsat == BDS3ref)
			{
				ref_bds3 = k;
				break;
			}
		}

		int tp = 0;  //记录单系统方程数
		int tw = 0;  //记录多系统方程数
		double DDerror = 0;  //初始化误差改正项
		//BDS-2
		for (k = 0; k < bds2_data.sat_num; k++)
		{
			//bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = bds2_data.sat1[k].numofsat;
			string nameofsat = "C" + to_string(PRN) + "-" + to_string(BDS2ref);
			//参考卫星跳过
			if (PRN == BDS2ref)
			{
				continue;
			}

			if (ref_bds2 < 0)
			{
				assert(0);
			}

			Sat tempSat1 = bds2_data.sat1[k];
			Sat RefSat1 = bds2_data.sat1[ref_bds2];
			Sat tempSat2 = bds2_data.sat2[k];
			Sat RefSat2 = bds2_data.sat2[ref_bds2];

			//误差改正
			DDerror = Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2);

			//计算卫星到测站距离（4个距离）
			double length_O1, length_R1, length_O2, length_R2;
			length_O1 = sqrt(pow((tempSat1.POS_X - POSITION_X1), 2) + pow((tempSat1.POS_Y - POSITION_Y1), 2) + pow((tempSat1.POS_Z - POSITION_Z1), 2));
			length_R1 = sqrt(pow((RefSat1.POS_X - POSITION_X1), 2) + pow((RefSat1.POS_Y - POSITION_Y1), 2) + pow((RefSat1.POS_Z - POSITION_Z1), 2));

			length_O2 = sqrt(pow((tempSat2.POS_X - POSITION_X2), 2) + pow((tempSat2.POS_Y - POSITION_Y2), 2) + pow((tempSat2.POS_Z - POSITION_Z2), 2));
			length_R2 = sqrt(pow((RefSat2.POS_X - POSITION_X2), 2) + pow((RefSat2.POS_Y - POSITION_Y2), 2) + pow((RefSat2.POS_Z - POSITION_Z2), 2));

			//三个坐标系数
			double Cof1 = (tempSat2.POS_X - POSITION_X2) / length_O2 - (RefSat2.POS_X - POSITION_X2) / length_R2;
			double Cof2 = (tempSat2.POS_Y - POSITION_Y2) / length_O2 - (RefSat2.POS_Y - POSITION_Y2) / length_R2;
			double Cof3 = (tempSat2.POS_Z - POSITION_Z2) / length_O2 - (RefSat2.POS_Z - POSITION_Z2) / length_R2;

			//系数矩阵
			//P1
			B[tw][0] = Cof1; B[tw][1] = Cof2; B[tw][2] = Cof3;
			//P2
			B[tw + num_s][0] = Cof1; B[tw + num_s][1] = Cof2; B[tw + num_s][2] = Cof3;
			//P3
			B[tw + num_s * 2][0] = Cof1; B[tw + num_s * 2][1] = Cof2; B[tw + num_s * 2][2] = Cof3;
			//EWL
			B[tw + num_s * 3][0] = Cof1; B[tw + num_s * 3][1] = Cof2; B[tw + num_s * 3][2] = Cof3; B[tw + num_s * 3][3 + 2 * num_s + tw] = lambda_L1_C;
			//WL
			B[tw + num_s * 4][0] = Cof1; B[tw + num_s * 4][1] = Cof2; B[tw + num_s * 4][2] = Cof3;
			B[tw + num_s * 4][3 + num_s + tw] = -lambda_L2_C; B[tw + num_s * 4][3 + 2 * num_s + tw] = lambda_L2_C;
			//L
			B[tw + num_s * 5][0] = Cof1; B[tw + num_s * 5][1] = Cof2; B[tw + num_s * 5][2] = Cof3;
			B[tw + num_s * 5][3 + tw] = lambda_L3_C; B[tw + num_s * 5][3 + num_s + tw] = -lambda_L3_C; B[tw + num_s * 5][3 + 2 * num_s + tw] = lambda_L3_C;

			//站星距
			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;

			//观测值矩阵
			double deltP1 = oper_TriP[1](tempSat2, 1, 0, 0) - oper_TriP[1](RefSat2, 1, 0, 0) - (oper_TriP[1](tempSat1, 1, 0, 0) - oper_TriP[1](RefSat1, 1, 0, 0));
			double deltP2 = oper_TriP[1](tempSat2, 0, 1, 0) - oper_TriP[1](RefSat2, 0, 1, 0) - (oper_TriP[1](tempSat1, 0, 1, 0) - oper_TriP[1](RefSat1, 0, 1, 0));
			double deltP3 = oper_TriP[1](tempSat2, 0, 0, 1) - oper_TriP[1](RefSat2, 0, 0, 1) - (oper_TriP[1](tempSat1, 0, 0, 1) - oper_TriP[1](RefSat1, 0, 0, 1));
			double deltL1 = oper_TriPhi[1](tempSat2, 1, 0, 0) - oper_TriPhi[1](RefSat2, 1, 0, 0) - (oper_TriPhi[1](tempSat1, 1, 0, 0) - oper_TriPhi[1](RefSat1, 1, 0, 0));
			double deltL2 = oper_TriPhi[1](tempSat2, 0, 1, 0) - oper_TriPhi[1](RefSat2, 0, 1, 0) - (oper_TriPhi[1](tempSat1, 0, 1, 0) - oper_TriPhi[1](RefSat1, 0, 1, 0));
			double deltL3 = oper_TriPhi[1](tempSat2, 0, 0, 1) - oper_TriPhi[1](RefSat2, 0, 0, 1) - (oper_TriPhi[1](tempSat1, 0, 0, 1) - oper_TriPhi[1](RefSat1, 0, 0, 1));
			L[tw][0] = deltP1 - DDgeo; L[tw + num_s][0] = deltP2 - DDgeo; L[tw + 2 * num_s][0] = deltP3 - DDgeo;
			L[tw + 3 * num_s][0] = deltL1 * lambda_L1_C - DDgeo; L[tw + 4 * num_s][0] = deltL2 * lambda_L2_C - DDgeo; L[tw + 5 * num_s][0] = deltL3 * lambda_L3_C - DDgeo;

			//判断旧参考星是否被剔除
			/*for (int ii = 0; ii < SatNum[0]; ii++)
			{
				if (OldRefPRN[0] == SatPRN[0][ii])
				{
					Init = true;
					break;
				}
				else
					Init = false;
			}*/

			/*if (i >= 1)
			{
				if (fabs(X_WL[tw + 3][0] - preN_WL[nameofsat]) > 0.2)
				{
					double var = (X_WL[3 + tw][0] + preN_WL[nameofsat]) / 2;
					X_WL[3 + tw][0] = var;
				}
			}*/

			//初始化模糊度
			if (fabs(X[3 + tw][0] - 0.0) < 1.0 || !Initbds2)
			{
				//初始化超宽巷(0,-1,1)
				const double deltLWN = -deltL2 + deltL3 - (deltP2 + deltP3) / lambda_longwideC1;
				X[3 + tw][0] = deltLWN;
				//初始化宽巷(1,-1,0)
				const double deltP = oper_TriP[1](tempSat2, 0.4972, 0.2187, 0.2841) - oper_TriP[1](RefSat2, 0.4972, 0.2187, 0.2841) -
					(oper_TriP[1](tempSat1, 0.4972, 0.2187, 0.2841) - oper_TriP[1](RefSat1, 0.4972, 0.2187, 0.2841));
				const double deltWN = deltL1 - deltL2 - deltP / lambda_wideC;
				X[3 + num_s + tw][0] = deltWN;
				//初始化基础模糊度(1,0,0)
				X[3 + 2 * num_s + tw][0] = deltL1 - deltP1 / lambda_L1_C;
				//初始化超宽巷协方差
				Q[3 + tw][3 + tw] = 10.0 * 10.0;
				//初始化宽巷协方差
				Q[3 + tw + num_s][3 + tw + num_s] = 30.0 * 30.0;
				//初始化基础模糊度协方差
				Q[3 + tw + 2 * num_s][3 + tw + 2 * num_s] = 50.0 * 50.0;
			}

			if (!bInit)
			{
				//初始化超宽巷(0,-1,1)
				const double deltLWN = -deltL2 + deltL3 - (deltP2 + deltP3) / lambda_longwideC1;
				X[3 + tw][0] = deltLWN;
				//初始化宽巷(1,-1,0)
				const double deltP = oper_TriP[1](tempSat2, 0.4972, 0.2187, 0.2841) - oper_TriP[1](RefSat2, 0.4972, 0.2187, 0.2841) -
					(oper_TriP[1](tempSat1, 0.4972, 0.2187, 0.2841) - oper_TriP[1](RefSat1, 0.4972, 0.2187, 0.2841));
				const double deltWN = deltL1 - deltL2 - deltP / lambda_wideC;
				X[3 + num_s + tw][0] = deltWN;
				//初始化基础模糊度(1,0,0)
				X[3 + 2 * num_s + tw][0] = deltL1 - deltP1 / lambda_L1_C;
			}
			else
			{
				for (int ii = 0; ii < NewSatNum[0]; ii++)
				{
					if (NewSatPRN[0][ii] == PRN)
					{
						//初始化超宽巷(0,-1,1)
						const double deltLWN = -deltL2 + deltL3 - (deltP2 + deltP3) / lambda_longwideC1;
						X[3 + tw][0] = deltLWN;
						//初始化宽巷(1,-1,0)
						const double deltP = oper_TriP[1](tempSat2, 0.4972, 0.2187, 0.2841) - oper_TriP[1](RefSat2, 0.4972, 0.2187, 0.2841) -
							(oper_TriP[1](tempSat1, 0.4972, 0.2187, 0.2841) - oper_TriP[1](RefSat1, 0.4972, 0.2187, 0.2841));
						const double deltWN = deltL1 - deltL2 - deltP / lambda_wideC;
						X[3 + num_s + tw][0] = deltWN;
						//初始化基础模糊度(1,0,0)
						X[3 + 2 * num_s + tw][0] = deltL1 - deltP1 / lambda_L1_C;
						//初始化超宽巷协方差
						Q[3 + tw][3 + tw] = 10.0 * 10.0;
						//初始化宽巷协方差
						Q[3 + tw + num_s][3 + tw + num_s] = 30.0 * 30.0;
						//初始化基础模糊度协方差
						Q[3 + tw + 2 * num_s][3 + tw + 2 * num_s] = 50.0 * 50.0;
						break;
					}
				}
			}
			tw++;
		}

		tp = 0;           //单系统方程数初始化
		DDerror = 0;
		//BDS-3
		for (k = 0; k < bds3_data.sat_num; k++)
		{
			//bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = bds3_data.sat1[k].numofsat;
			string nameofsat = "B" + to_string(PRN) + "-" + to_string(BDS3ref);
			//参考卫星跳过
			if (PRN == BDS3ref)
			{
				continue;
			}

			if (ref_bds3 < 0)
			{
				assert(0);
			}

			Sat tempSat1 = bds3_data.sat1[k];
			Sat RefSat1 = bds3_data.sat1[ref_bds3];
			Sat tempSat2 = bds3_data.sat2[k];
			Sat RefSat2 = bds3_data.sat2[ref_bds3];

			//误差改正
			DDerror = Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2);

			//计算卫星到测站距离（4个距离）
			double length_O1, length_R1, length_O2, length_R2;
			length_O1 = sqrt(pow((tempSat1.POS_X - POSITION_X1), 2) + pow((tempSat1.POS_Y - POSITION_Y1), 2) + pow((tempSat1.POS_Z - POSITION_Z1), 2));
			length_R1 = sqrt(pow((RefSat1.POS_X - POSITION_X1), 2) + pow((RefSat1.POS_Y - POSITION_Y1), 2) + pow((RefSat1.POS_Z - POSITION_Z1), 2));

			length_O2 = sqrt(pow((tempSat2.POS_X - POSITION_X2), 2) + pow((tempSat2.POS_Y - POSITION_Y2), 2) + pow((tempSat2.POS_Z - POSITION_Z2), 2));
			length_R2 = sqrt(pow((RefSat2.POS_X - POSITION_X2), 2) + pow((RefSat2.POS_Y - POSITION_Y2), 2) + pow((RefSat2.POS_Z - POSITION_Z2), 2));

			//三个坐标系数
			double Cof1 = (tempSat2.POS_X - POSITION_X2) / length_O2 - (RefSat2.POS_X - POSITION_X2) / length_R2;
			double Cof2 = (tempSat2.POS_Y - POSITION_Y2) / length_O2 - (RefSat2.POS_Y - POSITION_Y2) / length_R2;
			double Cof3 = (tempSat2.POS_Z - POSITION_Z2) / length_O2 - (RefSat2.POS_Z - POSITION_Z2) / length_R2;

			//系数矩阵
			//P1
			B[tw][0] = Cof1; B[tw][1] = Cof2; B[tw][2] = Cof3;
			//P2
			B[tw + num_s][0] = Cof1; B[tw + num_s][1] = Cof2; B[tw + num_s][2] = Cof3;
			//P3
			B[tw + num_s * 2][0] = Cof1; B[tw + num_s * 2][1] = Cof2; B[tw + num_s * 2][2] = Cof3;
			//EWL
			B[tw + num_s * 3][0] = Cof1; B[tw + num_s * 3][1] = Cof2; B[tw + num_s * 3][2] = Cof3; B[tw + num_s * 3][3 + 2 * num_s + tw] = lambda_L1_B;
			//WL
			B[tw + num_s * 4][0] = Cof1; B[tw + num_s * 4][1] = Cof2; B[tw + num_s * 4][2] = Cof3;
			B[tw + num_s * 4][3 + tw] = -lambda_L2_B; B[tw + num_s * 4][3 + num_s + tw] = -lambda_L2_B; B[tw + num_s * 4][3 + 2 * num_s + tw] = lambda_L2_B;
			//L
			B[tw + num_s * 5][0] = Cof1; B[tw + num_s * 5][1] = Cof2; B[tw + num_s * 5][2] = Cof3;
			B[tw + num_s * 5][3 + num_s + tw] = -lambda_L3_B; B[tw + num_s * 5][3 + 2 * num_s + tw] = lambda_L3_B;

			//站星距
			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;        

			//观测值矩阵
			double deltP1 = oper_TriP[2](tempSat2, 1, 0, 0) - oper_TriP[2](RefSat2, 1, 0, 0) - (oper_TriP[2](tempSat1, 1, 0, 0) - oper_TriP[2](RefSat1, 1, 0, 0));
			double deltP2 = oper_TriP[2](tempSat2, 0, 1, 0) - oper_TriP[2](RefSat2, 0, 1, 0) - (oper_TriP[2](tempSat1, 0, 1, 0) - oper_TriP[2](RefSat1, 0, 1, 0));
			double deltP3 = oper_TriP[2](tempSat2, 0, 0, 1) - oper_TriP[2](RefSat2, 0, 0, 1) - (oper_TriP[2](tempSat1, 0, 0, 1) - oper_TriP[2](RefSat1, 0, 0, 1));
			double deltL1 = oper_TriPhi[2](tempSat2, 1, 0, 0) - oper_TriPhi[2](RefSat2, 1, 0, 0) - (oper_TriPhi[2](tempSat1, 1, 0, 0) - oper_TriPhi[2](RefSat1, 1, 0, 0));
			double deltL2 = oper_TriPhi[2](tempSat2, 0, 1, 0) - oper_TriPhi[2](RefSat2, 0, 1, 0) - (oper_TriPhi[2](tempSat1, 0, 1, 0) - oper_TriPhi[2](RefSat1, 0, 1, 0));
			double deltL3 = oper_TriPhi[2](tempSat2, 0, 0, 1) - oper_TriPhi[2](RefSat2, 0, 0, 1) - (oper_TriPhi[2](tempSat1, 0, 0, 1) - oper_TriPhi[2](RefSat1, 0, 0, 1));
			L[tw][0] = deltP1 - DDgeo; L[tw + num_s][0] = deltP2 - DDgeo; L[tw + 2 * num_s][0] = deltP3 - DDgeo;
			L[tw + 3 * num_s][0] = deltL1 * lambda_L1_B - DDgeo; L[tw + 4 * num_s][0] = deltL2 * lambda_L2_B - DDgeo; L[tw + 5 * num_s][0] = deltL3 * lambda_L3_B - DDgeo;

			//判断旧参考星是否被剔除
			/*for (int ii = 0; ii < SatNum[1]; ii++)
			{
				if (OldRefPRN[1] == SatPRN[1][ii])
				{
					Init1 = true;
					break;
				}
				else
					Init1 = false;
			}*/

			/*if (i >= 1)
			{
				if (fabs(X_WL[tw + 3][0] - preN_WL[nameofsat]) > 0.2)
				{
					double var = (X_WL[3 + tw][0] + preN_WL[nameofsat]) / 2;
					X_WL[3 + tw][0] = var;
				}
			}*/

			//初始化模糊度
			if (fabs(X[3 + tw][0] - 0.0) < 1.0 || !Initbds3)
			{
				//初始化超宽巷模糊度(0,-1,1)
				const double deltP_EWL = oper_TriP[2](tempSat2, 0.0401, 0.5649, 0.3950) - oper_TriP[2](RefSat2, 0.0401, 0.5649, 0.3950) -
					(oper_TriP[2](tempSat1, 0.0401, 0.5649, 0.3950) - oper_TriP[2](RefSat1, 0.0401, 0.5649, 0.3950));
				const double deltLWN = deltL3 - deltL2 - deltP_EWL / lambda_longwideB1;
				X[3 + tw][0] = deltLWN;
				//初始化宽巷模糊度(1,0,-1)
				const double deltP_WL = oper_TriP[2](tempSat2, 0.6076, 0.1167, 0.2757) - oper_TriP[2](RefSat2, 0.6076, 0.1167, 0.2757) -
					(oper_TriP[2](tempSat1, 0.6076, 0.1167, 0.2757) - oper_TriP[2](RefSat1, 0.6076, 0.1167, 0.2757));
				const double deltWN = deltL1 - deltL3 - deltP_WL / lambda_wideB;
				X[3 + num_s + tw][0] = deltWN;
				//初始化基础模糊度(1,0,0)
				X[3 + 2 * num_s + tw][0] = deltL1 - deltP1 / lambda_L1_B;
				//初始化超宽巷协方差
				Q[3 + tw][3 + tw] = 10.0 * 10.0;
				//初始化宽巷协方差
				Q[3 + tw + num_s][3 + tw + num_s] = 30.0 * 30.0;
				//初始化基础模糊度协方差
				Q[3 + tw + 2 * num_s][3 + tw + 2 * num_s] = 50.0 * 50.0;
			}


			if (!bInit)
			{
				//初始化超宽巷模糊度(0,-1,1)
				const double deltP_EWL = oper_TriP[2](tempSat2, 0.0401, 0.5649, 0.3950) - oper_TriP[2](RefSat2, 0.0401, 0.5649, 0.3950) -
					(oper_TriP[2](tempSat1, 0.0401, 0.5649, 0.3950) - oper_TriP[2](RefSat1, 0.0401, 0.5649, 0.3950));
				const double deltLWN = deltL3 - deltL2 - deltP_EWL / lambda_longwideB1;
				X[3 + tw][0] = deltLWN;
				//初始化宽巷模糊度(1,0,-1)
				const double deltP_WL = oper_TriP[2](tempSat2, 0.6076, 0.1167, 0.2757) - oper_TriP[2](RefSat2, 0.6076, 0.1167, 0.2757) -
					(oper_TriP[2](tempSat1, 0.6076, 0.1167, 0.2757) - oper_TriP[2](RefSat1, 0.6076, 0.1167, 0.2757));
				const double deltWN = deltL1 - deltL3 - deltP_WL / lambda_wideB;
				X[3 + num_s + tw][0] = deltWN;
				//初始化基础模糊度(1,0,0)
				X[3 + 2 * num_s + tw][0] = deltL1 - deltP1 / lambda_L1_B;
			}
			else
			{
				for (int ii = 0; ii < NewSatNum[1]; ii++)
				{
					if (NewSatPRN[1][ii] == PRN)
					{
						//初始化超宽巷模糊度(0,-1,1)
						const double deltP_EWL = oper_TriP[2](tempSat2, 0.0401, 0.5649, 0.3950) - oper_TriP[2](RefSat2, 0.0401, 0.5649, 0.3950) -
							(oper_TriP[2](tempSat1, 0.0401, 0.5649, 0.3950) - oper_TriP[2](RefSat1, 0.0401, 0.5649, 0.3950));
						const double deltLWN = deltL3 - deltL2 - deltP_EWL / lambda_longwideB1;
						X[3 + tw][0] = deltLWN;
						//初始化宽巷模糊度(1,0,-1)
						const double deltP_WL = oper_TriP[2](tempSat2, 0.6076, 0.1167, 0.2757) - oper_TriP[2](RefSat2, 0.6076, 0.1167, 0.2757) -
							(oper_TriP[2](tempSat1, 0.6076, 0.1167, 0.2757) - oper_TriP[2](RefSat1, 0.6076, 0.1167, 0.2757));
						const double deltWN = deltL1 - deltL3 - deltP_WL / lambda_wideB;
						X[3 + num_s + tw][0] = deltWN;
						//初始化基础模糊度(1,0,0)
						X[3 + 2 * num_s + tw][0] = deltL1 - deltP1 / lambda_L1_B;
						//初始化超宽巷协方差
						Q[3 + tw][3 + tw] = 10.0 * 10.0;
						//初始化宽巷协方差
						Q[3 + tw + num_s][3 + tw + num_s] = 30.0 * 30.0;
						//初始化基础模糊度协方差
						Q[3 + tw + 2 * num_s][3 + tw + 2 * num_s] = 50.0 * 50.0;
						break;
					}
				}
			}
			tw++;
		}


		if (bds2_data.sat_num + bds3_data.sat_num >= 5)
		{
			//B系数矩阵  L观测值矩阵  R误差矩阵  Qw动态噪声矩阵  Q协方差矩阵  X位置参数和模糊度

			double ratio[3] = { 1.0 };
			CMatrix aXb(3, 1);
			CMatrix fix(num_n, 1);

			//X.MyTRACE();
			//L.MyTRACE();
			//R.MyTRACE();
			//Qw.MyTRACE();
			//Q.MyTRACE();
			//B.MyTRACE();

			//kalman滤波
			KalmanFilter(B, L, R, Qw, Q, X);

			//lambda算法固定模糊度
			//Resamb_LAMBDA(Q, X, aXb, ratio, fix);
			ParAR_LAMBDA(X, Q, aXb, num_s, ratio, fix);
			
			//恢复出估计的坐标值
			CMatrix X_POS(3, 1);
			if (ratio[2] > 2)
			{
				X_POS[0][0] = POSITION_X2 - aXb[0][0];
				X_POS[1][0] = POSITION_Y2 - aXb[1][0];
				X_POS[2][0] = POSITION_Z2 - aXb[2][0];
			}
			else
			{
				X_POS[0][0] = POSITION_X2 - X[0][0];
				X_POS[1][0] = POSITION_Y2 - X[1][0];
				X_POS[2][0] = POSITION_Z2 - X[2][0];
			}

			//真实坐标写成矩阵
			CMatrix XT(3, 1);
			XT[0][0] = P_X2;
			XT[1][0] = P_Y2;
			XT[2][0] = P_Z2;

			//与准确值作差得到XYZ坐标系下的误差
			//CMatrix ErrorX = X_POS - XT;

			if (i == 0)
				fprintf(fp, "%% GPST                        x-ecef(m)       y-ecef(m)       z-ecef(m)   Q  ns  ratio\n");
			int Q = 2;
			if (ratio[2] > 2)
				Q = 1;
			//fprintf(fp, "EpNum: %2d, %2d, %2d, %10.6f, %2d, %2d,%15.8g,%15.8g, %15.8g, %15.8g\n", j + 1, tempData2.hour, tempData2.minute, tempData2.second, bds2_data.sat_num, bds3_data.sat_num, ratio_WL,
				//ErrorX_NEU[0][0], ErrorX_NEU[1][0], ErrorX_NEU[2][0]);
			/*fprintf(fp, "%4d/%2d/%2d %6.2f  %10.4f  %10.4f  %10.4f  %3d %3d  %4.3f\n", tempData2.year, tempData2.month, tempData2.day, tempData2.GPSTIME,
				X_POS[0][0], X_POS[1][0], X_POS[2][0], Q, bds2num + bds3num, ratio);*/
			fprintf(fp, "%4d %6.2f  %10.4f  %10.4f  %10.4f  %3d %3d  %4.3f\n", 2128, tempData2.GPSTIME, X_POS[0][0], X_POS[1][0], X_POS[2][0], Q, bds2num + bds3num, ratio[2]);

			//记录上历元的模糊度
			int a = 0;
			for (int i1 = 0; i1 < bds2_data.sat_num; i1++)
			{
				int PRN = bds2_data.sat2[i1].numofsat;
				//参考星跳过
				if (PRN == BDS2ref)
					continue;
				string nameofsat = "C" + to_string(PRN) + "-" + to_string(BDS2ref);
				preN_EWL[nameofsat] = X[a + 3][0];
				preN_WL[nameofsat] = X[a + 3 + num_s][0];
				preN[nameofsat] = X[a + 3 + 2 * num_s][0];
				a++;
				fprintf(fp1, "%2d %4s  %10.4f  %10.4f  %10.4f\n", i + 1, nameofsat, preN_EWL[nameofsat], preN_WL[nameofsat], preN[nameofsat]);
			}

			for (int i2 = 0; i2 < bds3_data.sat_num; i2++)
			{
				int PRN = bds3_data.sat2[i2].numofsat;
				//参考星跳过
				if (PRN == BDS3ref)
					continue;
				string nameofsat = "B" + to_string(PRN) + "-" + to_string(BDS3ref);
				preN_EWL[nameofsat] = X[a + 3][0];
				preN_WL[nameofsat] = X[a + 3 + num_s][0];
				preN[nameofsat] = X[a + 3 + 2 * num_s][0];
				a++;
				fprintf(fp1, "%2d %4s  %10.4f  %10.4f  %10.4f\n", i + 1, nameofsat, preN_EWL[nameofsat], preN_WL[nameofsat], preN[nameofsat]);
			}

			bInit = true;  //初始化成功
			Initbds2 = true;
			Initbds3 = true;
		}

		for (int idx = 0; idx < 2; idx++)
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
	}
	fclose(fp);
	fclose(fp1);
}