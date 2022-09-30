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

//声明外部GetGPSTime函数，在ReadPfile.cpp中定义
double GetGPSTime(int year, int month, int day, int hour, int minute, double second, int& dayofy);

//声明EpochCommonPro函数，对基站和流动站卫星观测数据进行共视处理
void EpochCommonPro(Obs_epoch tempData1, Obs_epoch tempData2, Obs_epoch& useData);

//声明SelectRef函数，选择参考卫星
void SelectRef(Obs_epoch tempData, int& RefNo);

//声明ReadSPP_KIN函数，读取spp文件（针对本人自定义的spp文件格式，可修改）
bool ReadSPP_KIN(string filename, observe_spp& sppfile);

//声明GetDDR_Peo函数，依据误差传播定律获取噪声矩阵（伪距）
void GetDDR_Peo(Obs_epoch GPSData, vector<int> useprn, int GPSRef, CMatrix& R);

// 函数功能：依据误差传播定律获取噪声矩阵（载波）
void GetDDR_Leo(Obs_epoch GPSData, vector<int> useprn, int GPSRef, CMatrix& R) {
	//R矩阵设定
	int numgps = GPSData.sat_num;
	if (numgps != useprn.size()) {
		assert(0);
	}
	int numgps_x = numgps - 1;
	CMatrix TTR(numgps_x, numgps);
	CMatrix RT(numgps, numgps);
	int t = 0;
	int k = 0;
	for (k = 0; k < numgps; k++) {
		int match_O = -1;
		int PRN = useprn[k];
		for (int m = 0; m < GPSData.sat_num; m++) {
			if (GPSData.sat2[m].numofsat == PRN) {
				match_O = m;
				break;
			}
		}

		RT[k][k] = pow((0.003 + 0.003 * exp(-GPSData.sat2[match_O].E / 10.0)), 2) + pow((0.003 + 0.003 * exp(-GPSData.sat1[match_O].E / 10.0)), 2);

		if (PRN != GPSRef) {
			TTR[t][k] = 1.0;
			t++;
		}
		else {
			for (int m = 0; m < numgps_x; m++) {
				TTR[m][k] = -1.0;
			}
		}
	}

	R.first(numgps_x, numgps_x);
	//TTR矩阵的目的是根据误差传播定律，导出星间差分模式下的误差矩阵形式
	R = TTR * RT * TTR.T();

}

//函数功能：卡尔曼滤波
void KalmanFilter(CMatrix& B, CMatrix& L, CMatrix& R, CMatrix& Qw, CMatrix& Q, CMatrix& X) {
	int r = X.Row;
	CMatrix I, F;
	I.SetSize(r, r);
	F.SetSize(r, r);

	for (int i = 0; i < r; i++) {
		I[i][i] = 1.0;
		F[i][i] = 1.0;
	}
	CMatrix M = F * Q * F.T() + Qw;
	CMatrix K = M * B.T() * (B * M * B.T() + R).InvertGaussJordan();
	Q = (I - K * B) * M * (I - B.T() * K.T()) + K * R * K.T();
	CMatrix V = L - B * X;
	X = X + K * V;
}

//函数功能：周跳探测
void OnCycleDetect(Obs_epoch& tempEpoch, Obs_epoch& preEpoch) {
	double slip1 = 0;
	double slip2 = 0;
	for (int i = 0; i < tempEpoch.sat_num; i++) {
		Sat* sat = &tempEpoch.sat[i];
		if (sat->judge_use)
			continue;
		if ((sat->sattype != "G") && (sat->sattype != "C") && (sat->sattype != "B") && (sat->sattype != "E"))
			continue;
		for (int j = 0; j < preEpoch.sat_num; j++) {
			if ((tempEpoch.sat[i].sattype == preEpoch.sat[j].sattype) && (tempEpoch.sat[i].numofsat == preEpoch.sat[j].numofsat)) {
				Sat* presat = &preEpoch.sat[j];
				if (presat->judge_use)
					continue;
				if ((presat->sattype != "G") && (presat->sattype != "C") && (presat->sattype != "B") && (presat->sattype != "E"))
					continue;
				//GPS
				if (presat->sattype == "G") {
					slip1 = sat->data[2] - presat->data[2] - (FREQ1 / FREQ2) * (sat->data[3] - presat->data[3]);//电离层残差
					slip2 = (-(sat->data[2] - sat->data[3]) + 1 / lambda_L_wide * lambda_L_narrow * (sat->data[0] / lambda_L1 + sat->data[1] / lambda_L2))
						- (-(presat->data[2] - presat->data[3]) + 1 / lambda_L_wide * lambda_L_narrow * (presat->data[0] / lambda_L1 + presat->data[1] / lambda_L2)); //MW
				}
				//BDS-2
				else if (presat->sattype == "C") {
					slip1 = sat->data[2] - presat->data[2] - (FREQ1_BDS / FREQ2_BDS) * (sat->data[3] - presat->data[3]);//电离层残差
					slip2 = (-(sat->data[2] - sat->data[3]) + 1 / lambda_L_wide_C * lambda_L_narrow_C * (sat->data[0] / lambda_L1_C + sat->data[1] / lambda_L2_C))
						- (-(presat->data[2] - presat->data[3]) + 1 / lambda_L_wide_C * lambda_L_narrow_C * (presat->data[0] / lambda_L1_C + presat->data[1] / lambda_L2_C)); //MW
				}
				//BDS-3
				else if (presat->sattype == "B") {
					slip1 = sat->data[2] - presat->data[2] - (f1B / f2B) * (sat->data[3] - presat->data[3]);//电离层残差
					slip2 = (-(sat->data[2] - sat->data[3]) + 1 / lambda_L_wide_B * lambda_L_narrow_B * (sat->data[0] / lambda_L1_B + sat->data[1] / lambda_L2_B))
						- (-(presat->data[2] - presat->data[3]) + 1 / lambda_L_wide_B * lambda_L_narrow_B * (presat->data[0] / lambda_L1_B + presat->data[1] / lambda_L2_B));//MW
				}
				//Galileo
				else if (presat->sattype == "E") {
					slip1 = sat->data[2] - presat->data[2] - (f1E / f2E) * (sat->data[3] - presat->data[3]);//电离层残差
					slip2 = (-(sat->data[2] - sat->data[3]) + 1 / lambda_L_wide_E * lambda_L_narrow_E * (sat->data[0] / lambda_L1_E + sat->data[1] / lambda_L2_E))
						- (-(presat->data[2] - presat->data[3]) + 1 / lambda_L_wide_E * lambda_L_narrow_E * (presat->data[0] / lambda_L1_E + presat->data[1] / lambda_L2_E)); //MW
				}
				if (fabs(slip1) > 0.05 || fabs(slip2) > 10) {
					sat->judge_use = 1;
				}
				break;
			}
		}
	}
}

//函数功能：生成转换矩阵
bool GetGRCQValue(int PRN_X[], int PRN_old[], int num_sat, int num_sat_old, int Ref, int Ref_old, CMatrix& TT_Matrix)
{
	int i, j;
	int M_PRN[40];

	int num_same = 0;
	for (i = 0; i < num_sat_old; i++)
	{
		for (j = 0; j < num_sat; j++)
		{
			if (PRN_old[i] == PRN_X[j])
			{
				M_PRN[num_same] = PRN_old[i];
				num_same++;
				break;
			}
		}
	}

	//先进行换星检查
	if (Ref != Ref_old)
	{
		bool iffind_old = false;
		bool iffind_new = false;
		for (i = 0; i < num_sat_old; i++)
		{
			if (Ref_old == PRN_old[i])
			{
				iffind_old = true;
			}
			if (Ref == PRN_old[i])
			{
				iffind_new = true;
			}
		}
		if (iffind_old == false || iffind_new == false)
		{
			return false;
		}
	}

	CMatrix TMatrix1;
	TMatrix1.SetSize(num_sat_old - 1, num_sat_old - 1);
	for (int tx = 0; tx < num_sat_old - 1; tx++)
	{
		TMatrix1[tx][tx] = 1;
	}

	//先进行参考星更换
	if (Ref != Ref_old)
	{
		int match_old = -1;
		int match_new = -1;
		bool iffind_old = false;
		bool iffind_new = false;
		for (i = 0; i < num_sat_old; i++)
		{
			if (Ref_old == PRN_old[i])
			{
				match_old = i;
				iffind_old = true;
			}
			if (Ref == PRN_old[i])
			{
				match_new = i;
				iffind_new = true;
			}
		}
		if (iffind_old == false)
		{
			assert(false);
		}
		if (iffind_new == false)
		{
			Ref = Ref_old;    //如果新的参考卫星在上一个历元找不到，则维持参考卫星不变
		}
		else
		{
			CMatrix TT(num_sat_old - 1, num_sat_old - 1);
			CMatrix TC(num_sat_old - 1, num_sat_old - 1);
			for (i = 0; i < num_sat_old - 1; i++)
			{
				TT[i][i] = 1.0;
			}
			if (match_new < match_old)
			{
				for (i = 0; i < num_sat_old - 1; i++)
				{
					TT[i][match_new] = -1.0;
				}
				for (i = 0; i < match_new; i++)
				{
					TC[i][i] = 1.0;
				}
				for (i = match_new; i < match_old - 1; i++)
				{
					TC[i][i + 1] = 1.0;
				}
				TC[match_old - 1][match_new] = 1.0;
				for (i = match_old; i < num_sat_old - 1; i++)
				{
					TC[i][i] = 1.0;
				}
			}
			else if (match_new > match_old)
			{
				for (i = 0; i < num_sat_old - 1; i++)
				{
					TT[i][match_new - 1] = -1.0;
				}
				for (i = 0; i < match_old; i++)
				{
					TC[i][i] = 1.0;
				}
				TC[match_old][match_new - 1] = 1.0;
				for (i = match_old + 1; i < match_new; i++)
				{
					TC[i][i - 1] = 1.0;
				}
				for (i = match_new; i < num_sat_old - 1; i++)
				{
					TC[i][i] = 1.0;
				}
			}
			TMatrix1 = TC * TT;
		}
	}

	CMatrix TMatrix2;
	TMatrix2.SetSize(num_same - 1, num_sat_old - 1);
	if (num_same < num_sat_old)//卫星消失
	{
		CMatrix TT(num_same - 1, num_sat_old - 1);
		int t_col = 0;
		int t_raw = 0;
		for (i = 0; i < num_sat_old; i++)
		{
			if (PRN_old[i] == Ref)
			{
				continue;
			}
			bool iffind = false;
			for (j = 0; j < num_same; j++)
			{
				if (M_PRN[j] == PRN_old[i])
				{
					iffind = true;
					break;
				}
			}
			if (iffind == true)
			{
				TT[t_col][t_raw] = 1.0;
				t_col++;
				t_raw++;
			}
			else if (iffind == false)
			{
				t_raw++;
			}
		}
		if (t_col != num_same - 1)
		{
			assert(false);
		}
		if (t_raw != num_sat_old - 1)
		{
			assert(false);
		}

		TMatrix2 = TT * TMatrix1;
	}
	else if (num_same == num_sat_old)
	{
		TMatrix2 = TMatrix1;
	}
	else
	{
		assert(false);
	}

	CMatrix TMatrix3;
	TMatrix3.SetSize(num_sat - 1, num_sat_old - 1);
	if (num_same < num_sat)//卫星增加
	{
		CMatrix TT(num_sat - 1, num_same - 1);
		CMatrix MAX_Q(num_sat - 1, num_sat - 1);

		int t_col = 0;
		int t_raw = 0;
		for (i = 0; i < num_sat; i++)
		{
			if (PRN_X[i] == Ref)
			{
				continue;
			}
			bool iffind = false;
			for (j = 0; j < num_same; j++)
			{
				if (M_PRN[j] == PRN_X[i])
				{
					iffind = true;
					break;
				}
			}
			if (iffind == true)
			{
				TT[t_col][t_raw] = 1.0;

				t_col++;
				t_raw++;
			}
			else if (iffind == false)
			{
				for (int i = 0; i < num_sat - 1; i++)
				{
					if (i == t_col)
					{
						MAX_Q[t_col][t_col] = 1.0E10;
					}
				}
				t_col++;
			}
		}

		if (t_col != num_sat - 1)
		{
			assert(false);
		}
		if (t_raw != num_same - 1)
		{
			assert(false);
		}
		TMatrix3 = TT * TMatrix2;
	}
	else if (num_same == num_sat)
	{
		TMatrix3 = TMatrix2;
	}
	else
	{
		assert(false);
	}

	//返回转换矩阵+转换矩阵的行和列+新卫星协因数矩阵
	TT_Matrix = TMatrix3;
	return true;
}

//函数功能：拼接矩阵
void ProGRCTrans(CMatrix& TT_G, CMatrix& TT_C, CMatrix& TT_ALL)
{
	int Row_G = TT_G.Row;
	int Col_G = TT_G.Col;
	int Row_C = TT_C.Row;
	int Col_C = TT_C.Col;


	int Tot_Row = Row_G + Row_C;
	int Tot_Col = Col_G + Col_C;
	TT_ALL.SetSize(Tot_Row, Tot_Col);

	int i = 0;
	int j = 0;

	//GPS部分赋值给新的转换矩阵
	for (i = 0; i < Row_G; i++)
		for (j = 0; j < Col_G; j++)
			TT_ALL[i][j] = TT_G[i][j];

	//BDS部分赋值给新的转换矩阵
	for (i = 0; i < Row_C; i++)
		for (j = 0; j < Col_C; j++)
			TT_ALL[Row_G + i][Col_G + j] = TT_C[i][j];
}

//函数功能：初始化状态方程
void InitState(CMatrix& QX)
{
	int num_x = QX.Row;
	for (int idx = 0; idx < num_x; idx++)
	{
		if (idx <= 2)
		{
			QX[idx][idx] = 10.0 * 10.0;    //初始化位置
		}
		else
		{
			QX[idx][idx] = 10.0 * 10.0;  //初始化模糊度
		}
	}
}

//函数功能：更新协方差
void UpdataState(CMatrix& QX, CMatrix& X, int* SatNum, int* oldSatNum, int SatPRN[2][20], int oldSatPRN[2][20], int* RefSat, int* OldRef)
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

//函数功能：检测新卫星
void CheckNewSats(int* SatNum, int* oldSatNum, int SatPRN[2][20], int oldSatPRN[2][20], int NewSatPRN[2][20], int* NewSatNum)
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

//函数功能：Lambda算法固定模糊度
void Resamb_LAMBDA(CMatrix& Q, CMatrix& X, CMatrix& aXb, double& ratio, CMatrix& fixX)
{
	CMatrix Qbb, Qaa, Qab, Xa, Xb, fixXa, Q_LAM, X_LAM;
	int Num_X = X.Row - 3;
	Qbb.SetSize(3, 3);
	Qaa.SetSize(Num_X, Num_X);
	Qab.SetSize(3, Num_X);
	Q_LAM.SetSize(Num_X, Num_X);
	X_LAM.SetSize(Num_X, 1);
	fixXa.SetSize(Num_X, 1);
	Xa.SetSize(Num_X, 1);
	Xb.SetSize(3, 1);

	for (int i = 0; i < Num_X; i++)
	{
		for (int j = 0; j < Num_X; j++)
		{
			Q_LAM[i][j] = Q[i + 3][j + 3];
		}
		X_LAM[i][0] = X[i + 3][0];
		Xa[i][0] = X[i + 3][0];
	}
	for (int ii = 0; ii < 3; ii++)
	{
		Xb[ii][0] = X[ii][0];
	}

	for (int ii = 0; ii < 3; ii++)
	{
		for (int jj = 0; jj < 3; jj++)
		{
			Qbb[ii][jj] = Q[ii][jj];
		}
	}

	for (int ii = 0; ii < 3; ii++)
	{
		for (int jj = 0; jj < Num_X; jj++)
		{
			Qab[ii][jj] = Q[ii][jj + 3];
		}
	}

	Qaa = Q_LAM;
	CLambda lambda;
	lambda.afloat.SetSize(Num_X, 1);
	lambda.Qahat.SetSize(Num_X, Num_X);
	lambda.Qahat = Q_LAM;
	lambda.afloat = X_LAM;
	lambda.lambda2(lambda.afloat, lambda.Qahat);

	//固定解
	for (int ii = 0; ii < Num_X; ii++)
	{
		fixXa[ii][0] = lambda.afixed[ii][0];
	}

	aXb = Xb - Qab * Qaa.InvertGaussJordan() * (Xa - fixXa);
	Qbb = Qbb - Qab * Qaa.InvertGaussJordan() * Qab.T();
	ratio = lambda.ratio;
	fixX = fixXa;
}

//函数功能：误差改正
double Error_Correction(Sat& tempSat1, Sat& RefSat1, Sat& tempSat2, Sat& RefSat2) {
	//基站-非参考卫星
	double error_satclock1 = tempSat1.Sat_clock;
	double error_trop1 = tempSat1.Trop_Delay;
	double error_relat1 = tempSat1.Relat;
	double error_sagnac1 = tempSat1.Sagnac;
	double error_tgd1 = tempSat1.TGD;
	double error_antenna_height1 = tempSat1.Antenna_Height;
	double error_sat_antenna1 = tempSat1.Sat_Antenna;
	double test_geo1 = -error_satclock1 + error_trop1 - error_relat1 + error_sagnac1 + error_tgd1 - error_sat_antenna1;

	//基站-参考卫星
	double Rerror_satclock1 = RefSat1.Sat_clock;
	double Rerror_trop1 = RefSat1.Trop_Delay;
	double Rerror_relat1 = RefSat1.Relat;
	double Rerror_sagnac1 = RefSat1.Sagnac;
	double Rerror_tgd1 = RefSat1.TGD;
	double Rerror_antenna_height1 = RefSat1.Antenna_Height;
	double Rerror_sat_antenna1 = RefSat1.Sat_Antenna;
	double Rtest_geo1 = -Rerror_satclock1 + Rerror_trop1 - Rerror_relat1 + Rerror_sagnac1 + Rerror_tgd1 - Rerror_sat_antenna1;

	//参考站2-非参考卫星
	double error_satclock2 = tempSat2.Sat_clock;
	double error_trop2 = tempSat2.Trop_Delay;
	double error_relat2 = tempSat2.Relat;
	double error_sagnac2 = tempSat2.Sagnac;
	double error_tgd2 = tempSat2.TGD;
	double error_antenna_height2 = tempSat2.Antenna_Height;
	double error_sat_antenna2 = tempSat2.Sat_Antenna;
	double test_geo2 = -error_satclock2 + error_trop2 - error_relat2 + error_sagnac2 + error_tgd2 - error_sat_antenna2;

	//参考站2-参考卫星
	double Rerror_satclock2 = RefSat2.Sat_clock;
	double Rerror_trop2 = RefSat2.Trop_Delay;
	double Rerror_relat2 = RefSat2.Relat;
	double Rerror_sagnac2 = RefSat2.Sagnac;
	double Rerror_tgd2 = RefSat2.TGD;
	double Rerror_antenna_height2 = RefSat2.Antenna_Height;
	double Rerror_sat_antenna2 = RefSat2.Sat_Antenna;
	double Rtest_geo2 = -Rerror_satclock2 + Rerror_trop2 - Rerror_relat2 + Rerror_sagnac2 + Rerror_tgd2 - Rerror_sat_antenna2;

	//双差形式的几何项和误差
	const double DDerror = (test_geo2 - Rtest_geo2) - (test_geo1 - Rtest_geo1);
	return DDerror;
}

//函数功能：计算B,L矩阵
/*void Get_BL(Obs_epoch& Data, int k, int ref_gps, int tw, int num_n, double POSITION_X1, double POSITION_Y1, double POSITION_Z1, double POSITION_X2, double POSITION_Y2, double POSITION_Z2, CMatrix& B, CMatrix& L) {
	Sat tempSat1 = Data.sat1[k];
	Sat RefSat1 = Data.sat1[ref_gps];
	Sat tempSat2 = Data.sat2[k];
	Sat RefSat2 = Data.sat2[ref_gps];

	double length_O1, length_R1, length_O2, length_R2;
	double DDerror = 0;

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
	if(tempSat1.sattype=="G" && tempSat2.sattype=="G")
		B[tw + num_n][3 + tw] = lambda_L1;
	else
		B[tw + num_n][3 + tw] = lambda_L1_C;

	//误差改正
	Error_Correction(tempSat1, RefSat1, tempSat2, RefSat2, DDerror);

	//双差形式的几何项和误差
	double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + DDerror;
	double deltp1 = (tempSat2.data[0] - RefSat2.data[0]) - (tempSat1.data[0] - RefSat1.data[0]);  //伪距
	double deltL1 = 0;
	if (tempSat1.sattype == "G" && tempSat2.sattype == "G")
		double deltL1 = ((tempSat2.data[2] - RefSat2.data[2]) - (tempSat1.data[2] - RefSat1.data[2])) * lambda_L1; //载波
	else
		double deltL1 = ((tempSat2.data[2] - RefSat2.data[2]) - (tempSat1.data[2] - RefSat1.data[2])) * lambda_L1_C; //载波
	L[tw][0] = deltp1 - DDgeo;         //伪距在前
	L[tw + num_n][0] = deltL1 - DDgeo; //载波在后
}*/

//函数功能：通过读取spp文件，进行RTK解算
bool SPP_Kinematic_Pro(string filename1, string filename2) {
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
	string name = filename2.substr(0, 3) + "_RTK.txt";
	const char* output_filename = name.c_str();
	FILE* fp;
	fopen_s(&fp, output_filename, "w");

	int OldSatPRN[2][20] = { 0 };
	int OldRefPRN[2] = { -1,-1 };
	int OldSatNum[2] = { 0 };

	CMatrix Q, X;

	//根据流动站匹配基准站对应历元数据
	int posk = 0;   //用以记录基站数据的匹配
	for (int j = 0; j < sppfile2.liyuan_num; j++)   //流动站逐历元解算
	{
		Obs_epoch tempData;

		//流动站当前时刻
		Obs_epoch tempData2 = sppfile2.epoch[j];

		if (j >= 1) {
			Obs_epoch preData2 = sppfile2.epoch[j - 1];
			OnCycleDetect(tempData2, preData2);     //判断此此历元与上一个历元中的卫星载波测量值是否有周跳,如果有周跳就将此历元中的卫星记为不可用
		}

		//基准站数据相应时刻必须在前
		int i = 0;
		for (i = posk; i < sppfile1.liyuan_num; i++)
		{
			if (sppfile1.epoch[i].GPSTIME <= tempData2.GPSTIME)
				continue;
			else
				break;
		}

		//回溯一个历元：posk即为要找的最近的历元且满足posk的时刻小于tempData2的时间
		if (i > 0) posk = --i;
		else
		{
			posk = i;
			continue;
		}

		//判断时间差：两个是否相隔超过1s
		if (fabs(sppfile1.epoch[i].GPSTIME - tempData2.GPSTIME) >= 1)
		{
			printf("未匹配到基站数据\n");
			continue;
		}

		//如果基准站数据滞后
		if (sppfile1.epoch[i].GPSTIME > tempData2.GPSTIME)
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

		//分离GPS和BDS卫星，便于后续单系统、双系统处理
		Obs_epoch GPSData;
		Obs_epoch BDSData;

		int SatPRN[2][20] = { 0 };		//当前历元卫星号 依次存储
		int SatNum[2] = { 0 };	//当前历元卫星数目
		int RefPRN[2] = { 0 };  //参考卫星号

		int NewSatNum[2] = { 0 };
		int NewSatPRN[2][20] = { 0 };

		int k = 0;
		int m = 0;
		int gpsnum, bdsnum;
		gpsnum = bdsnum = 0;
		for (m = 0; m < tempData.sat_num; m++)
		{
			if (tempData.sat1[m].sattype == "G")
			{
				GPSData.sat1.push_back(tempData.sat1[m]);
				GPSData.sat2.push_back(tempData.sat2[m]);
				SatPRN[0][gpsnum] = tempData.sat1[m].numofsat;
				gpsnum++;
			}

			if (tempData.sat1[m].sattype == "C")
			{
				BDSData.sat1.push_back(tempData.sat1[m]);
				BDSData.sat2.push_back(tempData.sat2[m]);
				SatPRN[1][bdsnum] = tempData.sat1[m].numofsat;
				bdsnum++;
			}
		}

		SatNum[0] = gpsnum;		//当前历元的卫星数目 GPS/BDS
		SatNum[1] = bdsnum;

		GPSData.sat_num = GPSData.sat1.size();
		BDSData.sat_num = BDSData.sat1.size();
		tempData.GPSTIME = tempData2.GPSTIME;

		//暂时处理GPS和BDS都大于4颗卫星的数据，真正应用时一般两者相加大于5颗即可
		if (GPSData.sat_num < 4 || BDSData.sat_num < 4)
		{
			printf("GPS或BDS卫星数小于4颗\n");
			bInit = false;
			continue;
		}

		//采用站间星间二次差分模型，选择参考卫星
		int GPSRef = -1;
		int BDSRef = -2;
		GPSRef = OldRefPRN[0];
		BDSRef = OldRefPRN[1];
		SelectRef(GPSData, GPSRef);  //选择最大高度角卫星作为参考卫星 
		SelectRef(BDSData, BDSRef);  //卫星分开看

		RefPRN[0] = GPSRef;  //GPS参考卫星序号
		RefPRN[1] = BDSRef;	 //BDS参考卫星序号

		//当前历元GPS和BDS卫星序列
		int PRN_X_G[32];
		int PRN_X_B[20];

		memset(PRN_X_G, 0, sizeof(PRN_X_G));	//数列初始化
		memset(PRN_X_B, 0, sizeof(PRN_X_B));

		for (k = 0; k < GPSData.sat_num; k++)
			PRN_X_G[k] = GPSData.sat1[k].numofsat;
		for (k = 0; k < BDSData.sat_num; k++)
			PRN_X_B[k] = BDSData.sat1[k].numofsat;

		int num_X_G = GPSData.sat_num;	          //GPS卫星个数
		int num_X_B = BDSData.sat_num;	          //BDS卫星个数

		int numgps_x = num_X_G - 1;               //GPS双差方程个数
		int numbds_x = num_X_B - 1;               //BDS双差方程个数

		int num_x_p = (numgps_x + numbds_x) * 2;  //双系统融合双差方程个数  P1 + L1		载波观测方程加伪距观测方程

		int num_x = numgps_x + numbds_x + 3;      //待估参数个数：ΔX，ΔY，ΔZ，每颗卫星的模糊度
		int num_n = numgps_x + numbds_x;          //模糊度个数

		//定位解算相关的矩阵	
		//1-双系统
		CMatrix B(num_x_p, num_x);                //系数矩阵  
		CMatrix L(num_x_p, 1);                    //观测值矩阵
		CMatrix R(num_x_p, num_x_p);              //误差矩阵（逆阵即为权矩阵）
		CMatrix Qw(num_x, num_x);                 //动态噪声阵

		//位置的动态噪声阵，如果没有周跳，模糊度不变所以噪声为0
		for (int idx = 0; idx < 3; idx++) {
			Qw[idx][idx] = 30.0 * 30.0;
		}

		//2-单GPS
		CMatrix RG_P(numgps_x, numgps_x); //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RG_L(numgps_x, numgps_x); //误差矩阵（逆阵即为权矩阵） 载波

		//3-单BDS
		CMatrix RC_P(numbds_x, numbds_x); //误差矩阵（逆阵即为权矩阵） 伪距
		CMatrix RC_L(numbds_x, numbds_x); //误差矩阵（逆阵即为权矩阵） 载波

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

			//更新协方差阵
			UpdataState(Q, X, SatNum, OldSatNum, SatPRN, OldSatPRN, RefPRN, OldRefPRN);
		}

		vector<int> GPS_PRN;	//卫星序号存储容器
		vector<int> BDS_PRN;

		for (k = 0; k < GPSData.sat_num; k++)
			GPS_PRN.push_back(PRN_X_G[k]);
		for (k = 0; k < BDSData.sat_num; k++)
			BDS_PRN.push_back(PRN_X_B[k]);

		//伪距
		GetDDR_Peo(GPSData, GPS_PRN, GPSRef, RG_P);
		GetDDR_Peo(BDSData, BDS_PRN, BDSRef, RC_P);
		//载波
		GetDDR_Leo(GPSData, GPS_PRN, GPSRef, RG_L);
		GetDDR_Leo(BDSData, BDS_PRN, BDSRef, RC_L);

		//融合系统的R矩阵，分块对角矩阵
		//伪距
		for (int k1 = 0; k1 < numgps_x; k1++)
			for (int k2 = 0; k2 < numgps_x; k2++)
			{
				R[k1][k2] = RG_P[k1][k2];
			}

		for (int k1 = 0; k1 < numbds_x; k1++)
			for (int k2 = 0; k2 < numbds_x; k2++)
			{
				R[k1 + numgps_x][k2 + numgps_x] = RC_P[k1][k2];
			}

		//载波
		for (int k1 = 0; k1 < numgps_x; k1++)
			for (int k2 = 0; k2 < numgps_x; k2++)
			{
				R[k1 + numgps_x + numbds_x][k2 + numgps_x + numbds_x] = RG_L[k1][k2];
			}
		for (int k1 = 0; k1 < numbds_x; k1++)
			for (int k2 = 0; k2 < numbds_x; k2++)
			{
				R[k1 + numgps_x * 2 + numbds_x][k2 + numgps_x * 2 + numbds_x] = RC_L[k1][k2];
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
		int ref_gps = -1;
		int ref_bds = -1;
		for (k = 0; k < GPSData.sat_num; k++)
		{
			if (GPSData.sat1[k].numofsat == GPSRef)
			{
				ref_gps = k;
				break;
			}
		}

		for (k = 0; k < BDSData.sat_num; k++)
		{
			if (BDSData.sat1[k].numofsat == BDSRef)
			{
				ref_bds = k;
				break;
			}
		}

		//GPS
		int tp = 0;  //记录单系统方程数
		int tw = 0;  //记录双系统方程数
		double DDerror = 0;  //初始化误差改正项
		for (k = 0; k < GPSData.sat_num; k++)
		{
			bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = GPSData.sat1[k].numofsat;

			//参考卫星跳过
			if (PRN == GPSRef)
			{
				continue;
			}

			if (ref_gps < 0)
			{
				assert(0);
			}

			Sat tempSat1 = GPSData.sat1[k];
			Sat RefSat1 = GPSData.sat1[ref_gps];
			Sat tempSat2 = GPSData.sat2[k];
			Sat RefSat2 = GPSData.sat2[ref_gps];

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

		//BDS
		tp = 0;  //记录单系统方程数
		DDerror = 0;
		for (k = 0; k < BDSData.sat_num; k++)
		{
			bool Init = true;                      //记录是否有旧参考星被剔除
			int PRN = BDSData.sat1[k].numofsat;

			//参考卫星跳过
			if (PRN == BDSRef)
			{
				continue;
			}

			if (ref_bds < 0)
			{
				assert(0);
			}

			Sat tempSat1 = BDSData.sat1[k];
			Sat RefSat1 = BDSData.sat1[ref_bds];
			Sat tempSat2 = BDSData.sat2[k];
			Sat RefSat2 = BDSData.sat2[ref_bds];

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
				X[3 + tw][0] = (L[tw + num_n][0] - L[tw][0]) / lambda_L1_C;   //初始化模糊度
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

		if (GPSData.sat_num >= 4 && BDSData.sat_num >= 4)
		{
			//权矩阵PG、法方程矩阵NG(逆矩阵即为参数方差)、参数求解(相对概率坐标的改正数)

			//B系数矩阵  L观测值矩阵  R误差矩阵  Qw动态噪声矩阵  Q协方差矩阵  X位置和整周模糊度

			KalmanFilter(B, L, R, Qw, Q, X);

			double ratio = 1.0;

			CMatrix aXb(3, 1);
			CMatrix fix(num_n, 1);
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

			fprintf(fp, "EpNum: %2d, %2d, %2d, %10.6f, %2d, %2d,%15.8g,%15.8g, %15.8g, %15.8g\n", j + 1, tempData2.hour, tempData2.minute, tempData2.second, GPSData.sat_num, BDSData.sat_num, ratio,
				ErrorX_NEU[0][0], ErrorX_NEU[1][0], ErrorX_NEU[2][0]);
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
		bInit = true;  //初始化成功
	}
	fclose(fp);
	return true;
}

