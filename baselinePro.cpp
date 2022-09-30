#define _CRT_SECURE_NO_WARNINGS

//#include"stdafx.h"
#include"common.h"
#include"lambda.h"
//#include"Select_Sys.h"
#include<iostream>
#include "GPT2_1w_World.h"

#include<string>
using namespace std;
common coo;
/*
void GNSS::BDS23(vector<observe_ppp>baseline)//基于BDS-2和BDS-3的考虑系统间偏差的紧组合载波双差定位
{
	if (baseline.size() < 1)
	{
		cout << "共视文件为空";

		return;
	}

	size_t filenum = baseline.size();
	for (int a = 0; a < filenum; a++)
	{
		double BLH[3] = { 0 };
		coo.OnXYZtoBLH(baseline[a].rX, baseline[a].rY, baseline[a].rZ, BLH);
		double B = BLH[0];
		double L = BLH[1];
		double H = BLH[2];
		CMatrix TT(3, 3);
		TT[0][0] = -sin(B)*cos(L);    TT[1][0] = -sin(L);   TT[2][0] = cos(B)*cos(L);
		TT[0][1] = -sin(B)*sin(L);    TT[1][1] = cos(L);    TT[2][1] = cos(B)*sin(L);
		TT[0][2] = cos(B);            TT[1][2] = 0;         TT[2][2] = sin(B);
		FILE *p, *p1, *p2, *p3, *p4;
		CString s = "基于BDS-2和BDS-3的考虑系统间偏差的紧组合载波双差定位.txt";
		CString s1 = "模糊度N(0,-1,1).txt";
		CString s2 = "模糊度N(1,-1,0).txt";
		p = fopen((LPCSTR)s, "w");
		p1 = fopen((LPCSTR)s1, "w");
		p2 = fopen((LPCSTR)s2, "w");
		size_t epochnum = baseline[a].m_epoch.size();
		int refposkbds;
		pppsat refsatmS;//BDS-2
		pppsat refsatrS;
		pppsat refsatm3;//BDS-3
		pppsat refsatr3;
		ppp_epoch tempepoch;//存储上个历元
		ppp_epoch Bm;
		ppp_epoch Br;
		double posr[3]; posr[0] = baseline[a].rX; posr[1] = baseline[a].rY; posr[2] = baseline[a].rZ;
		double posm[3]; posm[0] = baseline[a].mX; posm[1] = baseline[a].mY; posm[2] = baseline[a].mZ;
		coo.check(baseline[a].m_epoch, baseline[a].r_epoch);
		for (int b = 0; b < epochnum; b++)//历元
		{
			//cout << "第"<<b+1<<"个历元" << endl;
			ppp_epoch BDSm;
			ppp_epoch BDSr;
			ppp_epoch CHm;
			ppp_epoch CHr;
			int bdsnum = 0;
			int chnum = 0;//bds-3

			for (int i = 0; i < baseline[a].m_epoch[b].sat.size(); i++)
			{
				if (baseline[a].m_epoch[b].sat[i].judge_use == 0)
				{
					
					if (baseline[a].m_epoch[b].sat[i].sattype == "C")
					{
						BDSm.sat.push_back(baseline[a].m_epoch[b].sat[i]);
						BDSr.sat.push_back(baseline[a].r_epoch[b].sat[i]);
						bdsnum++;
					}
					if (baseline[a].m_epoch[b].sat[i].sattype == "B")
					{
						CHm.sat.push_back(baseline[a].m_epoch[b].sat[i]);
						CHr.sat.push_back(baseline[a].r_epoch[b].sat[i]);
						chnum++;
					}
				}
			} 
			//bds-2
			CMatrix P(bdsnum - 1, 1);//伪距
			CMatrix fai0_11(bdsnum - 1, 1);
			CMatrix fai1_10(bdsnum - 1, 1);
			CMatrix N0_11(bdsnum - 1, 1);
			CMatrix N1_10(bdsnum - 1, 1);
			CMatrix BS(bdsnum - 1, 4);
			CMatrix L1(bdsnum - 1, 1);
			CMatrix LS(bdsnum - 1, 1);
			CMatrix PS(bdsnum - 1, bdsnum - 1);

			CMatrix IS(bdsnum - 1, bdsnum - 1);
			CMatrix FS(bdsnum - 1, bdsnum - 1);
			IS.ones(bdsnum - 1);
			FS.ones(bdsnum - 1);
			CMatrix RS(bdsnum - 1, bdsnum - 1);
			CMatrix QwS(bdsnum - 1, bdsnum - 1);

			pppsat tempsatmS = BDSm.sat[0];//r,m的sattype和prn号和卫星位置都相等  所以用sattype或prn或APP_XYZ判断时用r或m都可
			pppsat tempsatrS = BDSr.sat[0];
			pppsat tempsatm3 = CHm.sat[0];//r,m的sattype和prn号和卫星位置都相等  所以用sattype或prn或APP_XYZ判断时用r或m都可
			pppsat tempsatr3 = CHr.sat[0];
			double fbds1 = 1561.098*1E+6;
			double fbds2 = 1207.14*1E+6;
			double fbds3 = 1268.52*1E+6;
			double fch1 = 1575.42*1E+6;
			double fch2 = 1561.098*1E+6;
			//double fch3 = 1268.52*1E+6;
			double lambda011 = C_light / (fbds2 + fbds3);
			double lambda0_11 = C_light / (-fbds2 + fbds3);
			double lambda1_10 = C_light / (fbds1 - fbds2);
			double lambdaC = C_light / (fbds1 - fbds2);
			double lambdaB = C_light / (fch1 - fch2);
			double yita011 = fbds1 * fbds1*(1.0 / fbds2 + 1.0 / fbds3) / (fbds2 + fbds3);
			double yita0_11 = fbds1 * fbds1*(-1.0 / fbds2 + 1.0 / fbds3) / (-fbds2 + fbds3);
			double yita1_10 = fbds1 * fbds1*(1.0 / fbds1 + (-1.0) / fbds2) / (fbds1 - fbds2);
			//bds-2参考星选择
			for (int i = 0; i < bdsnum; i++)
			{
				if (BDSm.sat[i].E > tempsatmS.E)
				{
					tempsatmS = BDSm.sat[i];
					tempsatrS = BDSr.sat[i];
				}
			}


			if (b == 0)//第一个历元
			{
				refsatmS = tempsatmS;
				refsatrS = tempsatrS;
			}
			else//其他历元
			{
				bool find = false;
				for (int i = 0; i < tempepoch.sat.size(); i++)
				{
					if (tempepoch.sat[i].sattype == tempsatmS.sattype&&tempepoch.sat[i].PRN == tempsatmS.PRN)//若参考星在上一历元中未出现，继续使用上一历元的参考星
						find = true;
				}
				if (find == true)
				{
					refsatmS = tempsatmS;
					refsatrS = tempsatrS;
				}

			}//参考星选择完毕
			
			//bds-3参考星选择
			for (int i = 0; i < chnum; i++)
			{
				if (CHm.sat[i].E > tempsatm3.E)
				{
					tempsatm3 = CHm.sat[i];
					tempsatr3 = CHr.sat[i];
				}
			}


			if (b == 0)//第一个历元
			{
				refsatm3 = tempsatm3;
				refsatr3 = tempsatr3;
			}
			else//其他历元
			{
				bool find = false;
				for (int i = 0; i < tempepoch.sat.size(); i++)
				{
					if (tempepoch.sat[i].sattype == tempsatm3.sattype&&tempepoch.sat[i].PRN == tempsatm3.PRN)//若参考星在上一历元中未出现，继续使用上一历元的参考星
						find = true;
				}
				if (find == true)
				{
					refsatm3 = tempsatm3;
					refsatr3 = tempsatr3;
				}

			}//参考星选择完毕

			int tw = 0;
			for (int i = 0; i < bdsnum; i++)
			{
				if (refsatmS.PRN == BDSm.sat[i].PRN)
					continue;
				P[tw][0] = (fbds2* ((BDSr.sat[i].data[1] - BDSm.sat[i].data[1]) - (refsatrS.data[1] - refsatmS.data[1])) + fbds3 * ((BDSr.sat[i].data[4] - BDSm.sat[i].data[4]) - (refsatrS.data[4] - refsatmS.data[4]))) / (fbds2 + fbds3);
				fai0_11[tw][0] = (-C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3])) + C_light * ((BDSr.sat[i].data[5] - BDSm.sat[i].data[5]) - (refsatrS.data[5] - refsatmS.data[5]))) / (-fbds2 + fbds3);
				fai1_10[tw][0] = (C_light * ((BDSr.sat[i].data[2] - BDSm.sat[i].data[2]) - (refsatrS.data[2] - refsatmS.data[2])) - C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3]))) / (fbds1 - fbds2);
				N0_11[tw][0] = round((P[tw][0] - fai0_11[tw][0]) / lambda0_11);
				double TGD = (BDSr.sat[i].TGD - BDSm.sat[i].TGD) - (refsatrS.TGD - refsatmS.TGD);
				double TROP = (BDSr.sat[i].Trop_Delay - BDSm.sat[i].Trop_Delay) - (refsatrS.Trop_Delay - refsatmS.Trop_Delay);
				N1_10[tw][0] = round(-1.0 / lambda1_10 * (fai1_10[tw][0] - fai0_11[tw][0] - (yita0_11 - yita1_10)*TGD - lambda0_11 * N0_11[tw][0]));
				double fai = lambdaC * ((BDSr.sat[i].data[2] - BDSm.sat[i].data[2]) - (BDSr.sat[i].data[3] - BDSm.sat[i].data[3])) - lambdaB * ((refsatr3.data[7] - refsatm3.data[7]) - (refsatr3.data[2] - refsatm3.data[2]));

				double length0 = sqrt(pow(refsatr3.POS_X - baseline[a].rX, 2) + pow(refsatr3.POS_Y - baseline[a].rY, 2) + pow(refsatr3.POS_Z - baseline[a].rZ, 2));//3代参考卫星到流动站的距离
				double length1 = sqrt(pow(BDSr.sat[i].POS_X - baseline[a].rX, 2) + pow(BDSr.sat[i].POS_Y - baseline[a].rY, 2) + pow(BDSr.sat[i].POS_Z - baseline[a].rZ, 2));//2代某颗卫星到流动站的距离
				double length2 = sqrt(pow(refsatm3.POS_X - baseline[a].mX, 2) + pow(refsatm3.POS_Y - baseline[a].mY, 2) + pow(refsatm3.POS_Z - baseline[a].mZ, 2));//3代参考卫星到基站的距离
				double length3 = sqrt(pow(BDSm.sat[i].POS_X - baseline[a].mX, 2) + pow(BDSm.sat[i].POS_Y - baseline[a].mY, 2) + pow(BDSm.sat[i].POS_Z - baseline[a].mZ, 2));//2代某颗卫星到基站的距离
				BS[tw][0] = (BDSr.sat[i].POS_X - baseline[a].rX) / length1 - (refsatr3.POS_X - baseline[a].rX) / length0;
				BS[tw][1] = (BDSr.sat[i].POS_Y - baseline[a].rY) / length1 - (refsatr3.POS_Y - baseline[a].rY) / length0;
				BS[tw][2] = (BDSr.sat[i].POS_Z - baseline[a].rZ) / length1 - (refsatr3.POS_Z - baseline[a].rZ) / length0;
				BS[tw][3] = lambdaB;
				LS[tw][0] = fai - lambdaB * N1_10[tw][0] - ((length1 - length3) - (length0 - length2)) ;
				PS[tw][tw] = 0.09 + 0.09 / pow(sin(refsatr3.E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(BDSm.sat[i].E*PI / 180), 2);

				tw++;
			}

			CMatrix NN0_11(12, 1);
			for (int i = 0; i < bdsnum - 1; i++)
			{
				NN0_11[i][0] = N0_11[i][0];
			}
			fprintf(p1, "%14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f,\n", NN0_11[0][0], NN0_11[1][0], NN0_11[2][0], NN0_11[3][0], NN0_11[4][0], NN0_11[5][0], NN0_11[6][0], NN0_11[7][0], NN0_11[8][0], NN0_11[9][0], NN0_11[10][0], NN0_11[11][0]);

			CMatrix NN1_10(12, 1);
			for (int i = 0; i < bdsnum - 1; i++)
			{
				NN1_10[i][0] = N1_10[i][0];
			}
			fprintf(p2, "%14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f\n", NN1_10[0][0], NN1_10[1][0], NN1_10[2][0], NN1_10[3][0], NN1_10[4][0], NN1_10[5][0], NN1_10[6][0], NN1_10[7][0], NN1_10[8][0], NN1_10[9][0], NN1_10[10][0], NN1_10[11][0]);
			//最小二乘
			CMatrix XS = (BS.T()*PS*BS).InvertGaussJordan()*BS.T()*PS*LS;
			CMatrix NEUS = TT * XS;
			fprintf(p, " XYZ: %14.4f, %14.4f, %14.4f,  NEU: %14.4f, %14.4f, %14.4f, 系统间偏差：%14.4f,\n", XS[0][0], XS[1][0], XS[2][0],  NEUS[0][0], NEUS[1][0], NEUS[2][0], XS[3][0]);
			tempepoch = baseline[a].m_epoch[b];//存储上个历元
			Br = BDSr;
			Bm = BDSm;

			TRACE("第%d个历元生成结束\n", b + 1);
		}
		fclose(p);
		fclose(p1);
		fclose(p2);
	}
}
*/
/*******************************************************************************************************************************************************/
/*******************************************************************************************************************************************************/
/*******************************************************************************************************************************************************/
void GNSS::BDS23(vector<observe_ppp>baseline)//基于BDS-2和BDS-3的考虑系统间偏差的紧组合载波双差定位
{
	size_t filenum = baseline.size();
	for (int a = 0; a < filenum; a++)
	{
		double BLH[3] = { 0 };
		coo.OnXYZtoBLH(baseline[a].rX, baseline[a].rY, baseline[a].rZ, BLH);
		double B = BLH[0];
		double L = BLH[1];
		double H = BLH[2];
		CMatrix TT(3, 3);
		TT[0][0] = -sin(B)*cos(L);    TT[1][0] = -sin(L);   TT[2][0] = cos(B)*cos(L);
		TT[0][1] = -sin(B)*sin(L);    TT[1][1] = cos(L);    TT[2][1] = cos(B)*sin(L);
		TT[0][2] = cos(B);            TT[1][2] = 0;         TT[2][2] = sin(B);

		pppsat refsatmC;
		pppsat refsatrC;
		pppsat refsatmS;
		pppsat refsatrS;
		ppp_epoch tempepoch;//存储上个历元
		ppp_epoch Cm;
		ppp_epoch Cr;
		ppp_epoch Bm;
		ppp_epoch Br;

		FILE *p, *p1, *p2, *p3, *p4,*C1, *C2, *C3, *C4,*C5, *C6, *C7, *C8, *C9, *C10, *C11, *C12, *C13,*C14, *B19, *B20, *B21, *B22, *B23, *B24, *B25, *B26, *B27, *B28, *B29, *B30, *B32, *B33, *B34, *B35, *B36, *B37;;
		CString ss1 = baseline[a].marker_name;
		ss1 = ss1.Left(4);
		CString s1;
		s1.Format("短基线紧组合a%s.txt", ss1);
		

		CString s = "短基线模糊度NWC.txt";
		
		CString s2 = "短基线模糊度NW.txt";
		CString s3 = "xyz.txt";
		CString s4 = "xyz1.txt";
		CString c1 = "短基线C01.txt";
		CString c2 = "短基线C02.txt";
		CString c3 = "短基线C03.txt";
		CString c4 = "短基线C04.txt";
		CString c5 = "短基线C05.txt";
		CString c6 = "短基线C06.txt";
		CString c7 = "短基线C07.txt";
		CString c8 = "短基线C08.txt";
		CString c9 = "短基线C09.txt";
		CString c10 = "短基线C10.txt";
		CString c11 = "短基线C11.txt";
		CString c12 = "短基线C12.txt";
		CString c13 = "短基线C13.txt";
		CString c14 = "短基线C14.txt";

		CString b19 = "短基线B19.txt";
		CString b20 = "短基线B20.txt";
		CString b21 = "短基线B21.txt";
		CString b22 = "短基线B22.txt";
		CString b23 = "短基线B23.txt";
		CString b24 = "短基线B24.txt";
		CString b25 = "短基线B25.txt";
		CString b26 = "短基线B26.txt";
		CString b27 = "短基线B27.txt";
		CString b28 = "短基线B28.txt";
		CString b29 = "短基线B29.txt";
		CString b30 = "短基线B30.txt";
		CString b32 = "短基线B32.txt";
		CString b33 = "短基线B33.txt";
		CString b34 = "短基线B34.txt";
		CString b35 = "短基线B35.txt";
		CString b36 = "短基线B36.txt";
		CString b37 = "短基线B37.txt";
		p = fopen((LPCSTR)s, "w");
		p1 = fopen((LPCSTR)s1, "w");
		p2 = fopen((LPCSTR)s2, "w");
		p3 = fopen((LPCSTR)s3, "w");
		p4 = fopen((LPCSTR)s4, "w");

		C1 = fopen((LPCSTR)c1, "w");     C2 = fopen((LPCSTR)c2, "w");      C3 = fopen((LPCSTR)c3, "w");
		C4 = fopen((LPCSTR)c4, "w");     C6 = fopen((LPCSTR)c6, "w");      C7 = fopen((LPCSTR)c7, "w");
		C8 = fopen((LPCSTR)c8, "w");     C9 = fopen((LPCSTR)c9, "w");      C10 = fopen((LPCSTR)c10, "w");
		C11 = fopen((LPCSTR)c11, "w");   C12 = fopen((LPCSTR)c12, "w");    C13 = fopen((LPCSTR)c13, "w");
		C5= fopen((LPCSTR)c5, "w");      C14 = fopen((LPCSTR)c14, "w");


		B19 = fopen((LPCSTR)b19, "w");    B20 = fopen((LPCSTR)b20, "w");     B21 = fopen((LPCSTR)b21, "w");
		B22 = fopen((LPCSTR)b22, "w");    B23 = fopen((LPCSTR)b23, "w");     B24 = fopen((LPCSTR)b24, "w");
		B25 = fopen((LPCSTR)b25, "w");    B26 = fopen((LPCSTR)b26, "w");     B27 = fopen((LPCSTR)b27, "w");
		B28 = fopen((LPCSTR)b28, "w");    B29 = fopen((LPCSTR)b29, "w");     B30 = fopen((LPCSTR)b30, "w");
		B32 = fopen((LPCSTR)b32, "w");    B33 = fopen((LPCSTR)b33, "w");     B34 = fopen((LPCSTR)b34, "w");
		B35 = fopen((LPCSTR)b35, "w");    B36 = fopen((LPCSTR)b36, "w");     B37 = fopen((LPCSTR)b37, "w");

		CMatrix X, XWC;
		CMatrix Q, QWC;
		CMatrix  XWS;
		CMatrix  QWS;
		size_t epochnum = baseline[a].m_epoch.size();
		int refposkchs;//存储上历元参考星的序号
		int refposkbds;
		double posr[3]; posr[0] = baseline[a].rX; posr[1] = baseline[a].rY; posr[2] = baseline[a].rZ;
		double posm[3]; posm[0] = baseline[a].mX; posm[1] = baseline[a].mY; posm[2] = baseline[a].mZ;
		coo.check(baseline[a].m_epoch, baseline[a].r_epoch);
		double RMSN = 0;
		double RMSE = 0;
		double RMSU = 0;
		for (int b = 0; b < epochnum; b++)//历元
		{
			//cout << "第"<<b+1<<"个历元" << endl;
			ppp_epoch BDSm;
			ppp_epoch BDSr;
			ppp_epoch CHm;
			ppp_epoch CHr;

			int bdsnum = 0;
			int chnum = 0;//bds-3

			for (int i = 0; i < baseline[a].m_epoch[b].sat.size(); i++)
			{
				if (baseline[a].m_epoch[b].sat[i].judge_use == 0)
				{
					
					if (baseline[a].m_epoch[b].sat[i].sattype == "C"&&baseline[a].m_epoch[b].sat[i].E>15)
					{
						BDSm.sat.push_back(baseline[a].m_epoch[b].sat[i]);
						BDSr.sat.push_back(baseline[a].r_epoch[b].sat[i]);
						bdsnum++;
					}
					if (baseline[a].m_epoch[b].sat[i].sattype == "B"&&baseline[a].m_epoch[b].sat[i].E >15)
					{
						CHm.sat.push_back(baseline[a].m_epoch[b].sat[i]);
						CHr.sat.push_back(baseline[a].r_epoch[b].sat[i]);
						chnum++;
					}
					
				}
			} 

			if (chnum < 3)
				continue;
			//BDS-3
			
			CMatrix PC(chnum - 1, chnum - 1);
			CMatrix IC(chnum + 2, chnum + 2);
			CMatrix FC(chnum + 2, chnum + 2);
			IC.ones(chnum + 2);
			FC.ones(chnum + 2);
			CMatrix RC(chnum - 1, chnum - 1);
			CMatrix QwC(chnum + 2, chnum + 2);
			CMatrix QWwC(chnum + 2, chnum + 2);

			//BDS-23
			CMatrix B(bdsnum+chnum - 1, 4);
			CMatrix L(bdsnum + chnum - 1, 1);
			CMatrix P(bdsnum + chnum - 1, bdsnum + chnum - 1);
			CMatrix I(4, 4);
			CMatrix F(4, 4);
			I.ones(4);
			F.ones(4);
			CMatrix R(bdsnum + chnum - 1, bdsnum + chnum - 1);
			//CMatrix B(bdsnum -1, 3);
			//CMatrix L(bdsnum - 1, 1);
			//CMatrix P(bdsnum -1, bdsnum -1);




			//计算BDS-3的系统内双差模糊度NWC
			pppsat tempsatmC = CHm.sat[0];//r,m的sattype和prn号和卫星位置都相等  所以用sattype或prn或APP_XYZ判断时用r或m都可
			pppsat tempsatrC = CHr.sat[0];
			double fch1 = 1575.42*1E+6;
			double fch2 = 1561.098*1E+6;
			for (int i = 0; i < chnum; i++)
			{
				if (CHm.sat[i].E > tempsatmC.E)
				{
					tempsatmC = CHm.sat[i];
					tempsatrC = CHr.sat[i];
				}
			}


			if (b == 0)//第一个历元
			{
				refsatmC = tempsatmC;
				refsatrC = tempsatrC;
				//卡尔曼滤波参数初始化
				X.SetSize(4, 1);
				Q.SetSize(4, 4);
				for (int ii = 0; ii < 3; ii++)
				{
					Q[ii][ii] = 10.0*10.0;
				}
				Q[4][4] = 1000 * 1000;
				QWC.SetSize(chnum + 2, chnum + 2);
				XWC.SetSize(chnum + 2, 1);

				for (int j = 0; j < 3; j++)
				{
					QWC[j][j] = 30.0 * 30.0;
				}

			}
			else//其他历元
			{
				CMatrix TT_matrixC(chnum - 1, Cr.sat.size() - 1);
				CMatrix matrixC(chnum + 2, Cr.sat.size() + 2);
				QwC.SetSize(chnum + 2, chnum + 2);
				bool find = false;
				for (int i = 0; i < tempepoch.sat.size(); i++)
				{
					if (tempepoch.sat[i].sattype == tempsatmC.sattype&&tempepoch.sat[i].PRN == tempsatmC.PRN)//若参考星在上一历元中未出现，继续使用上一历元的参考星
						find = true;
				}
				if (find == true)
				{
					refsatmC = tempsatmC;
					refsatrC = tempsatrC;
				}
				coo.GetGRCQValue(CHm, Cm, chnum, Cm.sat.size(), refsatmC.PRN, Cm.sat[refposkchs].PRN, TT_matrixC);


				matrixC[0][0] = 1.0;  matrixC[1][1] = 1.0;  matrixC[2][2] = 1.0;
				for (int ii = 0; ii < chnum - 1; ii++)
				{
					for (int jj = 0; jj < Cr.sat.size() - 1; jj++)
					{
						matrixC[ii + 3][jj + 3] = TT_matrixC[ii][jj];
					}
				}
				
				QWC = matrixC * QWC*matrixC.T();
				XWC = matrixC * XWC;

			}//参考星选择完毕

			//存储上历元的参考星的序号
			for (int i = 0; i < chnum; i++)
			{
				if (CHm.sat[i].PRN == refsatmC.PRN)
				{
					refposkchs = i;
					break;
				}
			}
			CMatrix LWC(chnum - 1, 1);
			CMatrix BWC(chnum - 1, chnum + 2);

			coo.MW(b, posm, posr, CHm, CHr, refsatmC, refsatrC, Cr, XWC, LWC, BWC, QWC, PC);
			RC = PC.InvertGaussJordan();
			coo.Kalman(IC, FC, BWC, LWC, RC, QWC, QWwC, XWC);

			CLambda lamNC;
			lamNC.numLiyuan = b + 1;//当前历元个数
			lamNC.N = chnum - 1;//模糊度个数
			CMatrix X_LAM(chnum - 1, 1);
			CMatrix Q_LAM(chnum - 1, chnum - 1);
			for (int i = 0; i < chnum - 1; i++)
			{
				X_LAM[i][0] = XWC[i + 3][0];
				Q_LAM[i][i] = QWC[i + 3][i + 3];
			}
			lamNC.a = X_LAM;
			lamNC.Qa = Q_LAM;
			lamNC.OnFixAmbigultyMain();
			double ratioNWC = lamNC.Ratio;
			CMatrix NWC(chnum - 1, 1);
			NWC = lamNC.afix;
			CMatrix NNWC(12, 1);
			for (int i = 0; i < chnum - 1; i++)
				NNWC[i][0] = NWC[i][0];

			fprintf(p, "固定解NW：chnum-1：%4d ratio:%7.3f  %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f  \n", chnum - 1, ratioNWC, NNWC[0][0], NNWC[1][0], NNWC[2][0], NNWC[3][0], NNWC[4][0], NNWC[5][0], NNWC[6][0], NNWC[7][0], NNWC[8][0], NNWC[9][0], NNWC[10][0], NNWC[11][0]);
			

			/**********************************统计三代每一颗卫星的模糊度及偏差*********************************/
			int e = 0;
			for (int i = 0; i < chnum; i++)
			{
				if (refsatmC.PRN == CHm.sat[i].PRN)
					continue;
				if (CHr.sat[i].PRN == 19)
					fprintf(B19, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 20)
					fprintf(B20, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 21)
					fprintf(B21, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 22)
					fprintf(B22, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 23)
					fprintf(B23, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 24)
					fprintf(B24, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 25)
					fprintf(B25, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 26)
					fprintf(B26, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 27)
					fprintf(B27, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 28)
					fprintf(B28, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 29)
					fprintf(B29, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 30)
					fprintf(B30, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 32)
					fprintf(B32, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 33)
					fprintf(B33, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 34)
					fprintf(B34, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 35)
					fprintf(B35, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 36)
					fprintf(B36, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				if (CHr.sat[i].PRN == 37)
					fprintf(B37, " %14.4f  %14.4f\n ", NWC[e][0], X_LAM[e][0]);
				e++;
			}


			/********************************************************************************************************/


			CMatrix PP(bdsnum - 1, 1);//伪距
			CMatrix fai0_11(bdsnum - 1, 1);
			CMatrix fai1_10(bdsnum - 1, 1);
			CMatrix N0_11(bdsnum - 1, 1);
			CMatrix N1_10(bdsnum - 1, 1);
			CMatrix doubleN0_11(bdsnum - 1, 1);
			CMatrix doubleN1_10(bdsnum - 1, 1);

			pppsat tempsatmS = BDSm.sat[0];//r,m的sattype和prn号和卫星位置都相等  所以用sattype或prn或APP_XYZ判断时用r或m都可
			pppsat tempsatrS = BDSr.sat[0];

			double fbds1 = 1561.098*1E+6;
			double fbds2 = 1207.14*1E+6;
			double fbds3 = 1268.52*1E+6;
			
			double lambda011 = C_light / (fbds2 + fbds3);
			double lambda0_11 = C_light / (-fbds2 + fbds3);
			double lambda1_10 = C_light / (fbds1 - fbds2);
			double lambdaC = C_light / (fbds1 - fbds2);
			double lambdaB = C_light / (fch1 - fch2);
			double yita011 = fbds1 * fbds1*(1.0 / fbds2 + 1.0 / fbds3) / (fbds2 + fbds3);
			double yita0_11 = fbds1 * fbds1*(-1.0 / fbds2 + 1.0 / fbds3) / (-fbds2 + fbds3);
			double yita1_10 = fbds1 * fbds1*(1.0 / fbds1 + (-1.0) / fbds2) / (fbds1 - fbds2);
			double yitaC = yita1_10;//BDS2
			double yitaB= fch1 * fch1*(1.0 / fch1 + (-1.0) / fch2) / (fch1 - fch2);
			for (int i = 0; i < bdsnum; i++)
			{
				if (BDSm.sat[i].E > tempsatmS.E)
				{
					tempsatmS = BDSm.sat[i];
					tempsatrS = BDSr.sat[i];
				}
			}


			if (b == 0)//第一个历元
			{
				refsatmS = tempsatmS;
				refsatrS = tempsatrS;
			}
			else//其他历元
			{
				bool find = false;
				for (int i = 0; i < tempepoch.sat.size(); i++)
				{
					if (tempepoch.sat[i].sattype == tempsatmS.sattype&&tempepoch.sat[i].PRN == tempsatmS.PRN)//若参考星在上一历元中未出现，继续使用上一历元的参考星
						find = true;
				}
				if (find == true)
				{
					refsatmS = tempsatmS;
					refsatrS = tempsatrS;
				}

			}//参考星选择完毕
			//存储上历元的参考星的序号
			for (int i = 0; i < bdsnum; i++)
			{
				if (BDSm.sat[i].PRN == refsatmS.PRN)
				{
					refposkbds = i;
					break;
				}
			}
			int tw = 0;
			for (int i = 0; i < bdsnum; i++)
			{
				if (refsatmS.PRN == BDSm.sat[i].PRN)
					continue;
				PP[tw][0] = (fbds2* ((BDSr.sat[i].data[1] - BDSm.sat[i].data[1]) - (refsatrS.data[1] - refsatmS.data[1])) + fbds3 * ((BDSr.sat[i].data[4] - BDSm.sat[i].data[4]) - (refsatrS.data[4] - refsatmS.data[4]))) / (fbds2 + fbds3);
				fai0_11[tw][0] = (-C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3])) + C_light * ((BDSr.sat[i].data[5] - BDSm.sat[i].data[5]) - (refsatrS.data[5] - refsatmS.data[5]))) / (-fbds2 + fbds3);
				fai1_10[tw][0] = (C_light * ((BDSr.sat[i].data[2] - BDSm.sat[i].data[2]) - (refsatrS.data[2] - refsatmS.data[2])) - C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3]))) / (fbds1 - fbds2);
				N0_11[tw][0] = round((PP[tw][0] - fai0_11[tw][0]) / lambda0_11);
				doubleN0_11[tw][0] = (PP[tw][0] - fai0_11[tw][0]) / lambda0_11 - N0_11[tw][0];
				double TGD = (BDSr.sat[i].TGD - BDSm.sat[i].TGD) - (refsatrS.TGD - refsatmS.TGD);
				double TROP = (BDSr.sat[i].Trop_Delay - BDSm.sat[i].Trop_Delay) - (refsatrS.Trop_Delay - refsatmS.Trop_Delay);
				N1_10[tw][0] = round(-1.0 / lambda1_10 * (fai1_10[tw][0] - fai0_11[tw][0] - (yita0_11 - yita1_10)*TGD - lambda0_11 * N0_11[tw][0]));
				doubleN1_10[tw][0] = -1.0 / lambda1_10 * (fai1_10[tw][0] - fai0_11[tw][0] - (yita0_11 - yita1_10)*TGD - lambda0_11 * N0_11[tw][0]) - N1_10[tw][0];

				double length0 = sqrt(pow(refsatrS.POS_X - baseline[a].rX, 2) + pow(refsatrS.POS_Y - baseline[a].rY, 2) + pow(refsatrS.POS_Z - baseline[a].rZ, 2));//参考卫星到流动站的距离
				double length1 = sqrt(pow(BDSr.sat[i].POS_X - baseline[a].rX, 2) + pow(BDSr.sat[i].POS_Y - baseline[a].rY, 2) + pow(BDSr.sat[i].POS_Z - baseline[a].rZ, 2));//某颗卫星到流动站的距离
				double length2 = sqrt(pow(refsatmS.POS_X - baseline[a].mX, 2) + pow(refsatmS.POS_Y - baseline[a].mY, 2) + pow(refsatmS.POS_Z - baseline[a].mZ, 2));//参考卫星到基站的距离
				double length3 = sqrt(pow(BDSm.sat[i].POS_X - baseline[a].mX, 2) + pow(BDSm.sat[i].POS_Y - baseline[a].mY, 2) + pow(BDSm.sat[i].POS_Z - baseline[a].mZ, 2));//某颗卫星到基站的距离
				B[tw][0] = (BDSr.sat[i].POS_X - baseline[a].rX) / length1 - (refsatrS.POS_X - baseline[a].rX) / length0;
				B[tw][1] = (BDSr.sat[i].POS_Y - baseline[a].rY) / length1 - (refsatrS.POS_Y - baseline[a].rY) / length0;
				B[tw][2] = (BDSr.sat[i].POS_Z - baseline[a].rZ) / length1 - (refsatrS.POS_Z - baseline[a].rZ) / length0;
				L[tw][0] = fai1_10[tw][0] + lambda1_10 * N1_10[tw][0] - ((length1 - length3) - (length0 - length2)) + yita1_10 * TGD - TROP;
				P[tw][tw] = 0.09 + 0.09 / pow(sin(refsatrS.E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(BDSm.sat[i].E*PI / 180), 2);

				tw++;
			}
			
			/**********************************统计二代每一颗卫星的模糊度及偏差*********************************/
			int n = 0;
			for (int i = 0; i < bdsnum; i++)
			{
				if (refsatmS.PRN == BDSm.sat[i].PRN)
					continue;
				if (BDSr.sat[i].PRN == 1)
					fprintf(C1, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 2)
					fprintf(C2, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 3)
					fprintf(C3, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 4)
					fprintf(C4, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 5)
					fprintf(C5, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 6)
					fprintf(C6, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 7)
					fprintf(C7, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 8)
					fprintf(C8, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 9)
					fprintf(C9, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 10)
					fprintf(C10, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 11)
					fprintf(C11, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 12)
					fprintf(C12, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 13)
					fprintf(C13, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				if (BDSr.sat[i].PRN == 14)
					fprintf(C14, "N0-11  %14.4f  %14.4f  N1_10  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N1_10[n][0], doubleN1_10[n][0]);
				n++;
			}


			/********************************************************************************************************/

			for (int i = 0; i < chnum; i++)
			{
				
				if (refsatmC.PRN == CHm.sat[i].PRN)
					continue;
				double fai = lambdaB * ((CHr.sat[i].data[7] - CHm.sat[i].data[7]) - (CHr.sat[i].data[2] - CHm.sat[i].data[2])) - lambdaC * ((refsatrS.data[2] - refsatmS.data[2]) - (refsatrS.data[3] - refsatmS.data[3]));

				double length0 = sqrt(pow(refsatrS.POS_X - baseline[a].rX, 2) + pow(refsatrS.POS_Y - baseline[a].rY, 2) + pow(refsatrS.POS_Z - baseline[a].rZ, 2));//2代参考卫星到流动站的距离
				double length1 = sqrt(pow(CHr.sat[i].POS_X - baseline[a].rX, 2) + pow(CHr.sat[i].POS_Y - baseline[a].rY, 2) + pow(CHr.sat[i].POS_Z - baseline[a].rZ, 2));//3代某颗卫星到流动站的距离
				double length2 = sqrt(pow(refsatmS.POS_X - baseline[a].mX, 2) + pow(refsatmS.POS_Y - baseline[a].mY, 2) + pow(refsatmS.POS_Z - baseline[a].mZ, 2));//2代参考卫星到基站的距离
				double length3 = sqrt(pow(CHm.sat[i].POS_X - baseline[a].mX, 2) + pow(CHm.sat[i].POS_Y - baseline[a].mY, 2) + pow(CHm.sat[i].POS_Z - baseline[a].mZ, 2));//3代某颗卫星到基站的距离
				B[tw][0] = (CHr.sat[i].POS_X - baseline[a].rX) / length1 - (refsatrS.POS_X - baseline[a].rX) / length0;
				B[tw][1] = (CHr.sat[i].POS_Y - baseline[a].rY) / length1 - (refsatrS.POS_Y - baseline[a].rY) / length0;
				B[tw][2] = (CHr.sat[i].POS_Z - baseline[a].rZ) / length1 - (refsatrS.POS_Z - baseline[a].rZ) / length0;
				B[tw][3] = lambdaB;
				double TGD = (yitaB*(CHr.sat[i].TGD - CHm.sat[i].TGD) - yitaC * (refsatrS.TGD - refsatmS.TGD));
				double TROP = ((CHr.sat[i].Trop_Delay - CHm.sat[i].Trop_Delay) - (refsatrS.Trop_Delay - refsatmS.Trop_Delay));
				L[tw][0] = fai + lambdaB * NWC[tw-bdsnum+1][0] - ((length1 - length3) - (length0 - length2))+TGD-TROP;
				P[tw][tw] = 0.09 + 0.09 / pow(sin(refsatrS.E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(CHm.sat[i].E*PI / 180), 2);
				
				tw++;

			}
			double fai = lambdaB * ((refsatrC.data[7] - refsatmC.data[7]) - (refsatrC.data[2] - refsatmC.data[2])) - lambdaC * ((refsatrS.data[2] - refsatmS.data[2]) - (refsatrS.data[3] - refsatmS.data[3]));

			double length0 = sqrt(pow(refsatrS.POS_X - baseline[a].rX, 2) + pow(refsatrS.POS_Y - baseline[a].rY, 2) + pow(refsatrS.POS_Z - baseline[a].rZ, 2));//2代参考卫星到流动站的距离
			double length1 = sqrt(pow(refsatrC.POS_X - baseline[a].rX, 2) + pow(refsatrC.POS_Y - baseline[a].rY, 2) + pow(refsatrC.POS_Z - baseline[a].rZ, 2));//3代参考卫星到流动站的距离
			double length2 = sqrt(pow(refsatmS.POS_X - baseline[a].mX, 2) + pow(refsatmS.POS_Y - baseline[a].mY, 2) + pow(refsatmS.POS_Z - baseline[a].mZ, 2));//2代参考卫星到基站的距离
			double length3 = sqrt(pow(refsatmC.POS_X - baseline[a].mX, 2) + pow(refsatmC.POS_Y - baseline[a].mY, 2) + pow(refsatmC.POS_Z - baseline[a].mZ, 2));//3代参考卫星到基站的距离
			B[tw][0] = (refsatrC.POS_X - baseline[a].rX) / length1 - (refsatrS.POS_X - baseline[a].rX) / length0;
			B[tw][1] = (refsatrC.POS_Y - baseline[a].rY) / length1 - (refsatrS.POS_Y - baseline[a].rY) / length0;
			B[tw][2] = (refsatrC.POS_Z - baseline[a].rZ) / length1 - (refsatrS.POS_Z - baseline[a].rZ) / length0;
			B[tw][3] = lambdaB;
			double TGD = (yitaB*(refsatrC.TGD - refsatmC.TGD) - yitaC * (refsatrS.TGD - refsatmS.TGD));
			double TROP= ((refsatrC.Trop_Delay - refsatmC.Trop_Delay) - (refsatrS.Trop_Delay - refsatmS.Trop_Delay));
			L[tw][0] = fai  - ((length1 - length3) - (length0 - length2)) + TGD-TROP;
			P[tw][tw] = 0.09 + 0.09 / pow(sin(refsatrS.E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(refsatrC.E*PI / 180), 2);
			if (b == 0)
			{
				double sum=0;
				for (int i = 0; i < chnum ; i++)
				{
					sum += L[bdsnum-1+ i][0];
				}
				sum = sum / chnum ;
				X[3][0] = sum / lambdaB;
			}
			else
			{
				if (tempepoch.sat[refposkbds].PRN != refsatmS.PRN || tempepoch.sat[refposkchs].PRN != refsatmC.PRN)//若bds-2或bds-3参考星改变
				{
					double sum = 0;
					for (int i = 0; i < chnum; i++)
					{
						sum += L[bdsnum-1 + i][0];
					}
					sum = sum / chnum ;
					X[3][0] = sum / lambdaB;
				}
			}
			R = P.InvertGaussJordan();
			CMatrix Qw(4,4);
			Qw[0][0] = 100.0; Qw[1][1] = 100.0; Qw[2][2] = 100.0; Qw[3][3] = 0.05*0.05;
			coo.Kalman(I,F,B,L,R,Q,Qw,X);
			CMatrix NEU =  TT*X;
			
			fprintf(p1, "XYZ: %14.4f, %14.4f, %14.4f,NEU: %14.4f, %14.4f, %14.4f, 系统间偏差：%14.4f, \n", X[0][0], X[1][0], X[2][0],NEU[0][0], NEU[1][0], NEU[2][0],X[3][0] );
			
			RMSN += NEU[0][0] * NEU[0][0];
			RMSE += NEU[1][0] * NEU[1][0];
			RMSU += NEU[2][0] * NEU[2][0];
			tempepoch = baseline[a].m_epoch[b];//存储上个历元
			Cr = CHr;
			Cm = CHm;
			Br = BDSr;
			Bm = BDSm;

			TRACE("第%d个历元生成结束\n", b + 1);
		}
		RMSN = sqrt(RMSN / epochnum);
		RMSE = sqrt(RMSE / epochnum);
		RMSU = sqrt(RMSU / epochnum);
		fprintf(p1, "RMSN: %13.4f,RMSE: %13.4f, RMSU:%13.4f, \n", RMSN, RMSE, RMSU);
		fclose(p);
		fclose(p1);
		fclose(p2);
		fclose(p3);
		fclose(p4);
		fclose(C1);     fclose(C2);    fclose(C3);     fclose(C4); fclose(C5);    fclose(C6);   fclose(C7);    fclose(C8);   fclose(C9);   fclose(C10);    fclose(C11);    fclose(C12);    fclose(C13); fclose(C14);
		fclose(B19);  fclose(B20);  fclose(B21);   fclose(B22); fclose(B23); fclose(B24);  fclose(B25);  fclose(B26);  fclose(B27); fclose(B28);  fclose(B29);  fclose(B30); fclose(B32); fclose(B33); fclose(B34); fclose(B35); fclose(B36); fclose(B37);
	}
}



void GNSS::BDS2_BDS3(vector<observe_ppp>baseline)//短基线松组合
{
	size_t filenum = baseline.size();
	for (int a = 0; a < filenum; a++)
	{
		double BLH[3] = { 0 };
		coo.OnXYZtoBLH(baseline[a].rX, baseline[a].rY, baseline[a].rZ, BLH);
		double B = BLH[0];
		double L = BLH[1];
		double H = BLH[2];
		CMatrix TT(3, 3);
		TT[0][0] = -sin(B)*cos(L);    TT[1][0] = -sin(L);   TT[2][0] = cos(B)*cos(L);
		TT[0][1] = -sin(B)*sin(L);    TT[1][1] = cos(L);    TT[2][1] = cos(B)*sin(L);
		TT[0][2] = cos(B);            TT[1][2] = 0;         TT[2][2] = sin(B);

		pppsat refsatmC;
		pppsat refsatrC;
		pppsat refsatmS;
		pppsat refsatrS;
		ppp_epoch tempepoch;//存储上个历元
		ppp_epoch Cm;
		ppp_epoch Cr;
		ppp_epoch Bm;
		ppp_epoch Br;

		FILE *p, *p1, *p2, *p3, *p4;
		CString ss1 = baseline[a].marker_name;
		ss1 = ss1.Left(4);
		CString s1;
		
		s1.Format("短基线松组合%s.txt", ss1);
		CString s = "短基线松组合模糊度NWC.txt";
		
		CString s2 = "短基线松组合模糊度NW.txt";
		CString s3 = "xyz.txt";
		CString s4 = "xyz1.txt";
		p = fopen((LPCSTR)s, "w");
		p1 = fopen((LPCSTR)s1, "w");
		p2 = fopen((LPCSTR)s2, "w");
		p3 = fopen((LPCSTR)s3, "w");
		p4 = fopen((LPCSTR)s4, "w");
		CMatrix X, XWC;
		CMatrix Q, QWC;
		CMatrix  XWS;
		CMatrix  QWS;
		size_t epochnum = baseline[a].m_epoch.size();
		int refposkchs;//存储上历元参考星的序号
		int refposkbds;
		double posr[3]; posr[0] = baseline[a].rX; posr[1] = baseline[a].rY; posr[2] = baseline[a].rZ;
		double posm[3]; posm[0] = baseline[a].mX; posm[1] = baseline[a].mY; posm[2] = baseline[a].mZ;
		coo.check(baseline[a].m_epoch, baseline[a].r_epoch);
		double RMSN = 0;
		double RMSE = 0;
		double RMSU = 0;
		for (int b = 0; b < epochnum; b++)//历元
		{
			//cout << "第"<<b+1<<"个历元" << endl;
			ppp_epoch BDSm;
			ppp_epoch BDSr;
			ppp_epoch CHm;
			ppp_epoch CHr;

			int bdsnum = 0;
			int chnum = 0;//bds-3

			for (int i = 0; i < baseline[a].m_epoch[b].sat.size(); i++)
			{
				if (baseline[a].m_epoch[b].sat[i].judge_use == 0)
				{

					if (baseline[a].m_epoch[b].sat[i].sattype == "C"&&baseline[a].m_epoch[b].sat[i].E > 15 && baseline[a].m_epoch[b].sat[i].PRN != 5)
					{
						BDSm.sat.push_back(baseline[a].m_epoch[b].sat[i]);
						BDSr.sat.push_back(baseline[a].r_epoch[b].sat[i]);
						bdsnum++;
					}
					if (baseline[a].m_epoch[b].sat[i].sattype == "B"&&baseline[a].m_epoch[b].sat[i].E > 15)
					{
						CHm.sat.push_back(baseline[a].m_epoch[b].sat[i]);
						CHr.sat.push_back(baseline[a].r_epoch[b].sat[i]);
						chnum++;
					}

				}
			}
			if (chnum < 3)
				continue;
			//BDS-3

			CMatrix PC(chnum - 1, chnum - 1);
			CMatrix IC(chnum + 2, chnum + 2);
			CMatrix FC(chnum + 2, chnum + 2);
			IC.ones(chnum + 2);
			FC.ones(chnum + 2);
			CMatrix RC(chnum - 1, chnum - 1);
			CMatrix QwC(chnum + 2, chnum + 2);
			CMatrix QWwC(chnum + 2, chnum + 2);

			//BDS-23
			CMatrix B(bdsnum + chnum - 2, 3);
			CMatrix L(bdsnum + chnum - 2, 1);
			CMatrix P(bdsnum + chnum - 2, bdsnum + chnum - 2);
			CMatrix I(3, 3);
			CMatrix F(3, 3);
			I.ones(3);
			F.ones(3);
			CMatrix R(bdsnum + chnum - 2, bdsnum + chnum - 2);
			//CMatrix B(bdsnum -1, 3);
			//CMatrix L(bdsnum - 1, 1);
			//CMatrix P(bdsnum -1, bdsnum -1);




			//计算BDS-3的系统内双差模糊度NWC
			pppsat tempsatmC = CHm.sat[0];//r,m的sattype和prn号和卫星位置都相等  所以用sattype或prn或APP_XYZ判断时用r或m都可
			pppsat tempsatrC = CHr.sat[0];
			double fch1 = 1575.42*1E+6;
			double fch2 = 1561.098*1E+6;
			for (int i = 0; i < chnum; i++)
			{
				if (CHm.sat[i].E > tempsatmC.E)
				{
					tempsatmC = CHm.sat[i];
					tempsatrC = CHr.sat[i];
				}
			}


			if (b == 0)//第一个历元
			{
				refsatmC = tempsatmC;
				refsatrC = tempsatrC;
				//卡尔曼滤波参数初始化
				X.SetSize(3, 1);
				Q.SetSize(3, 3);
				for (int ii = 0; ii < 3; ii++)
				{
					Q[ii][ii] = 10.0*10.0;
				}

				QWC.SetSize(chnum + 2, chnum + 2);
				XWC.SetSize(chnum + 2, 1);

				for (int j = 0; j < 3; j++)
				{
					QWC[j][j] = 30.0 * 30.0;
				}

			}
			else//其他历元
			{
				CMatrix TT_matrixC(chnum - 1, Cr.sat.size() - 1);
				CMatrix matrixC(chnum + 2, Cr.sat.size() + 2);
				QwC.SetSize(chnum + 2, chnum + 2);
				bool find = false;
				for (int i = 0; i < tempepoch.sat.size(); i++)
				{
					if (tempepoch.sat[i].sattype == tempsatmC.sattype&&tempepoch.sat[i].PRN == tempsatmC.PRN)//若参考星在上一历元中未出现，继续使用上一历元的参考星
						find = true;
				}
				if (find == true)
				{
					refsatmC = tempsatmC;
					refsatrC = tempsatrC;
				}
				coo.GetGRCQValue(CHm, Cm, chnum, Cm.sat.size(), refsatmC.PRN, Cm.sat[refposkchs].PRN, TT_matrixC);


				matrixC[0][0] = 1.0;  matrixC[1][1] = 1.0;  matrixC[2][2] = 1.0;
				for (int ii = 0; ii < chnum - 1; ii++)
				{
					for (int jj = 0; jj < Cr.sat.size() - 1; jj++)
					{
						matrixC[ii + 3][jj + 3] = TT_matrixC[ii][jj];
					}
				}

				QWC = matrixC * QWC*matrixC.T();
				XWC = matrixC * XWC;

			}//参考星选择完毕

			//存储上历元的参考星的序号
			for (int i = 0; i < chnum; i++)
			{
				if (CHm.sat[i].PRN == refsatmC.PRN)
				{
					refposkchs = i;
					break;
				}
			}
			CMatrix LWC(chnum - 1, 1);
			CMatrix BWC(chnum - 1, chnum + 2);

			coo.MW(b, posm, posr, CHm, CHr, refsatmC, refsatrC, Cr, XWC, LWC, BWC, QWC, PC);
			RC = PC.InvertGaussJordan();
			coo.Kalman(IC, FC, BWC, LWC, RC, QWC, QWwC, XWC);

			CLambda lamNC;
			lamNC.numLiyuan = b + 1;//当前历元个数
			lamNC.N = chnum - 1;//模糊度个数
			CMatrix X_LAM(chnum - 1, 1);
			CMatrix Q_LAM(chnum - 1, chnum - 1);
			for (int i = 0; i < chnum - 1; i++)
			{
				X_LAM[i][0] = XWC[i + 3][0];
				Q_LAM[i][i] = QWC[i + 3][i + 3];
			}
			lamNC.a = X_LAM;
			lamNC.Qa = Q_LAM;
			lamNC.OnFixAmbigultyMain();
			double ratioNWC = lamNC.Ratio;
			CMatrix NWC(chnum - 1, 1);
			NWC = lamNC.afix;
			CMatrix NNWC(12, 1);
			for (int i = 0; i < chnum - 1; i++)
				NNWC[i][0] = NWC[i][0];

			fprintf(p, "固定解NW：chnum-1：%4d ratio:%7.3f  %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f  \n", chnum - 1, ratioNWC, NNWC[0][0], NNWC[1][0], NNWC[2][0], NNWC[3][0], NNWC[4][0], NNWC[5][0], NNWC[6][0], NNWC[7][0], NNWC[8][0], NNWC[9][0], NNWC[10][0], NNWC[11][0]);

			CMatrix PP(bdsnum - 1, 1);//伪距
			CMatrix fai0_11(bdsnum - 1, 1);
			CMatrix fai1_10(bdsnum - 1, 1);
			CMatrix N0_11(bdsnum - 1, 1);
			CMatrix N1_10(bdsnum - 1, 1);


			pppsat tempsatmS = BDSm.sat[0];//r,m的sattype和prn号和卫星位置都相等  所以用sattype或prn或APP_XYZ判断时用r或m都可
			pppsat tempsatrS = BDSr.sat[0];

			double fbds1 = 1561.098*1E+6;
			double fbds2 = 1207.14*1E+6;
			double fbds3 = 1268.52*1E+6;

			double lambda011 = C_light / (fbds2 + fbds3);
			double lambda0_11 = C_light / (-fbds2 + fbds3);
			double lambda1_10 = C_light / (fbds1 - fbds2);
			double lambdaC = C_light / (fbds1 - fbds2);
			double lambdaB = C_light / (fch1 - fch2);
			double yita011 = fbds1 * fbds1*(1.0 / fbds2 + 1.0 / fbds3) / (fbds2 + fbds3);
			double yita0_11 = fbds1 * fbds1*(-1.0 / fbds2 + 1.0 / fbds3) / (-fbds2 + fbds3);
			double yita1_10 = fbds1 * fbds1*(1.0 / fbds1 + (-1.0) / fbds2) / (fbds1 - fbds2);
			double yitaC = yita1_10;//BDS2
			double yitaB = fch1 * fch1*(1.0 / fch1 + (-1.0) / fch2) / (fch1 - fch2);
			for (int i = 0; i < bdsnum; i++)
			{
				if (BDSm.sat[i].E > tempsatmS.E)
				{
					tempsatmS = BDSm.sat[i];
					tempsatrS = BDSr.sat[i];
				}
			}


			if (b == 0)//第一个历元
			{
				refsatmS = tempsatmS;
				refsatrS = tempsatrS;
			}
			else//其他历元
			{
				bool find = false;
				for (int i = 0; i < tempepoch.sat.size(); i++)
				{
					if (tempepoch.sat[i].sattype == tempsatmS.sattype&&tempepoch.sat[i].PRN == tempsatmS.PRN)//若参考星在上一历元中未出现，继续使用上一历元的参考星
						find = true;
				}
				if (find == true)
				{
					refsatmS = tempsatmS;
					refsatrS = tempsatrS;
				}

			}//参考星选择完毕
			//存储上历元的参考星的序号
			for (int i = 0; i < bdsnum; i++)
			{
				if (BDSm.sat[i].PRN == refsatmS.PRN)
				{
					refposkbds = i;
					break;
				}
			}
			int tw = 0;
			for (int i = 0; i < bdsnum; i++)
			{
				if (refsatmS.PRN == BDSm.sat[i].PRN)
					continue;
				PP[tw][0] = (fbds2* ((BDSr.sat[i].data[1] - BDSm.sat[i].data[1]) - (refsatrS.data[1] - refsatmS.data[1])) + fbds3 * ((BDSr.sat[i].data[4] - BDSm.sat[i].data[4]) - (refsatrS.data[4] - refsatmS.data[4]))) / (fbds2 + fbds3);
				fai0_11[tw][0] = (-C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3])) + C_light * ((BDSr.sat[i].data[5] - BDSm.sat[i].data[5]) - (refsatrS.data[5] - refsatmS.data[5]))) / (-fbds2 + fbds3);
				fai1_10[tw][0] = (C_light * ((BDSr.sat[i].data[2] - BDSm.sat[i].data[2]) - (refsatrS.data[2] - refsatmS.data[2])) - C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3]))) / (fbds1 - fbds2);
				N0_11[tw][0] = round((PP[tw][0] - fai0_11[tw][0]) / lambda0_11);
				double TGD = (BDSr.sat[i].TGD - BDSm.sat[i].TGD) - (refsatrS.TGD - refsatmS.TGD);
				double TROP = (BDSr.sat[i].Trop_Delay - BDSm.sat[i].Trop_Delay) - (refsatrS.Trop_Delay - refsatmS.Trop_Delay);
				N1_10[tw][0] = round(-1.0 / lambda1_10 * (fai1_10[tw][0] - fai0_11[tw][0] - (yita0_11 - yita1_10)*TGD - lambda0_11 * N0_11[tw][0]));

				double length0 = sqrt(pow(refsatrS.POS_X - baseline[a].rX, 2) + pow(refsatrS.POS_Y - baseline[a].rY, 2) + pow(refsatrS.POS_Z - baseline[a].rZ, 2));//参考卫星到流动站的距离
				double length1 = sqrt(pow(BDSr.sat[i].POS_X - baseline[a].rX, 2) + pow(BDSr.sat[i].POS_Y - baseline[a].rY, 2) + pow(BDSr.sat[i].POS_Z - baseline[a].rZ, 2));//某颗卫星到流动站的距离
				double length2 = sqrt(pow(refsatmS.POS_X - baseline[a].mX, 2) + pow(refsatmS.POS_Y - baseline[a].mY, 2) + pow(refsatmS.POS_Z - baseline[a].mZ, 2));//参考卫星到基站的距离
				double length3 = sqrt(pow(BDSm.sat[i].POS_X - baseline[a].mX, 2) + pow(BDSm.sat[i].POS_Y - baseline[a].mY, 2) + pow(BDSm.sat[i].POS_Z - baseline[a].mZ, 2));//某颗卫星到基站的距离
				B[tw][0] = (BDSr.sat[i].POS_X - baseline[a].rX) / length1 - (refsatrS.POS_X - baseline[a].rX) / length0;
				B[tw][1] = (BDSr.sat[i].POS_Y - baseline[a].rY) / length1 - (refsatrS.POS_Y - baseline[a].rY) / length0;
				B[tw][2] = (BDSr.sat[i].POS_Z - baseline[a].rZ) / length1 - (refsatrS.POS_Z - baseline[a].rZ) / length0;
				L[tw][0] = fai1_10[tw][0] + lambda1_10 * N1_10[tw][0] - ((length1 - length3) - (length0 - length2)) + yita1_10 * TGD - TROP;
				P[tw][tw] = 0.09 + 0.09 / pow(sin(refsatrS.E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(BDSm.sat[i].E*PI / 180), 2);

				tw++;
			}

			for (int i = 0; i < chnum; i++)
			{
				if (refsatmC.PRN == CHm.sat[i].PRN)
					continue;

				double fai = lambdaB * ((CHr.sat[i].data[7] - CHm.sat[i].data[7]) - (CHr.sat[i].data[2] - CHm.sat[i].data[2])) - lambdaB * ((refsatrC.data[7] - refsatmC.data[7]) - (refsatrC.data[2] - refsatmC.data[2]));

				double length0 = sqrt(pow(refsatrC.POS_X - baseline[a].rX, 2) + pow(refsatrC.POS_Y - baseline[a].rY, 2) + pow(refsatrC.POS_Z - baseline[a].rZ, 2));//3代参考卫星到流动站的距离
				double length1 = sqrt(pow(CHr.sat[i].POS_X - baseline[a].rX, 2) + pow(CHr.sat[i].POS_Y - baseline[a].rY, 2) + pow(CHr.sat[i].POS_Z - baseline[a].rZ, 2));//3代某颗卫星到流动站的距离
				double length2 = sqrt(pow(refsatmC.POS_X - baseline[a].mX, 2) + pow(refsatmC.POS_Y - baseline[a].mY, 2) + pow(refsatmC.POS_Z - baseline[a].mZ, 2));//3代参考卫星到基站的距离
				double length3 = sqrt(pow(CHm.sat[i].POS_X - baseline[a].mX, 2) + pow(CHm.sat[i].POS_Y - baseline[a].mY, 2) + pow(CHm.sat[i].POS_Z - baseline[a].mZ, 2));//3代某颗卫星到基站的距离
				B[tw][0] = (CHr.sat[i].POS_X - baseline[a].rX) / length1 - (refsatrC.POS_X - baseline[a].rX) / length0;
				B[tw][1] = (CHr.sat[i].POS_Y - baseline[a].rY) / length1 - (refsatrC.POS_Y - baseline[a].rY) / length0;
				B[tw][2] = (CHr.sat[i].POS_Z - baseline[a].rZ) / length1 - (refsatrC.POS_Z - baseline[a].rZ) / length0;
				
				L[tw][0] = fai + lambdaB * NWC[tw - bdsnum + 1][0] - ((length1 - length3) - (length0 - length2)) + (yitaB*(CHr.sat[i].TGD - CHm.sat[i].TGD) - yitaB * (refsatrC.TGD - refsatmC.TGD)) - ((CHr.sat[i].Trop_Delay - CHm.sat[i].Trop_Delay) - (refsatrC.Trop_Delay - refsatmC.Trop_Delay));
				P[tw][tw] = 0.09 + 0.09 / pow(sin(refsatrC.E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(CHm.sat[i].E*PI / 180), 2);

				tw++;

			}
			
			R = P.InvertGaussJordan();
			CMatrix Qw(3, 3);
			coo.Kalman(I, F, B, L, R, Q, Qw, X);
			CMatrix NEU = TT * X;

			fprintf(p1, "XYZ: %14.4f, %14.4f, %14.4f,NEU: %14.4f, %14.4f, %14.4f,  \n", X[0][0], X[1][0], X[2][0], NEU[0][0], NEU[1][0], NEU[2][0]);

			RMSN += NEU[0][0] * NEU[0][0];
			RMSE += NEU[1][0] * NEU[1][0];
			RMSU += NEU[2][0] * NEU[2][0];
			tempepoch = baseline[a].m_epoch[b];//存储上个历元
			Cr = CHr;
			Cm = CHm;
			Br = BDSr;
			Bm = BDSm;

			TRACE("第%d个历元生成结束\n", b + 1);
		}
		RMSN = sqrt(RMSN / epochnum);
		RMSE = sqrt(RMSE / epochnum);
		RMSU = sqrt(RMSU / epochnum);
		fprintf(p1, "RMSN: %13.4f,RMSE: %13.4f, RMSU:%13.4f, \n", RMSN, RMSE, RMSU);
		fclose(p);
		fclose(p1);
		fclose(p2);
		fclose(p3);
		fclose(p4);
			}
}

void GNSS::BDS2_BDS3_long(vector<observe_ppp>baseline)//长基线松组合
{
	size_t filenum = baseline.size();
	for (int a = 0; a < filenum; a++)
	{
		double BLH[3] = { 0 };
		coo.OnXYZtoBLH(baseline[a].rX, baseline[a].rY, baseline[a].rZ, BLH);
		double B = BLH[0];
		double L = BLH[1];
		double H = BLH[2];
		CMatrix TT(3, 3);
		TT[0][0] = -sin(B)*cos(L);    TT[1][0] = -sin(L);   TT[2][0] = cos(B)*cos(L);
		TT[0][1] = -sin(B)*sin(L);    TT[1][1] = cos(L);    TT[2][1] = cos(B)*sin(L);
		TT[0][2] = cos(B);            TT[1][2] = 0;         TT[2][2] = sin(B);

		pppsat refsatmC;
		pppsat refsatrC;
		pppsat refsatmS;
		pppsat refsatrS;
		ppp_epoch tempepoch;//存储上个历元
		ppp_epoch Cm;
		ppp_epoch Cr;
		ppp_epoch Bm;
		ppp_epoch Br;

		FILE *p1, *p2, *p3, *p4;
		CString ss1 = baseline[a].marker_name;
		ss1 = ss1.Left(4);
		CString s1;

		s1.Format("长基线松组合%s.txt", ss1);
		//CString s = "长基线松组合模糊度NWC01_1.txt";

		CString s2 = "长基线松组合模糊度NW.txt";
		CString s3 = "长基线松组合模糊度NWC01_1.txt";
		CString s4 = "长基线松组合模糊度NWC01_1.txt";
		
		p1 = fopen((LPCSTR)s1, "w");
		p2 = fopen((LPCSTR)s2, "w");
		p3 = fopen((LPCSTR)s3, "w");
		p4 = fopen((LPCSTR)s4, "w");
		CMatrix X;
		CMatrix Q;
		CMatrix filter1;//NWC的多历元平滑值
		CMatrix filter2;//N1-10的多历元平滑值
		
		size_t epochnum = baseline[a].m_epoch.size();
		int refposkchs;//存储上历元参考星的序号
		int refposkbds;
		double posr[3]; posr[0] = baseline[a].rX; posr[1] = baseline[a].rY; posr[2] = baseline[a].rZ;
		double posm[3]; posm[0] = baseline[a].mX; posm[1] = baseline[a].mY; posm[2] = baseline[a].mZ;
		coo.check(baseline[a].m_epoch, baseline[a].r_epoch);
		double RMSN = 0;
		double RMSE = 0;
		double RMSU = 0;
		int NUM = 0;
		double satnum2[14];//0-13对应C01-C14
		for (int i = 0; i < 14; i++)
			satnum2[i] = 1.0;
		double satnum3[19];//0-18对应B19-B37
		for (int i = 0; i < 19; i++)
			satnum3[i] = 1.0;
		for (int b = 0; b < epochnum; b++)//历元
		{
			//cout << "第"<<b+1<<"个历元" << endl;
			ppp_epoch BDSm;
			ppp_epoch BDSr;
			ppp_epoch CHm;
			ppp_epoch CHr;
			int junknum2 = 0;
			int junknum3 = 0;
			int bdsnum = 0;
			int chnum = 0;//bds-3
			int flag = 0;
			for (int i = 0; i < baseline[a].m_epoch[b].sat.size(); i++)
			{
				if (baseline[a].m_epoch[b].sat[i].judge_use == 0)
				{

					if (baseline[a].m_epoch[b].sat[i].sattype == "C"&&baseline[a].m_epoch[b].sat[i].E > 15/* && baseline[a].m_epoch[b].sat[i].PRN != 5*/)
					{
						BDSm.sat.push_back(baseline[a].m_epoch[b].sat[i]);
						BDSr.sat.push_back(baseline[a].r_epoch[b].sat[i]);
						bdsnum++;
					}
					if (baseline[a].m_epoch[b].sat[i].sattype == "B"&&baseline[a].m_epoch[b].sat[i].E > 15 )
					{
						CHm.sat.push_back(baseline[a].m_epoch[b].sat[i]);
						CHr.sat.push_back(baseline[a].r_epoch[b].sat[i]);
						chnum++;
					}

				}
			}

			if (b > 0)
			{
				//防止出现上历元参考星本历元未出现的情况,如果出现这种情况，转换矩阵会报错
				for (int i = 0; i < chnum; i++)
				{

					if (Cm.sat[refposkchs].PRN == CHm.sat[i].PRN)
					{
						flag++; break;
					}

				}
				for (int i = 0; i < bdsnum; i++)
				{

					if (Bm.sat[refposkbds].PRN == BDSm.sat[i].PRN)
					{
						flag++; break;
					}

				}
				
			}
			//BDS-3

			//BDS-3
			CMatrix PC01_1(chnum - 1, 1);//伪距
			CMatrix faiC01_1(chnum - 1, 1);
			CMatrix faiC10_1(chnum - 1, 1);
			CMatrix NWC01_1(chnum - 1, 1);
			CMatrix NWC10_1(chnum - 1, 1);





			//计算BDS-3的系统内双差模糊度NWC
			pppsat tempsatmC = CHm.sat[0];//r,m的sattype和prn号和卫星位置都相等  所以用sattype或prn或APP_XYZ判断时用r或m都可
			pppsat tempsatrC = CHr.sat[0];
			double fch1 = 1561.098*1E+6;
			double fch2 = 1268.52*1E+6;
			double fch3 = 1176.45*1E+6;
			double lambdaB01_1 = C_light / (fch2 - fch3);
			double lambdaB10_1 = C_light / (fch1 - fch3);
			
			for (int i = 0; i < chnum; i++)
			{
				if (CHr.sat[i].E > tempsatrC.E)
				{
					tempsatmC = CHm.sat[i];
					tempsatrC = CHr.sat[i];
				}
			}


			if (b == 0)//第一个历元
			{
				refsatmC = tempsatmC;
				refsatrC = tempsatrC;
				X.SetSize(3, 1);
				Q.SetSize(3, 3);
				for (int ii = 0; ii < 3; ii++)
				{
					Q[ii][ii] = 10.0*10.0;
				}
				filter1.SetSize(chnum-1,1);
				filter2.SetSize(bdsnum - 1, 1);
			}
			else//其他历元
			{
				bool find = false;
				for (int i = 0; i < tempepoch.sat.size(); i++)
				{
					if (tempepoch.sat[i].sattype == tempsatmC.sattype&&tempepoch.sat[i].PRN == tempsatmC.PRN)//若参考星在上一历元中未出现，继续使用上一历元的参考星
						find = true;
				}
				if (find == true)
				{
					refsatmC = tempsatmC;
					refsatrC = tempsatrC;
				}
				
					CMatrix TT_matrix33(Cr.sat.size() - 1, Cr.sat.size() - 1);
					coo.GetGRCQValue(Cm, Cm, Cm.sat.size(), Cm.sat.size(), refsatmC.PRN, Cm.sat[refposkchs].PRN, TT_matrix33);//都带入上历元序列
					filter1 = TT_matrix33 * filter1;//上历元序列的参考星转换成本历元序列的参考星
					CMatrix TT_matrix3(chnum - 1, Cr.sat.size() - 1);
					coo.GetGRCQValue(CHm, Cm, chnum, Cm.sat.size(), refsatmC.PRN, refsatmC.PRN, TT_matrix3);//都带入本历元参考星
					filter1 = TT_matrix3 * filter1;
			
			}//参考星选择完毕

			
			int w = 0;
			for (int i = 0; i < chnum; i++)
			{

				if (refsatmC.PRN == CHm.sat[i].PRN)
					continue;
				PC01_1[w][0] = 0.0402* ((CHr.sat[i].data[0] - CHm.sat[i].data[0]) - (refsatrC.data[0] - refsatmC.data[0])) + 0.3951* ((CHr.sat[i].data[4] - CHm.sat[i].data[4]) - (refsatrC.data[4] - refsatmC.data[4])) + 0.5647 * ((CHr.sat[i].data[1] - CHm.sat[i].data[1]) - (refsatrC.data[1] - refsatmC.data[1]));
				faiC01_1[w][0] = (C_light * ((CHr.sat[i].data[6] - CHm.sat[i].data[6]) - (refsatrC.data[6] - refsatmC.data[6])) - C_light * ((CHr.sat[i].data[3] - CHm.sat[i].data[3]) - (refsatrC.data[3] - refsatmC.data[3]))) / (fch2 - fch3);
				faiC10_1[w][0] = (C_light * ((CHr.sat[i].data[2] - CHm.sat[i].data[2]) - (refsatrC.data[2] - refsatmC.data[2])) - C_light * ((CHr.sat[i].data[3] - CHm.sat[i].data[3]) - (refsatrC.data[3] - refsatmC.data[3]))) / (fch1 - fch3);
				NWC01_1[w][0] = round((PC01_1[w][0] - faiC01_1[w][0]) / (C_light / (fch2 - fch3)));

				if (b == 0)
				{
					filter1[w][0] = -1.0 / lambdaB10_1 * (faiC10_1[w][0] - faiC01_1[w][0] - lambdaB01_1 * NWC01_1[w][0]);
				}
				else
				{
					double k = satnum3[CHm.sat[i].PRN - 19];
					double temp = -1.0 / lambdaB10_1 * (faiC10_1[w][0] - faiC01_1[w][0] - lambdaB01_1 * NWC01_1[w][0]);
					if (abs(round(temp) - temp) > 0.3 || ((abs(filter1[w][0] - temp) > 0.4) && (abs(filter1[w][0] - temp) < 1.5)))
					{
						junknum3++;
						CHm.sat[i].judge_use = 1;
						w++; continue;
					}
					if (abs(filter1[w][0] - temp) < 0.4)
					{
						filter1[w][0] = filter1[w][0] * k / (k + 1) + temp / (k + 1);
						satnum3[CHm.sat[i].PRN - 19]++;
					}
					else
					{
						filter1[w][0] = temp;
						satnum3[CHm.sat[i].PRN - 19] = 1;
					}
				}
				NWC10_1[w][0] = round(filter1[w][0]);
				//doubleNWC[w][0] = (PPC[w][0] - faiC[w][0]) / (C_light / (fch2 - fch3)) - NWC[w][0];
				if (abs(filter1[w][0] - NWC10_1[w][0]) > 0.15 || (b > 100 && satnum3[CHm.sat[i].PRN - 19] < 5)) 
				{
					junknum3++;
					CHm.sat[i].judge_use = 1;
				}
				w++;
			}
			for (int i = 19; i < 38; i++)//B19-B37中在本历元未出现的置1
			{
				if (coo.trip(CHm, i) == false)
					satnum3[i - 19] = 1.0;
			}
			CMatrix NNWC01_1(12, 1);
			for (int i = 0; i < chnum - 1; i++)
				NNWC01_1[i][0] = NWC01_1[i][0];
			fprintf(p3, "%4d 固定解NW：chnum-1：%4d  %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f  \n", b, chnum - 1, NNWC01_1[0][0], NNWC01_1[1][0], NNWC01_1[2][0], NNWC01_1[3][0], NNWC01_1[4][0], NNWC01_1[5][0], NNWC01_1[6][0], NNWC01_1[7][0], NNWC01_1[8][0], NNWC01_1[9][0], NNWC01_1[10][0], NNWC01_1[11][0]);
			CMatrix NNWC10_1(12, 1);

			for (int i = 0; i < chnum - 1; i++)
				NNWC10_1[i][0] = NWC10_1[i][0];
			fprintf(p4, "%4d 固定解NW：chnum-1：%4d  %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f  \n", b, chnum - 1, NNWC10_1[0][0], NNWC10_1[1][0], NNWC10_1[2][0], NNWC10_1[3][0], NNWC10_1[4][0], NNWC10_1[5][0], NNWC10_1[6][0], NNWC10_1[7][0], NNWC10_1[8][0], NNWC10_1[9][0], NNWC10_1[10][0], NNWC10_1[11][0]);
			
			//存储上历元的参考星的序号
			for (int i = 0; i < chnum; i++)
			{
				if (CHm.sat[i].PRN == refsatmC.PRN)
				{
					refposkchs = i;
					break;
				}
			}




			CMatrix PP(bdsnum - 1, 1);//伪距
			CMatrix fai0_11(bdsnum - 1, 1);
			CMatrix fai14_5(bdsnum - 1, 1);
			CMatrix fai1_10(bdsnum - 1, 1);
			CMatrix N0_11(bdsnum - 1, 1);
			CMatrix N14_5(bdsnum - 1, 1);
			CMatrix N1_10(bdsnum - 1, 1);


			pppsat tempsatmS = BDSm.sat[0];//r,m的sattype和prn号和卫星位置都相等  所以用sattype或prn或APP_XYZ判断时用r或m都可
			pppsat tempsatrS = BDSr.sat[0];

			double fbds1 = 1561.098*1E+6;
			double fbds2 = 1207.14*1E+6;
			double fbds3 = 1268.52*1E+6;

			double lambda1_10 = C_light / (fbds1 - fbds2);
			double lambda0_11 = C_light / (-fbds2 + fbds3);
			double lambda14_5 = C_light / (fbds1 + 4 * fbds2 - 5 * fbds3);
			double lambdaC = C_light / (fbds1 - fbds2);
			
			
			for (int i = 0; i < bdsnum; i++)
			{
				if (BDSr.sat[i].E > tempsatrS.E)
				{
					tempsatmS = BDSm.sat[i];
					tempsatrS = BDSr.sat[i];
				}
			}


			if (b == 0)//第一个历元
			{
				refsatmS = tempsatmS;
				refsatrS = tempsatrS;
			}
			else//其他历元
			{
				bool find = false;
				for (int i = 0; i < tempepoch.sat.size(); i++)
				{
					if (tempepoch.sat[i].sattype == tempsatmS.sattype&&tempepoch.sat[i].PRN == tempsatmS.PRN)//若参考星在上一历元中未出现，继续使用上一历元的参考星
						find = true;
				}
				if (find == true)
				{
					refsatmS = tempsatmS;
					refsatrS = tempsatrS;
				}
				
					CMatrix TT_matrix(Br.sat.size() - 1, Br.sat.size() - 1);
					coo.GetGRCQValue(Bm, Bm, Bm.sat.size(), Bm.sat.size(), refsatmS.PRN, Bm.sat[refposkbds].PRN, TT_matrix);//都带入上历元序列
					filter2 = TT_matrix * filter2;//上历元序列的参考星转换成本历元序列的参考星
					CMatrix TT_matrix2(bdsnum - 1, Br.sat.size() - 1);
					coo.GetGRCQValue(BDSm, Bm, bdsnum, Bm.sat.size(), refsatmS.PRN, refsatmS.PRN, TT_matrix2);//都带入本历元参考星
					filter2 = TT_matrix2 * filter2;
				
				
			}//参考星选择完毕
		
			int tw = 0; 
			for (int i = 0; i < bdsnum; i++)
			{
				if (refsatmS.PRN == BDSm.sat[i].PRN)
					continue;
				PP[tw][0] = 0.0199*((BDSr.sat[i].data[0] - BDSm.sat[i].data[0]) - (refsatrS.data[0] - refsatmS.data[0])) + 0.5526* ((BDSr.sat[i].data[1] - BDSm.sat[i].data[1]) - (refsatrS.data[1] - refsatmS.data[1])) + 0.4275 * ((BDSr.sat[i].data[4] - BDSm.sat[i].data[4]) - (refsatrS.data[4] - refsatmS.data[4]));
				double P110 = (fbds1* ((BDSr.sat[i].data[0] - BDSm.sat[i].data[0]) - (refsatrS.data[0] - refsatmS.data[0])) + fbds2 * ((BDSr.sat[i].data[1] - BDSm.sat[i].data[1]) - (refsatrS.data[1] - refsatmS.data[1]))) / (fbds1 + fbds2);
				fai0_11[tw][0] = (-C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3])) + C_light * ((BDSr.sat[i].data[5] - BDSm.sat[i].data[5]) - (refsatrS.data[5] - refsatmS.data[5]))) / (-fbds2 + fbds3);
				fai1_10[tw][0] = (C_light * ((BDSr.sat[i].data[2] - BDSm.sat[i].data[2]) - (refsatrS.data[2] - refsatmS.data[2])) - C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3]))) / (fbds1 - fbds2);
				fai14_5[tw][0] = (C_light * ((BDSr.sat[i].data[2] - BDSm.sat[i].data[2]) - (refsatrS.data[2] - refsatmS.data[2])) + 4 * C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3])) - 5 * C_light * ((BDSr.sat[i].data[5] - BDSm.sat[i].data[5]) - (refsatrS.data[5] - refsatmS.data[5]))) / (fbds1 + 4 * fbds2 - 5 * fbds3);
				N0_11[tw][0] = round((PP[tw][0] - fai0_11[tw][0]) / lambda0_11);
				double n0_11 = (PP[tw][0] - fai0_11[tw][0]) / lambda0_11 - N0_11[tw][0];
				if (b == 0)
				{
					filter2[tw][0] = -1.0 / lambda1_10 * (fai1_10[tw][0] - fai0_11[tw][0] - lambda0_11 * N0_11[tw][0]);
				}
				else
				{
					double k = satnum2[BDSm.sat[i].PRN - 1];

					double temp = -1.0 / lambda1_10 * (fai1_10[tw][0] - fai0_11[tw][0] - lambda0_11 * N0_11[tw][0]);
					if (abs(round(temp) - temp) > 0.3 || ((abs(filter2[tw][0] - temp) > 0.4) && (abs(filter2[tw][0] - temp) < 1.5)))
					{
						junknum2++;
						BDSm.sat[i].judge_use = 1;
						tw++; continue;
					}
					if (abs(filter2[tw][0] - temp) < 0.4)
					{
						
						filter2[tw][0] = filter2[tw][0] * k / (k + 1) + temp / (k + 1);
						satnum2[BDSm.sat[i].PRN - 1]++;

					}
					else
					{

						filter2[tw][0] = -1.0 / lambda1_10 * (fai1_10[tw][0] - fai0_11[tw][0] - lambda0_11 * N0_11[tw][0]);
						satnum2[BDSm.sat[i].PRN - 1] = 1;
					}
				}

				//N14_5[tw][0] = round(filter2[tw][0]);
				//N1_10[tw][0] = 5 * N0_11[tw][0] + N14_5[tw][0];
				N1_10[tw][0] = round(filter2[tw][0]);
				if (abs(filter2[tw][0] - N1_10[tw][0]) > 0.2 || abs(n0_11) > 0.2 || (b > 100 && satnum2[BDSm.sat[i].PRN - 1] < 5))
				{
					junknum2++;
					BDSm.sat[i].judge_use = 1;
				}
				tw++;
			}
			for (int i = 1; i < 15; i++)//C01-C14中在本历元未出现的置1
			{
				if (coo.trip(BDSm, i) == false)
					satnum2[i - 1] = 1.0;
			}





			//BDS-23
			CMatrix B(chnum + bdsnum - 2 - junknum2 - junknum3, 3);
			CMatrix L(chnum + bdsnum - 2 - junknum2 - junknum3, 1);
			CMatrix P(chnum + bdsnum - 2 - junknum2 - junknum3, chnum + bdsnum - 2 - junknum2 - junknum3);
			CMatrix I(3, 3);
			CMatrix F(3, 3);
			I.ones(3);
			F.ones(3);
			CMatrix R(chnum + bdsnum - 2 - junknum2 - junknum3, chnum + bdsnum - 2 - junknum2 - junknum3);
			//chnum + bdsnum - 2 - junknum2 - junknum3
			tw = 0; int grow = 0;
			for (int i = 0; i < bdsnum; i++)
			{

				if (refsatmS.PRN == BDSm.sat[i].PRN)
					continue;
				if (BDSm.sat[i].judge_use == 1)
				{
					grow++;
					continue;
				}
				double TGD = (BDSr.sat[i].TGD - BDSm.sat[i].TGD) - (refsatrS.TGD - refsatmS.TGD);
				double TROP = (BDSr.sat[i].Trop_Delay - BDSm.sat[i].Trop_Delay) - (refsatrS.Trop_Delay - refsatmS.Trop_Delay);
				double relate = (BDSr.sat[i].Relat - BDSm.sat[i].Relat) - (refsatrS.Relat - refsatmS.Relat);
				double satclock = (BDSr.sat[i].Sat_clock - BDSm.sat[i].Sat_clock) - (refsatrS.Sat_clock - refsatmS.Sat_clock);
				double sagnac = (BDSr.sat[i].Sagnac - BDSm.sat[i].Sagnac) - (refsatrS.Sagnac - refsatmS.Sagnac);
				//double antenna_height1= (BDSr.sat[i].Antenna_Height - BDSm.sat[i].Antenna_Height) - (refsatrS.Antenna_Height - refsatmS.Antenna_Height);
				double sat_antenna1 = (BDSr.sat[i].Sat_Antenna - BDSm.sat[i].Sat_Antenna) - (refsatrS.Sat_Antenna - refsatmS.Sat_Antenna);
				double error = - TGD - TROP + relate + satclock - sagnac + sat_antenna1;

				double length0 = sqrt(pow(refsatrS.POS_X - baseline[a].rX, 2) + pow(refsatrS.POS_Y - baseline[a].rY, 2) + pow(refsatrS.POS_Z - baseline[a].rZ, 2));//参考卫星到流动站的距离
				double length1 = sqrt(pow(BDSr.sat[i].POS_X - baseline[a].rX, 2) + pow(BDSr.sat[i].POS_Y - baseline[a].rY, 2) + pow(BDSr.sat[i].POS_Z - baseline[a].rZ, 2));//某颗卫星到流动站的距离
				double length2 = sqrt(pow(refsatmS.POS_X - baseline[a].mX, 2) + pow(refsatmS.POS_Y - baseline[a].mY, 2) + pow(refsatmS.POS_Z - baseline[a].mZ, 2));//参考卫星到基站的距离
				double length3 = sqrt(pow(BDSm.sat[i].POS_X - baseline[a].mX, 2) + pow(BDSm.sat[i].POS_Y - baseline[a].mY, 2) + pow(BDSm.sat[i].POS_Z - baseline[a].mZ, 2));//某颗卫星到基站的距离
				B[tw][0] = (BDSr.sat[i].POS_X - baseline[a].rX) / length1 - (refsatrS.POS_X - baseline[a].rX) / length0;
				B[tw][1] = (BDSr.sat[i].POS_Y - baseline[a].rY) / length1 - (refsatrS.POS_Y - baseline[a].rY) / length0;
				B[tw][2] = (BDSr.sat[i].POS_Z - baseline[a].rZ) / length1 - (refsatrS.POS_Z - baseline[a].rZ) / length0;
				//double N = N1_10[tw + grow][0];
				L[tw][0] = fai1_10[tw + grow][0] + lambda1_10 * N1_10[tw + grow][0] - ((length1 - length3) - (length0 - length2)) + error;
				P[tw][tw] = 0.09 + 0.09 / pow(sin(BDSm.sat[i].E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(BDSr.sat[i].E*PI / 180), 2);
				tw++;
			}
			
			grow = 0;
			for (int i = 0; i < chnum; i++)
			{
				if (refsatmC.PRN == CHm.sat[i].PRN)
					continue;
				if (CHm.sat[i].judge_use == 1)
				{
					grow++;
					continue;
				}
				double TGD = (CHr.sat[i].TGD - CHm.sat[i].TGD) - (refsatrC.TGD - refsatmC.TGD);
				double TROP = (CHr.sat[i].Trop_Delay - CHm.sat[i].Trop_Delay) - (refsatrC.Trop_Delay - refsatmC.Trop_Delay);
				double relate = (CHr.sat[i].Relat - CHm.sat[i].Relat) - (refsatrC.Relat - refsatmC.Relat);
				double satclock = (CHr.sat[i].Sat_clock - CHm.sat[i].Sat_clock) - (refsatrC.Sat_clock - refsatmC.Sat_clock);
				double sagnac = (CHr.sat[i].Sagnac - CHm.sat[i].Sagnac) - (refsatrC.Sagnac - refsatmC.Sagnac);
				//double antenna_height1= (BDSr.sat[i].Antenna_Height - BDSm.sat[i].Antenna_Height) - (refsatrS.Antenna_Height - refsatmS.Antenna_Height);
				double sat_antenna1 = (CHr.sat[i].Sat_Antenna - CHm.sat[i].Sat_Antenna) - (refsatrC.Sat_Antenna - refsatmC.Sat_Antenna);
				double error = -TGD - TROP + relate + satclock - sagnac + sat_antenna1;

				double fai = lambdaB10_1 * ((CHr.sat[i].data[2] - CHm.sat[i].data[2]) - (CHr.sat[i].data[3] - CHm.sat[i].data[3])) - lambdaB10_1 * ((refsatrC.data[2] - refsatmC.data[2]) - (refsatrC.data[3] - refsatmC.data[3]));
				double length0 = sqrt(pow(refsatrC.POS_X - baseline[a].rX, 2) + pow(refsatrC.POS_Y - baseline[a].rY, 2) + pow(refsatrC.POS_Z - baseline[a].rZ, 2));//3代参考卫星到流动站的距离
				double length1 = sqrt(pow(CHr.sat[i].POS_X - baseline[a].rX, 2) + pow(CHr.sat[i].POS_Y - baseline[a].rY, 2) + pow(CHr.sat[i].POS_Z - baseline[a].rZ, 2));//3代某颗卫星到流动站的距离
				double length2 = sqrt(pow(refsatmC.POS_X - baseline[a].mX, 2) + pow(refsatmC.POS_Y - baseline[a].mY, 2) + pow(refsatmC.POS_Z - baseline[a].mZ, 2));//3代参考卫星到基站的距离
				double length3 = sqrt(pow(CHm.sat[i].POS_X - baseline[a].mX, 2) + pow(CHm.sat[i].POS_Y - baseline[a].mY, 2) + pow(CHm.sat[i].POS_Z - baseline[a].mZ, 2));//3代某颗卫星到基站的距离
				B[tw][0] = (CHr.sat[i].POS_X - baseline[a].rX) / length1 - (refsatrC.POS_X - baseline[a].rX) / length0;
				B[tw][1] = (CHr.sat[i].POS_Y - baseline[a].rY) / length1 - (refsatrC.POS_Y - baseline[a].rY) / length0;
				B[tw][2] = (CHr.sat[i].POS_Z - baseline[a].rZ) / length1 - (refsatrC.POS_Z - baseline[a].rZ) / length0;
				L[tw][0] = fai + lambdaB10_1 * NWC10_1[tw - bdsnum + 1 + junknum2 + grow][0] - ((length1 - length3) - (length0 - length2)) + error;
				P[tw][tw] = 0.09 + 0.09 / pow(sin(CHm.sat[i].E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(CHr.sat[i].E*PI / 180), 2);

				tw++;

			}
			//存储上历元的参考星的序号
			for (int i = 0; i < bdsnum; i++)
			{
				if (BDSm.sat[i].PRN == refsatmS.PRN)
				{
					refposkbds = i;
					break;
				}
			}
			R = P.InvertGaussJordan();
			CMatrix Qw(3, 3);
			Qw[0][0] = 100.0; Qw[1][1] = 100.0; Qw[2][2] = 100.0;
			coo.Kalman(I, F, B, L, R, Q, Qw, X);
			CMatrix NEU = TT * X;

			fprintf(p1, "%4d bdsnum:%4d chnum:%4d XYZ: %14.4f, %14.4f, %14.4f,NEU: %14.4f, %14.4f, %14.4f,  \n", b,bdsnum,chnum,X[0][0], X[1][0], X[2][0], NEU[0][0], NEU[1][0], NEU[2][0]);

			RMSN += NEU[0][0] * NEU[0][0];
			RMSE += NEU[1][0] * NEU[1][0];
			RMSU += NEU[2][0] * NEU[2][0];
			NUM++;
			tempepoch = baseline[a].m_epoch[b];//存储上个历元
			Cr = CHr;
			Cm = CHm;
			Br = BDSr;
			Bm = BDSm;

			TRACE("第%d个历元生成结束\n", b + 1);
		}
		RMSN = sqrt(RMSN / NUM);
		RMSE = sqrt(RMSE / NUM);
		RMSU = sqrt(RMSU / NUM);
		fprintf(p1, "RMSN: %13.4f,RMSE: %13.4f, RMSU:%13.4f, \n", RMSN, RMSE, RMSU);
		fclose(p1);
		fclose(p2);
		fclose(p3);
		fclose(p4);
	}
}

void GNSS::BDS23_long(vector<observe_ppp>baseline)//基于BDS-2和BDS-3的考虑系统间偏差的紧组合载波双差定位
{
	size_t filenum = baseline.size();
	for (int a = 0; a < filenum; a++)
	{
		double BLH[3] = { 0 };
		coo.OnXYZtoBLH(baseline[a].rX, baseline[a].rY, baseline[a].rZ, BLH);
		double B = BLH[0];
		double L = BLH[1];
		double H = BLH[2];
		CMatrix TT(3, 3);
		TT[0][0] = -sin(B)*cos(L);    TT[1][0] = -sin(L);   TT[2][0] = cos(B)*cos(L);
		TT[0][1] = -sin(B)*sin(L);    TT[1][1] = cos(L);    TT[2][1] = cos(B)*sin(L);
		TT[0][2] = cos(B);            TT[1][2] = 0;         TT[2][2] = sin(B);

		pppsat refsatmC;
		pppsat refsatrC;
		pppsat refsatmS;
		pppsat refsatrS;
		ppp_epoch tempepoch;//存储上个历元
		tempepoch = baseline[a].m_epoch[0];
		ppp_epoch Cm;
		ppp_epoch Cr;
		ppp_epoch Bm;
		ppp_epoch Br;

		

		FILE *p, *p1, *p2, *p3, *p4,*p5,*p6,*p7,*C1,*C2,*C3,*C4,*C6,*C7,*C8,*C9,*C10,*C11,*C12,*C13,*B19,*B20,*B21,*B22,*B23,*B24,*B25,*B26,*B27,*B28,*B29,*B30,*B32,*B33,*B34,*B35,*B36,*B37;
		CString ss1 = baseline[a].marker_name;
		ss1 = ss1.Left(4);
		CString s1;
		s1.Format("长基线紧组合a%s.txt", ss1);


		CString s = "长基线模糊度NWC.txt";

		CString s2 = "长基线北斗二代N14-5.txt";
		CString s3 = "长基线北斗二代N0-11.txt";
		CString s4 = "长基线北斗二代N1-10.txt";
		CString s5 = "长基线参考卫星.txt";
		CString s6 = "长基线北斗三代NWC01-1.txt";
		CString s7 = "长基线北斗三代NWC10-1.txt";

		CString c1 = "长基线C01.txt";
		CString c2 = "长基线C02.txt";
		CString c3 = "长基线C03.txt";
		CString c4 = "长基线C04.txt";
		CString c6 = "长基线C06.txt";
		CString c7=  "长基线C07.txt";
		CString c8 = "长基线C08.txt";
		CString c9 = "长基线C09.txt";
		CString c10 = "长基线C10.txt";
		CString c11 = "长基线C11.txt";
		CString c12 = "长基线C12.txt";
		CString c13 = "长基线C13.txt";
		
		CString b19 = "长基线B19.txt";
		CString b20 = "长基线B20.txt";
		CString b21= "长基线B21.txt";
		CString b22= "长基线B22.txt";
		CString b23= "长基线B23.txt";
		CString b24= "长基线B24.txt";
		CString b25= "长基线B25.txt";
		CString b26= "长基线B26.txt";
		CString b27 = "长基线B27.txt";
		CString b28 = "长基线B28.txt";
		CString b29 = "长基线B29.txt";
		CString b30 = "长基线B30.txt";
		CString b32 = "长基线B32.txt";
		CString b33 = "长基线B33.txt";
		CString b34 = "长基线B34.txt";
		CString b35 = "长基线B35.txt";
		CString b36 = "长基线B36.txt";
		CString b37 = "长基线B37.txt";

		p = fopen((LPCSTR)s, "w");
		p1 = fopen((LPCSTR)s1, "w");
		p2 = fopen((LPCSTR)s2, "w");
		p3 = fopen((LPCSTR)s3, "w");
		p4 = fopen((LPCSTR)s4, "w");
		p5 = fopen((LPCSTR)s5, "w");
		p6 = fopen((LPCSTR)s6, "w");
		p7 = fopen((LPCSTR)s7, "w");

		C1 = fopen((LPCSTR)c1, "w");     C2 = fopen((LPCSTR)c2, "w");      C3 = fopen((LPCSTR)c3, "w");
		C4 = fopen((LPCSTR)c4, "w");     C6 = fopen((LPCSTR)c6, "w");      C7 = fopen((LPCSTR)c7, "w");
		C8 = fopen((LPCSTR)c8, "w");     C9 = fopen((LPCSTR)c9, "w");      C10 = fopen((LPCSTR)c10, "w");
		C11 = fopen((LPCSTR)c11, "w");   C12 = fopen((LPCSTR)c12, "w");    C13 = fopen((LPCSTR)c13, "w");
		
		

		B19 = fopen((LPCSTR)b19, "w");    B20 = fopen((LPCSTR)b20, "w");     B21 = fopen((LPCSTR)b21, "w");
		B22 = fopen((LPCSTR)b22, "w");    B23 = fopen((LPCSTR)b23, "w");     B24 = fopen((LPCSTR)b24, "w");
		B25 = fopen((LPCSTR)b25, "w");    B26 = fopen((LPCSTR)b26, "w");     B27 = fopen((LPCSTR)b27, "w");
		B28 = fopen((LPCSTR)b28, "w");    B29 = fopen((LPCSTR)b29, "w");     B30 = fopen((LPCSTR)b30, "w");
		B32 = fopen((LPCSTR)b32, "w");    B33 = fopen((LPCSTR)b33, "w");     B34 = fopen((LPCSTR)b34, "w");
		B35 = fopen((LPCSTR)b35, "w");    B36 = fopen((LPCSTR)b36, "w");     B37 = fopen((LPCSTR)b37, "w");
		CMatrix X;
		CMatrix Q;
		CMatrix filter1;
		CMatrix filter2;
		size_t epochnum = baseline[a].m_epoch.size();
		int refposkchs;//存储上历元参考星的序号
		int refposkbds;
		double posr[3]; posr[0] = baseline[a].rX; posr[1] = baseline[a].rY; posr[2] = baseline[a].rZ;
		double posm[3]; posm[0] = baseline[a].mX; posm[1] = baseline[a].mY; posm[2] = baseline[a].mZ;
		coo.check(baseline[a].m_epoch, baseline[a].r_epoch);
		double RMSN = 0;
		double RMSE = 0;
		double RMSU = 0;
		int NUM = 0;
		double satnum2[14] ;//0-13对应C01-C14  //这一部分不能紧接着放在CMatirix xxx的下面，会报错
		for (int i = 0; i < 14; i++)
			satnum2[i] = 1.0;
		double satnum3[19];//0-18对应B19-B37
		for (int i = 0; i < 19; i++)
			satnum3[i] = 1.0;
		for (int b = 0; b < epochnum; b++)//历元
		{
			//cout << "第"<<b+1<<"个历元" << endl;
			ppp_epoch BDSm;
			ppp_epoch BDSr;
			ppp_epoch CHm;
			ppp_epoch CHr;
			int flag = 0;
			int junknum2 = 0;
			int junknum3 = 0;
			int bdsnum = 0;
			int chnum = 0;//bds-3
			

			for (int i = 0; i < baseline[a].m_epoch[b].sat.size(); i++)
			{
				if (baseline[a].m_epoch[b].sat[i].judge_use == 0)
				{

					if (baseline[a].m_epoch[b].sat[i].sattype == "C"&&baseline[a].m_epoch[b].sat[i].E > 15 /*&& baseline[a].m_epoch[b].sat[i].PRN != 5*/)
					{
						BDSm.sat.push_back(baseline[a].m_epoch[b].sat[i]);
						BDSr.sat.push_back(baseline[a].r_epoch[b].sat[i]);
						bdsnum++;
					}
					if (baseline[a].m_epoch[b].sat[i].sattype == "B"&&baseline[a].m_epoch[b].sat[i].E > 15)
					{
						CHm.sat.push_back(baseline[a].m_epoch[b].sat[i]);
						CHr.sat.push_back(baseline[a].r_epoch[b].sat[i]);
						chnum++;
					}

				}
			}
			
			
			if (b > 0)
			{
				//防止出现上历元参考星本历元未出现的情况,如果出现这种情况，转换矩阵会报错
				for (int i = 0; i < chnum; i++)
				{
				
					if (Cm.sat[refposkchs].PRN == CHm.sat[i].PRN)
					{
						flag++; break;
					}

				}
				for (int i = 0; i < bdsnum; i++)
				{
					
					if (Bm.sat[refposkbds].PRN == BDSm.sat[i].PRN)
					{
						flag++; break;
					}

				}
				
			}
			/////////////////


			//BDS-3
			CMatrix PC01_1(chnum - 1, 1);//伪距
			CMatrix faiC01_1(chnum - 1, 1);
			CMatrix faiC10_1(chnum - 1, 1);
			CMatrix NWC01_1(chnum - 1, 1);
			CMatrix NWC10_1(chnum - 1, 1);
			
			
	




			//计算BDS-3的系统内双差模糊度NWC
			pppsat tempsatmC = CHm.sat[0];//r,m的sattype和prn号和卫星位置都相等  所以用sattype或prn或APP_XYZ判断时用r或m都可
			pppsat tempsatrC = CHr.sat[0];
			double fch1=1561.098*1E+6;   //B1I C2
			double fch2 = 1268.52*1E+6;  //B3I C6
			double fch3 = 1176.45*1E+6;  //B2a C5
			double lambdaB01_1 = C_light / (fch2 - fch3);
			double lambdaB10_1 = C_light / (fch1 - fch3);
			
			for (int i = 0; i < chnum; i++)
			{
				if (CHr.sat[i].E > tempsatrC.E)
				{
					tempsatmC = CHm.sat[i];
					tempsatrC = CHr.sat[i];
				}
			}


			if (b == 0)//第一个历元
			{
				refsatmC = tempsatmC;
				refsatrC = tempsatrC;
				X.SetSize(4, 1);
				Q.SetSize(4, 4);
				for (int ii = 0; ii < 3; ii++)
				{
					Q[ii][ii] = 10.0*10.0;
				}
				Q[4][4] = 10000.0 * 10000.0;
				filter1.SetSize(chnum - 1, 1);
				filter2.SetSize(bdsnum - 1, 1);
			}
			else//其他历元
			{
				
				bool find = false;
				for (int i = 0; i < tempepoch.sat.size(); i++)
				{
					if (tempepoch.sat[i].sattype == tempsatmC.sattype&&tempepoch.sat[i].PRN == tempsatmC.PRN)//若参考星在上一历元中未出现，继续使用上一历元的参考星
					{
						find = true; break;
					}
				}
				if (find == true)
				{
					refsatmC = tempsatmC;
					refsatrC = tempsatrC;
				}
				//if (flag != 2)//若发生上历元BDS2,3参考星本历元未出现的情况
				//{
					CMatrix TT_matrix33(Cr.sat.size() - 1, Cr.sat.size() - 1);
					coo.GetGRCQValue(Cm, Cm, Cm.sat.size(), Cm.sat.size(), refsatmC.PRN, Cm.sat[refposkchs].PRN, TT_matrix33);//都带入上历元序列
					filter1 = TT_matrix33 * filter1;//上历元序列的参考星转换成本历元序列的参考星
					CMatrix TT_matrix3(chnum - 1, Cr.sat.size() - 1);
					coo.GetGRCQValue(CHm, Cm, chnum, Cm.sat.size(), refsatmC.PRN, refsatmC.PRN, TT_matrix3);//都带入本历元参考星
					filter1 = TT_matrix3 * filter1;
				//}
				/*else
				{
					CMatrix TT_matrix3(chnum - 1, Cr.sat.size() - 1);
					coo.GetGRCQValue(CHm, Cm, chnum, Cm.sat.size(), refsatmC.PRN, Cm.sat[refposkchs].PRN, TT_matrix3);
					filter1 = TT_matrix3 * filter1;
					
				}
				*/
			}//参考星选择完毕

			
			int w = 0;
			for (int i = 0; i < chnum; i++)
			{

				if (refsatmC.PRN == CHm.sat[i].PRN)
					continue;
				PC01_1[w][0] = 0.0402* ((CHr.sat[i].data[0] - CHm.sat[i].data[0]) - (refsatrC.data[0] - refsatmC.data[0])) + 0.3951* ((CHr.sat[i].data[4] - CHm.sat[i].data[4]) - (refsatrC.data[4] - refsatmC.data[4])) + 0.5647 * ((CHr.sat[i].data[1] - CHm.sat[i].data[1]) - (refsatrC.data[1] - refsatmC.data[1]));
				faiC01_1[w][0] = (C_light * ((CHr.sat[i].data[6] - CHm.sat[i].data[6]) - (refsatrC.data[6] - refsatmC.data[6])) - C_light * ((CHr.sat[i].data[3] - CHm.sat[i].data[3]) - (refsatrC.data[3] - refsatmC.data[3]))) / (fch2 - fch3);
				faiC10_1[w][0] = (C_light * ((CHr.sat[i].data[2] - CHm.sat[i].data[2]) - (refsatrC.data[2] - refsatmC.data[2])) - C_light * ((CHr.sat[i].data[3] - CHm.sat[i].data[3]) - (refsatrC.data[3] - refsatmC.data[3]))) / (fch1 - fch3);
				NWC01_1[w][0] = round((PC01_1[w][0] - faiC01_1[w][0]) / (C_light / (fch2 - fch3)));

				if (b == 0)
				{
					filter1[w][0] = -1.0 / lambdaB10_1 * (faiC10_1[w][0] - faiC01_1[w][0] - lambdaB01_1 * NWC01_1[w][0]);
				}
				else
				{
					double k = satnum3[CHm.sat[i].PRN - 19];
					double temp = -1.0 / lambdaB10_1 * (faiC10_1[w][0] - faiC01_1[w][0] - lambdaB01_1 * NWC01_1[w][0]);
					if (abs(round(temp) - temp) > 0.3 || ((abs(filter1[w][0] - temp) > 0.4) && (abs(filter1[w][0] - temp) < 1.5)))
					{
						junknum3++;
						CHm.sat[i].judge_use = 1;
						w++; continue;
					}
					if (abs(filter1[w][0] - temp) < 0.4)
					{
						filter1[w][0] = filter1[w][0] * k / (k + 1) + temp / (k + 1);
						satnum3[CHm.sat[i].PRN - 19]++;
					}
					else
					{
						filter1[w][0] = temp;
						satnum3[CHm.sat[i].PRN - 19] = 1;
					}
				}
				NWC10_1[w][0] = round(filter1[w][0]);
				//doubleNWC[w][0] = (PPC[w][0] - faiC[w][0]) / (C_light / (fch2 - fch3)) - NWC[w][0];
				if (abs(filter1[w][0] - NWC10_1[w][0]) > 0.15 || (b > 100 && satnum3[CHm.sat[i].PRN - 19] < 5))
				{
					junknum3++;
					CHm.sat[i].judge_use = 1;
				}
				w++;
			}
			for (int i = 19; i < 38; i++)//B19-B37中在本历元未出现的置1
			{
				if (coo.trip(CHm, i) == false)
					satnum3[i - 19] = 1.0;
			}
			CMatrix NNWC01_1(12, 1);
			for (int i = 0; i < chnum - 1; i++)
				NNWC01_1[i][0] = NWC01_1[i][0];
			fprintf(p6, "%4d 固定解NW：chnum-1：%4d  %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f  \n", b, chnum - 1, NNWC01_1[0][0], NNWC01_1[1][0], NNWC01_1[2][0], NNWC01_1[3][0], NNWC01_1[4][0], NNWC01_1[5][0], NNWC01_1[6][0], NNWC01_1[7][0], NNWC01_1[8][0], NNWC01_1[9][0], NNWC01_1[10][0], NNWC01_1[11][0]);
			CMatrix NNWC10_1(12, 1);

			for (int i = 0; i < chnum - 1; i++)
				NNWC10_1[i][0] = NWC10_1[i][0];
			fprintf(p7, "%4d 固定解NW：chnum-1：%4d  %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f  \n", b, chnum - 1, NNWC10_1[0][0], NNWC10_1[1][0], NNWC10_1[2][0], NNWC10_1[3][0], NNWC10_1[4][0], NNWC10_1[5][0], NNWC10_1[6][0], NNWC10_1[7][0], NNWC10_1[8][0], NNWC10_1[9][0], NNWC10_1[10][0], NNWC10_1[11][0]);

			
			
			/**********************************统计三代每一颗卫星的模糊度及偏差*********************************/
			/*
			int e = 0;
			for (int i = 0; i < chnum; i++)
			{
				if (refsatmC.PRN == CHm.sat[i].PRN)
					continue;
				if (CHr.sat[i].PRN == 19)
					fprintf(B19, " %14.4f  %14.4f\n ", NWC10_1[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 20)
					fprintf(B20, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 21)
					fprintf(B21, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 22)
					fprintf(B22, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 23)
					fprintf(B23, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 24)
					fprintf(B24, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 25)
					fprintf(B25, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 26)
					fprintf(B26, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 27)
					fprintf(B27, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 28)
					fprintf(B28, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 29)
					fprintf(B29, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 30)
					fprintf(B30, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 32)
					fprintf(B32, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 33)
					fprintf(B33, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 34)
					fprintf(B34, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 35)
					fprintf(B35, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 36)
					fprintf(B36, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				if (CHr.sat[i].PRN == 37)
					fprintf(B37, " %14.4f  %14.4f\n ", NWC[e][0], doubleNWC[e][0]);
				e++;
			}
			*/

			/********************************************************************************************************/






			CMatrix PP(bdsnum - 1, 1);//伪距
			CMatrix fai0_11(bdsnum - 1, 1);
			CMatrix fai14_5(bdsnum - 1, 1);
			CMatrix fai1_10(bdsnum - 1, 1);
			CMatrix N0_11(bdsnum - 1, 1);
			CMatrix N14_5(bdsnum - 1, 1);
			CMatrix N1_10(bdsnum - 1, 1);
			CMatrix doubleN0_11(bdsnum - 1, 1);
			CMatrix doubleN14_5(bdsnum - 1, 1);
			

			pppsat tempsatmS = BDSm.sat[0];//r,m的sattype和prn号和卫星位置都相等  所以用sattype或prn或APP_XYZ判断时用r或m都可
			pppsat tempsatrS = BDSr.sat[0];

			double fbds1 = 1561.098*1E+6;
			double fbds2 = 1207.14*1E+6;
			double fbds3 = 1268.52*1E+6;
			//double fch11=1575.42*1E+6;
			//double fch22= 1561.098*1E+6;
			double lambda1_10 = C_light / (fbds1 - fbds2);
			double lambda0_11 = C_light / (-fbds2 + fbds3);
			double lambda14_5 = C_light / (fbds1 +4* fbds2-5*fbds3);
			double lambdaC = C_light / (fbds1 - fbds2);
			
			
		
			
			for (int i = 0; i < bdsnum; i++)
			{
				if (BDSr.sat[i].E > tempsatrS.E)
				{
					tempsatmS = BDSm.sat[i];
					tempsatrS = BDSr.sat[i];
				}
			}


			if (b == 0)//第一个历元
			{
				refsatmS = tempsatmS;
				refsatrS = tempsatrS;
			}
			else//其他历元
			{
				bool find = false;
				for (int i = 0; i < tempepoch.sat.size(); i++)
				{
					if (tempepoch.sat[i].sattype == tempsatmS.sattype&&tempepoch.sat[i].PRN == tempsatmS.PRN)//若参考星在上一历元中未出现，继续使用上一历元的参考星
					{
						find = true; break;
					}
				}
				if (find == true)
				{
					refsatmS = tempsatmS;
					refsatrS = tempsatrS;
				}
				//if (flag != 2)
				//{
					CMatrix TT_matrix(Br.sat.size() - 1, Br.sat.size() - 1);
					coo.GetGRCQValue(Bm, Bm, Bm.sat.size(), Bm.sat.size(), refsatmS.PRN, Bm.sat[refposkbds].PRN, TT_matrix);//都带入上历元序列
					filter2 = TT_matrix * filter2;//上历元序列的参考星转换成本历元序列的参考星
					CMatrix TT_matrix2(bdsnum - 1, Br.sat.size() - 1);
					coo.GetGRCQValue(BDSm, Bm, bdsnum, Bm.sat.size(), refsatmS.PRN, refsatmS.PRN, TT_matrix2);//都带入本历元参考星
					filter2 = TT_matrix2 * filter2;
				//}
				/*else
				{
					CMatrix TT_matrix2(bdsnum - 1, Br.sat.size() - 1);
					coo.GetGRCQValue(BDSm, Bm, bdsnum, Bm.sat.size(), refsatmS.PRN, Bm.sat[refposkbds].PRN, TT_matrix2);
					filter2 = TT_matrix2 * filter2;
				}
				*/
			}//参考星选择完毕
			
			int tw = 0; 
			for (int i = 0; i < bdsnum; i++)
			{
				if (refsatmS.PRN == BDSm.sat[i].PRN)
					continue;
				PP[tw][0] = 0.0199*((BDSr.sat[i].data[0] - BDSm.sat[i].data[0]) - (refsatrS.data[0] - refsatmS.data[0])) + 0.5526* ((BDSr.sat[i].data[1] - BDSm.sat[i].data[1]) - (refsatrS.data[1] - refsatmS.data[1])) + 0.4275 * ((BDSr.sat[i].data[4] - BDSm.sat[i].data[4]) - (refsatrS.data[4] - refsatmS.data[4]));
				double P110 = (fbds1* ((BDSr.sat[i].data[0] - BDSm.sat[i].data[0]) - (refsatrS.data[0] - refsatmS.data[0])) + fbds2 * ((BDSr.sat[i].data[1] - BDSm.sat[i].data[1]) - (refsatrS.data[1] - refsatmS.data[1]))) / (fbds1 + fbds2);
				fai0_11[tw][0] = (-C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3])) + C_light * ((BDSr.sat[i].data[5] - BDSm.sat[i].data[5]) - (refsatrS.data[5] - refsatmS.data[5]))) / (-fbds2 + fbds3);
				fai1_10[tw][0] = (C_light * ((BDSr.sat[i].data[2] - BDSm.sat[i].data[2]) - (refsatrS.data[2] - refsatmS.data[2])) - C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3]))) / (fbds1 - fbds2);
				fai14_5[tw][0] = (C_light * ((BDSr.sat[i].data[2] - BDSm.sat[i].data[2]) - (refsatrS.data[2] - refsatmS.data[2])) + 4 * C_light * ((BDSr.sat[i].data[3] - BDSm.sat[i].data[3]) - (refsatrS.data[3] - refsatmS.data[3])) - 5 * C_light * ((BDSr.sat[i].data[5] - BDSm.sat[i].data[5]) - (refsatrS.data[5] - refsatmS.data[5]))) / (fbds1 + 4 * fbds2 - 5 * fbds3);
				N0_11[tw][0] = round((PP[tw][0] - fai0_11[tw][0]) / lambda0_11);
				double n0_11 = (PP[tw][0] - fai0_11[tw][0]) / lambda0_11 - N0_11[tw][0];
				if (b == 0)
				{
					filter2[tw][0] = -1.0 / lambda1_10 * (fai1_10[tw][0] - fai0_11[tw][0] - lambda0_11 * N0_11[tw][0]);
				}
				else
				{
					double k = satnum2[BDSm.sat[i].PRN - 1];

					double temp = -1.0 / lambda1_10 * (fai1_10[tw][0] - fai0_11[tw][0] - lambda0_11 * N0_11[tw][0]);
					if (abs(round(temp) - temp) > 0.3 || ((abs(filter2[tw][0] - temp) > 0.4) && (abs(filter2[tw][0] - temp) < 1.5)))
					{
						junknum2++;
						BDSm.sat[i].judge_use = 1;
						tw++; continue;
					}
					if (abs(filter2[tw][0] - temp) < 0.4)
					{
					
						filter2[tw][0] = filter2[tw][0] * k / (k + 1) + temp / (k + 1);
						satnum2[BDSm.sat[i].PRN - 1]++;

					}
					else
					{
						filter2[tw][0] = -1.0 / lambda1_10 * (fai1_10[tw][0] - fai0_11[tw][0] - lambda0_11 * N0_11[tw][0]);
						satnum2[BDSm.sat[i].PRN - 1] = 1;
					}
				}

				//N14_5[tw][0] = round(filter2[tw][0]);
				//N1_10[tw][0] = 5 * N0_11[tw][0] + N14_5[tw][0];
				N1_10[tw][0] = round(filter2[tw][0]);
				if (abs(filter2[tw][0] - N1_10[tw][0]) > 0.2 || abs(n0_11) > 0.2 || (b > 100 && satnum2[BDSm.sat[i].PRN - 1] < 5))
				{
					junknum2++;
					BDSm.sat[i].judge_use = 1;
				}
				tw++;
	
			}
			for (int i = 1; i < 15; i++)//C01-C14中在本历元未出现的置1
			{
				if (coo.trip(BDSm, i) == false)
					satnum2[i - 1] = 1.0;
			}
			
			CMatrix NN0_11(12, 1);
			for (int i = 0; i < bdsnum - 1; i++)
			{
				NN0_11[i][0] = N0_11[i][0];
			}
			fprintf(p3, "%14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f,\n", NN0_11[0][0], NN0_11[1][0], NN0_11[2][0], NN0_11[3][0], NN0_11[4][0], NN0_11[5][0], NN0_11[6][0], NN0_11[7][0], NN0_11[8][0], NN0_11[9][0], NN0_11[10][0], NN0_11[11][0]);

			CMatrix NN14_5(12, 1);
			for (int i = 0; i < bdsnum - 1; i++)
			{
				NN14_5[i][0] = N14_5[i][0];
			}
			fprintf(p2, "%14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f,\n", NN14_5[0][0], NN14_5[1][0], NN14_5[2][0], NN14_5[3][0], NN14_5[4][0], NN14_5[5][0], NN14_5[6][0], NN14_5[7][0], NN14_5[8][0], NN14_5[9][0], NN14_5[10][0], NN14_5[11][0]);
			
			CMatrix NN1_10(12, 1);
			for (int i = 0; i < bdsnum - 1; i++)
			{
				NN1_10[i][0] = N1_10[i][0];
			}
			fprintf(p4, "%14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f, %14.1f\n", NN1_10[0][0], NN1_10[1][0], NN1_10[2][0], NN1_10[3][0], NN1_10[4][0], NN1_10[5][0], NN1_10[6][0], NN1_10[7][0], NN1_10[8][0], NN1_10[9][0], NN1_10[10][0], NN1_10[11][0]);
			
			
			/**********************************统计二代每一颗卫星的模糊度及偏差*********************************/
			int n=0;
			for (int i = 0; i < bdsnum; i++)
			{
				if (refsatmS.PRN == BDSm.sat[i].PRN)
					continue;
				if (BDSr.sat[i].PRN == 1)
					fprintf(C1, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 2)
					fprintf(C2, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 3)
					fprintf(C3, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 4)
					fprintf(C4, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 6)
					fprintf(C6, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 7)
					fprintf(C7, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 8)
					fprintf(C8, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 9)
					fprintf(C9, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 10)
					fprintf(C10, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 11)
					fprintf(C11, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 12)
					fprintf(C12, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				if (BDSr.sat[i].PRN == 13)
					fprintf(C13, "N0-11  %14.4f  %14.4f  N14-5  %14.4f  %14.4f\n", N0_11[n][0], doubleN0_11[n][0], N14_5[n][0], doubleN14_5[n][0]);
				n++;
			}
			
			
			/********************************************************************************************************/
			//BDS-23
			CMatrix B(chnum + bdsnum - 1-junknum2-junknum3, 4);
			CMatrix L(chnum + bdsnum - 1 - junknum2 - junknum3, 1);
			CMatrix P(chnum + bdsnum - 1 - junknum2 - junknum3, chnum + bdsnum - 1 - junknum2 - junknum3);
			CMatrix I(4, 4);
			CMatrix F(4, 4);
			I.ones(4);
			F.ones(4);
			CMatrix R(chnum + bdsnum - 1 - junknum2 - junknum3, chnum + bdsnum - 1 - junknum2 - junknum3);
			
			 tw = 0; int grow = 0;
			for (int i = 0; i < bdsnum; i++)//bds2系统内双差
			{
				
				if (refsatmS.PRN == BDSm.sat[i].PRN)
					continue;
				if (BDSm.sat[i].judge_use == 1)
				{
					grow++;
					continue;
				}
				double TGD = (BDSr.sat[i].TGD - BDSm.sat[i].TGD) - (refsatrS.TGD - refsatmS.TGD);
				double TROP = (BDSr.sat[i].Trop_Delay - BDSm.sat[i].Trop_Delay) - (refsatrS.Trop_Delay - refsatmS.Trop_Delay);
				double relate = (BDSr.sat[i].Relat - BDSm.sat[i].Relat) - (refsatrS.Relat - refsatmS.Relat);
				double satclock = (BDSr.sat[i].Sat_clock - BDSm.sat[i].Sat_clock) - (refsatrS.Sat_clock - refsatmS.Sat_clock);
				double sagnac = (BDSr.sat[i].Sagnac - BDSm.sat[i].Sagnac) - (refsatrS.Sagnac - refsatmS.Sagnac);
				//double antenna_height1= (BDSr.sat[i].Antenna_Height - BDSm.sat[i].Antenna_Height) - (refsatrS.Antenna_Height - refsatmS.Antenna_Height);
				double sat_antenna1 = (BDSr.sat[i].Sat_Antenna - BDSm.sat[i].Sat_Antenna) - (refsatrS.Sat_Antenna - refsatmS.Sat_Antenna);
				double error = - TGD - TROP + relate + satclock - sagnac + sat_antenna1;

				double length0 = sqrt(pow(refsatrS.POS_X - baseline[a].rX, 2) + pow(refsatrS.POS_Y - baseline[a].rY, 2) + pow(refsatrS.POS_Z - baseline[a].rZ, 2));//参考卫星到流动站的距离
				double length1 = sqrt(pow(BDSr.sat[i].POS_X - baseline[a].rX, 2) + pow(BDSr.sat[i].POS_Y - baseline[a].rY, 2) + pow(BDSr.sat[i].POS_Z - baseline[a].rZ, 2));//某颗卫星到流动站的距离
				double length2 = sqrt(pow(refsatmS.POS_X - baseline[a].mX, 2) + pow(refsatmS.POS_Y - baseline[a].mY, 2) + pow(refsatmS.POS_Z - baseline[a].mZ, 2));//参考卫星到基站的距离
				double length3 = sqrt(pow(BDSm.sat[i].POS_X - baseline[a].mX, 2) + pow(BDSm.sat[i].POS_Y - baseline[a].mY, 2) + pow(BDSm.sat[i].POS_Z - baseline[a].mZ, 2));//某颗卫星到基站的距离
				B[tw][0] = (BDSr.sat[i].POS_X - baseline[a].rX) / length1 - (refsatrS.POS_X - baseline[a].rX) / length0;
				B[tw][1] = (BDSr.sat[i].POS_Y - baseline[a].rY) / length1 - (refsatrS.POS_Y - baseline[a].rY) / length0;
				B[tw][2] = (BDSr.sat[i].POS_Z - baseline[a].rZ) / length1 - (refsatrS.POS_Z - baseline[a].rZ) / length0;
				//double N = N1_10[tw + grow][0];
				L[tw][0] = fai1_10[tw+grow][0] + lambda1_10 * N1_10[tw+grow][0] - ((length1 - length3) - (length0 - length2)) + error;
				P[tw][tw] = 0.09 + 0.09 / pow(sin(BDSr.sat[i].E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(BDSm.sat[i].E*PI / 180), 2);
				tw++;
			}
			grow = 0;
			for (int i = 0; i < chnum; i++)
			{
				
				if (refsatmC.PRN == CHm.sat[i].PRN)
					continue;
				if (CHm.sat[i].judge_use == 1)
				{
					grow++;
					continue;
				}
				double fai = lambdaB10_1 * ( (CHr.sat[i].data[2] - CHm.sat[i].data[2])- (CHr.sat[i].data[3] - CHm.sat[i].data[3])) - lambdaC * ((refsatrS.data[2] - refsatmS.data[2]) - (refsatrS.data[3] - refsatmS.data[3]));
	

				double length0 = sqrt(pow(refsatrS.POS_X - baseline[a].rX, 2) + pow(refsatrS.POS_Y - baseline[a].rY, 2) + pow(refsatrS.POS_Z - baseline[a].rZ, 2));//2代参考卫星到流动站的距离
				double length1 = sqrt(pow(CHr.sat[i].POS_X - baseline[a].rX, 2) + pow(CHr.sat[i].POS_Y - baseline[a].rY, 2) + pow(CHr.sat[i].POS_Z - baseline[a].rZ, 2));//3代某颗卫星到流动站的距离
				double length2 = sqrt(pow(refsatmS.POS_X - baseline[a].mX, 2) + pow(refsatmS.POS_Y - baseline[a].mY, 2) + pow(refsatmS.POS_Z - baseline[a].mZ, 2));//2代参考卫星到基站的距离
				double length3 = sqrt(pow(CHm.sat[i].POS_X - baseline[a].mX, 2) + pow(CHm.sat[i].POS_Y - baseline[a].mY, 2) + pow(CHm.sat[i].POS_Z - baseline[a].mZ, 2));//3代某颗卫星到基站的距离
				B[tw][0] = (CHr.sat[i].POS_X - baseline[a].rX) / length1 - (refsatrS.POS_X - baseline[a].rX) / length0;
				B[tw][1] = (CHr.sat[i].POS_Y - baseline[a].rY) / length1 - (refsatrS.POS_Y - baseline[a].rY) / length0;
				B[tw][2] = (CHr.sat[i].POS_Z - baseline[a].rZ) / length1 - (refsatrS.POS_Z - baseline[a].rZ) / length0;
				B[tw][3] = lambdaB10_1;
				double TGD = (CHr.sat[i].TGD - CHm.sat[i].TGD) -  (refsatrS.TGD - refsatmS.TGD);
				double TROP = (CHr.sat[i].Trop_Delay - CHm.sat[i].Trop_Delay) - (refsatrS.Trop_Delay - refsatmS.Trop_Delay);
				double relate = (CHr.sat[i].Relat - CHm.sat[i].Relat) - (refsatrS.Relat - refsatmS.Relat);
				double satclock = (CHr.sat[i].Sat_clock - CHm.sat[i].Sat_clock) - (refsatrS.Sat_clock - refsatmS.Sat_clock);
				double sagnac = (CHr.sat[i].Sagnac - CHm.sat[i].Sagnac) - (refsatrS.Sagnac - refsatmS.Sagnac);
				//double antenna_height1= (BDSr.sat[i].Antenna_Height - BDSm.sat[i].Antenna_Height) - (refsatrS.Antenna_Height - refsatmS.Antenna_Height);
				double sat_antenna1 = (CHr.sat[i].Sat_Antenna - CHm.sat[i].Sat_Antenna) - (refsatrS.Sat_Antenna - refsatmS.Sat_Antenna);
				double error = -TGD - TROP + relate + satclock - sagnac + sat_antenna1;

				double length = ((length1 - length3) - (length0 - length2));
				//double N = NWC10_1[tw - bdsnum + 1 + junknum2 + grow][0];
				L[tw][0] = fai + lambdaB10_1 * NWC10_1[tw-bdsnum+1+junknum2+grow][0] - length +error;
				P[tw][tw] = 0.09 + 0.09 / pow(sin(CHr.sat[i].E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(CHm.sat[i].E*PI / 180), 2);

				tw++;

			}
			
			//double fai = lambdaB1_10 * ((refsatrC.data[2] - refsatmC.data[2]) -3*(refsatrC.data[6] - refsatmC.data[6])+2* (refsatrC.data[3] - refsatmC.data[3])) - lambdaC * ((refsatrS.data[2] - refsatmS.data[2]) - (refsatrS.data[3] - refsatmS.data[3]));
			double faii = lambdaB10_1 * ((refsatrC.data[2] - refsatmC.data[2])-(refsatrC.data[3] - refsatmC.data[3])) - lambdaC * ((refsatrS.data[2] - refsatmS.data[2]) - (refsatrS.data[3] - refsatmS.data[3]));
			double length00 = sqrt(pow(refsatrS.POS_X - baseline[a].rX, 2) + pow(refsatrS.POS_Y - baseline[a].rY, 2) + pow(refsatrS.POS_Z - baseline[a].rZ, 2));//2代参考卫星到流动站的距离
			double length11 = sqrt(pow(refsatrC.POS_X - baseline[a].rX, 2) + pow(refsatrC.POS_Y - baseline[a].rY, 2) + pow(refsatrC.POS_Z - baseline[a].rZ, 2));//3代参考卫星到流动站的距离
			double length22 = sqrt(pow(refsatmS.POS_X - baseline[a].mX, 2) + pow(refsatmS.POS_Y - baseline[a].mY, 2) + pow(refsatmS.POS_Z - baseline[a].mZ, 2));//2代参考卫星到基站的距离
			double length33 = sqrt(pow(refsatmC.POS_X - baseline[a].mX, 2) + pow(refsatmC.POS_Y - baseline[a].mY, 2) + pow(refsatmC.POS_Z - baseline[a].mZ, 2));//3代参考卫星到基站的距离
			B[tw][0] = (refsatrC.POS_X - baseline[a].rX) / length11 - (refsatrS.POS_X - baseline[a].rX) / length00;
			B[tw][1] = (refsatrC.POS_Y - baseline[a].rY) / length11 - (refsatrS.POS_Y - baseline[a].rY) / length00;
			B[tw][2] = (refsatrC.POS_Z - baseline[a].rZ) / length11 - (refsatrS.POS_Z - baseline[a].rZ) / length00;
			B[tw][3] = lambdaB10_1;
			double TGD = (refsatrC.TGD - refsatmC.TGD) - (refsatrS.TGD - refsatmS.TGD);
			double TROP = (refsatrC.Trop_Delay - refsatmC.Trop_Delay) - (refsatrS.Trop_Delay - refsatmS.Trop_Delay);
			double relate = (refsatrC.Relat - refsatmC.Relat) - (refsatrS.Relat - refsatmS.Relat);
			double satclock = (refsatrC.Sat_clock - refsatmC.Sat_clock) - (refsatrS.Sat_clock - refsatmS.Sat_clock);
			double sagnac = (refsatrC.Sagnac - refsatmC.Sagnac) - (refsatrS.Sagnac - refsatmS.Sagnac);
			//double antenna_height1= (BDSr.sat[i].Antenna_Height - BDSm.sat[i].Antenna_Height) - (refsatrS.Antenna_Height - refsatmS.Antenna_Height);
			double sat_antenna1 = (refsatrC.Sat_Antenna - refsatmC.Sat_Antenna) - (refsatrS.Sat_Antenna - refsatmS.Sat_Antenna);
			double error = -TGD - TROP + relate + satclock - sagnac + sat_antenna1;
			L[tw][0] = faii - ((length11 - length33) - (length00- length22)) + error;
			P[tw][tw] = 0.09 + 0.09 / pow(sin(refsatmC.E*PI / 180), 2) + 0.09 + 0.09 / pow(sin(refsatrC.E*PI / 180), 2);
			if (b == 0)
			{
				X[0][0] = 0;
				X[1][0] = 0;
				X[2][0] = 0;
				X[3][0] = L[ bdsnum - 1 - junknum2][0] / lambdaB10_1;
			}
			else
			{
				int x = 0; double CHN=0, BDSN=0;
				for (int i = 0; i < chnum; i++)
				{
					if (refsatmC.PRN == CHm.sat[i].PRN)
						continue;
					if (Cm.sat[refposkchs].PRN == CHm.sat[i].PRN)
						CHN = NWC10_1[x][0];
					x++;
				}
				x = 0;
				for (int i = 0; i < bdsnum; i++)
				{
					if (refsatmS.PRN == BDSm.sat[i].PRN)
						continue;
					if (Bm.sat[refposkbds].PRN == BDSm.sat[i].PRN)
						BDSN = N1_10[x][0];
					x++;
				}
				if (Br.sat[refposkbds].PRN != refsatmS.PRN&&flag==2)//若bds-2或bds-3参考星改变
				{
					X[3][0] = X[3][0] - lambda1_10 / lambdaB10_1 * BDSN;
				}
				if (Cm.sat[refposkchs].PRN != refsatmC.PRN&&flag == 2)
				{
					X[3][0] = X[3][0] + CHN;
				}
				if (flag != 2)//flag不为2代表着上历元参考本历元未出现  
				{
					double sum = 0;
					for (int i = 0; i < chnum - junknum3; i++)
						sum += L[bdsnum - 1 - junknum2+i][0];
					sum = sum / (chnum - junknum3);
					X[3][0] = sum / lambdaB10_1;
				}
			}
			R = P.InvertGaussJordan();
			CMatrix Qw(4, 4);
			Qw[0][0] = 100.0; Qw[1][1] = 100.0; Qw[2][2] = 100.0; Qw[3][3] = 0.05*0.05;
			coo.Kalman(I, F, B, L, R, Q, Qw, X);
			CMatrix NEU = TT * X;

			fprintf(p1, "XYZ: %14.4f, %14.4f, %14.4f,NEU: %14.4f, %14.4f, %14.4f, 系统间偏差：%14.4f, \n", X[0][0], X[1][0], X[2][0], NEU[0][0], NEU[1][0], NEU[2][0], X[3][0]);
			fprintf(p5," %5d  bdsnum：%4d ,chnum：%4d  C%4d,  B%4d,\n",b,bdsnum,chnum,refsatmS.PRN, refsatmC.PRN);

			RMSN += NEU[0][0] * NEU[0][0];
			RMSE += NEU[1][0] * NEU[1][0];
			RMSU += NEU[2][0] * NEU[2][0];
			NUM++;
			tempepoch = baseline[a].m_epoch[b];//存储上个历元
			Cr = CHr;
			Cm = CHm;
			Br = BDSr;
			Bm = BDSm;
			//存储上历元的参考星的序号
			for (int i = 0; i < bdsnum; i++)
			{
				if (BDSm.sat[i].PRN == refsatmS.PRN)
				{
					refposkbds = i;
					break;
				}
			}
			//存储上历元的参考星的序号
			for (int i = 0; i < chnum; i++)
			{
				if (CHm.sat[i].PRN == refsatmC.PRN)
				{
					refposkchs = i;
					break;
				}
			}
			TRACE("第%d个历元生成结束\n", b + 1);

		}
		RMSN = sqrt(RMSN / NUM);
		RMSE = sqrt(RMSE / NUM);
		RMSU = sqrt(RMSU / NUM);
		
		fprintf(p1, "RMSN: %13.4f,RMSE: %13.4f, RMSU:%13.4f, \n", RMSN, RMSE, RMSU);
		fclose(p);     fclose(p1);    fclose(p2);      fclose(p3);     fclose(p4);  fclose(p5);
		fclose(C1);     fclose(C2);    fclose(C3);     fclose(C4);      fclose(C6);   fclose(C7);    fclose(C8);   fclose(C9);   fclose(C10);    fclose(C11);    fclose(C12);    fclose(C13);
		fclose(B19);  fclose(B20);  fclose(B21);   fclose(B22); fclose(B23); fclose(B24);  fclose(B25);  fclose(B26);  fclose(B27); fclose(B28);  fclose(B29);  fclose(B30); fclose(B32); fclose(B33); fclose(B34); fclose(B35); fclose(B36); fclose(B37);
		
		
		
		
	}
}