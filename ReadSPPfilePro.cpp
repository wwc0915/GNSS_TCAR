#include "stdafx.h"
#include <fstream>     //文件流 
#include <string>      //字符串
#include <vector>      //容器
#include "Myheader.h"        //数据结构
#include "Matrix.h"          //矩阵类
#include <assert.h>          //断言，触发报错用
using namespace std;

#define READ_TXT 0

//声明外部函数声明，在ReadOfilePro.cpp中定义
void OnXYZtoBLH(double X, double Y, double Z, double BLH[]);

//声明外部GetGPSTime函数，在ReadPfile.cpp中定义
double GetGPSTime(int year, int month, int day, int hour, int minute, double second, int& dayofy);

//函数功能：针对本人自定义的spp文件格式，可修改
bool ReadSPP_KIN(string filename, observe_spp& sppfile)
{
	fstream infile;
	infile.open(filename, ios::in);
	if (!infile)
	{
		printf("文件打开失败");
		return false;
	}
#if READ_TXT
	fstream infileTxt;
	infileTxt.open("snfx_SPP.txt", ios::in);
	if (!infileTxt)
	{
		printf("文件打开失败");
		return false;
	}
#endif
	string strs;
	string str2;

	double pi = 3.14159265;
	double aa = 6378137.000, bb = 6356752.3142, ee = 0.006694379990;

	do
	{
		getline(infile, strs);
		if (strs.substr(4, 11) == "MARKER_NAME")
		{
			str2 = strs.substr(17, 4);
			sppfile.marker_name = str2;
			continue;
		}
		
		if (strs.substr(7, 8) == "INTERVAL")
		{
			str2 = strs.substr(17, 3);
			sppfile.INTERVAL = atof(str2.c_str());
			continue;
		}

		if (strs.substr(0, 15) == "APPROX_POSITION")
		{
			str2 = strs.substr(17, 15);
			sppfile.APP_X = atof(str2.c_str());
			str2 = strs.substr(33, 15);
			sppfile.APP_Y = atof(str2.c_str());
			str2 = strs.substr(49, 15);
			sppfile.APP_Z = atof(str2.c_str());

			//****计算测站大地坐标**********
			double BLH[3] = {0};
			OnXYZtoBLH(sppfile.APP_X, sppfile.APP_Y, sppfile.APP_Z, BLH);

			sppfile.B = BLH[0];
			sppfile.L = BLH[1];
			sppfile.DH = BLH[2];

			continue;
		}
	} while (strs.substr(5, 10) != "Clock type");                           //头文件读取结束

	getline(infile, strs);//读取GLONASS卫星标志号
	
	//连续读取3行，下面直接读取到观测数据
	getline(infile, strs);
	getline(infile, strs);
	getline(infile, strs);

	//开始读取每个历元的数据
	int i = 0;
	while (getline(infile, strs))
	{
		Obs_epoch temp_epochp;
		temp_epochp.liyuan_num = i;

		str2 = strs.substr(17, 3);
		temp_epochp.sat_num = atoi(str2.c_str());

		//时间信息
		int ye, mo, da, ho, mi, se;
		double second;
		str2 = strs.substr(43, 4);
		ye = atoi(str2.c_str());
		str2 = strs.substr(48, 2);
		mo = atoi(str2.c_str());
		str2 = strs.substr(51, 2);
		da = atoi(str2.c_str());
		str2 = strs.substr(54, 2);
		ho = atoi(str2.c_str());
		str2 = strs.substr(57, 2);
		mi = atoi(str2.c_str());
		str2 = strs.substr(60, 10);
		second = atoi(str2.c_str());
		se = floor(second + 0.5);

		temp_epochp.year = ye;
		temp_epochp.month = mo;
		temp_epochp.day = da;
		temp_epochp.hour = ho;
		temp_epochp.minute = mi;
		temp_epochp.second = (double)se;

		//年积日
		int dayofy;

		double gpstime = GetGPSTime(ye, mo, da, ho, mi, se, dayofy);
		temp_epochp.GPSTIME = gpstime;

		//读取历元基本信息
		str2 = strs.substr(75, 8);
		temp_epochp.ZTD = atof(str2.c_str());

#if	READ_TXT
		string tmpStr, str3;
		getline(infileTxt, tmpStr);
		str3 = tmpStr.substr(5, 14);
		temp_epochp.posX = atof(str3.c_str());

		str3 = tmpStr.substr(20, 12);
		temp_epochp.posY = atof(str3.c_str());

		str3 = tmpStr.substr(33, 13);
		temp_epochp.posZ = atof(str3.c_str());
#else
		str2 = strs.substr(85, 14);
		temp_epochp.posX = atof(str2.c_str());

		str2 = strs.substr(101, 14);
		temp_epochp.posY = atof(str2.c_str());

		str2 = strs.substr(117, 14);
		temp_epochp.posZ = atof(str2.c_str());
#endif


		for (int k = 0; k<temp_epochp.sat_num; k++)
		{
			Sat temp_sat;
			getline(infile, strs);

			str2 = strs.substr(0, 1);
			temp_sat.sattype = str2.c_str();

			str2 = strs.substr(1, 2);
			temp_sat.numofsat = atoi(str2.c_str());

			str2 = strs.substr(4, 1);
			temp_sat.judge_use = atoi(str2.c_str());

			//卫星的位置
			str2 = strs.substr(6, 15);
			temp_sat.POS_X = atof(str2.c_str());
			str2 = strs.substr(22, 15);
			temp_sat.POS_Y = atof(str2.c_str());
			str2 = strs.substr(38, 15);
			temp_sat.POS_Z = atof(str2.c_str());

			//卫星钟、高度角、方位角
			str2 = strs.substr(54, 15);
			temp_sat.Sat_clock = atof(str2.c_str());
			str2 = strs.substr(70, 15);
			temp_sat.E = atof(str2.c_str());
			str2 = strs.substr(86, 15);
			temp_sat.A = atof(str2.c_str());

			//伪距和载波观测值（双频4个）
			str2 = strs.substr(102, 15);
			temp_sat.data[0] = atof(str2.c_str());
			str2 = strs.substr(118, 15);
			temp_sat.data[1] = atof(str2.c_str());
			str2 = strs.substr(134, 15);
			temp_sat.data[2] = atof(str2.c_str());
			str2 = strs.substr(150, 15);
			temp_sat.data[3] = atof(str2.c_str());

			//BDS-2,GPS三频观测值
			if (temp_sat.sattype == "C" || temp_sat.sattype == "G")
			{
				str2 = strs.substr(326, 15);
				temp_sat.data[4] = atof(str2.c_str());
				str2 = strs.substr(342, 15);
				temp_sat.data[5] = atof(str2.c_str());
			}

			//GLONASS卫星的FDMA频率号
			if (temp_sat.sattype == "R")
			{
				str2 = strs.substr(326, 15);
				temp_sat.fre_glo = atof(str2.c_str());
			}

			//Galileo，BDS-3的三频，四频观测值，如有
			if ((temp_sat.sattype == "E" || temp_sat.sattype == "B") && strs.length() > 380)
			{
				str2 = strs.substr(326, 15);
				temp_sat.data[4] = atof(str2.c_str());
				str2 = strs.substr(342, 15);
				temp_sat.data[5] = atof(str2.c_str());
				str2 = strs.substr(358, 15);
				temp_sat.data[6] = atof(str2.c_str());
				str2 = strs.substr(374, 15);
				temp_sat.data[7] = atof(str2.c_str());
			}

			//各类误差改正
			str2 = strs.substr(166, 15);
			temp_sat.Trop_Delay = atof(str2.c_str());  //总延迟（倾斜，已投影）
			str2 = strs.substr(182, 15);
			temp_sat.Trop_Map = atof(str2.c_str());
			str2 = strs.substr(198, 15);
			temp_sat.Relat = atof(str2.c_str());
			str2 = strs.substr(214, 15);
			temp_sat.Sagnac = atof(str2.c_str());
			str2 = strs.substr(230, 15);
			temp_sat.TGD = atof(str2.c_str());
			str2 = strs.substr(246, 15);
			temp_sat.Antenna_Height = atof(str2.c_str());
			str2 = strs.substr(262, 15);
			temp_sat.Sat_Antenna = atof(str2.c_str());
			str2 = strs.substr(278, 15);
			temp_sat.OffsetL1 = atof(str2.c_str());
			str2 = strs.substr(294, 15);
			temp_sat.OffsetL2 = atof(str2.c_str());
			str2 = strs.substr(310, 15);
			temp_sat.Windup = atof(str2.c_str());

			temp_epochp.sat.push_back(temp_sat);
		}

		temp_epochp.sat_num = temp_epochp.sat.size();
		sppfile.epoch.push_back(temp_epochp);
		temp_epochp.sat.clear();
		printf("历元%d读入\n", i+1);
		i++;
	}//while
	sppfile.liyuan_num = i;
	infile.close();	
#if READ_TXT
	infileTxt.close();
#endif
	printf("spp站：%s读取结束, 总历元数：%5d\n", sppfile.marker_name.c_str(), sppfile.liyuan_num);
}

//函数功能：对基站和流动站卫星观测数据进行共视处理
void EpochCommonPro(Obs_epoch tempData1, Obs_epoch tempData2, Obs_epoch &useData)
{
	double CutEle = 10;
	vector<int> SatPRN_GPS;
	vector<int> SatPRN_BDS2;
	vector<int> SatPRN_BDS3;
	vector<int> SatPRN_Galileo;

	int PRNGPS[32] = { 0 };
	int PRNBDS2[20] = { 0 };
	int PRNBDS3[40] = { 0 };
	int PRNGalileo[40] = { 0 };

	for (auto& i : tempData1.sat)
	{
		if (i.sattype == "G" && i.E > CutEle && i.judge_use == 0 && i.data[0] > 0 && i.data[1] > 0 && i.data[2] > 0 && i
			.data[3] > 0 && i.data[4] > 0 && i.data[5] > 0)
			PRNGPS[i.numofsat - 1] += 1;
		if (i.sattype == "C" && i.E > CutEle && i.judge_use == 0 && i.data[0] > 0 && i.data[1] > 0 && i.data[2] > 0 && i
			.data[3] > 0 && i.data[4] > 0 && i.data[5] > 0)
			PRNBDS2[i.numofsat - 1] += 1;
		if (i.sattype == "B" && i.E > CutEle && i.judge_use == 0 && i.data[0] > 0 && i.data[1] > 0 && i.data[2] > 0 && i
			.data[3] > 0 && i.data[4] > 0 && i.data[5] > 0)
			PRNBDS3[i.numofsat - 1] += 1;
		if (i.sattype == "E" && i.E > CutEle && i.judge_use == 0 && i.data[0] > 0 && i.data[1] > 0 && i.data[2] > 0 && i
			.data[3] > 0 && i.data[4] > 0 && i.data[5] > 0)
			PRNGalileo[i.numofsat - 1] += 1;
	}
	for (auto& i : tempData2.sat)
	{
		if (i.sattype == "G" && i.E > CutEle && i.judge_use == 0 && i.data[0] > 0 && i.data[1] > 0 && i.data[2] > 0 && i
			.data[3] > 0 && i.data[4] > 0 && i.data[5] > 0)
			PRNGPS[i.numofsat - 1] += 1;
		if (i.sattype == "C" && i.E > CutEle && i.judge_use == 0 && i.data[0] > 0 && i.data[1] > 0 && i.data[2] > 0 && i
			.data[3] > 0 && i.data[4] > 0 && i.data[5] > 0)
			PRNBDS2[i.numofsat - 1] += 1;
		if (i.sattype == "B" && i.E > CutEle && i.judge_use == 0 && i.data[0] > 0 && i.data[1] > 0 && i.data[2] > 0 && i
			.data[3] > 0 && i.data[4] > 0 && i.data[5] > 0)
			PRNBDS3[i.numofsat - 1] += 1;
		if (i.sattype == "E" && i.E > CutEle && i.judge_use == 0 && i.data[0] > 0 && i.data[1] > 0 && i.data[2] > 0 && i
			.data[3] > 0 && i.data[4] > 0 && i.data[5] > 0)
			PRNGalileo[i.numofsat - 1] += 1;
	}

	int i = 0;
	for (i = 0; i < 32; i++)
		if (PRNGPS[i] == 2)
			SatPRN_GPS.push_back(i + 1);
	for (i = 0; i < 20; i++)
		if (PRNBDS2[i] == 2)
			SatPRN_BDS2.push_back(i + 1);
	for (i = 0; i < 40; i++)
		if (PRNBDS3[i] == 2)
			SatPRN_BDS3.push_back(i + 1);
	for (i = 0; i < 40; i++)
		if (PRNGalileo[i] == 2)
			SatPRN_Galileo.push_back(i + 1);


	for (i = 0; i < 32; i++)
	{
		if (PRNGPS[i] == 2)
		{
			for (auto& j : tempData1.sat)
			{
				if (j.numofsat == i + 1 && j.sattype == "G")
				{
					useData.sat1.push_back(j);
					break;
				}
			}
		}
	}

	for (i = 0; i < 20; i++)
	{
		if (PRNBDS2[i] == 2)
		{
			for (auto& j : tempData1.sat)
			{
				if (j.numofsat == i + 1 && j.sattype == "C")
				{
					useData.sat1.push_back(j);
					break;
				}
			}
		}
	}

	for (i = 0; i < 40; i++)
	{
		if (PRNBDS3[i] == 2)
		{
			for (auto& j : tempData1.sat)
			{
				if (j.numofsat == i + 1 && j.sattype == "B")
				{
					useData.sat1.push_back(j);
					break;
				}
			}
		}
	}

	for (i = 0; i < 40; i++)
	{
		if (PRNGalileo[i] == 2)
		{
			for (auto& j : tempData1.sat)
			{
				if (j.numofsat == i + 1 && j.sattype == "E")
				{
					useData.sat1.push_back(j);
					break;
				}
			}
		}
	}

	int TotSatNum = SatPRN_GPS.size() + SatPRN_BDS2.size() + SatPRN_BDS3.size() + SatPRN_Galileo.size();

	if (useData.sat1.size() != TotSatNum)
		assert(0);

	for (size_t k = 0; k < useData.sat1.size(); k++)
	{
		for (auto& j : tempData2.sat)
		{
			if (j.sattype == useData.sat1[k].sattype && j.numofsat == useData.sat1[k].numofsat)
			{
				useData.sat2.push_back(j);
				break;
			}
		}
	}

	//匹配检核
	if (useData.sat1.size() != useData.sat2.size())
		assert(0);

	for (size_t k = 0; k < useData.sat1.size(); k++)
	{
		if (useData.sat1[k].sattype != useData.sat2[k].sattype)
			assert(0);
		if (useData.sat1[k].numofsat != useData.sat2[k].numofsat)
			assert(0);
	}

	useData.liyuan_num = tempData2.liyuan_num;
	useData.GPSTIME = tempData2.GPSTIME;
	useData.sat_num = useData.sat1.size();
}

//函数功能：选择参考卫星
void SelectRef(Obs_epoch tempData, int &RefNo)
{
	bool iffind = false;
	//用于多历元连续定位判断
	if (RefNo > 0)
	{
		for (int i = 0; i < tempData.sat_num; i++)
		{
			if (tempData.sat1[i].numofsat == RefNo && tempData.sat1[i].E > 30 && tempData.sat1[i].judge_use == 0)
			{
				iffind = true;
				return;
			}
		}
	}
	//单历元定位模式会进下面if条件,选择最大高度角卫星作参考星
	if (iffind == false)
	{
		if (tempData.sat_num < 1)
		{
			printf("Error,SATNUM:%d\n", tempData.sat_num);
		}
		RefNo = -1;
		double max_ele = 10.0;
		for (int i = 0; i < tempData.sat_num; i++)
		{
			if (tempData.sat1[i].E >= max_ele && tempData.sat1[i].judge_use == 0)
			{
				RefNo = tempData.sat1[i].numofsat;
				max_ele = tempData.sat1[i].E;
			}
		}
	}

	//检核
	bool ifOk = false;
	for (int i = 0; i < tempData.sat_num; i++)
	{
		if (tempData.sat1[i].numofsat == RefNo && tempData.sat1[i].judge_use == 0)
		{
			ifOk = true;
			break;
		}
	}

	if (ifOk == false)
	{
		assert(0);
	}
}

//函数功能：依据误差传播定律获取噪声矩阵
void GetDDR_Peo(Obs_epoch GPSData, vector<int> useprn, int GPSRef, CMatrix &R)
{
	//R矩阵设定
	int numgps = GPSData.sat_num;
	if (numgps != useprn.size())
		assert(0);
	int numgps_x = numgps - 1;

	CMatrix TTR(numgps_x, numgps);
	CMatrix RT(numgps, numgps);
	int t = 0;
	int k = 0;
	for (k = 0; k < numgps; k++)
	{
		int match_O = -1;
		int PRN = useprn[k];
		for (int m = 0; m < GPSData.sat_num; m++)
		{
			if (GPSData.sat2[m].numofsat == PRN)
			{
				match_O = m;
				break;
			}
		}

		RT[k][k] = pow((0.3 + 0.3*exp(-GPSData.sat2[match_O].E / 10.0)), 2) + pow((0.3 + 0.3*exp(-GPSData.sat1[match_O].E / 10.0)), 2);

		if (PRN != GPSRef)
		{
			TTR[t][k] = 1.0;
			t++;
		}
		else
		{
			for (int m = 0; m < numgps_x; m++)
			{
				TTR[m][k] = -1.0;
			}
		}
	}

	R.first(numgps_x, numgps_x);
	//TTR矩阵的目的是根据误差传播定律，导出星间差分模式下的误差矩阵形式
	R = TTR * RT * TTR.T();
}

//函数功能：差分定位处理
bool SPP_DPOS_Pro(string filename1, string filename2)
{
	//默认第一个站为基准站，第二个站为流动站
	observe_spp sppfile1;
	observe_spp sppfile2;
	ReadSPP_KIN(filename1, sppfile1);
	ReadSPP_KIN(filename2, sppfile2);

	//投影矩阵：东北天坐标系（此处简化直接采用头文件坐标）
	CMatrix TT(3, 3);
	TT[0][0] = -sin(sppfile2.B)*cos(sppfile2.L);
	TT[0][1] = -sin(sppfile2.B)*sin(sppfile2.L);
	TT[0][2] = cos(sppfile2.B);

	TT[1][0] = -sin(sppfile2.L);
	TT[1][1] = cos(sppfile2.L);
	TT[1][2] = 0;

	TT[2][0] = cos(sppfile2.B)*cos(sppfile2.L);
	TT[2][1] = cos(sppfile2.B)*sin(sppfile2.L);
	TT[2][2] = sin(sppfile2.B);


	//输出差分定位结果用
	string name = filename2.substr(0, 4) + "_DPOS.txt";
	const char* output_filename = name.c_str();
	FILE *fp;
	fopen_s(&fp, output_filename, "w");

	//根据流动站匹配基准站对应历元数据
	int posk = 0;   //用以记录基站数据的匹配
	for (int j = 0; j < sppfile2.liyuan_num; j++)   //流动站逐历元解算
	{
		Obs_epoch tempData;

		//流动站当前时刻
		Obs_epoch tempData2 = sppfile2.epoch[j];

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

		//判断时间差：两个是否相隔超过20s
		if (fabs(sppfile1.epoch[i].GPSTIME - tempData2.GPSTIME) >= 20)
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

		//差分定位共视处理
		EpochCommonPro(tempData1, tempData2, tempData);

		//分离GPS和BDS卫星，便于后续单系统、双系统处理
		Obs_epoch GPSData;
		Obs_epoch BDSData;

		int k = 0;
		int m = 0;
		for (m = 0; m < tempData.sat_num; m++)
		{
			if (tempData.sat1[m].sattype == "G")
			{
				GPSData.sat1.push_back(tempData.sat1[m]);
				GPSData.sat2.push_back(tempData.sat2[m]);
			}

			if (tempData.sat1[m].sattype == "C")
			{
				BDSData.sat1.push_back(tempData.sat1[m]);
				BDSData.sat2.push_back(tempData.sat2[m]);
			}
		}

		GPSData.sat_num = GPSData.sat1.size();
		BDSData.sat_num = BDSData.sat1.size();
		tempData.GPSTIME = tempData2.GPSTIME;

		//暂时处理GPS和BDS都大于4颗卫星的数据，真正应用时一般两者相加大于5颗即可
		if (GPSData.sat_num < 4 || BDSData.sat_num < 4)
		{
			printf("GPS或BDS卫星数小于4颗\n");
			continue;
		}

		//采用站间星间二次差分模型，选择参考卫星
		int GPSRef = -1;
		int BDSRef = -2;
		SelectRef(GPSData, GPSRef);
		SelectRef(BDSData, BDSRef);

		//当前历元GPS和BDS卫星序列
		int PRN_X_G[32];
		int PRN_X_B[20];

		memset(PRN_X_G, 0, sizeof(PRN_X_G));
		memset(PRN_X_B, 0, sizeof(PRN_X_B));

		for (k = 0; k < GPSData.sat_num; k++)
			PRN_X_G[k] = GPSData.sat1[k].numofsat;
		for (k = 0; k < BDSData.sat_num; k++)
			PRN_X_B[k] = BDSData.sat1[k].numofsat;

		int num_X_G = GPSData.sat_num;
		int num_X_B = BDSData.sat_num;

		int numgps_x = num_X_G - 1;         //GPS双差方程个数
		int numbds_x = num_X_B - 1;         //BDS双差方程个数
		int num_x_p = numgps_x + numbds_x;  //双系统融合双差方程个数

		//定位解算相关的矩阵
		//1-双系统
		CMatrix B(num_x_p, 3);       //系数矩阵
		CMatrix L(num_x_p, 1);       //观测值矩阵
		CMatrix R(num_x_p, num_x_p); //误差矩阵（逆阵即为权矩阵）

		//2-单GPS
		CMatrix BG(numgps_x, 3);       //系数矩阵
		CMatrix LG(numgps_x, 1);       //观测值矩阵
		CMatrix RG(numgps_x, numgps_x); //误差矩阵（逆阵即为权矩阵）

		//3-单BDS
		CMatrix BC(numbds_x, 3);       //系数矩阵
		CMatrix LC(numbds_x, 1);       //观测值矩阵
		CMatrix RC(numbds_x, numbds_x); //误差矩阵（逆阵即为权矩阵）

		//获取GPS和BDS的误差矩阵

		vector<int> GPS_PRN;
		vector<int> BDS_PRN;

		for (k = 0; k < GPSData.sat_num; k++)
			GPS_PRN.push_back(PRN_X_G[k]);
		for (k = 0; k < BDSData.sat_num; k++)
			BDS_PRN.push_back(PRN_X_B[k]);

		GetDDR_Peo(GPSData, GPS_PRN, GPSRef, RG);
		GetDDR_Peo(BDSData, BDS_PRN, BDSRef, RC);

		//融合系统的R矩阵，分块对角矩阵
		for (int k1 = 0; k1 < numgps_x; k1++)
			for (int k2 = 0; k2 < numgps_x; k2++)
			{
				R[k1][k2] = RG[k1][k2];
			}

		for (int k1 = 0; k1 < numbds_x; k1++)
			for (int k2 = 0; k2 < numbds_x; k2++)
			{
				R[k1 + numgps_x][k2 + numgps_x] = RC[k1][k2];
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
		for (k = 0; k < GPSData.sat_num; k++)
		{
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

			//系数矩阵
			BG[tp][0] = Cof1;
			BG[tp][1] = Cof2;
			BG[tp][2] = Cof3;

			B[tw][0] = Cof1;
			B[tw][1] = Cof2;
			B[tw][2] = Cof3;

			//各类误差改正
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
			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + (test_geo2 - Rtest_geo2) - (test_geo1 - Rtest_geo1);
			double deltp1 = (tempSat2.data[0] - RefSat2.data[0]) - (tempSat1.data[0] - RefSat1.data[0]);

			L[tw][0] = deltp1 - DDgeo;
			LG[tp][0] = deltp1 - DDgeo;

			tp++;
			tw++;
		}

		//BDS
		tp = 0;  //记录单系统方程数
		for (k = 0; k < BDSData.sat_num; k++)
		{
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
			BC[tp][0] = Cof1;
			BC[tp][1] = Cof2;
			BC[tp][2] = Cof3;

			B[tw][0] = Cof1;
			B[tw][1] = Cof2;
			B[tw][2] = Cof3;

			//各类误差改正
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
			double DDgeo = (length_O2 - length_O1 - length_R2 + length_R1) + (test_geo2 - Rtest_geo2) - (test_geo1 - Rtest_geo1);
			double deltp1 = (tempSat2.data[0] - RefSat2.data[0]) - (tempSat1.data[0] - RefSat1.data[0]);

			L[tw][0] = deltp1 - DDgeo;
			LC[tp][0] = deltp1 - DDgeo;

			tp++;
			tw++;
		}

		//最小二乘解算
		//单GPS
		if (GPSData.sat_num >= 4 && BDSData.sat_num >= 4)
		{
			//权矩阵PG、法方程矩阵NG(逆矩阵即为参数方差)、参数求解(相对概率坐标的改正数)
			CMatrix PG = RG.InvertGaussJordan();
			CMatrix NG = BG.T()*PG*BG;
			CMatrix XG = NG.InvertGaussJordan() * BG.T()*PG*LG;

			CMatrix PC = RC.InvertGaussJordan();
			CMatrix NC = BC.T()*PC*BC;
			CMatrix XC = NC.InvertGaussJordan() * BC.T()*PC*LC;

			CMatrix P = R.InvertGaussJordan();
			CMatrix N = B.T()*P*B;
			CMatrix X = N.InvertGaussJordan() * B.T()*P*L;

			//恢复出估计的坐标值
			CMatrix XG_POS(3, 1); XG_POS[0][0] = POSITION_X2 - XG[0][0]; XG_POS[1][0] = POSITION_Y2 - XG[1][0]; XG_POS[2][0] = POSITION_Z2 - XG[2][0];
			CMatrix XC_POS(3, 1); XC_POS[0][0] = POSITION_X2 - XC[0][0]; XC_POS[1][0] = POSITION_Y2 - XC[1][0]; XC_POS[2][0] = POSITION_Z2 - XC[2][0];
			CMatrix X_POS(3, 1);  X_POS[0][0] = POSITION_X2 - X[0][0];   X_POS[1][0] = POSITION_Y2 - X[1][0];   X_POS[2][0] = POSITION_Z2 - X[2][0];

			//真实坐标写成矩阵
			CMatrix XT(3, 1);  XT[0][0] = P_X2; XT[1][0] = P_Y2; XT[2][0] = P_Z2;

			//与准确值作差得到XYZ坐标系下的误差
			CMatrix ErrorXG = XG_POS - XT;
			CMatrix ErrorXC = XC_POS - XT;
			CMatrix ErrorX = X_POS - XT;

			//将误差投影到NEU东北天坐标系
			CMatrix ErrorXG_NEU = TT * ErrorXG;
			CMatrix ErrorXC_NEU = TT * ErrorXC;
			CMatrix ErrorX_NEU = TT * ErrorX;

			fprintf(fp, "EpNum: %2d, %2d, %2d, %10.6f, %2d, %2d, %15.8g, %15.8g, %15.8g, %15.8g, %15.8g, %15.8g, %15.8g, %15.8g, %15.8g\n", j + 1, tempData2.hour, tempData2.minute, tempData2.second, GPSData.sat_num, BDSData.sat_num, 
				                                         ErrorXG_NEU[0][0], ErrorXG_NEU[1][0], ErrorXG_NEU[2][0], ErrorXC_NEU[0][0], ErrorXC_NEU[1][0], ErrorXC_NEU[2][0], ErrorX_NEU[0][0], ErrorX_NEU[1][0], ErrorX_NEU[2][0]);
		}
	}
	fclose(fp);
	return true;
}