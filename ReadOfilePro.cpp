#include "stdafx.h"
#include <fstream>     //文件流 
#include <string>      //字符串
#include <vector>      //容器
#include "Myheader.h"        //数据结构
#include "GPT2_1w_World.h"   //对流层改正用
#include "Matrix.h"          //矩阵类
using namespace std;

//声明外部GetGPSTime函数，在ReadPfile.cpp中定义
double GetGPSTime(int year, int month, int day, int hour, int minute, double second, int& dayofy);

//函数功能：得到GPS周
double getGPSweek(double& gpstime)
{
	int dayofy = 0;
	double t0 = GetGPSTime(1980, 1, 6, 0, 0, 0, dayofy);   //GPS起始日
	time_t sec = gpstime - t0;
	int week = (int)(sec / (86400 * 7));
	return week;
}

//函数功能：对卫星信息进行初始化
void inti_Sat(Sat &tempSat)
{
	tempSat.numofsat = 0;
	for (int i = 0; i < 20; i++)
	{
		tempSat.data[i] = 0;
	}

	tempSat.GPSTIME = 0;
	tempSat.ttlsec = 0;
	tempSat.TGD = 0;
	tempSat.a0 = 0;
	tempSat.a1 = 0;
	tempSat.a2 = 0;
	tempSat.tk = 0;
	tempSat.deltt = 0.0;
	tempSat.POS_X = 0.0;
	tempSat.POS_Y = 0.0;
	tempSat.POS_Z = 0.0;
	tempSat.r = 0;
	tempSat.A = 0;
	tempSat.E = 0;

	//误差改正
	tempSat.judge_use = 0;
	tempSat.Sat_clock = 0;
	tempSat.Trop_Delay = 0;
	tempSat.Trop_Map = 0;
	tempSat.Relat = 0;
	tempSat.Sagnac = 0;
	tempSat.Tide_Effect = 0;
	tempSat.Antenna_Height = 0;
	tempSat.Sat_Antenna = 0;
	tempSat.OffsetL1 = 0;
	tempSat.OffsetL2 = 0;
	tempSat.Windup = 0;
}

//函数功能：根据文件头中“SYS / # / OBS TYPES”信息，确定观测值顺序
void ChDateType_Auto(Sat& tempsat, vector<string> GType, vector<string> CType, vector<string> EType)
{
	double dt[8] = { 0 };
	//GPS
	if (tempsat.sattype == "G")
	{
		for (size_t k = 0; k < GType.size(); k++)
		{
			//频率L1_1575.42MHZ，编码C/A，伪距
			if (GType[k] == "C1C")
				dt[0] = tempsat.data[k];

			//频率L2_1227.60MHZ，编码Z，伪距
			if (GType[k] == "C2W")
				dt[1] = tempsat.data[k];

			//频率L2_1227.60MHZ，编码P,伪距
			if (GType[k] == "C2P")
				dt[6] = tempsat.data[k];

			//频率L5_1176.45MHZ，编码I+Q/Q/I,伪距
			if (GType[k] == "C5X" || GType[k] == "C5Q" || GType[k] == "C5I")
				dt[2] = tempsat.data[k];

			//频率L1_1575.42MHZ，编码C/A，载波相位
			if (GType[k] == "L1C")
				dt[3] = tempsat.data[k];

			//频率L2_1227.60MHZ，编码Z，载波相位
			if (GType[k] == "L2W")
				dt[4] = tempsat.data[k];

			//频率L2_1227.60MHZ，编码P,载波相位
			if (GType[k] == "L2P")
				dt[7] = tempsat.data[k];

			//频率L5_1176.45MHZ，编码I+Q/Q/I,载波相位
			if (GType[k] == "L5X" || GType[k] == "L5Q" || GType[k] == "L5I")
				dt[5] = tempsat.data[k];
		}

		tempsat.data[0] = dt[0];
		if (dt[1] > 100)
		{
			tempsat.data[1] = dt[1];
			tempsat.data[4] = dt[4];
		}
		else
		{
			tempsat.data[1] = dt[6];
			tempsat.data[4] = dt[7];
		}

		tempsat.data[2] = dt[2];
		tempsat.data[3] = dt[3];
		tempsat.data[5] = dt[5];
	}

	//BDS
	if (tempsat.sattype == "C")
	{
		for (size_t k = 0; k < CType.size(); k++)
		{
			if (CType[k] == "C1I" || CType[k] == "C2I")
				dt[0] = tempsat.data[k];
			if (CType[k] == "C7I")
				dt[1] = tempsat.data[k];
			if (CType[k] == "C6I")
				dt[2] = tempsat.data[k];

			if (CType[k] == "L1I" || CType[k] == "L2I")
				dt[3] = tempsat.data[k];
			if (CType[k] == "L7I")
				dt[4] = tempsat.data[k];
			if (CType[k] == "L6I")
				dt[5] = tempsat.data[k];
		}

		tempsat.data[0] = dt[0];
		tempsat.data[1] = dt[1];
		tempsat.data[2] = dt[2];
		tempsat.data[3] = dt[3];
		tempsat.data[4] = dt[4];
		tempsat.data[5] = dt[5];
	}

	//BDS-3
	if (tempsat.sattype == "B")
	{
		for (size_t k = 0; k < CType.size(); k++)
		{
			if (CType[k] == "C1I" || CType[k] == "C2I")
				dt[0] = tempsat.data[k];
			if (CType[k] == "C5X" || CType[k] == "C5P")
				dt[1] = tempsat.data[k];
			if (CType[k] == "C6I")
				dt[2] = tempsat.data[k];
			if (CType[k] == "C1X" || CType[k] == "C1P")
				dt[3] = tempsat.data[k];

			if (CType[k] == "L1I" || CType[k] == "L2I")
				dt[4] = tempsat.data[k];
			if (CType[k] == "L5X" || CType[k] == "L5P")
				dt[5] = tempsat.data[k];
			if (CType[k] == "L6I")
				dt[6] = tempsat.data[k];
			if (CType[k] == "L1X" || CType[k] == "L1P")
				dt[7] = tempsat.data[k];
		}

		tempsat.data[0] = dt[0];
		tempsat.data[1] = dt[1];
		tempsat.data[2] = dt[2];
		tempsat.data[3] = dt[3];
		tempsat.data[4] = dt[4];
		tempsat.data[5] = dt[5];
		tempsat.data[6] = dt[6];
		tempsat.data[7] = dt[7];
	}

	//Galileo
	if (tempsat.sattype == "E")
	{
		for (size_t k = 0; k < EType.size(); k++)
		{
			if (EType[k] == "C1X")
				dt[0] = tempsat.data[k];
			if (EType[k] == "C5X")
				dt[1] = tempsat.data[k];
			if (EType[k] == "C7X")
				dt[2] = tempsat.data[k];
			if (EType[k] == "C8X" || EType[k] == "C8Q")
				dt[3] = tempsat.data[k];

			if (EType[k] == "L1X")
				dt[4] = tempsat.data[k];
			if (EType[k] == "L5X")
				dt[5] = tempsat.data[k];
			if (EType[k] == "L7X")
				dt[6] = tempsat.data[k];
			if (EType[k] == "L8X" || EType[k] == "L8Q")
				dt[7] = tempsat.data[k];
		}

		tempsat.data[0] = dt[0];
		tempsat.data[1] = dt[1];
		tempsat.data[2] = dt[2];
		tempsat.data[3] = dt[3];
		tempsat.data[4] = dt[4];
		tempsat.data[5] = dt[5];
		tempsat.data[6] = dt[6];
		tempsat.data[7] = dt[7];
	}
}

//函数功能：根据卫星系统、PRN号、GPSTIME在对应的星历数据中进行匹配（时间最近原则）
bool Find_Eph(int system, int PRN, double GPSTIME, int &posk, vector<Ephemeris_CN_epoch> &EpochNG, vector<Ephemeris_CN_epoch> &EpochNC, vector<Ephemeris_CN_epoch> &EpochNB, vector < Ephemeris_CN_epoch>& EpochNE)
{
	posk = -1;

	double t1 = 0;
	double min = 1e4;

	if (system < 1 || system > 8)
	{
		printf("卫星系统参数错误！");
		_ASSERT(false);
	}

	switch (system)
	{
	case 1://GPS
	{
			  if (EpochNG.size() < 1)
			  {
				  printf("GPS星历为空！");
				  _ASSERT(false);
			  }

			  for (size_t k = 0; k < EpochNG.size(); k++)
			  {
				  if (PRN == EpochNG[k].PRN)
				  {
					  t1 = fabs(GPSTIME - EpochNG[k].GPSTIME);
					  min = t1;
					  posk = k;
					  break;
				  }
			  }

			  if (posk == -1)
			  {
				  return false;
				  break;
			  }

			  for (size_t k = posk; k < EpochNG.size(); k++)
			  {
				  if (PRN == EpochNG[k].PRN)
				  {
					  t1 = fabs(GPSTIME - EpochNG[k].GPSTIME);
					  if (t1 < min)
					  {
						  min = t1;
						  posk = k;
					  }
				  }
			  }

			  if (min > 7200)
			  {
				  return false;
				  break;
			  }

			  return true;
			  break;
	}

	case 2://BDS2
	{
			  if (EpochNC.size() < 1)
			  {
				  printf("BDS-2星历为空！");
				  _ASSERT(false);
			  }

			  for (size_t k = 0; k < EpochNC.size(); k++)
			  {
				  if (PRN == EpochNC[k].PRN)
				  {
					  t1 = fabs(GPSTIME - EpochNC[k].GPSTIME);
					  min = t1;
					  posk = k;
					  break;
				  }
			  }

			  if (posk == -1)
			  {
				  return false;
				  break;
			  }

			  for (size_t k = posk; k < EpochNC.size(); k++)
			  {
				  if (PRN == EpochNC[k].PRN)
				  {
					  t1 = fabs(GPSTIME - EpochNC[k].GPSTIME);
					  if (t1 < min)
					  {
						  min = t1;
						  posk = k;
					  }
				  }
			  }

			  if (min > 3600)
			  {
				  return false;
				  break;
			  }

			  return true;
			  break;
	}

	case 3://Galileo
	{
		if (EpochNE.size() < 1)
		{
			printf("Galileo星历为空！");
			_ASSERT(false);
		}

		for (size_t k = 0; k < EpochNE.size(); k++)
		{
			if (PRN == EpochNE[k].PRN)
			{
				t1 = fabs(GPSTIME - EpochNE[k].GPSTIME);
				min = t1;
				posk = k;
				break;
			}
		}

		if (posk == -1)
		{
			return false;
			break;
		}

		for (size_t k = posk; k < EpochNE.size(); k++)
		{
			if (PRN == EpochNE[k].PRN)
			{
				t1 = fabs(GPSTIME - EpochNE[k].GPSTIME);
				if (t1 < min)
				{
					min = t1;
					posk = k;
				}
			}
		}

		if (min > 600)
		{
			return false;
			break;
		}

		return true;
		break;
	}

	case 8://BDS3
	{
			  if (EpochNB.size() < 1)
			  {
				  printf("BDS-3星历为空！");
				  _ASSERT(false);
			  }

			  for (size_t k = 0; k < EpochNB.size(); k++)
			  {
				  if (PRN == EpochNB[k].PRN)
				  {
					  t1 = fabs(GPSTIME - EpochNB[k].GPSTIME);
					  min = t1;
					  posk = k;
					  break;
				  }
			  }

			  if (posk == -1)
			  {
				  return false;
				  break;
			  }

			  for (size_t k = posk; k < EpochNB.size(); k++)
			  {
				  if (PRN == EpochNB[k].PRN)
				  {
					  t1 = fabs(GPSTIME - EpochNB[k].GPSTIME);
					  if (t1 < min)
					  {
						  min = t1;
						  posk = k;
					  }
				  }
			  }

			  if (min > 3600)
			  {
				  return false;
				  break;
			  }

			  return true;
			  break;
	}
	return false;
	}
}

//函数功能：根据匹配到的星历进行卫星位置计算（GPS和BDS）,参考李征航教材《GPS测量与数据处理》3.5节GPS卫星位置的计算
void Calc_Eph_GCEJ(int system, int PRN, double GPSTIME, int posk, double rclk, Sat& sat, int CalTime, double XYZ[], vector<Ephemeris_CN_epoch> &EpochNG, vector<Ephemeris_CN_epoch> &EpochNC, vector<Ephemeris_CN_epoch> &EpochNB, vector< Ephemeris_CN_epoch > & EpochNE)
{
	if (sat.judge_use == 1)
		return;

	if (posk < 0)
	{
		sat.judge_use = 1;
		return;
	}

	Ephemeris_CN_epoch tempEpoch;

	switch (system)
	{
	case 1:
	{
			  tempEpoch = EpochNG[posk];
			  break;
	}

	case 2:
	{
			  tempEpoch = EpochNC[posk];
			  GPSTIME = GPSTIME - 14.000000;
			  break;
	}

	case 3:
	{
		      tempEpoch = EpochNE[posk];
		      break;
	}

	case 8:
	{
			  tempEpoch = EpochNB[posk];
			  GPSTIME = GPSTIME - 14.000000;
			  break;
	}
		break;
	}

	double c = 299792458.458;

	double we = 7.2921151467e-5;
	double GM = 3.9860050e14;

	if (system == 2 || system == 8)
	{
		we = 7.292115e-5;
		GM = 3.986004418e14;
	}

	if (system == 3)
	{
		we = 7.2921151467e-5;
		GM = 3.986004418e14;
	}

	double n0, n, tk, Mk, ek1, ek2, Ek, Vk, fk;
	double deltau, deltar, deltai, uk, rk, ik, xk, yk, Lk;

	//卫星运动的平均角速度n
	n0 = sqrt(GM) / pow(tempEpoch.sqrtA, 3);
	n = n0 + tempEpoch.delta_n;

	tk = GPSTIME - tempEpoch.GPSTIME;

	if (tk > 302400) tk -= 604800.0;
	else if (tk < -302400) tk += 604800.0;

	sat.Sat_clock = tempEpoch.a0 + tempEpoch.a1 * tk + tempEpoch.a2 * pow(tk, 2);
	sat.TGD = tempEpoch.TGD * C_light;

	tk = tk - sat.Sat_clock - rclk;

	if (CalTime == 2 && system == 2 && PRN == 3)
	{
		double xx = tk;
		double ax = 0.2;
	}

	//计算观测瞬间卫星的平近点角Mk
	Mk = tempEpoch.M0 + n*tk;

	//计算偏近点角Ek
	ek1 = Mk;
	do
	{
		ek2 = Mk + tempEpoch.e*sin(ek1);
		if (fabs(ek2 - ek1) <= 1.0e-13)  break;
		ek1 = ek2;
	} while (1);
	Ek = ek1;

	//计算真近点角Vk
	Vk = 2 * atan(sqrt((1.0 + tempEpoch.e) / (1.0 - tempEpoch.e))*tan(Ek / 2));

	//计算升交角距fk
	fk = Vk + tempEpoch.w;

	//计算摄动改正项deltau、deltar、deltai
	deltau = tempEpoch.Cuc * cos(2.0*fk) + tempEpoch.Cus * sin(2.0*fk);
	deltar = tempEpoch.Crc * cos(2.0*fk) + tempEpoch.Crs * sin(2.0*fk);
	deltai = tempEpoch.Cic * cos(2.0*fk) + tempEpoch.Cis * sin(2.0*fk);

	//摄动改正
	uk = fk + deltau;
	rk = tempEpoch.sqrtA * tempEpoch.sqrtA * (1.0 - tempEpoch.e*cos(Ek)) + deltar;
	ik = tempEpoch.i0 + deltai + tempEpoch.i_DOT * tk;

	//计算卫星在轨道面坐标系中的位置
	xk = rk * cos(uk);
	yk = rk * sin(uk);

	if (system == 2 && PRN <= 5)
	{
		//计算观测瞬间升交点的经度
		Lk = tempEpoch.OMEGA + tempEpoch.OMEGA_DOT * tk - we * tempEpoch.TOE;

		//计算卫星在瞬时地球坐标系中的位置
		double xx, yy, zz;
		xx = xk * cos(Lk) - yk * cos(ik) * sin(Lk);
		yy = xk * sin(Lk) + yk * cos(ik) * cos(Lk);
		zz = yk * sin(ik);

		double we_tk = we * tk;
		double an = -5 * PI / 180;

		sat.POS_X = cos(we_tk) * xx + sin(we_tk) * cos(an) * yy + sin(we_tk) * sin(an) * zz;
		sat.POS_Y = -sin(we_tk) * xx + cos(we_tk) * cos(an) * yy + cos(we_tk) * sin(an) * zz;
		sat.POS_Z = -sin(an) * yy + cos(an) * zz;
	}
	else
	{
		//计算观测瞬间升交点的经度
		Lk = tempEpoch.OMEGA + (tempEpoch.OMEGA_DOT - we)*tk - we*tempEpoch.TOE;

		//计算卫星在瞬时地球坐标系中的位置
		sat.POS_X = xk * cos(Lk) - yk * cos(ik) * sin(Lk);
		sat.POS_Y = xk * sin(Lk) + yk * cos(ik) * cos(Lk);
		sat.POS_Z = yk * sin(ik);
	}

	sat.xdl_t = -2 * sqrt(GM) * tempEpoch.e * tempEpoch.sqrtA * sin(Ek) / C_light;

	//计算地球自转效应:粗略解算
	sat.Sagnac = we * (sat.POS_X * XYZ[1] - sat.POS_Y * XYZ[0]) / C_light;
}

//函数功能：初步单点定位功能（不考虑高度角定权）
bool Get_SPP_POS(Obs_epoch &tepoch, CMatrix &XYZ)
{
	int gps_num = 0;
	int bds_num = 0;
	int gal_num = 0;

	bool GPS_Y = false;
	bool BDS_Y = false;
	bool GAL_Y = false;

	//存储BL矩阵信息，不超过55颗卫星
	CMatrix BL(55, 5);

	int tw = 0;

	for (int i = 0; i < tepoch.sat_num; i++)
	{
		Sat tempSat = tepoch.sat[i];

		//条件1：成功匹配卫星
		if (tempSat.posk < 0)
			continue;

		//条件2：GPS或BDS卫星
		if (tempSat.sattype != "G" && tempSat.sattype != "C" && tempSat.sattype != "E")
			continue;

		//条件3：伪距观测值不为0，且两个频率上的伪距之差不大于100
		double DiffPsedu = fabs(tempSat.data[0] - tempSat.data[1]);
		if (tempSat.data[0] < 100 || tempSat.data[1] < 100 || DiffPsedu > 100 || tempSat.health == 1)
			continue;

		//满足以上三个条件认为卫星伪距信息可用，开始进行单点定位

		if (tempSat.sattype == "C")
		{
			double IF_P = (f1C * f1C * tempSat.data[0] - f2C * f2C * tempSat.data[1]) / (f1C * f1C - f2C * f2C);

			if (IF_P >= 50000000)
			{
				tepoch.sat[i].judge_use = 1;
				continue;
			}

			if (tempSat.numofsat <= 10 && IF_P <= 30000000)
			{
				tepoch.sat[i].judge_use = 1;
				continue;
			}
		}

		if (tempSat.sattype == "G")
		{
			GPS_Y = true;
			gps_num++;
			BL[tw][4] = 1.0;
		}
		if (tempSat.sattype == "C")
		{
			BDS_Y = true;
			bds_num++;
			BL[tw][4] = 2.0;
		}
		if (tempSat.sattype == "E")
		{
			GAL_Y = true;
			gal_num++;
			BL[tw][4] = 3.0;
		}

		double Length = sqrt(pow(tempSat.POS_X - XYZ[0][0], 2) + pow(tempSat.POS_Y - XYZ[1][0], 2) + pow(tempSat.POS_Z - XYZ[2][0], 2));

		double Cof1 = (tempSat.POS_X - XYZ[0][0]) / Length;
		double Cof2 = (tempSat.POS_Y - XYZ[1][0]) / Length;
		double Cof3 = (tempSat.POS_Z - XYZ[2][0]) / Length;

		BL[tw][0] = -Cof1;
		BL[tw][1] = -Cof2;
		BL[tw][2] = -Cof3;

		//参考站1-非参考卫星
		double error_satclock1 = tempSat.Sat_clock * C_light;
		double error_relat1 = tempSat.xdl_t;
		double we = 7.2921151467e-5;

		if (tempSat.sattype == "C" || tempSat.sattype == "B")
		{
			we = 7.292115e-5;
		}
		tempSat.Sagnac = we * (tempSat.POS_X * XYZ[1][0] - tempSat.POS_Y * XYZ[0][0]) / C_light;

		//Sagnac效应也需要改正
		double error_sagnac1 = tempSat.Sagnac;
		double error_tgd = tempSat.TGD;
		double error_antenna_height1 = tempSat.Antenna_Height;
		double error_sat_antenna1 = tempSat.Sat_Antenna;
		double test_geo1 = -error_satclock1 - error_relat1 + error_sagnac1 + error_tgd - error_sat_antenna1;

		double test_P11;
		test_P11 = -(Length - tempSat.data[0] + test_geo1);

		BL[tw][3] = test_P11;

		tw++;
	}

	//大于5颗表示可以定位
	if (tw > 5)
	{
		//定义B矩阵
		int col_num = 4;

		if (GPS_Y == false && BDS_Y == false && GAL_Y == false)
			return false;

		if (GPS_Y == true && BDS_Y == true && GAL_Y == true)
			col_num = 6;
		else if ((GPS_Y == true && BDS_Y == true && GAL_Y == false) || (GPS_Y == true && BDS_Y == false && GAL_Y == true) || (GPS_Y == false && BDS_Y == true && GAL_Y == true))
			col_num = 5;
		else if ((GPS_Y == true && BDS_Y == false && GAL_Y == false) || (GPS_Y == false && BDS_Y == true && GAL_Y == false) || (GPS_Y == false && BDS_Y == false && GAL_Y == true))
			col_num = 4;

		CMatrix B(tw, col_num);
		CMatrix L(tw, 1);
		CMatrix X(col_num, 1);

		for (int i = 0; i < tw; i++)
		{
			B[i][0] = BL[i][0];
			B[i][1] = BL[i][1];
			B[i][2] = BL[i][2];

			L[i][0] = BL[i][3];

			//一般情况下三种系统都有，不需要考虑单系统的情况
			if (fabs(BL[i][4] - 1.0) < 0.1 && GPS_Y == true)
				B[i][3] = 1.0;
			if (fabs(BL[i][4] - 2.0) < 0.1 && BDS_Y == true)
				B[i][4] = 1.0;
			if (fabs(BL[i][4] - 3.0) < 0.1 && GAL_Y == true)
				B[i][5] = 1.0;
		}

		X = (B.T()*B).InvertGaussJordan() * B.T() * L;

		XYZ[0][0] += X[0][0];
		XYZ[1][0] += X[1][0];
		XYZ[2][0] += X[2][0];

		return true;
	}
	else
		return false;
}

//函数功能：XYZ转为BLH（空间直角坐标转为大地坐标）
void OnXYZtoBLH(double X, double Y, double Z, double BLH[])
{
	double a = 6378137.000, e2 = 0.00669437999014132;
	double R = sqrt(X*X + Y*Y);
	double B0 = atan2(Z, R);
	double N;
	double B, L, H;

	L = atan2(Y, X);

	while (1)
	{
		N = a / sqrt(1 - e2 * sin(B0) * sin(B0));
		B = atan2(Z + N * e2 * sin(B0), R);
		H = R / cos(B) - N;

		if (fabs(B - B0) < 1e-12)
			break;
		B0 = B;
	}

	BLH[0] = B;
	BLH[1] = L;
	BLH[2] = H;
}

//函数功能：获取卫星的高度角和方位角
void CalSatEA(CMatrix HH, double StaXYZ[], Sat& sat)
{
	if (sat.posk < 0)
		return;

	if (HH.Row != 3 || HH.Col != 3)
	{
		_ASSERT(false);
	}

	//计算高度角和方位角
	CMatrix DeltXYZ(3, 1);
	CMatrix Sta2Sat(3, 1);

	Sta2Sat[0][0] = sat.POS_X - StaXYZ[0];
	Sta2Sat[1][0] = sat.POS_Y - StaXYZ[1];
	Sta2Sat[2][0] = sat.POS_Z - StaXYZ[2];

	DeltXYZ = HH * Sta2Sat;
	double r = sqrt(pow(DeltXYZ[0][0], 2) + pow(DeltXYZ[1][0], 2) + pow(DeltXYZ[2][0], 2));
	sat.A = atan2(DeltXYZ[1][0], DeltXYZ[0][0]);
	sat.A = sat.A * 180 / PI;

	if (sat.A < 0)
	{
		sat.A = sat.A + 360;
	}

	sat.E = atan(DeltXYZ[2][0] / sqrt(pow(DeltXYZ[0][0], 2) + pow(DeltXYZ[1][0], 2)));
	sat.E = sat.E * 180 / PI;
}

//函数功能：单点定位功能（考虑高度角定权）
bool Get_SPP_POS_Ele(Obs_epoch &tepoch, CMatrix &XYZ)
{
	int gps_num = 0;
	int bds_num = 0;
	int gal_num = 0;

	bool GPS_Y = false;
	bool BDS_Y = false;
	bool GAL_Y = false;

	//存储BL矩阵信息，不超过55颗卫星
	CMatrix BL(55, 6);

	int tw = 0;

	for (int i = 0; i < tepoch.sat_num; i++)
	{
		Sat tempSat = tepoch.sat[i];

		if (tempSat.E < 15.0)
			continue;

		//条件1：成功匹配卫星
		if (tempSat.posk < 0)
			continue;

		//条件2：GPS、BDS或Galileo卫星
		if (tempSat.sattype != "G" && tempSat.sattype != "C" && tempSat.sattype != "E")
			continue;

		//条件3：伪距观测值不为0，且两个频率上的伪距之差不大于100
		double DiffPsedu = fabs(tempSat.data[0] - tempSat.data[1]);
		if (tempSat.data[0] < 100 || tempSat.data[1] < 100 || DiffPsedu > 100 || tempSat.health == 1)
			continue;

		//满足以上三个条件认为卫星伪距信息可用，开始进行单点定位

		if (tempSat.sattype == "C")
		{
			double IF_P = (f1C * f1C * tempSat.data[0] - f2C * f2C * tempSat.data[1]) / (f1C * f1C - f2C * f2C);

			if (IF_P >= 50000000)
			{
				tepoch.sat[i].judge_use = 1;
				continue;
			}

			if (tempSat.numofsat <= 10 && IF_P <= 30000000)
			{
				tepoch.sat[i].judge_use = 1;
				continue;
			}
		}

		if (tempSat.sattype == "G")
		{
			GPS_Y = true;
			gps_num++;
			BL[tw][4] = 1.0;
		}
		if (tempSat.sattype == "C")
		{
			BDS_Y = true;
			bds_num++;
			BL[tw][4] = 2.0;
		}
		if (tempSat.sattype == "E")
		{
			GAL_Y = true;
			gal_num++;
			BL[tw][4] = 3.0;
		}

		BL[tw][5] = tempSat.E;

		double Length = sqrt(pow(tempSat.POS_X - XYZ[0][0], 2) + pow(tempSat.POS_Y - XYZ[1][0], 2) + pow(tempSat.POS_Z - XYZ[2][0], 2));

		double Cof1 = (tempSat.POS_X - XYZ[0][0]) / Length;
		double Cof2 = (tempSat.POS_Y - XYZ[1][0]) / Length;
		double Cof3 = (tempSat.POS_Z - XYZ[2][0]) / Length;

		BL[tw][0] = -Cof1;
		BL[tw][1] = -Cof2;
		BL[tw][2] = -Cof3;

		//参考站1-非参考卫星
		double error_satclock1 = tempSat.Sat_clock * C_light;
		double error_relat1 = tempSat.xdl_t;
		double error_trop1 = tempSat.Trop_Delay;

		double we = 7.2921151467e-5;

		if (tempSat.sattype == "C" || tempSat.sattype == "B")
		{
			we = 7.292115e-5;
		}
		tempSat.Sagnac = we * (tempSat.POS_X * XYZ[1][0] - tempSat.POS_Y * XYZ[0][0]) / C_light;

		//Sagnac效应也需要改正
		double error_sagnac1 = tempSat.Sagnac;
		double error_tgd = tempSat.TGD;
		double error_antenna_height1 = tempSat.Antenna_Height;
		double error_sat_antenna1 = tempSat.Sat_Antenna;
		double test_geo1 = -error_satclock1 + error_trop1 - error_relat1 + error_sagnac1 + error_tgd - error_sat_antenna1;

		double test_P11;
		test_P11 = -(Length - tempSat.data[0] + test_geo1);

		BL[tw][3] = test_P11;

		tw++;
	}

	//大于5颗表示可以定位
	if (tw > 5)
	{
		//定义B矩阵
		int col_num = 4;

		if (GPS_Y == false && BDS_Y == false && GAL_Y == false)
			return false;

		if (GPS_Y == true && BDS_Y == true && GAL_Y == true)
			col_num = 6;
		else if ((GPS_Y == true && BDS_Y == true && GAL_Y == false) || (GPS_Y == true && BDS_Y == false && GAL_Y == true) || (GPS_Y == false && BDS_Y == true && GAL_Y == true))
			col_num = 5;
		else if ((GPS_Y == true && BDS_Y == false && GAL_Y == false) || (GPS_Y == false && BDS_Y == true && GAL_Y == false) || (GPS_Y == false && BDS_Y == false && GAL_Y == true))
			col_num = 4;

		CMatrix B(tw, col_num);
		CMatrix L(tw, 1);
		CMatrix R(tw, tw);
		CMatrix X(col_num, 1);

		for (int i = 0; i < tw; i++)
		{
			B[i][0] = BL[i][0];
			B[i][1] = BL[i][1];
			B[i][2] = BL[i][2];

			L[i][0] = BL[i][3];

			double elp = sin(BL[i][5] / 180.0*PI) * sin(BL[i][5] / 180.0*PI);
			R[i][i] = 0.09 + 0.09 / elp;

			//一般情况下三种系统都有，不需要考虑单系统的情况
			if (fabs(BL[i][4] - 1.0) < 0.1 && GPS_Y == true)
				B[i][3] = 1.0;
			if (fabs(BL[i][4] - 2.0) < 0.1 && BDS_Y == true)
				B[i][4] = 1.0;
			if (fabs(BL[i][4] - 3.0) < 0.1 && GAL_Y == true)
				B[i][5] = 1.0;
		}

		CMatrix P = R.InvertGaussJordan();
		X = (B.T()*P*B).InvertGaussJordan() * B.T() * P * L;

		XYZ[0][0] += X[0][0];
		XYZ[1][0] += X[1][0];
		XYZ[2][0] += X[2][0];

		return true;
	}
	else
		return false;
}

//函数功能：进行接收机天线改正
void AntOffsetCorrect(Obs_epoch& tempData, double PCO1[3], double PCO2[3])
{
	double c = 299792458.458;
	double off_L1, off_L2, lambda_l1, lambda_l2;
	off_L1 = off_L2 = lambda_l1 = lambda_l2 = 0;

	double f_bds1 = 1561.098*1E+6;
	double f_bds2 = 1207.14*1E+6;

	double f_gps1 = 154 * 10.23E6;
	double f_gps2 = 120 * 10.23E6;

	double f_geo1 = 1575.42E6;
	double f_geo2 = 1176.45E6;

	double f_glo1 = 1602E6;
	double f_glo2 = 1246E6;

	for (size_t i = 0; i < tempData.sat.size(); i++)
	{
		if (tempData.sat[i].judge_use == 2)
			continue;

		if (tempData.sat[i].sattype == "G")
		{
			lambda_l1 = c / f_gps1;
			lambda_l2 = c / f_gps2;
		}

		if (tempData.sat[i].sattype == "C")
		{
			lambda_l1 = c / f_bds1;
			lambda_l2 = c / f_bds2;
		}

		if (tempData.sat[i].sattype == "B")
		{
			lambda_l1 = c / f1B;
			lambda_l2 = c / f2B;
		}

		if (tempData.sat[i].sattype == "E")
		{
			lambda_l1 = c / f_geo1;
			lambda_l2 = c / f_geo2;
		}

		double elev = tempData.sat[i].E / 180.0 * PI;
		double A = tempData.sat[i].A / 180.0 * PI;

		//表明未能获得概略坐标
		if (elev < 0.5 && A < 0.5)
		{
			tempData.sat[i].OffsetL1 = 0;
			tempData.sat[i].OffsetL2 = 0;
			continue;
		}

		double cosel = cos(elev);
		double sinel = sin(elev);
		double cosaz = cos(A);
		double sinaz = sin(A);

		off_L1 = (PCO1[0] * cosel * cosaz + PCO1[1] * cosel * sinaz + PCO1[2] * sinel) / 1000.0;
		off_L2 = (PCO2[0] * cosel * cosaz + PCO2[1] * cosel * sinaz + PCO2[2] * sinel) / 1000.0;

		tempData.sat[i].OffsetL1 = off_L1 / lambda_l1;
		tempData.sat[i].OffsetL2 = off_L2 / lambda_l2;
	}
}

//函数功能：O文件读取、数据处理（单点定位解算）及误差文件输出
bool ReadOfileProcess(string Filename, vector <Ephemeris_CN_epoch> &EpochNG, vector <Ephemeris_CN_epoch> &EpochNC, vector <Ephemeris_CN_epoch> &EpochNB, vector <Ephemeris_CN_epoch>& EpochNE)
{
	fstream infile;
	infile.open(Filename, ios::in);
	if (!infile)
	{
		printf("文件打开失败");
		return false;
	}

	printf("开始观测数据读取！\n");

	string name1 = Filename.substr(0, 4) + "_SPP.txt"; //左边4位取出测站名称+"SPP.txt"构成文件名
	const char* filename1 = name1.c_str();

	string name2 = Filename.substr(0, 4) + "_error.txt"; //左边4位取出测站名称+"_error.txt"构成文件名
	const char* filename2 = name2.c_str();

	//输出单点定位结果用
	FILE *fp;
	fopen_s(&fp, filename1, "w");

	//输出中间过程（卫星位置、各类误差改正）
	FILE *fp1;
	fopen_s(&fp1, filename2, "w");

	string strs;
	string str2;

	//对流层改正初始化：细节可不看
	CGPT2_1w_World TROP;
	TROP.init_trop();

	//存储GPS和BDS观测数据类型：BDS-2与BDS-3一致
	vector <string> GType;
	vector <string> CType;
	vector <string> EType;
	double XYZ[3] = { 0 };

	//开始读取o文件中内容
	Ofileheader temp_o;

	//初始化测量的开始时刻和结束时刻
	for (int i = 0; i < 6; i++)
	{
		temp_o.FIRST_TIME[i] = 0;
		temp_o.LAST_TIME[i] = 0;
	}

	do
	{
		getline(infile, strs);
		if (strs.length() < 60) {
			continue;
		}

		if (strs.substr(60, 11) == "MARKER NAME")
		{
			str2 = strs.substr(0, 4);
			temp_o.marker_name = str2;
			continue;
		}

		//文件头概率坐标
		if (strs.substr(60, 19) == "APPROX POSITION XYZ")
		{
			str2 = strs.substr(0, 14);
			temp_o.APP_X = atof(str2.c_str());
			str2 = strs.substr(14, 14);
			temp_o.APP_Y = atof(str2.c_str());
			str2 = strs.substr(28, 14);
			temp_o.APP_Z = atof(str2.c_str());

			XYZ[0] = temp_o.APP_X;
			XYZ[1] = temp_o.APP_Y;
			XYZ[2] = temp_o.APP_Z;

			continue;
		}

		//天线类型
		if (strs.substr(60, 12) == "ANT # / TYPE")
		{
			temp_o.ANT_TYPE = strs.substr(20, 20);
			fstream infileAnt;
			infileAnt.open("antenna.txt", ios::in);
			if (!infileAnt)
			{
				printf("天线文件打开失败");
				_ASSERT(false);
			}

			temp_o.PCO1[0] = 0;
			temp_o.PCO1[1] = 0;
			temp_o.PCO1[2] = 0;

			temp_o.PCO2[0] = 0;
			temp_o.PCO2[1] = 0;
			temp_o.PCO2[2] = 0;

			string stra;
			while (getline(infileAnt, stra))
			{
				if (temp_o.ANT_TYPE == stra.substr(0, 20))
				{
					temp_o.PCO1[0] = atof(stra.substr(20,10).c_str());
					temp_o.PCO1[1] = atof(stra.substr(30,10).c_str());
					temp_o.PCO1[2] = atof(stra.substr(40,10).c_str());

					temp_o.PCO2[0] = atof(stra.substr(50,10).c_str());
					temp_o.PCO2[1] = atof(stra.substr(60,10).c_str());
					temp_o.PCO2[2] = atof(stra.substr(70,10).c_str());
					break;
				}
			}//如果找到相同的天线
			infileAnt.close();
			printf("天线文件读取完成！\n");
			continue;
		}

		//接收机类型（固定）
		if (strs.substr(60, 19) == "REC # / TYPE / VERS")
		{
			temp_o.REC_TYPE = "ComNav";
			continue;
		}

		//天线高
		if (strs.substr(60, 20) == "ANTENNA: DELTA H/E/N")
		{
			str2 = strs.substr(0, 14);
			temp_o.H = atof(str2.c_str());
			str2 = strs.substr(22, 6);
			temp_o.E = atof(str2.c_str());
			str2 = strs.substr(36, 6);
			temp_o.N = atof(str2.c_str());
			continue;
		}

		//观测值类型顺序
		if (strs.substr(60, 19) == "SYS / # / OBS TYPES")
		{
			//GPS
			if (strs.substr(0, 1) == "G")
			{
				str2 = strs.substr(4, 2);
				int type = atoi(str2.c_str());

				if (type <= 13)
				{
					for (int k = 0; k < type; k++)
					{
						str2 = strs.substr(7 + k * 4, 3);
						GType.push_back(str2);
					}
				}
				else
				{
					for (int k = 0; k < 13; k++)
					{
						str2 = strs.substr(7 + k * 4, 3);
						GType.push_back(str2);
					}

					getline(infile, strs);

					for (int k = 0; k < (type - 13); k++)
					{
						str2 = strs.substr(7 + k * 4, 3);
						GType.push_back(str2);
					}
				}
			}

			//BDS
			if (strs.substr(0, 1) == "C")
			{
				str2 = strs.substr(4, 2);
				int type = atoi(str2.c_str());

				if (type <= 13)
				{
					for (int k = 0; k < type; k++)
					{
						str2 = strs.substr(7 + k * 4, 3);
						CType.push_back(str2);
					}
				}
				else
				{
					for (int k = 0; k < 13; k++)
					{
						str2 = strs.substr(7 + k * 4, 3);
						CType.push_back(str2);
					}

					getline(infile, strs);

					for (int k = 0; k < (type - 13); k++)
					{
						str2 = strs.substr(7 + k * 4, 3);
						CType.push_back(str2);
					}
				}
			}

			//Galileo
			if (strs.substr(0, 1) == "E")
			{
				str2 = strs.substr(4, 2);
				int type = atoi(str2.c_str());

				if (type <= 13)
				{
					for (int k = 0; k < type; k++)
					{
						str2 = strs.substr(7 + k * 4, 3);
						EType.push_back(str2);
					}
				}
				else
				{
					for (int k = 0; k < 13; k++)
					{
						str2 = strs.substr(7 + k * 4, 3);
						EType.push_back(str2);
					}

					getline(infile, strs);

					for (int k = 0; k < (type - 13); k++)
					{
						str2 = strs.substr(7 + k * 4, 3);
						EType.push_back(str2);
					}
				}
			}

			continue;
		}

		//历元间隔，单位为秒
		if (strs.substr(60, 8) == "INTERVAL")
		{
			str2 = strs.substr(0, 10);
			temp_o.INTERVAL = atof(str2.c_str());
			continue;
		}

		//观测开始时刻
		if (strs.substr(60, 17) == "TIME OF FIRST OBS")
		{
			for (int i = 0; i < 5; i++)
			{
				str2 = strs.substr((0 + 6 * i), 6);
				temp_o.FIRST_TIME[i] = atoi(str2.c_str());
			}
			str2 = strs.substr(30, 13);
			temp_o.FIRST_TIME[5] = atoi(str2.c_str());
			continue;
		}

		//观测结束时刻
		if (strs.substr(60, 16) == "TIME OF LAST OBS")
		{
			for (int i = 0; i<5; i++)
			{
				str2 = strs.substr((0 + 6 * i), 6);
				temp_o.LAST_TIME[i] = atoi(str2.c_str());
			}
			str2 = strs.substr(30, 13);
			temp_o.LAST_TIME[5] = atoi(str2.c_str());
			continue;
		}

		//跳秒
		if (strs.substr(60, 12) == "LEAP SECONDS")
		{
			str2 = strs.substr(0, 6);
			temp_o.LEAPSEC = atoi(str2.c_str());
			continue;
		}
	} while (strs.substr(60, 3) != "END");   //文件头读取结束

	printf("文件头读取完成\n");

	//读取观测文件部分
	double X_POS = 0;
	double Y_POS = 0;
	double Z_POS = 0;
	int i = 0;

	while (getline(infile, strs))
	{
		Obs_epoch temp_oe;
		temp_oe.liyuan_num = i;

		//读取年月日时分秒
		str2 = strs.substr(2, 4);
		temp_oe.year = atoi(str2.c_str());
		str2 = strs.substr(6, 3);
		temp_oe.month = atoi(str2.c_str());
		str2 = strs.substr(9, 3);
		temp_oe.day = atoi(str2.c_str());
		str2 = strs.substr(12, 3);
		temp_oe.hour = atoi(str2.c_str());
		str2 = strs.substr(15, 3);
		temp_oe.minute = atoi(str2.c_str());
		str2 = strs.substr(18, 11);
		temp_oe.second = atof(str2.c_str());

		temp_oe.GPSTIME = GetGPSTime(temp_oe.year, temp_oe.month, temp_oe.day, temp_oe.hour, temp_oe.minute, temp_oe.second, temp_oe.dayofy);
		//对流层改正模型用
		double JD = TROP.Julday(temp_oe.year, temp_oe.month, temp_oe.day, temp_oe.hour, temp_oe.minute, temp_oe.second);

		str2 = strs.substr(29, 3);
		temp_oe.flag = atoi(str2.c_str());
		str2 = strs.substr(32, 3);
		temp_oe.sat_num = atoi(str2.c_str());

		for (int k = 0; k<temp_oe.sat_num; k++)
		{
			Sat tempsat;
			inti_Sat(tempsat);

			getline(infile, strs);

			str2 = strs.substr(0, 1);                //类型
			tempsat.sattype = str2.c_str();

			if (tempsat.sattype == " ")
				tempsat.sattype = "G";

			//只处理GPS、BDS和Galileo数据
			if (tempsat.sattype != "G" && tempsat.sattype != "C" && tempsat.sattype != "E")
				continue;

			str2 = strs.substr(1, 2);                //卫星号 
			tempsat.numofsat = atoi(str2.c_str());

			//区分BDS-3, 使用B表示
			if (tempsat.sattype == "C" && tempsat.numofsat > 16)
			{
				tempsat.sattype = "B";
			}

			//储存观测数据
			int datatype = (strs.length() - 1) / 16;
			for (int t = 0; t < datatype; t++)
			{
				str2 = strs.substr(3 + t * 16, 14);
				tempsat.data[t] = atof(str2.c_str());
			}

			//根据文件头中声明的观测顺序调整数据存储
			ChDateType_Auto(tempsat, GType, CType, EType);

			tempsat.system = -1;
			if (tempsat.sattype == "G")
				tempsat.system = 1;

			if (tempsat.sattype == "C")
				tempsat.system = 2;

			if (tempsat.sattype == "B")
				tempsat.system = 8;

			if (tempsat.sattype == "E")
				tempsat.system = 3;

			/*
			if (tempsat.sattype == "J")
			tempsat.system = 4;

			if (tempsat.sattype == "R")
			tempsat.system = 5;

			if (tempsat.sattype == "S")
			tempsat.system = 6;

			if (tempsat.sattype == "I")
			{
			tempsat.system = 7;
			continue;
			}
			*/

			if (tempsat.system < 0)
			{
				_ASSERT(false);
			}

			if (tempsat.sattype == "G" || tempsat.sattype == "C" || tempsat.sattype == "B" || tempsat.sattype == "E")
			{
				bool findYorN = Find_Eph(tempsat.system, tempsat.numofsat, temp_oe.GPSTIME, tempsat.posk, EpochNG, EpochNC, EpochNB, EpochNE);

				if (findYorN == false)
				{
					//未匹配到星历，标为不可用
					tempsat.judge_use = 1;
					continue;
				}

				double SpreadTime = tempsat.data[0] / C_light;

				if (tempsat.system <= 4 || tempsat.system == 8)
					Calc_Eph_GCEJ(tempsat.system, tempsat.numofsat, temp_oe.GPSTIME, tempsat.posk, tempsat.data[0] / C_light, tempsat, 1, XYZ, EpochNG, EpochNC, EpochNB, EpochNE);
				else
					_ASSERT(false);

				temp_oe.sat.push_back(tempsat);
			}
		}

		temp_oe.sat_num = temp_oe.sat.size();

		bool POS_IF = false;

		//初始历元进行伪距单点定位获得概略坐标，更新XYZ:
		double dCoordX1 = 0;
		double dCoordY1 = 0;
		double dCoordZ1 = 0;

		CMatrix Pos(3, 1);

		//从第二个历元开始，不用第一次迭代的结果
		if (i > 0)
		{
			Pos[0][0] = X_POS;
			Pos[1][0] = Y_POS;
			Pos[2][0] = Z_POS;

			dCoordX1 = Pos[0][0];
			dCoordY1 = Pos[1][0];
			dCoordZ1 = Pos[2][0];
		}

		//以下两句仅用于判断该历元是否满足SPP条件（获取POS_IF）
		CMatrix PosX = Pos;
		POS_IF = Get_SPP_POS(temp_oe, PosX);

		if (fabs(dCoordX1) < 1.0 && fabs(dCoordY1) < 1.0 && fabs(dCoordZ1) < 1.0)
		{
			do
			{
				dCoordX1 = Pos[0][0];
				dCoordY1 = Pos[1][0];
				dCoordZ1 = Pos[2][0];

				Get_SPP_POS(temp_oe, Pos);
				if (POS_IF == false)
					break;

			} while (fabs(dCoordX1 - Pos[0][0]) >= 0.1 || fabs(dCoordY1 - Pos[1][0]) >= 0.1 || fabs(dCoordZ1 - Pos[2][0]) >= 0.1);
		}

		//表示无法进行单点定位，也就无法获得初始坐标，暂时标为不可用
		if (POS_IF == false)
		{
			//将初始坐标置为0
			Pos[0][0] = 0;
			Pos[1][0] = 0;
			Pos[2][0] = 0;

			temp_oe.ZHD = 0;

			//用新的坐标更新高度角、方位角、对流层等信息
			for (int k = 0; k < temp_oe.sat_num; k++)
			{
				Sat tempsat;
				tempsat = temp_oe.sat[k];

				if (tempsat.sattype == "G" || tempsat.sattype == "C" || tempsat.sattype == "B" || tempsat.sattype == "E")
				{
					tempsat.E = 0;
					tempsat.A = 0;
					tempsat.Trop_Delay = 0;
					tempsat.Trop_Map = 0;
					tempsat.Sagnac = 0;
					temp_oe.sat[k] = tempsat;
				}
			}
		}

		else
		{
			//更新经纬度
			XYZ[0] = dCoordX1;
			XYZ[1] = dCoordY1;
			XYZ[2] = dCoordZ1;

			double BLH[3] = { 0 };
			OnXYZtoBLH(dCoordX1, dCoordY1, dCoordZ1, BLH);

			temp_o.B = BLH[0];
			temp_o.L = BLH[1];
			temp_o.DH = BLH[2];
			CMatrix HH(3, 3);

			HH[0][0] = -sin(temp_o.B)*cos(temp_o.L);
			HH[0][1] = -sin(temp_o.B)*sin(temp_o.L);
			HH[0][2] = cos(temp_o.B);
			HH[1][0] = -sin(temp_o.L);
			HH[1][1] = cos(temp_o.L);
			HH[1][2] = 0;
			HH[2][0] = cos(temp_o.B)*cos(temp_o.L);
			HH[2][1] = cos(temp_o.B)*sin(temp_o.L);
			HH[2][2] = sin(temp_o.B);

			//考虑对流层再进行一次坐标解算
			for (int k = 0; k < temp_oe.sat_num; k++)
			{
				Sat tempsat;
				tempsat = temp_oe.sat[k];

				if (tempsat.sattype == "G" || tempsat.sattype == "C" || tempsat.sattype == "B" || tempsat.sattype == "E")
				{
					//计算高度角和方位角
					CalSatEA(HH, XYZ, tempsat);

					//计算对流层
					double RTROP, ZHD, VHF, ZWD, VWF;
					TROP.gpt_vmftrop(JD, temp_o.B, temp_o.L, temp_o.DH, tempsat.E, &RTROP, &ZHD, &VHF, &ZWD, &VWF);
					tempsat.Trop_Delay = ZHD * VHF;
					tempsat.Trop_Map = VWF;
					tempsat.trop = ZHD * VHF + ZWD * VWF;
					temp_oe.ZHD = ZHD + ZWD;

					//计算Sagnac效应改正
					double we = 7.2921151467e-5;
					if (tempsat.sattype == "C" || tempsat.sattype == "B")
					{
						we = 7.292115e-5;
					}

					tempsat.Sagnac = we * (tempsat.POS_X * XYZ[1] - tempsat.POS_Y * XYZ[0]) / C_light;
					temp_oe.sat[k] = tempsat;
				}
			}

			//重新计算坐标
			do
			{
				dCoordX1 = Pos[0][0];
				dCoordY1 = Pos[1][0];
				dCoordZ1 = Pos[2][0];

				Get_SPP_POS_Ele(temp_oe, Pos);

			} while (fabs(dCoordX1 - Pos[0][0]) >= 0.1 || fabs(dCoordY1 - Pos[1][0]) >= 0.1 || fabs(dCoordZ1 - Pos[2][0]) >= 0.1);

			XYZ[0] = Pos[0][0];
			XYZ[1] = Pos[1][0];
			XYZ[2] = Pos[2][0];

			//用新的坐标更新高度角、方位角、对流层等信息
			for (int k = 0; k < temp_oe.sat_num; k++)
			{
				Sat tempsat;
				tempsat = temp_oe.sat[k];

				if (tempsat.sattype == "G" || tempsat.sattype == "C" || tempsat.sattype == "B" || tempsat.sattype == "E")
				{
					//计算高度角和方位角
					CalSatEA(HH, XYZ, tempsat);

					//计算对流层
					double RTROP, ZHD, VHF, ZWD, VWF;
					TROP.gpt_vmftrop(JD, temp_o.B, temp_o.L, temp_o.DH, tempsat.E, &RTROP, &ZHD, &VHF, &ZWD, &VWF);
					tempsat.Trop_Delay = ZHD * VHF;
					tempsat.Trop_Map = VWF;
					tempsat.trop = ZHD * VHF + ZWD * VWF;
					temp_oe.ZHD = ZHD + ZWD;

					//计算Sagnac效应改正（地球自转导致的坐标系变化等效改正）
					double we = 7.2921151467e-5;
					if (tempsat.sattype == "C" || tempsat.sattype == "B")
					{
						we = 7.292115e-5;
					}

					tempsat.Sagnac = we * (tempsat.POS_X * XYZ[1] - tempsat.POS_Y * XYZ[0]) / C_light;
					temp_oe.sat[k] = tempsat;
				}
			}
		}

		dCoordX1 = Pos[0][0];
		dCoordY1 = Pos[1][0];
		dCoordZ1 = Pos[2][0];

		//fprintf(fp, "EpNum: %d, %d, %d, %f, %15.8g, %15.8g, %15.8g, %10.5g, %10.5g, %10.5g\n", i+1, temp_oe.hour, temp_oe.minute, temp_oe.second, dCoordX1, dCoordY1, dCoordZ1, dCoordX1 - temp_o.APP_X, dCoordY1 - temp_o.APP_Y, dCoordZ1 - temp_o.APP_Z);
		
		//输出误差文件error file
		//(1)输出文件头
		if (i==0)
		{
			fprintf(fp, "%% GPST                        x-ecef(m)        y-ecef(m)        z-ecef(m)   Q  ns\n");
			fprintf(fp1, "    MARKER_NAME：%s\n", temp_o.marker_name.c_str());
			fprintf(fp1, "       INTERVAL：%3d second\n", (int)temp_o.INTERVAL);
			fprintf(fp1, "       ANT_TYPE：%s\n", temp_o.ANT_TYPE.c_str());
			fprintf(fp1, "       REC_TYPE：%s\n", temp_o.REC_TYPE.c_str());
			fprintf(fp1, "APPROX_POSITION： %14.4f, %14.4f, %14.4f  (m)\n", temp_o.APP_X, temp_o.APP_Y, temp_o.APP_Z);
			fprintf(fp1, "ANTENNA_DELTA_H：         0.0000  (m)\n");
			fprintf(fp1, "ANTENNA_DELTA_E：         0.0000  (m)\n");
			fprintf(fp1, "ANTENNA_DELTA_N：         0.0000  (m)\n");
			fprintf(fp1, "   ANTENNA_PCO1： %7.4f (m), %7.4f (m), %7.4f (m)\n", temp_o.PCO1[0] / 1000.0, temp_o.PCO1[1] / 1000.0, temp_o.PCO1[2] / 1000.0);
			fprintf(fp1, "   ANTENNA_PCO2： %7.4f (m), %7.4f (m), %7.4f (m)\n", temp_o.PCO2[0] / 1000.0, temp_o.PCO2[1] / 1000.0, temp_o.PCO2[2] / 1000.0);
			fprintf(fp1, "       GPS_TYPE：C1\n");
			fprintf(fp1, "   GLONASS_TYPE：C1\n");
			fprintf(fp1, "       BDS_TYPE：C1\n");
			fprintf(fp1, "   GALILEO_TYPE：\n");
			fprintf(fp1, "Tropospheric Delay           : GPT Model\n");
			fprintf(fp1, "Tropospheric Mapping Function: Global Mapping Function\n");
			fprintf(fp1, "Navigation type: broadcast\n");
			fprintf(fp1, "     Clock type: broadcast\n");
			fprintf(fp1, "GLONASS K:   1  -4   5   6   1  -4   5   6  -2  -7   0  -1  -2  -7   0  -1   4  -6   3   2   4  -3   3   2\n");
			fprintf(fp1, "Start Time: %4d   %2d   %2d   %2d   %2d   %2d\n", temp_o.FIRST_TIME[0], temp_o.FIRST_TIME[1], temp_o.FIRST_TIME[2], temp_o.FIRST_TIME[3], temp_o.FIRST_TIME[4], temp_o.FIRST_TIME[5]);
			fprintf(fp1, "  End Time: %4d   %2d   %2d   %2d   %2d   %2d\n", temp_o.LAST_TIME[0], temp_o.LAST_TIME[1], temp_o.LAST_TIME[2], temp_o.LAST_TIME[3], temp_o.LAST_TIME[4], temp_o.LAST_TIME[5]);
			fprintf(fp1, "PRN,s: Satposition(X), Satposition(Y), Satposition(Z),   Sat Clock(m),   Elevation(°),     Azimuth(°),          P1(m),          P2(m),     L1(cycles),     L2(cycles),     Trop Delay,       Trop Map,  Relativity(m),      Sagnac(m),  TGD Effect(m), Antenna Height, Sat Antenna(m),OffsetL1(cycles),OffsetL2(cycles),Windup(cycles),          P3(m)     L3(cycles),          P4(m)     L4(cycles)\n");
		}

		fprintf(fp, "%4d/%2d/%2d %6.2f %14.4f %14.4f %14.4f %3d %3d \n", temp_oe.year, temp_oe.month, temp_oe.day, temp_oe.GPSTIME, dCoordX1, dCoordY1, dCoordZ1, 5, temp_oe.sat_num);

		//进行天线PCO改正
		AntOffsetCorrect(temp_oe, temp_o.PCO1, temp_o.PCO2);
		fprintf(fp1, "Satellite Number: %2.2d,(yyyy-mm-dd-hh-mm-ss):%4d-%2d-%2d-%2d-%2d-%10.7f,zhd: %7.4f, %14.4f, %14.4f, %14.4f, %2d\n", temp_oe.sat_num, temp_oe.year, temp_oe.month, temp_oe.day, temp_oe.hour, temp_oe.minute, temp_oe.second, temp_oe.ZHD, XYZ[0], XYZ[1], XYZ[2], POS_IF);

		for (size_t j = 0; j < temp_oe.sat.size(); j++)
		{
			Sat tempsat = temp_oe.sat[j];

			if (tempsat.E < 10.0 && tempsat.judge_use == 0)
				tempsat.judge_use = 1;

			if (tempsat.sattype == "G")
			{
				fprintf(fp1, "%s%2.2d,%1.1d: %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f,         0.0000,         0.0000, %14.4f, %14.4f,         0.0000, %14.4f, %14.4f\n",
					tempsat.sattype.c_str(), tempsat.numofsat, tempsat.judge_use, tempsat.POS_X, tempsat.POS_Y, tempsat.POS_Z, tempsat.Sat_clock*C_light, tempsat.E, tempsat.A, tempsat.data[0], tempsat.data[1], tempsat.data[3], tempsat.data[4], tempsat.Trop_Delay, tempsat.Trop_Map, tempsat.xdl_t,
					tempsat.Sagnac, tempsat.TGD, tempsat.OffsetL1, tempsat.OffsetL2, tempsat.data[2], tempsat.data[5]);
			}

			if (tempsat.sattype == "C")
			{
				if (tempsat.numofsat == 15)
					tempsat.numofsat = 13;

				fprintf(fp1, "%s%2.2d,%1.1d: %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f,         0.0000,         0.0000, %14.4f, %14.4f,         0.0000, %14.4f, %14.4f\n",
					tempsat.sattype.c_str(), tempsat.numofsat, tempsat.judge_use, tempsat.POS_X, tempsat.POS_Y, tempsat.POS_Z, tempsat.Sat_clock*C_light, tempsat.E, tempsat.A, tempsat.data[0], tempsat.data[1], tempsat.data[3], tempsat.data[4], tempsat.Trop_Delay, tempsat.Trop_Map, tempsat.xdl_t,
					tempsat.Sagnac, tempsat.TGD, tempsat.OffsetL1, tempsat.OffsetL2, tempsat.data[2], tempsat.data[5]);
			}

			if (tempsat.sattype == "B")
			{
				fprintf(fp1, "%s%2.2d,%1.1d: %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f,         0.0000,         0.0000, %14.4f, %14.4f,         0.0000, %14.4f, %14.4f, %14.4f, %14.4f\n",
					tempsat.sattype.c_str(), tempsat.numofsat, tempsat.judge_use, tempsat.POS_X, tempsat.POS_Y, tempsat.POS_Z, tempsat.Sat_clock*C_light, tempsat.E, tempsat.A, tempsat.data[0], tempsat.data[1], tempsat.data[4], tempsat.data[5], tempsat.Trop_Delay, tempsat.Trop_Map, tempsat.xdl_t,
					tempsat.Sagnac, tempsat.TGD, tempsat.OffsetL1, tempsat.OffsetL2, tempsat.data[2], tempsat.data[6], tempsat.data[3], tempsat.data[7]);
			}

			if (tempsat.sattype == "E")
			{
				fprintf(fp1, "%s%2.2d,%1.1d: %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f, %14.4f,         0.0000,         0.0000, %14.4f, %14.4f,         0.0000, %14.4f, %14.4f, %14.4f, %14.4f\n",
					tempsat.sattype.c_str(), tempsat.numofsat, tempsat.judge_use, tempsat.POS_X, tempsat.POS_Y, tempsat.POS_Z, tempsat.Sat_clock*C_light, tempsat.E, tempsat.A, tempsat.data[0], tempsat.data[1], tempsat.data[4], tempsat.data[5], tempsat.Trop_Delay, tempsat.Trop_Map, tempsat.xdl_t,
					tempsat.Sagnac, tempsat.TGD, tempsat.OffsetL1, tempsat.OffsetL2, tempsat.data[2], tempsat.data[6], tempsat.data[3], tempsat.data[7]);
			}
		}

		i++;
		printf("第%5d个历元处理结束\n", i);
	}

	infile.close();
	fclose(fp);
	fclose(fp1);
	return true;
}
