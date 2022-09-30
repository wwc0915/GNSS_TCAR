#pragma once
#include <vector>      //容器
#include <string>      //字符串
#include "Matrix.h"

constexpr double C_light = 299792458.458;
constexpr double PI = 3.14159265359;

//GPS信号频率
#define  f1G 1575.42E+6
#define  f2G 1227.60E+6
#define  f5G 1176.45E+6

//BDS-2信号频率
#define  f1C 1561.098E+6
#define  f2C 1207.14E+6
#define  f3C 1268.52E+6

//BDS-3信号频率
#define  f1B 1561.098E+6
#define  f2B 1176.45E+6     //B2a频率
#define  f3B 1268.52E+6
#define  f4B 1575.42E+6     //B1c频率

//Galileo信号频率
#define f1E 1575.42E+6
#define f2E 1176.45E+6      //E5a频率
#define f3E 1207.14E+6      //E5b频率
#define f4E 1191.795E+6     //E5a+b频率

#define FREQ1       1.57542E9               /* L1/E1  frequency (Hz) */
#define FREQ2       1.22760E9               /* L2     frequency (Hz) */
#define FREQ5       1.17645E9               /* L5/E5a frequency (Hz) */
#define FREQ6       1.27875E9               /* E6/LEX frequency (Hz) */
#define FREQ7       1.20714E9               /* E5b    frequency (Hz) */
#define FREQ8       1.191795E9              /* E5a+b  frequency (Hz) */
#define FREQ9       2.492028E9              /* S      frequency (Hz) */
#define FREQ1_GLO   1.60200E9               /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO   0.56250E6               /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO   1.24600E9               /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO   0.43750E6               /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO   1.202025E9              /* GLONASS G3 frequency (Hz) */
#define FREQ1_BDS   1.561098E9              /* BeiDou B1 frequency (Hz) */
#define FREQ2_BDS   1.20714E9               /* BeiDou B2 frequency (Hz) */
#define FREQ3_BDS   1.26852E9               /* BeiDou B3 frequency (Hz) */

//用于中长基线解算
#define lambda_longwideC1 (C_light / (f3C-f2C))             //BDS2超宽巷（0，-1，1）波长
#define lambda_longwideC2 (C_light/ (f1C+4*f2C-5*f3C))      //BDS2超宽巷（1，4，-5）波长
#define lambda_wideC (C_light / (f1C-f2C))                  //BDS2宽巷（1，-1，0）波长

#define yita0_11 (f1C*f1C*(1/f3C-1/f2C)/(f3C-f2C))         //BDS2（0，-1，1）载波电离层延迟因子
#define yita1_10 (f1C*f1C*(1/f1C-1/f2C)/(f1C-f2C))         //BDS2（1，-1，0）载波电离层延迟因子

#define lambda_longwideB1 (C_light / (f3B-f2B))              //BDS3超宽巷（0，1，-1）（B1I,B3I,B2a）波长
#define lambda_longwideB2 (C_light / (f1B-4*f3B+3*f2B))      //BDS3超宽巷（1，-4，3）（B1I,B3I,B2a）波长
#define lambda_wideB (C_light / (f1B-f3B))                   //BDS3宽巷（1，-1，0）（B1I,B3I,B2a）波长

#define yita01_1 (f1B*f1B*(1/f3B-1/f2B)/(f3B-f2B))         //BDS3（0，1，-1）载波电离层延迟因子
#define yita10_1 (f1B*f1B*(1/f1B-1/f2B)/(f1B-f2B))         //BDS3（1，0，-1）载波电离层延迟因子

#define lambda_L1    (C_light / (FREQ1))		//GPS L1波长
#define lambda_L2    (C_light / (FREQ2))		//GPS L2波长

#define lambda_L1_C    (C_light / (FREQ1_BDS))   //BDS-2 B1波长
#define lambda_L2_C   (C_light / (FREQ2_BDS))    //BDS-2 B2波长
#define lambda_L3_C   (C_light / (FREQ3_BDS))    //BDS-2 B3波长

#define lambda_L1_B     (C_light / (f1B))    //BDS-3 B1波长
#define lambda_L2_B   (C_light / (f2B))      //BDS-3 B2a波长
#define lambda_L3_B   (C_light/ (f3B))       //BDS-3 B3波长

#define lambda_L1_E     (C_light / (f1E))    //Galileo E1波长
#define lambda_L2_E   (C_light / (f2E))     //Galileo E5a波长

#define lambda_L_wide     (C_light / (FREQ1 - FREQ2))   //GPS宽巷波长
#define lambda_L_narrow   (C_light / (FREQ1 + FREQ2))   //GPS窄巷波长

#define lambda_L_wide_C     (C_light / (FREQ1_BDS - FREQ2_BDS))   //BDS-2 宽巷波长
#define lambda_L13_wide_C    (C_light / (FREQ1_BDS - FREQ3_BDS))  //BDS-2 13频宽巷波长

#define lambda_L_narrow_C    (C_light / (FREQ1_BDS + FREQ2_BDS))  //BDS-2 窄巷波长
#define lambda_L13_narrow_C   (C_light / (FREQ1_BDS + FREQ3_BDS)) //BDS-2 13频窄巷波长

#define lambda_L_wide_B     (C_light / (f1B - f2B))     //BDS-3 宽巷波长
#define lambda_L_narrow_B     (C_light / (f1B + f2B))   //BDS-3 窄巷波长

#define lambda_L13_wide_B   (C_light/ (f1B-f3B))        //BDS-3 13频宽巷波长
#define lambda_L13_narrow_B   (C_light / (f1B+f3B))     //BDS-3 13频窄巷波长

#define lambda_L_wide_E     (C_light / (f1E - f2E))    //Galileo 宽巷波长
#define lambda_L_narrow_E     (C_light / (f1E + f2E))  //Galileo 窄巷波长

//GPS和BDS星历的数据结构（文件名后缀***.**p）
struct Ephemeris_CN_epoch
{
	int PRN;
	int year;
	int month;
	int day;
	int hour;
	int minute;
	double second;
	double GPSTIME;
	double TTLSEC;
	double a0;
	double a1;
	double a2;
	double IDOE;
	double Crs;
	double delta_n;
	double M0;
	double Cuc;
	double e;
	double Cus;
	double sqrtA;
	double TOE;
	double Cic;
	double OMEGA;
	double Cis;
	double i0;
	double Crc;
	double w;
	double OMEGA_DOT;
	double i_DOT;
	double code_L2;
	double gps_week;
	double mark_code_L2;
	double pre_sat;                   //精度
	double hel_sat;                   //健康状态
	double TGD;
	double IODC;
	double time_sig_send;             //电文发送时刻
	int doy;                          //年积日（备用）
};

//卫星数据结构
struct Sat
{
	std::string sattype;             //卫星的类型
	int numofsat;                    //卫星的序号（PRN）
	double data[20];                 //卫星的观测值
	double GPSTIME;
	double ttlsec;


	double TGD;                      //TGD改正
	double a0;                       //钟差改正系数1
	double a1;                       //钟差改正系数2
	double a2;                       //钟差改正系数3

	double tk;                       //距离星历节点的外推时间

	double deltt;
	double POS_X;                    //卫星位置X
	double POS_Y;                    //卫星位置Y
	double POS_Z;                    //卫星位置Z

	double r;                        //卫星与测站间距离
	double A;                        //方位角
	double E;                        //高度角

	int posk;                        //星历历元标志
	int health;                      //健康标志
	int system;                      //所属系统
	int judge_use;                   //可用性标志

	double xdl_t;                    //相对论效应的影响（时间 s）
	double trop;                     //对流层延迟1

	double fre_glo;

	//误差项
	double Sat_clock;                //卫星钟差
	double Trop_Delay;               //对流层延迟2
	double Trop_Map;                 //对流层湿延迟投影
	double Relat;                    //相对论
	double Sagnac;                   //地球自转
	double Tide_Effect;              //潮汐效应
	double Antenna_Height;
	double Sat_Antenna;              //卫星天线相位中心改正
	double OffsetL1;                 //L1相位偏差（PCO+PCV）
	double OffsetL2;                 //L2相位偏差（PCO+PCV）
	double Windup;                   //相位缠绕
};

//O文件历元数据结构
struct  Obs_epoch
{
	int liyuan_num;                  //历元序号
	int year;                        //年
	int month;                       //月
	int day;                         //日
	int hour;                        //时
	int minute;                      //分
	double second;                   //秒
	double GPSTIME;                  //GPS时间
	double ttlsec;                   //距离GPS时间原点的总时间
	int dayofy;                      //年积日，在对流层计算的时候用到
	int flag;                        //健康标志
	int sat_num;                     //此历元卫星数量
	std::vector<Sat> sat;
	std::vector<Sat> sat1;           //差分定位中，共视部分需用
	std::vector<Sat> sat2;           //差分定位中，共视部分需用
	double ZHD;                      //天顶对流层延迟（干）
	double ZTD;                      //天顶对流层延迟（总）
	double RClk[8];                  //接收机钟差，G/C/R/E/S/J
	int JClk[8];                     //判断接收机钟差，G/C/R/E/S/J

	double posX;
	double posY;
	double posZ;

	bool ifOK;                       //表征历元初始单点定位是否解算成功
};

//O文件头数据结构
struct Ofileheader
{
	std::string marker_name;         //测站名称
	double APP_X;                    //文件头中的概率坐标X
	double APP_Y;                    //文件头中的概率坐标Y
	double APP_Z;                    //文件头中的概率坐标Z

	//经纬度及高程
	double B;
	double L;
	double DH;

	std::string ANT_TYPE;            //接收机天线类型
	double PCO1[3];                  //根据接收机天线类型通过搜索天线列表匹配的POC改正（频率1）
	double PCO2[3];                  //根据接收机天线类型通过搜索天线列表匹配的POC改正（频率2）
	std::string REC_TYPE;            //接收机类型
	double H, E, N;                  //文件头中的天线偏置
	int FIRST_TIME[6];               //开始时刻
	int LAST_TIME[6];                //结束时刻
	double INTERVAL;                 //采样间隔
	int LEAPSEC;                     //跳秒
};

struct observe_spp
{
	std::string marker_name;         //测站名称
	double INTERVAL;                 //采样间隔
	double APP_X;                    //测站近似坐标
	double APP_Y;
	double APP_Z;
	double L;
	double B;
	double DH;
	int liyuan_num;
	std::vector<Obs_epoch> epoch;
};
