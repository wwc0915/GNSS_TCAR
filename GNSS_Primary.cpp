// GNSS_Primary.cpp : 定义控制台应用程序的入口点。
// GNSS_Primary：本程序将实现GNSS标准单点定位：GPS+BDS双系统单频模式

#include "stdafx.h"
#include <string>      //字符串
#include <vector>      //容器
#include "Function.h"
using namespace std;

//说明此程序为win32程序，仅供入门学习使用，后续请大家完善为MFC框架程序，可进行要处理的文件路径选择等
int _tmain(int argc, _TCHAR* argv[])
{
	//功能1：卫星基本处理和标准单点定位,执行功能1时将下述功能2下代码注释掉
	//定义存储星历的容器：EpochNG、EpochNC、EpochNB
	
	//vector <Ephemeris_CN_epoch> EpochNG; 
	//vector <Ephemeris_CN_epoch> EpochNC;
	//vector <Ephemeris_CN_epoch> EpochNB;
	//vector <Ephemeris_CN_epoch> EpochNE;

	//string pFilename = "20210709p.rnx";  //默认当前工程路径
	//bool pfileread = ReadPfile(pFilename, EpochNG, EpochNC, EpochNB, EpochNE);
	//if (!pfileread)
	//	exit(1);

	//string oFilename = "20210709-2.21O";  //默认当前工程路径
	//bool ofileread = ReadOfileProcess(oFilename, EpochNG, EpochNC, EpochNB, EpochNE);

	//功能2：读入两个spp文件进行伪距差分定位,执行功能2时将上述功能1下代码注释掉
	string sppFilename1 = "SINA_error.txt";  //默认当前工程路径
	string sppFilename2 = "mov7_error.txt";  //默认当前工程路径
	//SPP_DPOS_Pro(sppFilename1, sppFilename2);
	//SPP_Kinematic_Pro(sppFilename1, sppFilename2);
	//MSMF_RTK(sppFilename1, sppFilename2);
	//TCAR_fix(sppFilename1, sppFilename2);
	//TCAR_move(sppFilename1, sppFilename2);
	TCAR_test(sppFilename1, sppFilename2);

	return 0;
}