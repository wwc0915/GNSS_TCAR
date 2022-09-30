#include <fstream>     //文件流 
#include <string>      //字符串
#include <vector>      //容器
#include "Myheader.h"

//读入星历数据（p文件）并存储
bool ReadPfile(std::string Filename, std::vector <Ephemeris_CN_epoch> &EpochNG, std::vector <Ephemeris_CN_epoch> &EpochNC, std::vector <Ephemeris_CN_epoch> &EpochNB, std::vector < Ephemeris_CN_epoch> &EpochNE);

//读入观测数据（o文件）并进行相应的误差处理和单点定位s
bool ReadOfileProcess(std::string Filename, std::vector <Ephemeris_CN_epoch> &EpochNG, std::vector <Ephemeris_CN_epoch> &EpochNC, std::vector <Ephemeris_CN_epoch> &EpochNB, std::vector < Ephemeris_CN_epoch>& EpochNE);

//读入两个spp文件进行差分定位解算，后续组建自己程序时，建议读取两个o文件直接处理
bool SPP_DPOS_Pro(std::string filename1, std::string filename2);

//读入两个spp文件进行RTK解算
bool SPP_Kinematic_Pro(std::string filename1, std::string filename2);

//读入两个spp文件进行多系统RTK解算
bool MSMF_RTK(std::string filename1, std::string filename2);

//读入两个spp文件进行多系统RTK解算（除BDS3外）
bool trisys_RTK(std::string filename1, std::string filename2);

//函数功能，TCAR法固定模糊度,采用几何相关模型,基准站
void TCAR_fix(std::string& filename1, std::string& filename2);

//函数功能，TCAR法固定模糊度,采用几何相关模型,流动站
void TCAR_move(std::string& filename1, std::string& filename2);

//函数功能，TCAR法固定模糊度,采用几何相关模型,流动站,合起来算
void TCAR_test(std::string& filename1, std::string& filename2);