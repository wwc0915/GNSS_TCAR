// GPT2_1w_World.h: interface for the CGPT2_1w_World class.
//
//////////////////////////////////////////////////////////////////////
#include <vector>
#include <string>

#if !defined(AFX_GPT2_1W_WORLD_H__ED9DCFEA_C0EE_46E4_A30B_552F3C95964A__INCLUDED_)
#define AFX_GPT2_1W_WORLD_H__ED9DCFEA_C0EE_46E4_A30B_552F3C95964A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#define  PII 3.141592653590

#endif // _MSC_VER > 1000

class CGPT2_1w_World  
{
public:
	CGPT2_1w_World();
	~CGPT2_1w_World();

	double **pgrid;
	double **Tgrid;
	double **Qgrid;
	double **dTgrid;
	double *u;
	double *Hs;
	double **ahgrid;
	double **awgrid;
	double **lagrid;
	double **Tmgrid;
	
	void init_trop();
	int sign2(double x);
	void gpt_vmftrop(double jd,double latitude,double longitude,double height,double zenith,double *RTROP,double *ZHD,double *VHF,double *ZWD,double *VWF);
	void gpt2_1w (double dmjd,double dlat,double dlon,double hell,int nstat,int it,double *p,double *T,double *dT,double *Tm,double *e,double *ah,double *aw,double *la,double *undu);
	double saasthyd (double p,double dlat,double hell);
	double asknewet (double e,double Tm,double lambda);
	void vmf1_ht (double ah,double aw,double dmjd,double dlat,double ht,double zd,double *vmf1h,double *vmf1w);
	void free_trop();
	double Julday(int y,int m,int d,int h,int minute,double second);
	void split(std::string& str, std::string pattern, std::vector<std::string>& arrout);
};

#endif // !defined(AFX_GPT2_1W_WORLD_H__ED9DCFEA_C0EE_46E4_A30B_552F3C95964A__INCLUDED_)
