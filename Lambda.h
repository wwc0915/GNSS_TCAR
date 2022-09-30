#pragma once
#include <vector>
#include <iostream>
#include "Matrix.h"

class CLambda
{
public:
	CLambda(void);
	~CLambda(void);
	CMatrix afloat;//模糊度浮点解
	CMatrix Qahat;//协方差阵
	CMatrix sqnorm;
	CMatrix Z;
	CMatrix afixed;//固定解
	CMatrix D;
	CMatrix L;
	CMatrix z;
	double ratio;
	void lambda2(CMatrix afloat, CMatrix Qahat);
	void decorrel(CMatrix Q, CMatrix a);//去相关
	double chistart(CMatrix D, CMatrix L, CMatrix afloat, int ncands);
	void lsearch(CMatrix afloat, CMatrix L, CMatrix D, double Chi2, int ncands);

	double PsCaltor(CMatrix D, int num);
	double Normpdf(double u);
	double round(double r);
};


