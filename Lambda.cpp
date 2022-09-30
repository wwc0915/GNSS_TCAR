#include "stdafx.h"
#include "Lambda.h"
#include <vector>
#include <cmath>
#include <assert.h>

using namespace std;

CLambda::CLambda(void)
{
}


CLambda::~CLambda(void)
{
}

double CLambda::round(double r)
{
	return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

double gamma(int n)
{
	int i;
	double k;

	if (n % 2) // n为奇数
	{
		k = 1.772453850905516;
		i = 1;
	}
	else
	{
		k = 1.0;
		i = 2;
	}

	while (i <= n - 2)
	{
		k *= i / 2.0;
		i += 2;
	}

	return k;
}

int sign(double x)
{
	int b;
	if (x < 0)
		b = -1;
	if (x > 0)
		b = 1;
	if (fabs(x) < 1e-6)
		b = 0;
	return b;
}

void CLambda::lambda2(CMatrix afloat, CMatrix Qahat)
{
	int ncands = 2;
	int n = Qahat.Row;//获得协方差阵的维数
	afixed.first(n, ncands);//固定解的维数
	sqnorm.first(1, ncands);

	//中间少了一部分

	CMatrix incr(n, 1);
	for (int i = 0; i < n; i++)
	{
		incr[i][0] = afloat[i][0] - fmod(afloat[i][0], 1);
		afloat[i][0] = fmod(afloat[i][0], 1);
	}
	decorrel(Qahat, afloat);
	afloat = z;
	double Chi2 = chistart(D, L, afloat, ncands);
	lsearch(afloat, L, D, Chi2, ncands);
	afixed = (afixed.T() * Z.InvertGaussJordan()).T();
	CMatrix tt;
	tt.first(n, 2);
	for (int i = 0; i < n; i++)
	{
		tt[i][0] = incr[i][0];
		tt[i][1] = incr[i][0];
	}
	for (int i = 0; i < n; i++)
	{
		afixed[i][0] = floor(afixed[i][0] + tt[i][0] + 0.5);
		afixed[i][1] = floor(afixed[i][1] + tt[i][1] + 0.5);
	}
	ratio = sqnorm[0][1] / sqnorm[0][0];
}

void CLambda::decorrel(CMatrix Q, CMatrix a)
{
	int n = Q.Row;
	CMatrix Zti(n);//单位阵
	int i1 = n - 1;
	int sw = 1;

	//分解
	Q.LTDL(this->L, this->D);
	double mu, delta, lambda3, eta, Help;

	while (sw)
	{
		int i = n;
		sw = 0;
		while (!sw && i > 1)
		{
			i = i - 1;
			if (i <= i1)
			{
				for (int j = i + 1; j <= n; j++)
				{
					mu = round(L[j - 1][i - 1]);
					if (mu != 0)
					{
						for (int k = j; k <= n; k++)
						{
							L[k - 1][i - 1] = L[k - 1][i - 1] - mu * L[k - 1][j - 1];
						}//k
						for (int k = 1; k <= n; k++)
						{
							Zti[k - 1][j - 1] = Zti[k - 1][j - 1] + mu * Zti[k - 1][i - 1];
						}//k
					}//if mu.

				}//j
			}//if i<i1
			delta = D[i - 1][i - 1] + L[i][i - 1] * L[i][i - 1] * D[i][i];
			if (delta < D[i][i])
			{
				lambda3 = D[i][i] * L[i][i - 1] / delta;
				eta = D[i - 1][i - 1] / delta;
				D[i - 1][i - 1] = eta * D[i][i];
				D[i][i] = delta;
				for (int k = 1; k <= i - 1; k++)
				{
					Help = L[i][k - 1] - L[i][i - 1] * L[i - 1][k - 1];
					L[i][k - 1] = lambda3 * L[i][k - 1] + eta * L[i - 1][k - 1];
					L[i - 1][k - 1] = Help;
				}
				L[i][i - 1] = lambda3;
				for (int k = i + 2; k <= n; k++)
				{
					Help = L[k - 1][i - 1];
					L[k - 1][i - 1] = L[k - 1][i];
					L[k - 1][i] = Help;
				}
				for (int k = 1; k <= n; k++)
				{
					Help = Zti[k - 1][i - 1];
					Zti[k - 1][i - 1] = Zti[k - 1][i];
					Zti[k - 1][i] = Help;
				}
				i1 = i;
				sw = 1;
			}//if delta
		}//while !sw
	}

	Z = (Zti.T()).InvertGaussJordan();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Z[i][j] = round(Z[i][j]);

		}
	}
	z = Z.T() * a;
}

double  CLambda::chistart(CMatrix D, CMatrix L, CMatrix a, int ncands)
{
	ncands = 2;
	double factor = 1.5;
	int n = afloat.Row;
	double dw, Chi2;
	vector<double> Chi;
	CMatrix temp, temp1;
	CMatrix Linv, Dinv;
	Linv.first(n, n);
	Dinv.first(1, n);
	CMatrix tem2, aa;
	CMatrix temp_afixed(n, 1);
	aa = a;
	if (ncands <= n + 1)
	{
		for (int k = n; k >= 0; k--)
		{
			afloat = a;
			temp_afixed = a;

			for (int i = n; i >= 1; i--)
			{
				dw = 0;
				for (int j = n; j >= i; j--)
				{
					dw = dw + L[j - 1][i - 1] * (afloat[j - 1][0] - temp_afixed[j - 1][0]);
				}
				afloat[i - 1][0] = afloat[i - 1][0] - dw;
				if (i != k)
				{
					temp_afixed[i - 1][0] = round(afloat[i - 1][0]);
				}
				else
					if (fabs(afloat[i - 1][0] - temp_afixed[i - 1][0]) < 1E-6)
					{
						temp_afixed[i - 1][0] = round(temp_afixed[i - 1][0] + 1);
					}
					else
						temp_afixed[i - 1][0] = round(afloat[i - 1][0] + sign((afloat[i - 1][0] - temp_afixed[i - 1][0])));

			}
			//diag(D,temp);
			// afixed.MyTRACE();
			temp1 = (a - temp_afixed).T() * (L.T() * D * L).InvertGaussJordan() * (a - temp_afixed);
			Chi.push_back(temp1[0][0]);
		}//k
		int N = Chi.size();
		double tem;
		for (int i = 0; i < N; i++)
		{
			for (int j = i + 1; j < N; j++)
			{
				if (Chi[i] >= Chi[j])
				{

					tem = Chi[i];
					Chi[i] = Chi[j];
					Chi[j] = tem;
				}
			}
		}

		tem2.first(1, N);
		for (int i = 0; i < N; i++)
		{
			Chi[i] = Chi[i];
			tem2[0][i] = Chi[i];
		}
		Chi2 = Chi[ncands - 1] + 1e-6;
	}//if
	else
	{
		Linv = L.InvertGaussJordan();
		double produ = 1;
		for (int i = 0; i < n; i++)
		{
			Dinv[0][i] = 1 / D[0][i];
			produ *= Dinv[0][i];
		}
		double pi = 3.1415926535898;
		double Vn = (2.0 / n) * pow(pi, n / 2.0) / gamma(n / 2);
		Chi2 = factor * pow((ncands / sqrt(produ) * Vn), 2.0 / n);
	}
	afloat = aa;
	return Chi2;
}

void CLambda::lsearch(CMatrix afloat, CMatrix L, CMatrix D, double Chi2, int ncands)
{
	ncands = 2;
	CMatrix Linv, Dinv, right, left, dq;
	int n = afloat.Row;
	Linv.first(n, n);
	Dinv.first(1, n);
	right.first(n + 1, 1);
	left.first(n + 1, 1);
	dq.first(1, n);
	sqnorm.first(1, ncands);
	afixed.first(n, ncands);
	for (int i = 0; i < n; i++)
	{
		Dinv[0][i] = 1 / D[i][i];
	}
	right[n][0] = Chi2;
	Linv = L.InvertGaussJordan();
	for (int i = 0; i < n - 1; i++)
	{
		dq[0][i] = Dinv[0][i + 1] / Dinv[0][i];
	}
	dq[0][n - 1] = 1 / Dinv[0][n - 1];
	int True = 1;
	int False = 0;
	int cand_n = False;
	int c_stop = False;
	int endsearch = False;
	CMatrix lef, distl, endd;
	lef.first(1, n);
	distl.first(n, 1);
	endd.first(1, n);
	double reach, delta, maxnorm = 0;
	int ipos;
	int ncan = 0;
	int i = n + 1;
	int iold = i;
	int ierr = 0;
	while (!endsearch)
	{
		i = i - 1;
		if (iold <= i)
		{
			lef[0][i - 1] = lef[0][i - 1] + Linv[i][i - 1];
		}
		else
		{
			lef[0][i - 1] = 0;
			for (int j = i + 1; j <= n; j++)
			{
				lef[0][i - 1] = lef[0][i - 1] + Linv[j - 1][i - 1] * distl[j - 1][0];
			}

		}
		iold = i;
		right[i - 1][0] = (right[i][0] - left[i][0]) * dq[0][i - 1];////

		reach = sqrt(right[i - 1][0]);
		delta = afloat[i - 1][0] - reach - lef[0][i - 1];
		distl[i - 1][0] = ceil(delta) - afloat[i - 1][0];
		if (distl[i - 1][0] > (reach - lef[0][i - 1]))
		{
			cand_n = False;
			c_stop = False;
			while (!c_stop && (i < n))
			{
				i = i + 1;
				if (distl[i - 1][0] < endd[0][i - 1])
				{
					distl[i - 1][0] = distl[i - 1][0] + 1;
					left[i - 1][0] = pow(distl[i - 1][0] + lef[0][i - 1], 2);
					c_stop = True;
					if (i == n)
						cand_n = True;
				}
			}//while (! c_stop && i < n)
			if (i == n && !cand_n)
				endsearch = True;
		}//if(distl[i-1][0] >(reach - lef[i-1][0]))
		else
		{
			endd[0][i - 1] = reach - lef[0][i - 1] - 1;
			left[i - 1][0] = pow(distl[i - 1][0] + lef[0][i - 1], 2);
		}
		double t;
		if (i == 1)
		{
			t = Chi2 - (right[0][0] - left[0][0]) * Dinv[0][0];
			endd[0][0] = endd[0][0] + 1;
			while (distl[0][0] <= endd[0][0])
			{
				if (ncan < ncands)
				{
					ncan = ncan + 1;
					for (int j = 0; j < n; j++)
					{
						afixed[j][ncan - 1] = distl[j][0] + afloat[j][0];
					}
					sqnorm[0][ncan - 1] = t;
				}
				else
				{
					maxnorm = sqnorm[0][0];
					ipos = 1;
					for (int j = 1; j <= ncan; j++)
					{
						if (sqnorm[0][j - 1] > maxnorm)
						{
							maxnorm = sqnorm[0][j - 1];
							ipos = j;
						}
					}
				}
				if (t < maxnorm)
				{
					for (int j = 1; j <= n; j++)
					{
						afixed[j - 1][ipos - 1] = distl[j - 1][0] + afloat[j - 1][0];
					}
					sqnorm[0][ipos - 1] = t;
				}


				t = t + (2 * (distl[0][0] + lef[0][0]) + 1) * Dinv[0][0];
				distl[0][0] = distl[0][0] + 1;

			}//while(distl[0][0] <= endd[0][0])
			cand_n = False;
			c_stop = False;
			while (!c_stop && (i < n))
			{
				i = i + 1;
				if (distl[i - 1][0] < endd[0][i - 1])
				{
					distl[i - 1][0] = distl[i - 1][0] + 1;
					left[i - 1][0] = pow(distl[i - 1][0] + lef[0][i - 1], 2);
					c_stop = True;
					if (i == n)
						cand_n = True;
				}
			}
			if (i == n && !cand_n)
				endsearch = True;
		}//i==1
	}//!endsearchp

	CMatrix tem;
	tem.first(ncands, n + 1);
	sqnorm = sqnorm.T();//2x1
	afixed = afixed.T();//2*n
	tem[0][0] = sqnorm[0][0];
	tem[1][0] = sqnorm[1][0];

	for (int j = 0; j < n; j++)
	{
		tem[0][j + 1] = afixed[0][j];
		tem[1][j + 1] = afixed[1][j];
	}
	sqnorm.first(1, 2);
	if (tem[0][0] > tem[1][0])
	{

		sqnorm[0][0] = tem[1][0];
		sqnorm[0][1] = tem[0][0];
		afixed.first(n, 2);
		for (int i = 0; i < n; i++)
		{
			afixed[i][0] = tem[1][i + 1];
			afixed[i][1] = tem[0][i + 1];
		}
	}
	else
	{

		sqnorm[0][0] = tem[0][0];
		sqnorm[0][1] = tem[1][0];
		afixed.first(n, 2);
		for (int i = 0; i < n; i++)
		{
			afixed[i][0] = tem[0][i + 1];
			afixed[i][1] = tem[1][i + 1];
		}
	}
	//sqnorm.MyTRACE();
	ratio = sqnorm[0][1] / sqnorm[0][0];
}

double CLambda::PsCaltor(CMatrix D, int num)
{
	int AmbNum1 = D.Row;
	// 	int AmbNum2 = D.Col;

	if (AmbNum1 != num || AmbNum1 <= 2)
	{
		printf("协方差矩阵维数有问题or模糊度个数小于2");
		assert(false);
	}

	double ps = 1;

	for (int i = 0; i < AmbNum1; i++)
	{
		double sigma = 1 / (2 * sqrt(D[AmbNum1 - i - 1][0]));
		double ptemp = 2 * Normpdf(sigma) - 1;
		ps *= ptemp;
	}
	return ps;
}

double CLambda::Normpdf(double u)
{
	if (u < -5.0) return 0.0;
	if (u > 5.0) return 1.0;

	double y = fabs(u) / sqrt(2.0);

	double p = 1.0 + y * (0.0705230784 + y * (0.0422820123 + y * (0.0092705272 +
		y * (0.0001520143 + y * (0.0002765672 + y * 0.0000430638)))));

	double er = 1 - pow(p, -16.0);
	p = (u < 0.0) ? 0.5 - 0.5 * er : 0.5 + 0.5 * er;
	return p;
}