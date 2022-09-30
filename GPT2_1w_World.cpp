// GPT2_1w_World.cpp: implementation of the CGPT2_1w_World class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "GPT2_1w_World.h"
#include <cmath>
#include <string>
#include <fstream>
using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CGPT2_1w_World::CGPT2_1w_World()
{
	pgrid  = new double *[64800];
	Tgrid  = new double *[64800];
	Qgrid  = new double *[64800];
	dTgrid = new double *[64800];
	u      = new double[64800];  
	Hs     = new double[64800];
	ahgrid = new double *[64800];
	awgrid = new double *[64800];
	lagrid = new double *[64800];
	Tmgrid = new double *[64800];
	for(int i=0; i<64800; i++)
	{
		pgrid[i]  = new double [5];
		Tgrid[i]  = new double [5];
		Qgrid[i]  = new double [5];
		dTgrid[i] = new double [5];
		ahgrid[i] = new double [5];
		awgrid[i] = new double [5];
		lagrid[i] = new double [5];
		Tmgrid[i] = new double [5];
	}
}

CGPT2_1w_World::~CGPT2_1w_World()
{
	for(int i=0; i<64800; i++)        //列数
	{
		delete [] pgrid[i] ;
		delete [] Tgrid[i] ;
		delete [] Qgrid[i] ;
		delete [] dTgrid[i];
		delete [] ahgrid[i];
		delete [] awgrid[i];
		delete [] lagrid[i];
		delete [] Tmgrid[i];
	}
	
	delete [] pgrid ;
	delete [] Tgrid ;
	delete [] Qgrid ;
	delete [] dTgrid;
	delete [] u     ;  
	delete [] Hs    ;
	delete [] ahgrid;
	delete [] awgrid;
	delete [] lagrid;
	delete [] Tmgrid;
}

int CGPT2_1w_World::sign2(double x)
{
	if(x>0)
		return 1;
	else
		return -1;
}

void CGPT2_1w_World::split(std::string& str, std::string pattern, vector<std::string>& arrout)
{
	arrout.clear();
	size_t idx = 0, pos = 0;
	size_t size = str.size();
	for (idx = 0; idx < size; ++idx)
	{
		pos = str.find(pattern, idx);
		if (pos != std::string::npos)
		{
			std::string s = str.substr(idx, pos - idx);
			if (s.size() > 0) arrout.push_back(s);
			idx = pos + pattern.size() - 1;
		}
		else
		{
			if (idx < size)
			{
				std::string s = str.substr(idx, str.size() - idx);
				if (s.size() > 0) arrout.push_back(s);
				break;
			}
		}
	}
}

void CGPT2_1w_World::init_trop()
{
 	string filename ="gpt2_1wA_World.txt";
	fstream infile;
	infile.open(filename, ios::in);

	if (!infile)
	{
		printf("文件打开失败");
		_ASSERT(false);
	}
	string strLine;
	getline(infile, strLine);
	int n =0;

	printf("正在完成对流层信息初始化，请耐心等待！\n");

	while (getline(infile, strLine))
	{ 
		vector <string> arry;
		split(strLine, " ", arry);
		size_t len = arry.size();
		double vel[44];
		for (size_t i = 0; i < len; i++)
		{
			vel[i] = atof(arry[i].c_str());
		}

		for(int i=0;i<5;i++)
		{
			pgrid[n][i]  = vel[i+2];

			double temp = pgrid[n][i];
			Tgrid[n][i]  = vel[i+7];
			Qgrid[n][i]  = vel[i+12]/1000;
			dTgrid[n][i] = vel[i+17]/1000;
			ahgrid[n][i] = vel[i+24]/1000;
			awgrid[n][i] = vel[i+29]/1000;
			lagrid[n][i] = vel[i+34];
			Tmgrid[n][i] = vel[i+39];
		}
		u[n]  = vel[22];
		Hs[n] = vel[23];
		n++;
	} 
	infile.close();
	printf("对流层初始化完成！\n");
}

void CGPT2_1w_World::gpt_vmftrop(double jd,double latitude,double longitude,double height,double Elevation,double *RTROP,double *ZHD,double *VHF,double *ZWD,double *VWF)
{
	double zenith = (90 - Elevation)*PII/180;
	double p,T,dT,Tm,e,ah,aw,la,undu;
	gpt2_1w(jd-2400000.5,latitude,longitude,height,1,0,&p,&T,&dT,&Tm,&e,&ah,&aw,&la,&undu);
	*ZHD = saasthyd(p,latitude,height);
	*ZWD = asknewet(e,Tm,la);
	vmf1_ht(ah,aw,jd-2400000.5,latitude,height,zenith, VHF,VWF);
	*RTROP = (*ZHD)*(*VHF) + (*ZWD)*(*VWF);
}

void CGPT2_1w_World::gpt2_1w (double dmjd,double dlat,double dlon,double hell,int nstat,int it,double *p,double *T,double *dT,double *Tm,double *e,double *ah,double *aw,double *la,double *undu)
{
	//% change the reference epoch to January 1 2000
	double dmjd1 = dmjd-51544.5;

	//% mean gravity in m/s**2
	double gm = 9.80665;
	//% molar mass of dry air in kg/mol
	double dMtr = 28.965/1000;
	//% universal gas constant in J/K/mol
	double Rg = 8.3143;

	//% factors for amplitudes
	double cosfy,coshy,sinfy,sinhy ;
	if (it==1) //% then  constant parameters
	{
		cosfy = 0;
		coshy = 0;
		sinfy = 0;
		sinhy = 0;
	}
	else 
	{
		cosfy = cos(dmjd1/365.25*2*PII);
		coshy = cos(dmjd1/365.25*4*PII);
		sinfy = sin(dmjd1/365.25*2*PII);
		sinhy = sin(dmjd1/365.25*4*PII);
	}
	double plon;  
	//dlon = dlon - 53.5;
	if (dlon < 0) //正的经度
	{
		plon = (dlon + 2*PII)*180/PII;
	}
	else
	{
		plon = dlon*180/PII;
	}
	double ppod = (-dlat + PII/2)*180/PII; 
	int ipod = (int)floor((ppod+1)); 
	int ilon = (int)floor((plon+1)); 
	double diffpod = (ppod - (ipod - 0.5));
	double difflon = (plon - (ilon - 0.5));
	if (ipod == 181)
	{
		ipod = 180;
	}
	if (ilon == 361)
	{
		ilon = 1;
	}
	if( ilon == 0)
	{
		ilon = 360;
	}

	int indx[4];

	//indx[0] = (ipod - 1)*63 + ilon -36*63 - 73;
	indx[0] = (ipod - 1)*360 + ilon;

	int bilinear = 0;

	if (ppod > 0.5 && ppod < 179.5) 
	{
		bilinear = 1;          
	}    

	if (bilinear == 0)
	{
		int ix = indx[0];
		double undu = u[ix];
		double hgt = hell-undu;
	}
	else
	{
		int ipod1 = ipod + sign2(diffpod);
		int ilon1 = ilon + sign2(difflon);

		if (ilon1 == 361)
		{
			ilon1 = 1;
		}
		if (ilon1 == 0)
		{
			ilon1 = 360;
		}

        indx[1] = (ipod1 - 1)*360 + ilon;  //% along same longitude
		indx[2] = (ipod  - 1)*360 + ilon1; //% along same polar distance
        indx[3] = (ipod1 - 1)*360 + ilon1; //% diagonal

// 		indx[1] = (ipod1 - 1)*63+ ilon - 36*63 - 73;  //% along same longitude
// 		indx[2] = (ipod  - 1)*63+ ilon1 - 36*63 - 73; //% along same polar distance
// 		indx[3] = (ipod1 - 1)*63 + ilon1 - 36*63 - 73; //% diagonal
		
		double undul[4];
		double Ql[4];
		double dTl[4];
		double Tl[4];
		double pl[4];
		double ahl[4];
		double awl[4];
		double el[4];
		double lal[4];
		double Tml[4];
		for (int l = 0;l<4;l++)
		{
			undul[l] = u[indx[l]-1];
			//TRACE("%15f\n",undul[l]);
			double hgt = hell-undul[l];
			//TRACE("%15f\n",hgt);
			double T0 = Tgrid[indx[l]-1][0] + 
				Tgrid[indx[l]-1][1]*cosfy + Tgrid[indx[l]-1][2]*sinfy + 
				Tgrid[indx[l]-1][3]*coshy + Tgrid[indx[l]-1][4]*sinhy;
			//TRACE("%15f\n",T0);
			double p0 = pgrid[indx[l]-1][0] + 
				pgrid[indx[l]-1][1]*cosfy + pgrid[indx[l]-1][2]*sinfy + 
				pgrid[indx[l]-1][3]*coshy + pgrid[indx[l]-1][4]*sinhy;
			//TRACE("%15f\n",p0);
			Ql[l] = Qgrid[indx[l]-1][0] + 
				Qgrid[indx[l]-1][1]*cosfy + Qgrid[indx[l]-1][2]*sinfy + 
				Qgrid[indx[l]-1][3]*coshy + Qgrid[indx[l]-1][4]*sinhy;
			//TRACE("%15f\n",Ql[l]);
			//% reduction = stationheight - gridheight
			double Hs1 = Hs[indx[l]-1];
			//TRACE("%15f\n",Hs1);
			double redh = hgt - Hs1;
			//TRACE("%15f\n",redh);
			//% lapse rate of the temperature in degree / m
			dTl[l] = dTgrid[indx[l]-1][0] + 
				dTgrid[indx[l]-1][1]*cosfy + dTgrid[indx[l]-1][2]*sinfy + 
				dTgrid[indx[l]-1][3]*coshy + dTgrid[indx[l]-1][4]*sinhy; 

			//% temperature reduction to station height
			Tl[l] = T0 + dTl[l]*redh - 273.15;
			//TRACE("%15f\n",Tl[l]);
			//% virtual temperature
			double Tv = T0*(1+0.6077*Ql[l]);  
			//TRACE("%15f\n",Tv);
			double c = gm*dMtr/(Rg*Tv);
			//TRACE("%15f\n",c);
			//% pressure in hPa
			pl[l] = (p0*exp(-c*redh))/100;
			//TRACE("%15f\n",pl[l]);
			//% hydrostatic coefficient ah
			ahl[l] = ahgrid[indx[l]-1][0] + 
				ahgrid[indx[l]-1][1]*cosfy + ahgrid[indx[l]-1][2]*sinfy + 
				ahgrid[indx[l]-1][3]*coshy + ahgrid[indx[l]-1][4]*sinhy;
			//TRACE("%15f\n",ahl[l]);
			//% wet coefficient aw
			awl[l] = awgrid[indx[l]-1][0] + 
				awgrid[indx[l]-1][1]*cosfy + awgrid[indx[l]-1][2]*sinfy + 
				awgrid[indx[l]-1][3]*coshy + awgrid[indx[l]-1][4]*sinhy;
			
			//TRACE("%15f\n",awl[l]);
			//% water vapor decrease factor la - added by GP
			lal[l] = lagrid[indx[l]-1][0] + 
				lagrid[indx[l]-1][1]*cosfy + lagrid[indx[l]-1][2]*sinfy + 
				lagrid[indx[l]-1][3]*coshy + lagrid[indx[l]-1][4]*sinhy;
			//TRACE("%15f\n",lal[l]);
			//% mean temperature of the water vapor Tm - added by GP
			Tml[l] = Tmgrid[indx[l]-1][0] + 
				Tmgrid[indx[l]-1][1]*cosfy + Tmgrid[indx[l]-1][2]*sinfy + 
				Tmgrid[indx[l]-1][3]*coshy + Tmgrid[indx[l]-1][4]*sinhy;
			//TRACE("%15f\n",Tml[l]);
			//% water vapor pressure in hPa - changed by GP
			double e0 = Ql[l]*p0/(0.622+0.378*Ql[l])/100;// % on the grid
			//TRACE("%15f\n",e0);
			el[l] = e0*pow((100*pl[l]/p0),(lal[l]+1)); // % on the station height - (14) Askne and Nordius, 1987
			//TRACE("%15f\n",el[l]);
		}//l

		double dnpod1 = fabs(diffpod); //% distance nearer point
		double dnpod2 = 1 - dnpod1;  // % distance to distant point
		double dnlon1 = fabs(difflon);
		double dnlon2 = 1 - dnlon1;

		//% pressure
		double R1 = dnpod2*pl[0]+dnpod1*pl[1];
		double R2 = dnpod2*pl[2]+dnpod1*pl[3];
		*p = dnlon2*R1+dnlon1*R2;

		//% temperature
		R1 = dnpod2*Tl[0]+dnpod1*Tl[1];
		R2 = dnpod2*Tl[2]+dnpod1*Tl[3];
		*T = dnlon2*R1+dnlon1*R2;

		//% temperature in degree per km
		R1 = dnpod2*dTl[0]+dnpod1*dTl[1];
		R2 = dnpod2*dTl[2]+dnpod1*dTl[3];
		*dT = (dnlon2*R1+dnlon1*R2)*1000;

		//% water vapor pressure in hPa - changed by GP
		R1 = dnpod2*el[0]+dnpod1*el[1];
		R2 = dnpod2*el[2]+dnpod1*el[3];
		*e = dnlon2*R1+dnlon1*R2;

		//% hydrostatic
		R1 = dnpod2*ahl[0]+dnpod1*ahl[1];
		R2 = dnpod2*ahl[2]+dnpod1*ahl[3];
		*ah = dnlon2*R1+dnlon1*R2;

		//% wet
		R1 = dnpod2*awl[0]+dnpod1*awl[1];
		R2 = dnpod2*awl[2]+dnpod1*awl[3];
		*aw = dnlon2*R1+dnlon1*R2;

		//% undulation
		R1 = dnpod2*undul[0]+dnpod1*undul[1];
		R2 = dnpod2*undul[2]+dnpod1*undul[3];
		*undu = dnlon2*R1+dnlon1*R2;

		//% water vapor decrease factor la - added by GP
		R1 = dnpod2*lal[0]+dnpod1*lal[1];
		R2 = dnpod2*lal[2]+dnpod1*lal[3];
		*la = dnlon2*R1+dnlon1*R2;

		//% mean temperature of the water vapor Tm - added by GP
		R1 = dnpod2*Tml[0]+dnpod1*Tml[1];
		R2 = dnpod2*Tml[2]+dnpod1*Tml[3];
		*Tm = dnlon2*R1+dnlon1*R2;
	}//==1
}

double CGPT2_1w_World::saasthyd(double p,double dlat,double hell)
{
	double zhd;
	double f = 1-0.00266*cos(2*dlat) - 0.00000028*hell;

	//% calculate the zenith hydrostatic delay
	zhd = 0.0022768*p/f;

	return zhd;
}

double CGPT2_1w_World::asknewet(double e,double Tm,double lambda)
{
	double k1  = 77.604; 				   //% K/hPa
	double k2 = 64.79; 				 //  % K/hPa
	double k2p = k2 - k1*18.0152/28.9644;// % K/hPa
	double k3  = 377600; 				  // % KK/hPa

	//% mean gravity in m/s**2
	double gm = 9.80665;
	//% molar mass of dry air in kg/mol
	double dMtr = 28.965/1000;
	//% universal gas constant in J/K/mol
	double R = 8.3143;

	//% specific gas constant for dry consituents
	double Rd = R/dMtr ;  // % 

	double zwd = 1e-6*(k2p + k3/Tm)*Rd/(lambda + 1)/gm*e;
	return zwd;

}

void CGPT2_1w_World::vmf1_ht(double ah,double aw,double dmjd,double dlat,double ht,double zd,double *vmf1h,double *vmf1w)
{
	double doy = dmjd  - 44239.0 + 1 - 28;

	double bh = 0.0029;
	double c0h = 0.062;
	double phh,c11h,c10h;
	if (dlat<0)    //  %   ! southern hemisphere
	{
		phh  = PII;
		c11h = 0.007;
		c10h = 0.002;
	}
	else            // %   ! northern hemisphere
	{
		phh  = 0.0;
		c11h = 0.005;
		c10h = 0.001;
	}

	double ch = c0h + ((cos(doy/365.250*2.0*PII + phh)+1.0)*c11h/2.0  
		+ c10h)*(1.0-cos(dlat));

	double sine   = sin(PII/2.0 - zd);
	double beta   = bh/( sine + ch  );
	double gamma  = ah/( sine + beta);
	double topcon = (1.0 + ah/(1.0 + bh/(1.0 + ch)));
	*vmf1h   = topcon/(sine+gamma);

	//% C  height correction for hydrotatic part [Niell, 1996]     
	double a_ht = 2.53e-5;
	double b_ht = 5.49e-3;
	double c_ht = 1.14e-3;
	double hs_km  = ht/1000.0;
	beta         = b_ht/( sine + c_ht);
	gamma        = a_ht/( sine + beta);
	topcon       = (1.0 + a_ht/(1.0 + b_ht/(1.0 + c_ht)));
	double ht_corr_coef = 1.0/sine - topcon/(sine + gamma);
	double ht_corr      = ht_corr_coef * hs_km;
	*vmf1h        = *vmf1h + ht_corr;

	double bw = 0.00146;
	double cw = 0.04391;
	beta   = bw/( sine + cw );
	gamma  = aw/( sine + beta);
	topcon = (1.0 + aw/(1.0 + bw/(1.0 + cw)));
	*vmf1w   = topcon/(sine+gamma);
}

void CGPT2_1w_World::free_trop()
{
	delete []u;
	delete []Hs;
	for(int i=0;i<3213;i++)
	{
		delete []pgrid[i] ;
		delete []Tgrid[i] ;
		delete []Qgrid[i] ;
		delete []dTgrid[i] ;
		delete []ahgrid[i] ;
		delete []awgrid[i] ;
		delete []lagrid[i] ;
		delete []Tmgrid[i] ;
	}
	delete []pgrid ;
	delete []Tgrid ;
	delete []Qgrid ;
	delete []dTgrid ;
	delete []ahgrid ;
	delete []awgrid ;
	delete []lagrid ;
	delete []Tmgrid ;
}

double CGPT2_1w_World::Julday(int y,int m,int d,int h,int minute,double second)
{
	double jd ;
	if (m <= 2)
	{
		y = y-1; 
		m = m+12;
	}
	double hour = h + minute/60.0 + second/3600;
	jd = floor(365.25*y) + floor(30.6001 * (m + 1)) + d + hour/24.0 + 1720981.5;
	return jd;
}
