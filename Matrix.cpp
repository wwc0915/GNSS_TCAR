// Matrix.cpp: implementation of the CMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Matrix.h"
#include <math.h>
#include "math.h"
#include <cassert>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMatrix::CMatrix()
{
	p =NULL;
    Col = 0;
    Row = 0;
}

CMatrix::~CMatrix()
{
    delete []p;
}
int CMatrix::GetRows(void)
{
    return this->Row;
}

int CMatrix::GetCols(void)
{
    return this->Col;
}

void CMatrix::ones(int n)
{
	assert(this->Col == this->Row && this->Row>1);
	this->p = new double[n * n];
	for(int i = 0 ; i < n ; i++)
		for(int j = 0 ; j < n ; j++)
			*(p + n * i + j )  = 1;
}

CMatrix CMatrix::combine(CMatrix C)
{
	assert(this->Col==C.Col);
	CMatrix tt(this->Row + C.Row,this->Col);
	for(int i=0;i<this->Row;i++)
	{
		for(int j=0;j<this->Col;j++)
		{
			tt[i][j] = (*this)[i][j];
		}
	}
	
	for(int i=0;i<C.Row;i++)
	{
		for(int j=0;j<this->Col;j++)
		{
			tt[i+this->Row][j] = C[i][j];
		}
	}
	return tt;
}

void CMatrix::ones(int row,int col)
{
	//assert(this->Col==this->Row && this->Row>1);
	this->Row = row;
	this->Col = col;
	this->p = new double[row * col];
	for(int i = 0 ; i < row ; i++)
		for(int j = 0 ; j < col ; j++)
			*(p + row * i + j )  = 1;
}

CMatrix::CMatrix(int n)
{
    this->Row = this->Col = n;
    this->p = new double[n * n];
    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ;j < n; j++)
            if(i == j) *(p + n * i + j) = 1;
            else *(p + n * i + j) = 0;
}
CMatrix::CMatrix(int row ,int col)
{
    this->Row = row;
    this->Col = col;
    this->p = new double[row * col];
    for(int i = 0 ; i < row ; i++)
        for(int j = 0 ; j < col ; j ++)
            *(p +col * i + j) = 0;		
}
CMatrix::CMatrix(double *arrAddress ,int arrWidth)
{
    int arrHeight = 1;
    this->p = new double[arrWidth * arrHeight];
    for(int i = 0 ; i < arrHeight ; i++)
        for(int j = 0 ; j < arrWidth; j++)
            *(p + arrWidth * i + j ) = *(arrAddress + arrWidth * i + j);
		this->Col = arrWidth ; 
		this->Row = arrHeight;
}
CMatrix::CMatrix(double *arrAddress,int rows,int cols)
{
    this->p = new double[rows * cols];
    for(int i = 0 ; i < rows ; i++)
        for(int j = 0 ; j < cols ; j ++)
            *(p +cols * i + j) = *(arrAddress + cols * i + j );
		this->Col = cols;
		this->Row = rows;
}
CMatrix::CMatrix(const CMatrix &m)//拷贝构造函数
{
    this->Col = m.Col;
    this->Row = m.Row;
    p = new double[this->Col * this->Row];
    for(int  i = 0 ; i < this->Row ; i++)
        for(int j = 0 ; j < this->Col ; j ++)
            *(p + this->Col * i + j) = *(m.p + this->Col * i + j);
}
void CMatrix::first(int row ,int col)
{
    delete []this->p;
    this->Row = row;
    this->Col = col;
    this->p = new double[row * col];
    for(int i = 0 ; i < row ; i++)
        for(int j = 0 ; j < col ; j ++)
            *(p +col * i + j) = 0;	
}
bool CMatrix::IsVector(void)//判断矩阵是否为一个数
{
    return !(this->Col == 1 && this->Row ==1);
} 

//从方阵中提取左上角offset*offset的方阵
CMatrix CMatrix::SubMatrix(int offset)
{
    assert(this->Col==this->Row && offset<=this->Col && offset >= 0);
    double *t = new double[offset * offset] ;
    for(int i = 0 ; i < offset ; i++)
        for(int j = 0 ; j < offset ; j++)
            *(t + offset * i + j )  = *(p + this->Col * i + j);
	CMatrix m(t ,offset ,offset);
	delete []t;
	return m;
}
//从方阵中提取右下角offset*offset的方阵
CMatrix CMatrix::RudMatrix(int offset)
{
    assert(this->Col==this->Row && offset<=this->Col && offset >= 0);
    double *t = new double[offset * offset] ;
    for(int i =(this->Row)-offset ; i < this->Row ; i++)
        for(int j = (this->Col)-offset ; j < this->Col ; j++)
            *(t + offset * (i-((this->Row)-offset)) + (j-((this->Col)-offset)) )  = *(p + this->Col * i + j);
	CMatrix m(t ,offset ,offset);
	delete []t;
	return m;
}

void CMatrix::EyeMatrix(int offset,double m)
{
    delete []this->p;
    this->Row = offset;
    this->Col = offset;
    this->p = new double[offset * offset];
    for(int i = 0 ; i < offset ; i++)
        for(int j = 0 ; j < offset ; j ++)
			if(i==j)
				*(p +offset * i + j) = m;
			else
                *(p +offset * i + j) = 0;	
}

double CMatrix::Arg(void)
{
    assert(this->Row == this->Col);
    double result = 1;
    double zhk,max ;
	int i,j,k;
	int index;
	int n=this->GetCols();
    CMatrix temple = *this;
	
    for( i = 0 ; i < this->Row - 1 ; i++)
    {
		//列主元
		max=fabs(temple[i][i]);
		index=i;//最大值的行号
		for(k=i+1;k<n;k++)
		{
			if(fabs(temple[k][i])>max)
			{
				max=fabs(temple[k][i]);
				index=k;
			}
		}
		if(index!=i)//换行
		{
			for(j=i;j<n;j++)
			{
				zhk=temple[i][j];
				temple[i][j]=temple[index][j];
				temple[index][j]=zhk;
			}
		}
		
		
        for( j = i + 1; j < this->Row ; j++)
        {
            zhk = temple[j][i] / temple[i][i];
            temple[j][i] = 0 ;
            for(int nn = i + 1; nn < this->Col ; nn ++)
            {
                temple[j][nn] = temple[j][nn] - zhk * temple[i][nn];
            }
        }
    }
    for( i=0; i < this->Row ; i++)
    {
        for(int j = 0 ; j < this->Col ; j++)
        {
            if(i == j ) 
				result *= temple[i][j];
        }
    }
    return result;
}

//类的对象a可以直接通过a[m][n]来取其中数组的值,第一个[]为下面的重载，第二个[]为一般一维数组的取值
double * CMatrix::operator[](int heightPos)
{
    //assert(heightPos >= 0 && heightPos < this->Row);
    return this->p + this->Col * heightPos;
    return NULL;
}



CMatrix CMatrix::operator *(double alpha)
{
	CMatrix aa(*this);
    for(int i = 0 ; i < this->Row ; i++)
	{
        for(int j = 0 ; j < this->Col ; j++)
            aa[i][j]=aa[i][j]* alpha;
	}
	return aa;
} 

CMatrix CMatrix::operator / (double alpha)
{
	CMatrix aa(*this);
    for(int i = 0 ; i < this->Row ; i++)
	{
        for(int j = 0 ; j < this->Col ; j++)
            aa[i][j]=aa[i][j]/alpha;
	}
	return aa;
} 

/*
CMatrix CMatrix::operator /(CMatrix &m1)
{
CMatrix t = this->T();
t = t * *(this);
t = t.Inver();
t = t * this->T();
t = t * m1;
return t;
}
*/

CMatrix CMatrix::operator +(CMatrix &m1)
{
	assert(this->Row == m1.Row );
	assert(this->Col == m1.Col );

	CMatrix tt(m1.Row,m1.Col);
	
	for(int i=0;i<m1.Row ;i++)
	{
		for(int j=0;j<m1.Col;j++)
		{
			tt[i][j]=(*this)[i][j]+m1[i][j];
		}	
	}
	return tt;
}

CMatrix CMatrix::operator - (CMatrix &m1)
{
	CMatrix tt(m1.Row,m1.Col);
	for(int i=0;i<m1.Row ;i++)
	{
		for(int j=0;j<m1.Col;j++)
		{
			tt[i][j]=(*this)[i][j]-m1[i][j];
		}	
	}
	return tt;
}

CMatrix CMatrix::operator *(CMatrix &m1)
{
        CMatrix ttt(*this);
        int mr = this->Row;
        int mc = m1.Col;
        CMatrix tt(mr ,mc);
        for(int i = 0 ; i < mr ; i++)
        {
            for(int j = 0 ; j < mc; j++)
            {
                for(int ii = 0 ; ii < this->Col; ii++)
                {
                    
                    tt[i][j] += ttt[i][ii] * m1[ii][j];
                }
            }
        }
        return tt;	
}

CMatrix& CMatrix::operator =(const CMatrix &m)
{
    //if(&m == this)return *this;
    this->Row = m.Row;
    this->Col = m.Col;
    delete []this->p;
    p = new double[this->Row * this->Col];
    for(int i = 0 ; i < this->Row ; i++)
    {
        for(int j = 0 ; j < this->Col ; j++)
        {
            *(this->p + this->Col * i + j) = *(m.p + m.Col * i + j);
        }
    }
    return *this;
}

void CMatrix::SetArr(double *m,int p,int q)
{
	
    delete []this->p;
    this->p = new double[p*q];
    for(int i = 0 ; i<p*q; i++)
    {
		*(this->p + i) = *(m+i );
    }
}


bool CMatrix::IsPtv(void)//判断方阵是否可以求逆
{
    assert(this->Col == this->Row);//是方阵才能计算
    bool result = true;
    CMatrix m=*this;
    if(m.Arg() == 0)
    {
		printf("矩阵不可逆，不能求解!");
        result = false;
    }
    return result;
}

double CMatrix::Maxabs()//求矩阵元素绝对值的的最大值
{
	double max=0;
	for(int i=0 ;i<this->GetRows();i++)
	{
		for (int j=0;j<this->GetCols();j++)
		{
			if( fabs((*this)[i][j])>max)
				max=fabs((*this)[i][j]);
		}
	}
	return max;
}
			

CMatrix CMatrix::T(void)//矩阵转置
{
    double *t = new double[this->Col * this->Row] ;
    for(int i = 0 ; i< this->Row ;i++)
        for(int j = 0 ; j < this->Col ; j++)
            *(t + this->Row * j + i) = *(this->p + this->Col * i + j);
		CMatrix m(t ,this->Col ,this->Row);
		delete []t;
		return m;
}

CMatrix CMatrix::Sol(CMatrix &l)//列主元Gauss消去法
{
	assert(this->Row == this->Col && this->Col==l.Row && l.Col==1);
	int i , j, k, n;
	int index;
	double  sum=0;
	double max,zh;
	n=this->Col;
	CMatrix acopy(*this);
	CMatrix b;
	b=l;
	CMatrix c(n,1);
	for(k=0;k<n-1;k++) //k控制约化矩阵acopy时的的行
	{
		//比较大小，换主元
		max=fabs(acopy[k][k]);
		index=k;//最大值的行号
		for(i=k+1;i<n;i++)
		{
			if(fabs(acopy[i][k])>max)
			{
				max=fabs(acopy[i][k]);
				index=i;
			}
		}
		for(j=k;j<n;j++)
		{
			zh=acopy[k][j];
			acopy[k][j]=acopy[index][j];
			acopy[index][j]=zh;
		}
		zh=b[k][0];
		b[k][0]=b[index][0];
		b[index][0]=zh;
		//3个循环进行高斯消元
		for(i=k+1;i<n;i++)
		{
			//zh=acopy[i][k];
		
			for(j=k+1;j<n;j++)
			{
				acopy[i][j]=acopy[i][j]-acopy[k][j]*(acopy[i][k]/acopy[k][k]);
			}
			b[i][0]=b[i][0]-b[k][0]*(acopy[i][k]/acopy[k][k]);
			acopy[i][k]=0;
		}
	}
	c[n-1][0]=b[n-1][0]/acopy[n-1][n-1];
	for(i=n-2;i>=0;i--)
	{
		for(j=i+1;j<n;j++)
		{
			sum=sum+acopy[i][j]*c[j][0];
		}
		c[i][0]=(b[i][0]-sum)/acopy[i][i];
		sum=0;
	}
	return c;
	
}

CMatrix CMatrix::Sol2(CMatrix &l,double w,int &n)//逐次超松弛迭代法
{
	CMatrix a(*this);
	CMatrix b(l);
	CMatrix x(this->GetRows(),1);
	CMatrix xx;
	CMatrix xxx;

	int i,j;
	double sum=0;
	CMatrix f(this->GetRows(),1);
	CMatrix s(this->GetRows(),this->GetCols());
	n=0;
	//给S和f矩阵赋值
	for( i=0 ;i<this->GetRows();i++)
	{
		for (j=0;j<this->GetCols();j++)
		{ 
			if (i==j)
				s[i][j]=1-w;
			else
				s[i][j]=w*(-1)*a[i][j]/a[i][i];
		}
		f[i][0]=w*b[i][0]/a[i][i];
	}

	//迭带
	do
	{
		n++;
		xx=x;
		for(int i=0;i<this->GetRows();i++)
		{
			sum=0;
			for (int j=0;j<this->GetCols();j++)
			{
				sum=sum+s[i][j]*x[j][0];
			}
			x[i][0]=sum+f[i][0];
		}
		xxx=xx-x;
	}while(xxx.Maxabs()>=0.000005);

	return x;
}

CMatrix CMatrix::Inver()
{
    assert(this->Row == this->Col);
	assert(this->IsPtv());
	int i , j, k, n;
	double  sum=0;
	double zh,ch;
	n=this->Col;
	CMatrix e(n,n); //单位阵
	CMatrix acopy(*this);
	CMatrix c(n,n);
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==j)
				e[i][j]=1.0;
			else
				e[i][j]=0.0;
		}
	}
	
    
	for(int l=0;l<n;l++) //l控制单位阵e的列
	{
		for(k=0;k<n;k++) //k控制约化矩阵acopy时的的行
		{
			ch=acopy[k][k];
			e[k][l]=e[k][l]/ch;
			for(j=k;j<n;j++)
				acopy[k][j]=acopy[k][j]/ch;
			for(i=k+1;i<n;i++)
			{
				zh=acopy[i][k];
				for(j=k;j<n;j++)
				{
					acopy[i][j]=acopy[i][j]-zh*acopy[k][j];
				}
				e[i][l]=e[i][l]-zh*e[k][l];
				
			}
		}
		c[n-1][l]=e[n-1][l];
		for(i=n-2;i>=0;i--)
		{
			for(j=i+1;j<n;j++)
			{
				sum=sum+acopy[i][j]*c[j][l];
			}
			c[i][l]=e[i][l]-sum;
			sum=0;
		}
		acopy=*this;
	}
	return c;
}

void CMatrix::DeleteP()
{
	delete [] p;
	p=NULL;
}

void CMatrix::SetSize(int row, int col)
{
	int i,j;
    this->Row = row;
    this->Col = col;
	delete []p;
	p=NULL;
    this->p = new double[row * col];
    for( i = 0 ; i < row ; i++)
        for(j = 0 ; j < col ; j ++)
            *(p +col * i + j) = 0;

}

bool CMatrix::LDLT(CMatrix & LL,CMatrix &DD)
{
	CMatrix a=*(this);
	int i,j,k;
	double sum=0.0;
	int row=a.GetRows();
	int col=a.GetCols();
	if(row!=col)
	{
		printf("该矩阵不是对称正定矩阵，不能LDLT分解!!!");
		return false;
	}
	CMatrix L(row,row);
	CMatrix D(row,row);
	for(i=0;i<row;i++)
	{
		L[i][i]=1;
	}
	
	for(i=0;i<row;i++)
	{
		sum=0.0;
		for(j=0;j<i;j++)
		{
			sum+=L[i][j]*L[i][j]*D[j][j];
		}
		D[i][i]=a[i][i]-sum;
		for(j=i+1;j<row;j++)
		{
			sum=0.0;
			for(k=0;k<i;k++)
			{
				sum+=L[i][k]*L[j][k]*D[k][k];
			}
			L[j][i]=(a[i][j]-sum)/L[i][i]/D[i][i];
		}
	}
	LL=L;
	DD=D;
	return true;
}

bool CMatrix::LTDL(CMatrix & LL,CMatrix &DD)
{
	CMatrix a=*(this);
	int i,j,k;
	double sum=0.0;
	int row=a.GetRows();
	int col=a.GetCols();
	if(row!=col)
	{
		printf("该矩阵不是正定矩阵，不能LDLT分解!!!");
		return false;
	}

	CMatrix L(row,row);
	CMatrix D(row,row);

	for(i=0;i<row;i++)
	{
		L[i][i]=1;
	}
	
	for(i=row-1;i>=0;i--)
	{
		sum=0.0;
		for(j=row-1;j>i;j--)
		{
			sum+=L[j][i]*L[j][i]*D[j][j];
		}
		D[i][i]=a[i][i]-sum;
		for(j=i-1;j>=0;j--)
		{
			sum=0.0;
			for(k=row-1;k>i;k--)
			{
				sum+=L[k][i]*L[k][j]*D[k][k];
			}
			L[i][j]=(a[i][j]-sum)/L[i][i]/D[i][i];
		}
	}

	LL=L;
	DD=D;
    return true;
}

CMatrix CMatrix::Conver()
{
    assert(this->Row == this->Col);
 	assert(this->IsPtv());
	int i , j, k, n;
	int index;
	double max;
	double  sum=0;
	double zh,ch;
	n=this->Col;
	CMatrix e(n,n); //单位阵
	CMatrix acopy(*this);
	CMatrix c(n,n);
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==j)
				e[i][j]=1.0;
			else
				e[i][j]=0.0;
		}
	}
	
	//合并单位阵组成增广矩阵
	CMatrix temple(this->Row ,2*this->Col);
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			temple[i][j]=acopy[i][j];
		}
		temple[i][n+i]=1;
	}

	//列主元消去求逆
	for(k=0;k<n;k++) //k控制增广矩阵temple的行
	{
		//列主元
		max=fabs(temple[k][k]);
		index=k;//最大值的行号
		for(i=k+1;i<n;i++)
		{
			if(fabs(temple[i][k])>max)
			{
				max=fabs(temple[i][k]);
				index=i;
			}
		}
		if(index!=k)//换行
		{
			for(j=k;j<2*n;j++)
			{
				zh=temple[k][j];
				temple[k][j]=temple[index][j];
				temple[index][j]=zh;
			}
		}

		ch=temple[k][k];
		for(j=k;j<2*n;j++)//当行转换
			temple[k][j]=temple[k][j]/ch;
		for(i=0;i<n;i++)//其余行转换
		{
			if(i!=k)
			{
				zh=temple[i][k];
				for(j=k;j<2*n;j++)
				{
					temple[i][j]=temple[i][j]-zh*temple[k][j];
				}	
			}
		}
	}
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			c[i][j]=temple[i][j+n];
		}
	}
	return c;
}


double CMatrix::Tri()
{
    assert(this->Row == this->Col);
	double sum=0;
	CMatrix temp(*this);
	for(int i=0;i<this->Row;i++)
	{
		sum+=temp[i][i];
	}
	return sum;
}

double CMatrix::GetAt(int row, int col)
{
	double a=*(p + this->Col * row + col);
	return a;
}

CMatrix CMatrix::InvertGaussJordan()
{
	assert(this->Row == this->Col);
	assert(this->IsPtv());
	int *pnRow, *pnCol,i,j,k,l,u,v;
    double d = 0, p=0;
	int m_nNumColumns=this->Col;
	CMatrix temp(*this);
	double *m_pData;
	m_pData = new double[m_nNumColumns*m_nNumColumns];
	for(i=0;i<m_nNumColumns;i++)
	{
		for(j=0;j<m_nNumColumns;j++)
		{
			l=i*m_nNumColumns+j;
			m_pData[l]=temp[i][j];
		}
	}

	pnRow = new int[m_nNumColumns];
    pnCol = new int[m_nNumColumns];
	if (pnRow == NULL || pnCol == NULL)
		return false;

	// 消元
    for (k=0; k<=m_nNumColumns-1; k++)
    { 
		d=0.0;
        for (i=k; i<=m_nNumColumns-1; i++)
		{
			for (j=k; j<=m_nNumColumns-1; j++)
			{ 
				l=i*m_nNumColumns+j; p=fabs(m_pData[l]);
				if (p>d) 
				{ 
					d=p; 
					pnRow[k]=i; 
					pnCol[k]=j;
				}
			}
		}
        
		// 失败
		if (d == 0.0)
		{
			delete[] pnRow;
			delete[] pnCol;
			return false;
		}

        if (pnRow[k] != k)
		{
			for (j=0; j<=m_nNumColumns-1; j++)
			{ 
				u=k*m_nNumColumns+j; 
				v=pnRow[k]*m_nNumColumns+j;
				p=m_pData[u]; 
				m_pData[u]=m_pData[v]; 
				m_pData[v]=p;
			}
		}
        
		if (pnCol[k] != k)
		{
			for (i=0; i<=m_nNumColumns-1; i++)
            { 
				u=i*m_nNumColumns+k; 
				v=i*m_nNumColumns+pnCol[k];
				p=m_pData[u]; 
				m_pData[u]=m_pData[v]; 
				m_pData[v]=p;
            }
		}

        l=k*m_nNumColumns+k;
        m_pData[l]=1.0/m_pData[l];
        for (j=0; j<=m_nNumColumns-1; j++)
		{
			if (j != k)
            { 
				u=k*m_nNumColumns+j; 
				m_pData[u]=m_pData[u]*m_pData[l];
			}
		}

        for (i=0; i<=m_nNumColumns-1; i++)
		{
			if (i!=k)
			{
				for (j=0; j<=m_nNumColumns-1; j++)
				{
					if (j!=k)
					{ 
						u=i*m_nNumColumns+j;
						m_pData[u]=m_pData[u]-m_pData[i*m_nNumColumns+k]*m_pData[k*m_nNumColumns+j];
					}
                }
			}
		}

        for (i=0; i<=m_nNumColumns-1; i++)
		{
			if (i!=k)
            { 
				u=i*m_nNumColumns+k; 
				m_pData[u]=-m_pData[u]*m_pData[l];
			}
		}
    }

    // 调整恢复行列次序
    for (k=m_nNumColumns-1; k>=0; k--)
    { 
		if (pnCol[k]!=k)
		{
			for (j=0; j<=m_nNumColumns-1; j++)
            { 
				u=k*m_nNumColumns+j; 
				v=pnCol[k]*m_nNumColumns+j;
				p=m_pData[u]; 
				m_pData[u]=m_pData[v]; 
				m_pData[v]=p;
            }
		}

        if (pnRow[k]!=k)
		{
			for (i=0; i<=m_nNumColumns-1; i++)
            { 
				u=i*m_nNumColumns+k; 
				v=i*m_nNumColumns+pnRow[k];
				p=m_pData[u]; 
				m_pData[u]=m_pData[v]; 
				m_pData[v]=p;
            }
		}
    }


    for(i=0;i<m_nNumColumns;i++)
	{
		for(j=0;j<m_nNumColumns;j++)
		{
			l=i*m_nNumColumns+j;
			temp[i][j]=m_pData[l];
		}
	}
	// 清理内存
	delete[] pnRow;
	delete[] pnCol;
	delete[] m_pData;
	return temp;
   
}
void CMatrix::MyTRACE(void)
{	
    for(int i = 0 ; i < this->Row ; i ++)
    {
        for(int j = 0 ; j < this->Col ; j++)
        {
            printf("%15.7f  " ,*(this->p + this->Col * i + j));
        }
		printf("\n");
    }	
}
