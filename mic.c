#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>

#define Debug 5
#define LEVEL_1 6    /* the higest level of debug */
#define LEVEL_2 5    /* only important messages will show */
#define LEVEL_3 4
#define PI 3.1415926
#define LOWEST -FLT_MAX
#define CHECK_DATA 1 

#define EXACT(x) ((floor)(x*1000000))/1000000

#define MAX(a,b) (EXACT(a)>EXACT(b))?a:b
#define MIN(a,b) (EXACT(a)>EXACT(b))?b:a

#if (CHECK_DATA==1)
#undef Debug
#define Debug 0
#endif
/* the definition of a point */
typedef struct Point
{
	double val;
	int   pos;
}data_type;

typedef struct Label
{
	/* rn :row number */
	int rn;
	/* cn: clomn number */
	int cn;
}label_type;
/* the definition of the data */
typedef struct DataSet
{
	data_type * x;
	data_type * y;
	label_type * grid;
	int len;
}Points;

/* the definition of the parameters */
typedef struct Param
{
	float alpha;
	int	c;

}Para;


void gd ( Points * D);
void eqpY ( Points *D, data_type * Dy,int l,int y );
void q_sort( data_type * sd, int l );
void qs(data_type *sd ,int left , int right );
int p_qs(data_type *sd,int left,int right);
void swap(data_type * sd,int a,int b);
void vd(Points * D,data_type *sdx ,int l);
void vyp ( int y , int * yp );
int rnd( float x );
int eqpX (Points *D,data_type * Dx,data_type *Dy, int l, int x, int y, int f);
int clmX ( data_type * Dx, int l, int x, int y, int f, int *clm );
int smr(Points *D,data_type *sdx,data_type * sdy,int a,int b,int y);
int smc(data_type *sdx,int a ,int b,int clmn,int * clm);
void fu(int len ,data_type *SRC_DA,data_type *DES_DB );
void gdrs(Points *D ,int ** rs ); 
void bub_sort(data_type *sd, int l);
void hst(int clmn,int y,int ** rs,float ** Hpq,float **Hst,float **Ai,float *Hq);
void Optx(int x,int y, int clmn, int **rs, float **MI);
void cal_MI(Points *D,data_type * Dx,data_type * Dy,int x ,int y, int f,float **MI);

int main (int argc,char *argv[] )
{	
	Points * D=(Points *)malloc(sizeof(Points));
	D->len=250;
	D->x=(data_type *)malloc(sizeof(data_type)*D->len);
	D->y=(data_type *)malloc(sizeof(data_type)*D->len);
	D->grid=(label_type *)malloc(sizeof(label_type)*D->len);

	data_type *Dx=(data_type *)malloc(sizeof(data_type)*D->len);
	data_type *Dy=(data_type *)malloc(sizeof(data_type)*D->len);

	gd (D);
	fu(D->len,D->x,Dx);
	fu(D->len,D->y,Dy);

	q_sort(Dy,D->len);
	q_sort(Dx,D->len);
	if(Debug>LEVEL_1)
		vd ( D,Dy,D->len );

	float alpha=0.6;
	int f=15;
	int B=rnd((float)pow(D->len,alpha));
	int max_B=rnd((float)B/2);
	int y=0;
	int x=0;
	float **I=(float **)malloc(sizeof(float *)*(max_B+1));
	int i=0;
	for(i=0;i<=max_B;i++)
	{
		I[i]=(float *)malloc(sizeof(float)*(max_B+1));
	}

	float **MI_xy=(float **)malloc(sizeof(float *)*(max_B+1));
	float **MI_yx=(float **)malloc(sizeof(float *)*(max_B+1));
	for(i=0;i<=max_B;i++)
	{
		MI_xy[i]=(float *)malloc(sizeof(float)*(max_B+1));
		MI_yx[i]=(float *)malloc(sizeof(float)*(max_B+1));
	}
	if(Debug>LEVEL_3)
	{
		printf("B:%d\tB/2:%d\tf:%d\talpha:%f\n",B,rnd((float)B/2),f,alpha);
	}
	for(y=2;y<=rnd((float)B/2);y++)
	{
		/* MI will store the mutual information on a fixed y
		   and the clumps on x axis is begin from 2 to x */

		x=(B/y);
		if(Debug>LEVEL_3)
		{
			printf("%d\t%d\n",x,y);
		}
		cal_MI(D,Dx,Dy,x,y,f,MI_xy);
		cal_MI(D,Dy,Dx,x,y,f,MI_yx);
	}

	float **MIC=(float **)malloc(sizeof(float * )*(max_B+1));
	for(i=0;i<=max_B;i++)
	{
		MIC[i]=(float *)malloc(sizeof(float)*(max_B+1));
	}
	float maxMIC=LOWEST;
	for(y=2;y<=max_B;y++)
	{
		for(x=2;x<=max_B;x++)
		{
			if(x*y<=B)
			{
				I[x][y]=MAX(MI_xy[x][y],MI_yx[y][x]);
				MIC[x][y]=(float)I[x][y]/log(MIN(x,y));
				if(EXACT(MIC[x][y])>EXACT(maxMIC))
				{
					maxMIC=MIC[x][y];
				}
				if(Debug>LEVEL_3)
				{
					printf("%10f",MIC[x][y]);
				}
			}

		}
			printf("\n");
	}
	printf("MaxMIC:%f\n",maxMIC);

	free(Dx);
	free(Dy);
	free(D->x);
	free(D->y);
	free(D);

	return 0;
}

void cal_MI(Points *D,data_type * Dx,data_type * Dy,int x ,int y, int f,float **MI)
{
	/* equipartition y axis into y parts */	
	eqpY (D, Dx,D->len,y);

	/*get almost y superclumps on x axis */
	/* super clump number */
	int scln=0;
	scln=eqpX (D,Dy,Dx,D->len,x,y,f);

	/*fill the cell of current grid resolution*/
	int **rs=(int **)malloc(sizeof(int *)*y);
	int i=0;
	int j=0;
	for(i=0;i<y;i++)
	{
		rs[i]=(int *)malloc(sizeof(int)*scln);
		for(j=0;j<scln;j++)
		{
			rs[i][j]=0;
		}
	}
	gdrs(D,rs);

	/* calculate the f values without HQ */
	Optx(x,y,scln,rs,MI);

	/* free memory for cell record */
	for(i=0;i<y;i++)
	{
		free(rs[i]);
	}
	free(rs);

}

void Optx(int x,int y, int clmn, int **rs, float **MI)
{
	float **Hpq=(float **)malloc((clmn+1)*sizeof(float *));
	float **Hst=(float **)malloc((clmn+1)*sizeof(float *));
	float **Ai=(float **)malloc((clmn+1)*sizeof(float *));
	float *Hq=(float *)malloc(sizeof(float));

	int i=0;
	for(i=0;i<=clmn;i++)
	{
		Hpq[i]=(float *)malloc((clmn+1)*sizeof(float));
		Hst[i]=(float *)malloc((clmn+1)*sizeof(float));
		Ai[i]=(float *)malloc((clmn+1)*sizeof(float));
		int j=0;
		for(j=0;j<=clmn;j++)
		{
			Hpq[i][j]=0;
			Hst[i][j]=0;
			Ai[i][j]=0;
		}
	}

	if(Debug>LEVEL_3)
	{
		printf("%d\t%d\t%d\n",x,y,clmn);
	}

	hst(clmn,y,rs,Hpq,Hst,Ai,Hq);

	float *F=(float *)malloc(sizeof(float)*(clmn+1));
	for(i=0;i<=clmn;i++)
	{
		F[i]=0;
	}
	int t=0;
	int s=0;
	int k=clmn;
	float f=0.00;
	/* calculate all the f values of different l */
	float *tp=(float *)malloc(sizeof(float)*(k+1));
	int l=0;
	if(Debug>LEVEL_3)
	{
		printf("Normalized : [ ");
	}
	for(l=2;l<=x;l++)
	{
		for(t=l;t<=k;t++)
		{
			float fmax=LOWEST;
			for(s=l-1;s<t;s++)
			{
				if(l>2)
				{
					f=Ai[s][t]*F[s]+Hst[s][t];
				}
				else
				{
					f=Hpq[s][t];
				}

				if(EXACT(f)>EXACT(fmax))
				{
					fmax=f;
				}
			}
			tp[t]=fmax;

		}
		for(i=0;i<=k;i++)
		{
			F[i]=tp[i];
		}
		MI[l][y]=F[k]+(*Hq);
		if(Debug>LEVEL_3)
		{
			printf("%10f\t",MI[l][y]/log(MIN(l,y)));
		}
	}

	if(Debug>LEVEL_3)
	{
		printf(" ]\n \n");
	}

	/* memory free */

	for(i=0;i<=clmn;i++)
	{
		free(Hpq[i]);
		free(Hst[i]);
		free(Ai[i]);
	}

	free(Hpq);
	free(Hst);
	free(Ai);
	free(tp);
	free(F);
}

void hst(int clmn,int y,int **rs,float **Hpq,float **Hst,float **Ai,float *Hq)
{
	int i=0;
	int j=0;

	int *xq=(int *)malloc(sizeof(int)*clmn);
	for(i=0;i<clmn;i++)
	{
		xq[i]=0;
	}
	int c=0;
	int sum=0;
	for(i=0,c=0;i<clmn;i++)
	{
		for(j=0;j<y;j++)
		{
			/* the sum number of points */
			sum+=rs[j][i];
		}
		xq[c++]=sum;
	}

	int ** yp=(int **)malloc(sizeof(int *)*clmn);
	for(i=0;i<clmn;i++)
	{
		yp[i]=(int *)malloc(sizeof(int)*y);
		for(j=0;j<y;j++)
		{
			yp[i][j]=0;
		}
	}

	if(Debug>LEVEL_3)
	{
		printf("y axis partition: [ ");
	}
	for(j=0;j<y;j++)
	{

		int sum=0;
		for(i=0;i<clmn;i++)
		{

			sum+=rs[j][i];
			yp[i][j]=sum;

		}
		if(Debug>LEVEL_3)
		{

			printf("%d->%d\t",j,yp[clmn-1][j]);
		}

	}
	if(Debug>LEVEL_3)
	{
		printf(" ] \n");
	}

	*Hq=0.00;
	for(i=0;i<y;i++)
	{ 
		int sti=yp[clmn-1][i];
		int st=xq[clmn-1];
		float p1=(float)sti/st;
		if(p1)
		{
			(*Hq)-=(float)p1*log(p1);
		}
	}

	int t=0;
	int s=0;
	int k=clmn;
	for(t=2;t<=k;t++)
	{
		for(s=1;s<t;s++)
		{
			for(i=0;i<y;i++)
			{
				/*test F l=2*/
				int m=xq[t-1];
				int jl=xq[s-1];
				int jr=xq[t-1]-xq[s-1];
				int c1=yp[s-1][i];
				int c2=yp[t-1][i]-yp[s-1][i];

				float p1=(float)c1/m;
				float p2=(float)c2/m;
				float p3=(float)c1/jl;
				float p4=(float)c2/jr;
				if(c1>0)
				{
					Hpq[s][t]+=(float)p1*log(p3);
				}
				if(c2>0)
				{
					Hpq[s][t]+=(float)p2*log(p4);
				}

				int sl=xq[t-1]-xq[s-1];
				int il=yp[t-1][i]-yp[s-1][i];
				p1=(float)sl/m;
				p2=(float)il/sl;
				if(il>0)
				{
					Hst[s][t]+=p2*log(p2);
				}
			}
			int m=xq[t-1];
			int sl1=xq[s-1];
			int sl=xq[t-1]-xq[s-1];
			float p=(float)sl/m;
			Hst[s][t]=(float)Hst[s][t]*p;
			Ai[s][t]=(float)sl1/m;
		}
	}

	if(Debug>LEVEL_3)
	{
		printf("HQ:%f\n",*Hq);
	}
	for(i=0;i<clmn;i++)
	{
		free(yp[i]);
	}
	free(yp);
	free(xq);

}

/* grid resoultuion ,number of points in a cell */
void gdrs(Points *D ,int ** rs ) 
{
	int i=0;
	for(i=0;i< D->len;i++)
	{
		rs[D->grid[i].rn][D->grid[i].cn]+=1;
		if(Debug>LEVEL_1)
		{
			printf("i:%d\trn:%d\tcn:%d\trs:%d\n",i,D->grid[i].rn,D->grid[i].cn,rs[D->grid[i].rn][D->grid[i].cn]);
		}
	}

}

/* get super clumps  */
int eqpX (Points *D,data_type * Dx, data_type *Dy, int l, int x, int y, int f)
{

	int i=0;
	int clmn=f*x;
	int scz=rnd((float)l/clmn);
	int c=0;
	int sci=0;

	if(Debug>LEVEL_3)
	{
		printf("Desired size:%d\n",scz);
	}

	if(Debug>LEVEL_3)
	{
		printf("Clumps partition: [ ");
	}
	while(i<l)
	{	
		int t=0;
		int CAF=0;
		int can=1;
		int sc=1;	
		float T=EXACT(Dx[i].val);
		int flg=0;
		for(t=1;i+t<l;t++)
		{
			if(EXACT(T)==EXACT(Dx[i+t].val))
			{
				/*sc:the possible number of
				  points in current superclump*/
				sc+=1;
			}
			else
			{
				/*here the same row check should not 
				  include the i+t th point */
				if(!smr(D,Dx,Dy,i,t-1,y))
				{
					break;
				}
				else
				{
					flg=1;
					T=EXACT(Dx[i+t].val);
					sc+=1;
				}
			}
		}
		/* whether to cut down the current clump */
		if(flg)
		{
			int j=0;
			for(j=i+t-1;j>i;)
			{
				if(EXACT(Dx[i+t-1].val)==EXACT(Dx[j].val))
				{
					j-=1;
					sc-=1;
				}
				else
				{
					break;
				}
			}
		}

		if((abs(sci+sc-scz)>=abs(sci-scz))&&(sci))
		{
			sci=0;

			/*update grid's cloumn */
			D->grid[Dx[i-1].pos].cn=c;

			if(Debug>LEVEL_3)
			{
				printf("%-5d",i-1);	
			}

			if((i==l-1)||(c==clmn-1))
			{
				scz=l;
			}
			else
			{
				scz=(l-i-1)/(clmn-c);
			}
			c++;
		}
		int j=0;
		for(j=i;j<=i+sc-1;j++)
		{
			D->grid[Dx[j].pos].cn=c;
			if(Debug>LEVEL_1)
			{
				printf("sdx[%d]:%f\tpos:%d\tcn:%d\tc:%d\n",i,Dx[j].val,Dx[j].pos,D->grid[Dx[j].pos].cn,c);	
			}
		}
		i+=sc;
		sci+=sc;
	}
	if(Debug>LEVEL_3)
	{
		printf("%d ]\n",i-1);	
	}
	return (c+1);
}

/* check if all the points in one clump are in the same row */
int smr (Points *D,data_type *sdx,data_type *sdy,int a,int b,int y)
{
	int j=0;
	for(j=a;j+1<=a+b;j++)
	{
		if(D->grid[sdx[j].pos].rn!=D->grid[sdx[j+1].pos].rn)
		{
			return 0;
		}
	}
	return 1;

}

/* eqpY : equipartition the y axis 	*/
void eqpY (Points *D, data_type * Dy ,int l, int y )
{
	int i=0;
	int s=0;
	int c=0;
	int bi=0;
	int bsz=rnd((float)l/y);
	while(i<l)
	{
		/* calculate the number of 
		   points that have the same
		   y values				*/
		/* count begin */
		s=1;
		int j=0;
		for(j=1;j+i<l;j++)
		{
			if(Dy[i].val==Dy[i+j].val)
			{
				s+=1;
			}
			else
				break;
		}
		/* count end */
		if(abs(bi+s-bsz)>=abs(bi-bsz)&&bi!=0)
		{
			if(Debug>LEVEL_1)
			{
				printf("break point\n");
				printf("c:%d\ti:%d\tbi:%d\ts:%d\tbsz:%d\n",c,i,bi,s,bsz);
				printf("%f\t%f\t\n",Dy[i].val,Dy[i+1].val);
			}
			/* each new bin the 
			   bi should be reset */
			bi=0;

			/* update the grid */
			D->grid[Dy[i-1].pos].rn=c;
			if(Debug>LEVEL_2)
			{
				printf("sdy[%d]:%f\tpos:%d\trn:%d\tc:%d\n",i,Dy[i-1].val,Dy[i-1].pos,D->grid[Dy[i-1].pos].rn,c);	
			}
			if(c==y-2)
			{
				bsz=l;
			}
			else
			{
				bsz=rnd((float)(l-i)/(y-c-1));
			}
			c++;
		}

		/* update the grid */
		for(j=i;j<=i+s-1;j++)
		{
			D->grid[Dy[j].pos].rn=c;
			if(Debug>LEVEL_1)
			{
				printf("sdy[%d]:%f\tpos:%d\trn:%d\tc:%d\n",i,Dy[j].val,Dy[j].pos,D->grid[Dy[j].pos].rn,c);	
			}
		}
		bi+=s;
		i+=s;
	}

	if(Debug>LEVEL_2)
	{
		printf("sdy[%d]:%f\tpos:%d\trn:%d\tc:%d\n",i,Dy[l-1].val,Dy[l-1].pos,D->grid[Dy[l-1].pos].rn,c);	
	}
}

/* check the y axis partition */
void vyp ( int y , int * yp )
{
	int i=0;
	printf("Index\tPosition\tNumber\n");
	for(i=0;i<y;i++)
	{
		printf("%d\t%d\t",i,yp[i]);
		if(i>0)
		{
			printf("%d\n",yp[i]-yp[i-1]);		
		}
		else
		{
			printf("%d\n",yp[i]+1);
		}
	}

}

/* generate a test data set*/
void gd (Points * D )
{
	int i=0;
	int j=0;
	while(i < D->len)
	{
		D->x[i].pos=i;
		D->y[i].pos=i;

		D->x[i].val=i*PI/2.00;
		D->y[i].val=sin(D->x[i].val);

		D->x[i].val=EXACT(D->x[i].val);
		D->y[i].val=EXACT(D->y[i].val);
		if(Debug>LEVEL_1)
		{
			printf("%d:\t%f\t%f\n",i,D->x[i].val,D->y[i].val);
		}
		i++;
	}
	/*
	for(j=1,i=0;i<=23;i++,j++)
	{
		D->y[i].val=0;
		if(Debug>LEVEL_1)
		{
			printf("%d->%d:\t%f\t%f\n",j,i,D->x[i].val,D->y[i].val);
		}
	}
	*/
	if(CHECK_DATA==1)
	{
		i=0;
		printf("Index,Value\n");
		for(i=0;i<D->len;i++)
		{
			printf("%f,%f\n",D->x[i].val,D->y[i].val);
		}
	}
}

/*  fill DES_DB with SRC_DA*/
void fu (int len ,data_type *SRC_DA,data_type *DES_DB )
{
	int i=0;
	for(i=0;i<len;i++)
	{
		DES_DB[i].val=SRC_DA[i].val;
		DES_DB[i].pos=SRC_DA[i].pos;
	}
}

/* view the data x or y  */
void vd (Points *D,data_type * sdx, int l)
{
	int i=0;
	puts("\n");
	for(i=0;i<l;i++)
	{
		int flag=0;
		if(i+1<l)
		{

			if(sdx[i].val==sdx[i+1].val)
			{
				flag=1;
			}
		}
		printf("%d\t%lf\t%lf\t%d\n",i,sdx[i].val,D->x[sdx[i].pos].val,flag);
	}
}

/* round the float number  */
int rnd ( float x )
{
	return floor((x*10+4)/10.00);
}

/* quick sort an array 
   q_sort , qs , p_qs, swap */
void q_sort ( data_type * sd, int l )
{
	int	left=0;
	int right=l-1;
	qs(sd ,left , right );

}
void qs (data_type *sd ,int left , int right )
{
	if(left<right)
	{
		int m=p_qs(sd,left,right);
		qs(sd,left,m-1);
		qs(sd,m+1,right);
	}

}
int p_qs (data_type *sd,int left,int right)
{
	double t=sd[left].val;
	int b=left,i;
	for(i=left;i<=right;i++)
	{
		if(sd[i].val<t)
		{
			b++;
			swap(sd,b,i);

		}
	}
	swap(sd,b,left);
	return b;
}
void swap (data_type * sd,int a,int b)
{
	data_type temp=sd[a];
	sd[a]=sd[b];
	sd[b]=temp;
}

/* this is the very way to sort */
void bub_sort(data_type *sd, int l)
{
	int i=0;
	int j=0;
	for(i=0;i<l-1;i++)
	{
		for(j=0;j<l-1;j++)
		{
			if(sd[j].val>sd[j+1].val)
			{
				swap(sd,j,j+1);
			}
		}
	}
}
