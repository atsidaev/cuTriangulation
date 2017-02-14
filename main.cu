#include <GL/glew.h>
#include <GL/glu.h>
#include <GL/glut.h>
//#include <cutil_inline.h>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include < math.h > 
#include <malloc.h>
#define SIZE 50
typedef struct point 
{
    float x;
    float y;
} Point;

typedef struct line 
{
    Point a;
    Point b;
} Line;
__device__ __host__ void quickSort(Point *mas, int size,int sortType,int t)
{
	int i = 1;
    int j = 2;
	
	while (i<size+1)
	{
		bool m=false;
		if(sortType==1)
		{
			if(mas[i-1].x==mas[i].x)
				m=mas[i-1].y<mas[i].y;
			else if(mas[i-1].y==mas[i].y)
			{
				m=t*mas[i-1].x>t*mas[i].x;
			}
			if ((m)||(mas[i-1].y<mas[i].y)) { i=j; j++; }
			else
			{
				Point tmp=mas[i];
				mas[i]=mas[i-1];
				mas[i-1]=tmp;
				i--;
				if (i==0) { i=j; j++; }
			}
		}
		else
		{
			if(mas[i-1].x==mas[i].x)
					m=mas[i-1].y>mas[i].y;
			if ((m)||(mas[i-1].x<mas[i].x)) { i=j; j++; }
			else
			{
				Point tmp=mas[i];
				mas[i]=mas[i-1];
				mas[i-1]=tmp;
			i--;
			if (i==0) { i=j; j++; }
			}
			
		}
	}
}
__device__  __host__ bool isIntersect(Point ax,Point bx,Point cx,Point dx)
 {
	float x1,x2,x3,x4;
	float y1,y2,y3,y4;
	float k1,k2;
	float b1,b2;
	float x,y;
	x1=ax.x;
	x2=bx.x;
	y1=ax.y;
	y2=bx.y;

	x3=cx.x;
	x4=dx.x;
	y3=cx.y;
	y4=dx.y;


	if(y1==y2||x1==x2)k1=0;
	else k1=(y2-y1)/(x2-x1);
	if(y3==y4||x3==x4)k2=0;
	else k2=(y4-y3)/(x4-x3);
	if(k1==k2&&k1!=0)return false;
	b1=y1-k1*x1;
	b2=y3-k2*x3;
	if(x1==x2)
	{
		x=x1;
		y=k2*x+b2;
		if(k1==k2)
		y=k1*x+b1;
	}
	else if (x3==x4)
	{
		x=x3;
		y=k1*x+b1;
	}
	else 
	{
		x=(b2-b1)/(k1-k2);
		if(y1==y2)
			y=y1;
		else if(y3==y4)
			y=y3;
		else
		y=k1*x+b1;
	}
	if(x1>x2)
	{
		x1=bx.x;
		x2=ax.x;
	}
	if(y1>y2)
	{
		y1=bx.y;
		y2=ax.y;
	}
	if(x3>x4)
	{
		x3=dx.x;
		x4=cx.x;
		
	}
	if(y3>y4)
	{
		y3=dx.y;
		y4=cx.y;
	}
	if(((x1<=x&&x<=x2))&&((x3<=x&&x<=x4)))
	{	
		if(((y1<=y&&y<=y2))&&((y3<=y&&y<=y4)))return true;
		return false;
	}
	return false;

 }
  __host__ void drawLine(Point a,Point b)
 {
	glBegin(GL_LINES);
	glVertex2f(a.x,a.y);
	glVertex2f(b.x,b.y);
	glEnd();
 }
 __device__ __host__ bool isPointInCircle(Point *points,Point point)
 {

	 float a=points[0].y-points[1].y;
	 float b=points[1].x-points[0].x;
	 float c=points[2].x-points[0].x;
	 float d=points[2].y-points[0].y;

	 float x1,x2,x3;
	 float y1,y2,y3;
	 float cosin=((a*c+b*d)/(sqrt(a*a+b*b)*sqrt(c*c+d*d)));
	 float s;

	 if(cosin>0)
	 {
		 x1=points[1].x-point.x;
		 x2=points[0].x-point.x;
		 x3=points[2].x-point.x;

		 y1=points[1].y-point.y;
		 y2=points[0].y-point.y;
		 y3=points[2].y-point.y;
	 }
	 else if(cosin<0)
	 {
		 x1=points[0].x-point.x;
		 x2=points[1].x-point.x;
		 x3=points[2].x-point.x;

		 y1=points[0].y-point.y;
		 y2=points[1].y-point.y;
		 y3=points[2].y-point.y;
	 }
	 else
	 {
		return false;
	 }
	 s=(x1*x1 + y1*y1)*(y2*x3 - x2*y3) + (x2*x2 + y2*y2)*(x1*y3 - y1*x3) + (x3*x3 + y3*y3)*(y1*x2 - x1*y2) ;

	if( s <= 0)

		return true;
	return false;
 }
 __device__ __host__ void addLine(Line *lines,int *sizeL,Point a,Point b)
 {
	 lines[*sizeL].a=a;
	 lines[*sizeL].b=b;
	*sizeL=*sizeL+1;
	 
 }
 
 void connectWithFourthPoint(Point *points,Line *lines,int *sizeL)
 {
	float minDistance=sqrt((points[0].x-points[3].x)*(points[0].x-points[3].x)+(points[0].y-points[3].y)*(points[0].y-points[3].y));
	float distance1=sqrt((points[1].x-points[3].x)*(points[1].x-points[3].x)+(points[1].y-points[3].y)*(points[1].y-points[3].y));
	float distance2=sqrt((points[2].x-points[3].x)*(points[2].x-points[3].x)+(points[2].y-points[3].y)*(points[2].y-points[3].y));
	if(distance1<minDistance)
	{
		if(isIntersect(points[3],points[1],points[2],points[0])==false)
		{
		addLine(lines,sizeL,points[1],points[3]);
		}
		if(distance2<minDistance)
			{
				addLine(lines,sizeL,points[2],points[3]);
				if(isIntersect(points[3],points[0],points[2],points[1])==false)
				{
					addLine(lines,sizeL,points[3],points[0]);
				}
			}
		else if(distance2>=minDistance)
		{	
			addLine(lines,sizeL,points[0],points[3]);
			if(isIntersect(points[3],points[2],points[1],points[0])==false)
			{
				addLine(lines,sizeL,points[3],points[2]);
			}
		}
	}
	else if(distance1>=minDistance)
	{
		if(isIntersect(points[3],points[0],points[2],points[1])==false)
		{
		addLine(lines,sizeL,points[0],points[3]);
		}
		if(distance2<distance1)
		{
			addLine(lines,sizeL,points[2],points[3]);
			if(isIntersect(points[3],points[1],points[2],points[0])==false)
			{
				addLine(lines,sizeL,points[3],points[1]);
			}
		}
		else if(distance2>distance1)
		{
			addLine(lines,sizeL,points[1],points[3]);
			if(isIntersect(points[3],points[2],points[0],points[1])==false)
			{
				addLine(lines,sizeL,points[3],points[2]);
			}
		}
	}
 }
 __device__ __host__ bool isDelaunayCondition(Point *points,Point *mas,int length)
  {
	for(int i=0;i<length;++i)
	{
		if(isPointInCircle(points,mas[i])==false)
		return false;
	}
	return true;
  }
 
  __device__  float getAngle(Point Ax,Point Bx,Point Cx)
	{
		 float a=Ax.y-Bx.y;
		 float b=Ax.x-Bx.x;
		 float c=Cx.x-Bx.x;
		 float d=Cx.y-Bx.y;
		 if(sqrt(a*a+b*b)*sqrt(c*c+d*d)!=0)
		 return acos((a*d+b*c)/(sqrt(a*a+b*b)*sqrt(c*c+d*d)))*(180/3.14);
		 return 0.0;
	}
	__device__ void findBorderLine(int start1,int start2,int end1,int end2,short c,short dir,Point *mas1,Point *mas2,Point *SearchPoint1,Point *SearchPoint2,float Y1,float Y2)
{
	for(int i=start1;c*i<=c*end1;i+=c)
	{
		if(dir*mas1[i].x<dir*(*SearchPoint1).x)
		{
			Point tmp;
			tmp.x=(*SearchPoint1).x;
			tmp.y=Y1;
			if(!isIntersect(*SearchPoint2,mas1[i],*SearchPoint1,tmp))
			{
				*SearchPoint1=mas1[i];
			}
			else if((*SearchPoint1).x!=mas1[start1].x&&(*SearchPoint1).y!=mas1[start1].y)
			{
				break;
			}
			
		}
	}
	for(int i=start2;c*i<=c*end2;i+=c)
	{
		if(dir*mas2[i].x<dir*(*SearchPoint2).x)
		{		
			Point tmp;
			tmp.x=mas2[i].x;
			tmp.y=Y2;
			if(isIntersect(*SearchPoint2,*SearchPoint1,mas2[i],tmp))
			{
				*SearchPoint2=mas2[i];
				break;
			}
		}
	}


}

__device__ void buildTriangulation(Point *a,Point *b,int *sizeB,int *sizeA,Point *newA,Point *newB,Line *lines,int *sizeL,Line *newLines,int *sizeNL,Line downLine,int oldSizeA,int oldSizeB,bool *wasAddOldPoints)
 {

	 if((*newA).x!=downLine.a.x||(*newA).y!=downLine.a.y||(*newB).x!=downLine.b.x||(*newB).y!=downLine.b.y)
	 { 
		bool BA=false,AB=false;
		float min1=180.0;
		float min2=180.0;
	
		Point points[3];
		points[0]=(*newA);
		points[1]=(*newB);

		if((*newA).y==a[(*sizeA)-1].y&&(*newB).y==a[(*sizeA)-1].y)
		{
				(*sizeA)=(*sizeA)-1;
		}
		points[2]=a[(*sizeA)-1];
	
		BA=isDelaunayCondition(points, b,oldSizeB+1)&&isDelaunayCondition(points, a,oldSizeA+1);
		
		if((*newA).y==b[(*sizeB)-1].y&&(*newB).y==b[(*sizeB)-1].y)
		{
				(*sizeB)=(*sizeB)-1;
		}
		points[2]=b[(*sizeB)-1];
		AB=isDelaunayCondition(points, a,oldSizeA+1)&&isDelaunayCondition(points, b,oldSizeB+1);
	
		if(BA&&AB)
		{
			if(*wasAddOldPoints)
			{
				for(int i=0;i<*sizeNL;i++)
				{
					if((newLines[i].a.x== (*newA).x&&newLines[i].a.y== (*newA).y)&&( newLines[i].b.x==b[(*sizeB)-1].x&& newLines[i].b.y==b[(*sizeB)-1].y)||
						(newLines[i].a.x==b[(*sizeB)-1].x&&newLines[i].a.y== b[(*sizeB)-1].y)&&( newLines[i].b.x==(*newA).x&& newLines[i].b.y==(*newA).y))
						{
							AB=false;
							break;
						}
					if((newLines[i].a.x== (*newB).x&&newLines[i].a.y== (*newB).y)&&( newLines[i].b.x==a[(*sizeA)-1].x&& newLines[i].b.y==a[(*sizeA)-1].y)||
						(newLines[i].a.x==a[(*sizeA)-1].x&&newLines[i].a.y== a[(*sizeA)-1].y)&&( newLines[i].b.x==(*newB).x&& newLines[i].b.y==(*newB).y))
						{
							BA=false;
							break;
						}
				}
			}
			if(AB&&BA)
			{
				
				float alfa1=getAngle( (*newB), (*newA), a[(*sizeA)-1]);
				float alfa2=getAngle( (*newA), (*newB), b[(*sizeB)-1]);
				 min1=alfa1;
				 min2=alfa2;
				if(alfa1<180.0&&alfa2<180.0)
				{
 
					alfa1=getAngle( (*newA), (*newB), a[(*sizeA)-1]);
					if(180.0-(alfa1+min1)<min1)min1=180.0-(alfa1+min1);
					if(alfa1<min1)min1=alfa1;

 
					alfa2=getAngle( (*newB), (*newA), b[(*sizeB)-1]);
					if(180.0-(alfa2+min2)<min2)min2=180.0-(alfa2+min2);
					if(alfa2<min2)min2=alfa2; 
				}
				else if(alfa1>=180.0)min1=0;
				else if(alfa2>=180.0)min2=0;
			} 
		}
		if((min1>=min2&&min1!=180.0)||(BA&&!AB))
		{
			points[2]=a[(*sizeA)-1];
			if(*wasAddOldPoints)
			{
				for(int i=0;i<*sizeNL;i++)
				{
					
					if((newLines[i].a.x== (*newB).x&&newLines[i].a.y== (*newB).y)&&( newLines[i].b.x==a[(*sizeA)-1].x&& newLines[i].b.y==a[(*sizeA)-1].y)||
						(newLines[i].a.x==a[(*sizeA)-1].x&&newLines[i].a.y== a[(*sizeA)-1].y)&&( newLines[i].b.x==(*newB).x&& newLines[i].b.y==(*newB).y))
						{
							BA=false;
							(*sizeA)=(*sizeA)-1;
							break;
						}
				}
			}
			
	
			if(BA==true)
			{
				for(int i=0;i<*sizeL;i++)
				{
				if((lines[i].a.x!=(*newB).x||lines[i].a.y!= (*newB).y)&& (lines[i].b.x!=a[(*sizeA)-1].x|| lines[i].b.y!=a[(*sizeA)-1].y))
					if((lines[i].a.x!=a[(*sizeA)-1].x||lines[i].a.y!= a[(*sizeA)-1].y)&&( lines[i].b.x!=(*newB).x|| lines[i].b.y!=(*newB).y))
					if(isIntersect( lines[i].a, lines[i].b, (*newB), a[(*sizeA)-1]))
					{
						for (int j=i;j<*sizeL-1;j++)
						{
							lines[j]=lines[j+1];
						}
						*sizeL=*sizeL-1;
						i--;
					}
				}
				addLine(newLines,sizeNL,(*newB),a[(*sizeA)-1]);
				(*newA)=a[(*sizeA)-1];
				if((*sizeA)>1)
				{
					(*sizeA)=(*sizeA)-1;
				}
			}
			
		  
		}
		else if((min2>min1&&min2<180.0)||(!BA&&AB))
		{
			if(*wasAddOldPoints)
			{
				for(int i=0;i<*sizeNL;i++)
				{
					if((newLines[i].a.x== (*newA).x&&newLines[i].a.y== (*newA).y)&&( newLines[i].b.x==b[*sizeB-1].x&& newLines[i].b.y==b[*sizeB-1].y)||
						(newLines[i].a.x==b[*sizeB-1].x&&newLines[i].a.y== b[*sizeB-1].y)&&( newLines[i].b.x==(*newA).x&& newLines[i].b.y==(*newA).y))
						{
							AB=false;
							*sizeB=(*sizeB)-1;
							break;
						}
				}
			}
		
	
			if(AB==true)
			{
				for(int i=0;i<*sizeL;i++)
				{
				if((lines[i].a.x!= (*newA).x||lines[i].a.y!= (*newA).y)&&( lines[i].b.x!=b[*sizeB-1].x|| lines[i].b.y!=b[*sizeB-1].y))
					if((lines[i].a.x!=b[*sizeB-1].x||lines[i].a.y!= b[*sizeB-1].y)&&( lines[i].b.x!=(*newA).x|| lines[i].b.y!=(*newA).y))
					if(isIntersect( lines[i].a, lines[i].b, *newA, b[*sizeB-1]))
					{
						for (int j=i;j<*sizeL-1;j++)
						{
							lines[j]=lines[j+1];
						}
						*sizeL=*sizeL-1;
						i--;
					}
				}
				addLine(newLines,sizeNL,*newA,b[*sizeB-1]);
				*newB=b[*sizeB-1];
				if(*sizeB>1)
				{
					*sizeB=(*sizeB)-1;
				}
			}
		}
		else if(!BA&&!AB)
		{
			if(*sizeA>0&&a[*sizeA-1].y>(*newA).y&&a[*sizeA-1].y>(*newB).y)
			{
				*sizeA=*sizeA-1;
			}
			if(*sizeB>0&&b[*sizeB-1].y>(*newA).y&&b[*sizeB-1].y>(*newB).y)
			{
				*sizeB=*sizeB-1;
			}
		
			int end=0;
			Point tmp;
			if((*sizeA)>(*sizeB))end=(*sizeA);
			else end=(*sizeB);
			for (int i=0;i<end;i++)
			{
				if((*sizeA)-1-i>=0)
				{
					points[2]=a[(*sizeA)-1-i];
					if(isDelaunayCondition(points, a,oldSizeA+1)&&isDelaunayCondition(points, b,oldSizeB+1))
					{	tmp=a[(*sizeA)-1];
						a[(*sizeA)-1]=a[(*sizeA)-1-i];
						a[(*sizeA)-1-i]=tmp;
						break;
					}	
				}
				if((*sizeB)-1-i>=0)
				{
					points[2]=b[(*sizeB)-1-i];
					if(isDelaunayCondition(points, a,oldSizeA+1)&&isDelaunayCondition(points, b,oldSizeB+1))
					{
						tmp=b[(*sizeB)-1];
						b[(*sizeB)-1]=b[(*sizeB)-1-i];
						b[(*sizeB)-1-i]=tmp;
						break;
					}
				}
				if(i==end-1)
				{
						Point tmp;
						*wasAddOldPoints=true;
						for(int i=0;i<oldSizeA;i++)
						{
							if(a[oldSizeA-1-i].x==(*newA).x&&a[oldSizeA-1-i].y==(*newA).y)
							{
								tmp=a[oldSizeA-1-i];
								a[oldSizeA-1-i]=a[oldSizeA-1];
								a[oldSizeA-1]=tmp;
								(*sizeA)=oldSizeA+1;
								break;
							}	
						}
						for(int i=0;i<oldSizeB;i++)
						{
							if(b[oldSizeB-1-i].x==(*newB).x&&b[oldSizeB-1-i].y==(*newB).y)
							{
								tmp=b[oldSizeB-1-i];
								b[oldSizeB-1-i]=b[oldSizeB-1];
								b[oldSizeB-1]=tmp;
								(*sizeB)=oldSizeB+1;
								break;
							}
						}
				}
			}
		}
	 }
	 
 }
 
__device__ void makePolygonConvex(Point *a,Point *b,int *sizeA,int*sizeB,Line *newLines,int *sizeNL,Line *lines,int *sizeL)
 {
	Line topLine;
	Line downLine;
		
	quickSort(a,(*sizeA)-1,1,-1);
	quickSort(b, (*sizeB)-1,1,1);
	
	topLine.a=a[(*sizeA)-1];
	topLine.b=b[(*sizeB)-1];

	(*sizeA)=(*sizeA)-1;
	(*sizeB)=(*sizeB)-1;
	
	quickSort(a,(*sizeA)-1,1,1);
	quickSort(b, (*sizeB)-1,1,-1);
	
	
	downLine.a=a[0];
	downLine.b=b[0];
	
			
	if(topLine.a.y<topLine.b.y)
	{
		findBorderLine((*sizeA),(*sizeB),0,0,-1,1,a,b,&topLine.a,&topLine.b,downLine.a.y,downLine.b.y);
	}
	else if(topLine.a.y>topLine.b.y)
	{
	
		findBorderLine((*sizeB),(*sizeA),0,0,-1,-1,b,a,&topLine.b,&topLine.a,downLine.b.y,downLine.a.y);
	}
	for(int i=0;i<(*sizeA);i++)
	{
		if(a[(*sizeA)-1-i].x==topLine.a.x&&a[(*sizeA)-1-i].y==topLine.a.y)
		{
			Point tmp=a[(*sizeA)-1-i];
			a[(*sizeA)-1-i]=a[(*sizeA)];
			a[(*sizeA)]=tmp;
			break;
		}
	}
	for(int i=0;i<(*sizeB);i++)
	{
		if(b[(*sizeB)-1-i].x==topLine.b.x&&b[(*sizeB)-1-i].y==topLine.b.y)
		{
			Point  tmp=b[(*sizeB)-1-i];
			b[(*sizeB)-1-i]=b[(*sizeB)];
			b[(*sizeB)]=tmp;
			break;
		}
	}
				
	addLine(newLines,sizeNL,topLine.a,topLine.b);
	
	if(downLine.a.y>downLine.b.y)
	{
			findBorderLine(0,0,(*sizeA),(*sizeB),1,1,a,b,&downLine.a,&downLine.b,topLine.a.y,topLine.b.y);
			
	}
	else if(downLine.a.y<downLine.b.y)
	{
		findBorderLine(0,0,(*sizeB),(*sizeA),1,-1,b,a,&downLine.b,&downLine.a,topLine.b.y,topLine.a.y);
			
	}
	Point newA,newB;
	newA=topLine.a;
	newB=topLine.b;
	int oldsizeA=(*sizeA);
	int oldsizeB=(*sizeB);

	bool wasAddOldPoints=false;
	int c=0;
	while(newA.x!=downLine.a.x||newA.y!=downLine.a.y||newB.x!=downLine.b.x||newB.y!=downLine.b.y)
	{
	c++;
	buildTriangulation(a,b,sizeB,sizeA,&(newA),&(newB),lines,sizeL,newLines,sizeNL,downLine,oldsizeA,oldsizeB,&wasAddOldPoints);
	if(c>=50)break;
	}
 }


__global__ void connectTwoTriangulations(Point *points,int *counts,int sizeC,Line *lines,int *sizeL,int *countsSum,int *sizeLs,int *linesSum)
 {
	 int index=2*threadIdx.x;
	 int idx=threadIdx.x;
	 int sizeA=0;
		 int sizeB=0;
	
	
		__shared__ int cashe[SIZE/6+2];
		if(index+1<sizeC)
		{ 
			 sizeA=counts[index];
			 sizeB=counts[index+1];
		}
		 Line newLines[SIZE];
		int sizeNL=0;
		Point a[SIZE];
		Point b[SIZE]; 
		Line ls[4*SIZE];

		int end=sizeA;
		if(sizeB>sizeA&&sizeB>sizeLs[idx])
			end=sizeB;
		else if(sizeLs[idx]>sizeA&&sizeLs[idx]>sizeB)
			end=sizeLs[idx];
		for(int i=0;i<end;i++)
		{
			if(i<sizeA)
			{
				if(index-1>=0)
					a[i]=points[i+countsSum[index-1]];
				else 
					a[i]=points[i];
			}
			if(i<sizeB)
				b[i]=points[i+countsSum[index]];
			if(i<sizeLs[idx])
			{
				if(idx-1>=0)
					ls[i]=lines[i+linesSum[idx-1]];
				else 
					ls[i]=lines[i];
			}
				
		}
	
		if(index+1<sizeC)
		{
			makePolygonConvex(a,b,&sizeA,&sizeB,newLines,&sizeNL,ls,&sizeLs[idx]);//,indexes[idx],&sizeI[idx]
		}
		cashe[idx]=sizeNL+sizeLs[idx];
		__syncthreads();

		int s=0;
		for(int j=0;j<idx;j++)
		{
			s+=cashe[j];
			
		}

		for(int i=0;i<cashe[idx];i++)
		{
			if(i<sizeLs[idx])
				lines[s+i]=ls[i];
			else
				lines[s+i]=newLines[i-sizeLs[idx]];

		}
		sizeLs[idx]=cashe[idx];
		if(idx==blockDim.x-1)
		{

			*sizeL=s+sizeLs[idx];
		}

		linesSum[idx]=s+sizeLs[idx];
		__syncthreads();

		if(idx==0)
		{
			for(int i=0;i<blockDim.x;i++)
			{
				counts[i]=counts[2*i]+counts[2*i+1];
				countsSum[i]=countsSum[2*i+1];
				sizeLs[i]=sizeLs[2*i]+sizeLs[2*i+1];
				linesSum[i]=linesSum[2*i+1];
			}
			linesSum[blockDim.x]=linesSum[blockDim.x-1];
			sizeLs[blockDim.x]=0;
			counts[blockDim.x]=0;
		}
 } 
void divideIntoTriangles(int n,Point *points,Line *lines,int *sizeL,int *counts,int *sizeC,int *sizeLs);
 void splitPoints(int size1,int size2,Point *points,Line *lines,int *sizeL,int *counts,int *sizeC,int *sizeLs)
 {
	Point *p1= (Point *) malloc(size1*sizeof(Point));
	Point *p2= (Point *) malloc(size2*sizeof(Point));
	int size=size2;
	if(size1>size2)
	{
		size=size1;
	}
	 for(int i=0;i<size;++i)
	 {
			if(i<size1)
			{
				p1[i]=points[i];
			}
			if(i<size2)
			{
				p2[i]=points[i+size1];
			}
	 }
		 
	divideIntoTriangles(size1,p1,lines,sizeL,counts,sizeC,sizeLs);
	divideIntoTriangles(size2,p2,lines,sizeL,counts,sizeC,sizeLs);

	free(p1);
	free(p2);
	 
 }
  void divideIntoTriangles(int n,Point *points,Line *lines,int *sizeL,int *counts,int *sizeC,int *sizeLs)
 {
	 if(n==3)
	 { 
		counts[*sizeC]=n;
		sizeLs[*sizeC]=n;
		*sizeC=*sizeC+1;
		
		 addLine(lines,sizeL,points[0],points[1]);
		 addLine(lines,sizeL,points[1],points[2]);
		 addLine(lines,sizeL,points[2],points[0]);

	 }
	 else if(n==4)
	 {
		counts[*sizeC]=n;
		sizeLs[*sizeC]=*sizeL;
		if(isPointInCircle(points,points[3])==true)
		{
			 addLine(lines,sizeL,points[0],points[1]);
			 addLine(lines,sizeL,points[1],points[2]);
			 addLine(lines,sizeL,points[2],points[0]);
			 connectWithFourthPoint(points,lines,sizeL); 
		}
		else
		{
			for(int i=0;i<3;++i)
			{
				addLine(lines,sizeL,points[3],points[i]);
				if(i+1==3)
				{
					if(isIntersect(points[i],points[0],points[3],points[i-1])==false)
					{
						addLine(lines,sizeL,points[i],points[0]);
					}
				}
				else if(i-1<0)
				{
					if(isIntersect(points[i],points[i+1],points[3],points[2])==false)
					{
						addLine(lines,sizeL,points[i],points[i+1]);
					}
				}
				else
				{
					if(isIntersect(points[i],points[i+1],points[3],points[i-1])==false)
					{
						addLine(lines,sizeL,points[i],points[i+1]);
					}
				}
			}
		}
		 sizeLs[*sizeC]=*sizeL-sizeLs[*sizeC];
		*sizeC=*sizeC+1;

	 }
	 else if(n==8)
	 {
		splitPoints(4,4,points,lines,sizeL,counts,sizeC,sizeLs);
		 
	 }
	 else if(n<12&&n>2)
	 {
		 int t=n-3;
		 splitPoints(3,t,points,lines,sizeL,counts,sizeC,sizeLs);
		}
	 else if(n>=12)
	 {
		int t=n/2;
		splitPoints(t,n-t,points,lines,sizeL,counts,sizeC,sizeLs);
	 }
	 else if(n==2)
	 {
		counts[*sizeC]=n;
		sizeLs[*sizeC]=1;
		*sizeC=*sizeC+1;
		 addLine(lines,sizeL,points[0],points[1]);
	 }
	
 }
void mainFunction(void) 
{
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	int size=SIZE;
	float mas[SIZE][2]={{0.0,0.0},{-0.2,0.1},{0.32,-0.15},{0.17,-0.29},{-0.41,0.01},{-0.12,-0.2},{-0.48,-0.1},{0.54,0.29},{0.2,0.38},{-0.52,-0.09},
       {0.46,0.12},{0.49,-0.34},{-0.37,-0.02},{-0.09,0.05},{-0.22,0.2},{-0.55,-0.05},{0.02,-0.24},{0.3,-0.2},{-0.32,-0.04},{0.47,0.26},
       {-0.03,0.25},{-0.1,0.2},{0.16,0.4},{-0.16,0.15},{-0.55,0.15},{0.54,-0.01},{-0.28,-0.19},{0.28,-0.26},{0.48,0.02},{-0.46,0.01},
       {0.31,-0.06},{0.1,-0.21},{-0.41,-0.34},{0.0,-0.08},{0.21,-0.22},{-0.17,0.11},{-0.37,0.23},{-0.48,0.16},{0.45,0.38},{-0.39,0.32},
       {0.03,0.12},{0.29,0.24},{-0.36,0.41},{-0.54,0.19},{0.27,0.33},{-0.46,-0.33},{-0.11,-0.23},{-0.41,0.23},{-0.395,0.04},{0.34,0.15}};
	Point *points;
	Line *lines= (Line *) malloc(4*size*sizeof(Line));
	int *counts=(int *) malloc((size/3+2)*sizeof(int));
	int *countsSum=(int *) malloc((size/3+2)*sizeof(int));
	int *linesSum=(int *) malloc((size/3+2)*sizeof(int));
	int sizeL=0;
	int sizeC=0;
	points= (Point *) malloc(size*sizeof(Point));
	 int *sizeLs= (int *) malloc((size/3+2)*sizeof(int));

	  for (int i = 0; i < size; ++i)
	  {
		 points[i].x = mas[i][0];
		 points[i].y=mas[i][1];
	  }
	  quickSort(points, size-1,0,0);
	
	divideIntoTriangles(size, points,lines,&sizeL,counts,&sizeC,sizeLs);
	if(size>4)
	{
	counts[sizeC]=0;
	sizeLs[sizeC]=0;
	countsSum[0]=counts[0];
	linesSum[0]=sizeLs[0];
	for(int i=0;i<=sizeC;i++)
	{
		if(i-1>=0)
		{
			countsSum[i]=counts[i]+countsSum[i-1];
			linesSum[i]=sizeLs[i]+linesSum[i-1];
			if(i%2!=0)
				linesSum[(i-1)/2]=linesSum[i];

		}
		if(2*i+1<=sizeC)
			sizeLs[i]=sizeLs[2*i]+sizeLs[2*i+1];
	}
	
	Point *dev_points;
	int *dev_counts,*dev_sizeL,*dev_countsSum,*dev_sizeLs,*dev_linesSum;
	Line *dev_lines;
	cudaMalloc((void**)&dev_points,size*sizeof(Point));
	cudaMalloc((void**)&dev_counts,(size/3+2)*sizeof(int));
	cudaMalloc((void**)&dev_lines,4*size*sizeof(Line));
	cudaMalloc((void**)&dev_sizeL,sizeof(int));
	cudaMalloc((void**)&dev_countsSum,(size/3+2)*sizeof(int));
	cudaMalloc((void**)&dev_sizeLs,(size/3+2)*sizeof(int));
	cudaMalloc((void**)&dev_linesSum,(size/3+2)*sizeof(int));

	cudaMemcpy(dev_points,points,size*sizeof(Point),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_lines,lines,4*size*sizeof(Line),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_sizeL,&sizeL,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_counts,counts,(size/3+2)*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_sizeLs,sizeLs,(size/3+2)*sizeof(int),cudaMemcpyHostToDevice);
	
	cudaMemcpy(dev_countsSum,countsSum,(size/3+2)*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_linesSum,linesSum,(size/3+2)*sizeof(int),cudaMemcpyHostToDevice);
	int threads=0;
		while(sizeC!=1)
		{
			threads=sizeC/2;
			if(sizeC%2!=0)
				threads++;
			connectTwoTriangulations<<<1,threads>>>(dev_points, dev_counts, sizeC,dev_lines,dev_sizeL,dev_countsSum,dev_sizeLs,dev_linesSum);
			sizeC=sizeC-sizeC/2;
		}
	
cudaMemcpy(&sizeL,dev_sizeL,sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(lines,dev_lines,sizeL*sizeof(Line),cudaMemcpyDeviceToHost);
	
	cudaFree(dev_points);
	cudaFree(dev_counts);
	cudaFree(dev_lines);
	cudaFree(dev_sizeL);
	cudaFree(dev_sizeLs);
	cudaFree(dev_countsSum);
	cudaFree(dev_linesSum);
	
	}
cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	float time;
	cudaEventElapsedTime(&time,start,stop);
	printf("time: %f ms\n",time);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
printf("%d\n",sizeL);
	
for(int i=0;i<sizeL;i++)
{
	drawLine(lines[i].a,lines[i].b);
}
    glutSwapBuffers();
	free(lines);
	free(counts);
	free(points);
	free(countsSum);
	free(linesSum);
	free(sizeLs);
	
}
int main(int argc, char **argv)
{
	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(1000,800);
	glutCreateWindow("Triangulation");
	glutDisplayFunc(mainFunction);
	glutMainLoop();
	
	return 0;
}
