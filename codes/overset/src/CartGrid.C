//
// This file is part of the Tioga software library
//
// Tioga  is a tool for overset grid assembly on parallel distributed systems
// Copyright (C) 2015 Jay Sitaraman
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# include "codetypes.h"
# include "CartGrid.h"

void CartGrid::registerData(int nfin,int qstridein,double *qnodein,int *idata,
			    double *rdata,int ngridsin,int qnodesize)
{  //nfin:辅助网格文件个数
	//qstridein：三个坐标方向上辅助网格子块个数（Nx,Ny,Nz）
	//qnodein：辅助网格节点(1维数组，大小=节点个数)
	//ngridsin：辅助网格块数
	//qnodesize：辅助网格节点数
  int i,i3,i6,iloc,n;
  FILE *fp;  //file name
  ngrids = ngridsin; //辅助网格块个数
  global_id=(int *) malloc(sizeof(int)*ngrids); //辅助网格块全局id(1维数组，大小=辅助网格块个数)
  level_num=(int *) malloc(sizeof(int)*ngrids); //辅助网格块级数
  proc_id=(int *) malloc(sizeof(int)*ngrids); //进程id
  ilo=(int *) malloc(sizeof(int)*3*ngrids); //最小坐标
  ihi=(int *) malloc(sizeof(int)*3*ngrids); //最大坐标
  xlo=(double *) malloc(sizeof(double)*3*ngrids); //辅助网格块（全局）最小坐标xmin
  dx=(double *) malloc(sizeof(double)*3*ngrids);  //
  porder=(int *) malloc(sizeof(int)*ngrids); //辅助网格块顺序
  local_id=(int *)malloc(sizeof(int)*ngrids); //辅助网格块本地id（网格壁面边界被多个辅助网格块包围）
  qnode=(double *)malloc(sizeof(double)*qnodesize); //辅助网格节点（大小=节点个数）
  dims=(int *)malloc(sizeof(dims)*3*ngrids); //三坐标方向上子块个数(Nx,Ny,Nz)
  for(i=0;i<qnodesize;i++)  { qnode[i]=qnodein[i];} //辅助网格节点索引
  //                            if (myid==0) printf("qnode[%d]= %f\n",i,qnode[i]);}
  nf=nfin; //文件个数
  qstride=qstridein; //
  if (myid==0) fp=fopen("cartGrid.dat","w");
  for(i=0;i<ngrids;i++)
    {
      i3=3*i;
      i6=2*i3;
      iloc=11*i;

      global_id[i]=idata[iloc];    //全局id
      level_num[i]=idata[iloc+1]; //级数
      proc_id[i]=idata[iloc+2];  //处理器id
      porder[i]=idata[iloc+3];  //顺序
      local_id[i]=idata[iloc+4]; //本地id
      for(n=0;n<3;n++)
	{
	  ilo[i3+n]=idata[iloc+5+n]; //最小坐标
	  ihi[i3+n]=idata[iloc+8+n]; //最大坐标
	  dims[i3+n]=ihi[i3+n]-ilo[i3+n]+1; //三个坐标方向上的坐标差极值
	}
      xlo[i3]=rdata[i6]; //辅助网格三个坐标方向的最小值
      xlo[i3+1]=rdata[i6+1];
      xlo[i3+2]=rdata[i6+2];
      dx[i3]=rdata[i6+3]; //(xhi-xlo)/N
      dx[i3+1]=rdata[i6+4];
      dx[i3+2]=rdata[i6+5];
      if (myid==0) 
        fprintf(fp,"%d %d %d %d %d %f %f %f\n",global_id[i],level_num[i],proc_id[i],
                                   porder[i],local_id[i],dx[i3],dx[i3+1],
                                   dx[i3+2]);
    }
   if (myid==0) fclose(fp);
};

//
// Bare bone preprocessor now
// willl add more support data structures
// to promote efficient search once the concept
// works
//
void CartGrid::preprocess(void)
{
  int i,n;
  //
  // find the global minimum coord location 找到全局最小坐标的位置
  //
  xlosup[0]=xlosup[1]=xlosup[2]=BIGVALUE;
  maxlevel=-1;
  for (i=0;i<ngrids;i++)
    {
      for(n=0;n<3;n++)
	xlosup[n]=TIOGA_Min(xlosup[n],xlo[3*i+n]);
      maxlevel=TIOGA_Max(maxlevel,level_num[i]);
    }
    maxlevel++;
  lcount=(int *)malloc(sizeof(int)*maxlevel);
  dxlvl=(double *)malloc(sizeof(double)*3*maxlevel);
  for(i=0;i<maxlevel;i++) lcount[i]=0;
  for(i=0;i<ngrids;i++)
    {
      lcount[level_num[i]]++;
      for(n=0;n<3;n++)
	dxlvl[3*level_num[i]+n]=dx[3*i+n];
    }
}
//
// Basic search routine now
// will improve efficiency once it works
//
void CartGrid::search(double *x,int *donorid,int npts)
{
  int i,j,k,l,n,il[3];
  bool flag;
  int dcount;
  dcount=0;
  for(i=0;i<npts;i++)
    {
      flag=0;
      donorid[i]=-1;
      for(l=maxlevel-1;l>=0 && flag==0;l--)
	{
	  for(n=0;n<3;n++)
	    il[n]=floor((x[3*i+n]-xlosup[n])/dxlvl[3*l+n]);
	  for(j=0;j<ngrids && flag==0;j++)
	    {
	      if (level_num[j]==l) 
		{
		  flag=1;
		  // for(n=0;n<3;n++) flag=flag && (x[3*i+n] >=xlo[3*j+n]);
		  // for(n=0;n<3;n++) flag=flag && (x[3*i+n] <=xlo[3*j+n]+
		  //   			 dx[3*j+n]*(dims[3*j+n]));
		  //for(n=0;n<3;n++) flag = flag && (il[n] >=ilo[3*j+n]);
		  //for(n=0;n<3;n++) flag = flag && (il[n] <=ihi[3*j+n]);
          for(n=0;n<3;n++) flag=flag && ((x[3*i+n]-xlo[3*j+n]) > -TOL);
          for(n=0;n<3;n++) flag=flag && ((x[3*i+n]- (xlo[3*j+n]+
                                                     dx[3*j+n]*(dims[3*j+n]))) < TOL);
		  if (flag) { 
		    dcount++; 
		    donorid[i]=j; 
		}
	    }
	}
    }
  }
 //printf("CartGrid::search Processor %d located %d of %d points\n",myid,dcount,npts);
}
