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
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "codetypes.h"
#include "tioga.h"
#include "globals.h"
#include <string.h>
//
// All the interfaces that are 
// accessible to third party f90 and C
// flow solvers
//
//
// Jay Sitaraman
// 02/24/2014
//
extern "C" {

  void TIOGA_INIT_F90(int *scomm)
  {
    int id_proc,nprocs;  //某一个并行执行的进程的标识，所有参加计算的进程的个数
    MPI_Comm tcomm;  //通讯器
    //tcomm=(MPI_Comm) (*scomm);
    tcomm=MPI_Comm_f2c(*scomm); 
    //
    tg=new tioga[1];
    //
    //MPI_Comm_rank(MPI_COMM_WORLD,&id_proc);
    //MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(tcomm,&id_proc); //得到当前正在运行的进程的标识号
    MPI_Comm_size(tcomm,&nprocs);  //得到所有参加运算的进程的个数放在nprocs中
    //
    tg->setCommunicator(tcomm,id_proc,nprocs);  //设置id_proc的并行通讯
    for(int i=0;i<MAXBLOCKS;i++)
     {
      idata[i].nc=NULL;  //网格块中各不同单元的数量（三棱镜、金字塔、六角等）=0
      idata[i].nv=NULL;  //每种类型单元的顶点数=0
      idata[i].vconn=NULL; //每种单元的连通性=0
     }
  }

  void TIOGA_INIT(MPI_Comm tcomm)
  {
    int id_proc,nprocs;
    //MPI_Comm tcomm;
    //tcomm=(MPI_Comm) (*scomm);
    //tcomm=MPI_Comm_f2c(*scomm);
    //
    tg=new tioga[1];
    //
    //MPI_Comm_rank(MPI_COMM_WORLD,&id_proc);
    //MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(tcomm,&id_proc);
    MPI_Comm_size(tcomm,&nprocs);
    //
    tg->setCommunicator(tcomm,id_proc,nprocs);
    for (int i=0;i<MAXBLOCKS;i++)
     {
      idata[i].nc=NULL;
      idata[i].nv=NULL;
      idata[i].vconn=NULL;
     }
  }
  
  void TIOGA_REGISTERGRID_DATA(int *btag,int *nnodes,double *xyz,int *ibl,int *nwbc, int *nobc,int *wbcnode,
                               int *obcnode,int *ntypes,...)   //记录分区的网格数据
   {                 //*btag:此分区的网格标记（标量）
	                //*nnodes:此分区中的节点数（标量）
	               //*xyz:节点坐标（1维实数数组，大小=3*节点个数）
	              //*ibl:节点类型标志（整数数组大小=节点个数）
	             //*nwbc:重叠边界节点数
                //*nobc:壁面边界节点数
	           //*wbcnode:壁面边界节点索引
	           //*obcnode:重叠边界节点索引
	           //*ntypes:单元类型的个数（标量）
	           //......

    va_list arguments; //参数表（实参）
    int i;
    int iblk=0;
    va_start(arguments, ntypes);

    if(idata[iblk].nv) TIOGA_FREE(idata[iblk].nv); //nv，nc，vconn参数初始化设置为零
    if(idata[iblk].nc) TIOGA_FREE(idata[iblk].nc);
    if(idata[iblk].vconn) TIOGA_FREE(idata[iblk].vconn);
    idata[iblk].nv=(int *) malloc(sizeof(int)*(*ntypes)); //1维数列，大小=单元类型个数   
    idata[iblk].nc=(int *) malloc(sizeof(int)*(*ntypes));
    idata[iblk].vconn=(int **)malloc(sizeof(int *)*(*ntypes));
    for(i=0;i<*ntypes;i++)
     {
      idata[iblk].nv[i]=*(va_arg(arguments, int *)); //iblk网格块中单元类型i的单元顶点数
      idata[iblk].nc[i]=*(va_arg(arguments, int *)); //iblk网格块中单元类型i的单元个数
      idata[iblk].vconn[i]=va_arg(arguments, int *); //iblk网格块中单元类型i的连通性
     }
    tg->registerGridData(*btag,*nnodes,xyz,ibl,*nwbc,*nobc,wbcnode,obcnode,*ntypes,idata[iblk].nv,idata[iblk].nc,idata[iblk].vconn);
	//idata[iblk].nv：在此分区中iblk网格的各单元类型中每个单元的顶点数
	//idata[iblk].nc：在此分区中iblk网格的各单元类型的单元个数
	//idata[iblk].vconn：在此分区中iblk网格的各单元类型的连通性
  }


  void TIOGA_REGISTERGRID_DATA_MB(int *bid, int *btag,int *nnodes,double *xyz,int *ibl,int *nwbc, int *nobc,int *wbcnode, 
			       int *obcnode,int *ntypes,...)   //记录网格块网格数据
  {                 //*bid：网格块id
    va_list arguments;
    int i;
    int iblk=*bid-BASE;   //网格块索引

    va_start(arguments, ntypes);

    if(idata[iblk].nv) TIOGA_FREE(idata[iblk].nv);
    if(idata[iblk].nc) TIOGA_FREE(idata[iblk].nc);
    if(idata[iblk].vconn) TIOGA_FREE(idata[iblk].vconn);
    idata[iblk].nv=(int *) malloc(sizeof(int)*(*ntypes));    
    idata[iblk].nc=(int *) malloc(sizeof(int)*(*ntypes));
    idata[iblk].vconn=(int **)malloc(sizeof(int *)*(*ntypes));
    for(i=0;i<*ntypes;i++)
     {
      idata[iblk].nv[i]=*(va_arg(arguments, int *));
      idata[iblk].nc[i]=*(va_arg(arguments, int *));
      idata[iblk].vconn[i]=va_arg(arguments, int *);
     }
    tg->registerGridData(*btag,*nnodes,xyz,ibl,*nwbc,*nobc,wbcnode,obcnode,*ntypes,idata[iblk].nv,idata[iblk].nc,idata[iblk].vconn);
  }

  void TIOGA_REGISTER_AMR_GLOBAL_DATA(int *nf, int *qstride, double *qnodein,
				      int *idata,double *rdata,
				      int *ngridsin,int *qnodesize) //网格自适应（非必需）
  {
    tg->register_amr_global_data(*nf,*qstride,qnodein,idata,rdata,*ngridsin,*qnodesize);
  }


  void TIOGA_REGISTER_AMR_PATCH_COUNT(int *npatches)//（非必需）
  {
    tg->set_amr_patch_count(*npatches);
  }

  void TIOGA_REGISTER_AMR_LOCAL_DATA(int *ipatch,int *global_id,int *iblank,double *q) //（非必需）
	  //*ipatch：辅助网格的补丁（分区）
	  //*global_id：全局id
	  //*iblank：每个网格节点的Iblan k值
	  //*q:变量
  {
    tg->register_amr_local_data(*ipatch,*global_id,iblank,q);
  }

  void TIOGA_PREPROCESS_GRIDS(void) //网格预处理
  {
    tg->profile();
  }

  void TIOGA_PERFORMCONNECTIVITY(void)   //连通性
  {
	  tg->performConnectivity();
  }

  void TIOGA_PERFORMCONNECTIVITY_HIGHORDER(void)  //高阶连通性
  {
    tg->performConnectivityHighOrder();
  }

  void TIOGA_PERFORMCONNECTIVITY_AMR(void)
  {
   tg->performConnectivityAMR();
  }

  void TIOGA_REGISTERSOLUTION(int *bid,double *q)
  {      // bid：网格块id
	     //q：变量（比如密度、速度、压力等）
    tg->registerSolution(*bid,q);
  }

  void TIOGA_DATAUPDATE_MB(int *nvar,char *itype)
  {      //nvar：每个节点的场变量个数, 如果场是不同的数组，可以为每个场调用多次
	    //字符串结构类型，row(行)或column(列)
    int interptype;
    if (strstr(itype,"row")) 
      {
	interptype=0;
      }
    else if (strstr(itype,"column")) 
      {
	interptype=1;
      }
    else
      {
	printf("#tiogaInterface.C:dataupdate_:unknown data orientation\n");
	return;
      }

    //tg->dataUpdate(*nvar,interptype);

    if (tg->ihighGlobal==0) 
    {
      if (tg->iamrGlobal==0) 
      	{
     	    tg->dataUpdate(*nvar,interptype);
        }
     	else
     	  {
     	    tg->dataUpdate_AMR(*nvar,interptype);
     	  }
    }
    else
    {
     if (tg->iamrGlobal==0) 
       {
     	    tg->dataUpdate(*nvar,interptype,1);
       }
      else
      {
        printf("Data udpate between high-order near-body and AMR cartesian Not implemented yet\n");
      }
    }
  }

  void TIOGA_DATAUPDATE(double *q,int *nvar,char *itype)
  {
    int interptype;
    int bid=0;
    tg->registerSolution(bid,q);
    TIOGA_DATAUPDATE_MB(nvar,itype);
  }

  void TIOGA_WRITEOUTPUTFILES(int *nvar,char *itype) //编写输出文件
  {
    int interptype;
    if (strstr(itype,"row")) 
      {
	interptype=0;
      }
    else if (strstr(itype,"column")) 
      {
	interptype=1;
      }
    else
      {
	printf("#tiogaInterface.C:dataupdate_:unknown data orientation\n");
	return;
      }
    tg->writeData(*nvar,interptype);
  }    
  void TIOGA_GETDONORCOUNT(int *btag,int *dcount,int *fcount)
  {
    tg->getDonorCount(*btag,dcount,fcount); //计算捐献者人数
  }
  void TIOGA_GETDONORINFO(int *btag,int *receptors,int *indices,double *frac,int *dcount)
  {
    tg->getDonorInfo(*btag,receptors,indices,frac,dcount); //得到贡献单元信息
  }

  void TIOGA_SETSYMMETRY(int *isym)
  {
    tg->setSymmetry(*isym);
  }

  void TIOGA_SETRESOLUTIONS(double *nres,double *cres)
  {
    tg->setResolutions(nres,cres);
  }

  void TIOGA_SETRESOLUTIONS_MULTI(int *btag,double *nres,double *cres)
  {
    tg->setResolutions(*btag,nres,cres);
  }
  
  void TIOGA_SETCELLIBLANK(int *iblank_cell)
  {
    tg->set_cell_iblank(iblank_cell);
  }

  void TIOGA_SET_HIGHORDER_CALLBACK(void (*f1)(int*, int*),
				    void (*f2)(int *,int *,double *),
				    void (*f3)(int *,double *,int *,double *),
				    void (*f4)(int *,double *,int *,int *,double *,double *,int *),
				     void (*f5)(int *,int *,double *,int *,int *,double *))
  {
    tg->setcallback(f1,f2,f3,f4,f5);
    //    get_nodes_per_cell=f1;
    //get_receptor_nodes=f2;
    //donor_inclusion_test=f3;
    //donor_frac=f4;
    //convert_to_modal=f5;
  }
  
  void TIOGA_SET_P4EST(void)
  {
    tg->set_p4est();
  }
  void TIOGA_SET_AMR_CALLBACK(void (*f1)(int *,double *,int *,double *))
  {
    tg->set_amr_callback(f1);
  }
  void TIOGA_SET_P4EST_SEARCH_CALLBACK(void (*f1)(double *xsearch,int *process_id,int *cell_id,int *npts),
					void (*f2)(int *pid,int *iflag))
  {
    tg->setp4estcallback(f1,f2);
  //jayfixme  tg->set_p4est_search_callback(f1);
  }  

  void TIOGA_REDUCE_FRINGES(void)
  {
    tg->reduce_fringes();
  }

  void TIOGA_SETNFRINGE(int *nfringe)
  {
    tg->setNfringe(nfringe);
  }

  void TIOGA_SETMEXCLUDE(int *mexclude)
  {
   tg->setMexclude(mexclude);
  }

  void TIOGA_DELETE(void)
   {
    delete [] tg;
    for(int i=0;i<MAXBLOCKS;i++)
     {
       if (idata[i].nc) TIOGA_FREE(idata[i].nc);
       if (idata[i].nv) TIOGA_FREE(idata[i].nv);
       if (idata[i].vconn) TIOGA_FREE(idata[i].vconn);
     }
   }

}
