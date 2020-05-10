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

#ifndef CODETYPES_H
#define CODETYPES_H

#define MPICH_SKIP_MPICXX
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mpi.h"
/*====================================================================*/
/*  Floating point definition                                         */
/*====================================================================*/
# define REAL double
/*====================================================================*/
/*  Base for indexing (0 or 1) */
/*====================================================================*/
# define BASE 1

/*====================================================================*/
/*  Define arithmetic constants                                       */
/*====================================================================*/
// #define ZERO               0.0e+00
// #define ONE                1.0e+00
// #define TWO                2.0e+00
// #define THREE              3.0e+00
// #define FOUR               4.0e+00
// #define HALF               0.5e+00
// #define THIRD              0.333333333e+00
// #define FIFTH              0.2
// #define PI                 3.1415926535897932e+00
// #define RAD2DEG            (180.0/PI)
// #define DEG2RAD            (PI/180.0)
#define BIGVALUE           1.0e+15
#define BIGINT             2147483647
#define TOL                1.0e-10
#define HOLEMAPSIZE        192
// #define NFRINGE            3
// #define NVAR               6
/*==================================================================*/
/* inline debugging tools                                             */
/*==================================================================*/
# define TRACEI(x)  printf("#tioga:\t"#x" =%d\n",x);
# define TRACED(x)  printf("#tioga:\t"#x" =%.16e\n",x);
# define TIOGA_MIN(x,y)  (x) < (y) ? (x) : (y)
# define TIOGA_MAX(x,y)  (x) > (y) ? (x) : (y)
# define TIOGA_FREE(a1)  {free(a1);a1=NULL;}
// # define debug(x,y)  printf("#tioga:\t"#x"=%d,"#y"=%d\n",x,y);
// # define stdwrite(x) if (myid==0) printf("#tioga:\t"#x"\n");
// # define dstr(x) printf("#tioga:\t"#x"\n");
// # define ditch(x,y) {dstr(x); TRACEI(y); MPI_Abort(MPI_COMM_WORLD,ierr);}
/*====================================================================*/
/*  Numerical Tools                                                   */
/*====================================================================*/
// #define Sign(a1,a2) (((a2) < ZERO)? - fabs(a1): fabs(a1))
#define TIOGA_Max(a1,a2) (((a1) >= (a2))? (a1): (a2))
#define TIOGA_Min(a1,a2) (((a1) <= (a2))? (a1): (a2))
// #define Abs(aa) (((aa) >= 0)? (aa): -(aa))
// #define Round(aa) (int) ((fabs((aa) - floor(aa)) >= HALF)? ceil(aa): floor(aa))
// #define swap(a,b) { a=a+b;b=a-b;a=a-b;}
/*********************************************************************/
/* Code specific types */
/*********************************************************************/
typedef struct HOLEMAP  //构造结构辅助网格
{
  int existWall;
  int nx[3];//三个整数（Nx，Ny,Nz），辅助网格中三个方向的子块个数
  int *samLocal; //结构辅助网格
  int *sam;
  double extents[6];//6个实数，辅助网格最大坐标和最小坐标
} HOLEMAP;


typedef struct OBB //构造定向边界箱
{
  double xc[3];  //中心坐标
  double dxc[3]; //中心距（边框最小值-边框中心坐标）
  double vec[3][3]; //边框三个边的单位向量

  int comm_idx;    //* Index in comm map for this OBB                       */这个obb的comm映射中的索引
  int iblk_local;  //* Index of this mesh block                             */网格块的索引
  int iblk_remote;  //* Index of the remote mesh block (intersecting  pair)  */远程网格块的索引（相交对）
  int tag_remote;   //标签
  int send_tag;  //发送
  int recv_tag; //接收
}OBB;
typedef struct DONORLIST
{
  int donorData[4];
  double donorRes;  //贡献者分辨率
  double receptorRes; //接收者分辨率
  int cancel;
  struct DONORLIST *next;
} DONORLIST;

typedef struct PACKET
{
  int nints;
  int nreals;
  int *intData;
  REAL *realData;
} PACKET;

typedef struct INTERPLIST
{
  int cancel;
  int nweights;
  int receptorInfo[3];
  double xtmp[3];
  int *inode;
  double *weights;
} INTERPLIST;

typedef struct INTERPLIST2
{
  int cancel;
  int nweights;
  int receptorInfo[3];
  int *inode;
  double *weights;
  struct INTERPLIST2 *next;
} INTERPLIST2;

typedef struct INTEGERLIST
{
  int inode;
  struct INTEGERLIST *next;
} INTEGERLIST;

typedef struct INTEGERLIST2
{
  int intDataSize,realDataSize;
  int *intData;
  double *realData;
  struct INTEGERLIST2 *next;
}INTEGERLIST2; 

#endif /* CODETYPES_H */
