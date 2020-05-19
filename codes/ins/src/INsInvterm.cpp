/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

//#include "UINsCorrectSpeed.h"
#include "INsInvterm.h"
#include "INsVisterm.h"
#include "Iteration.h"
#include "UINsCom.h"
#include "Zone.h"
#include "DataBase.h"
#include "UCom.h"
#include "Com.h"
#include "INsCom.h"
#include "INsIDX.h"
#include "HXMath.h"
#include "Ctrl.h"
#include "Boundary.h"
#include "BcRecord.h"

BeginNameSpace( ONEFLOW )

INsInv iinv;



INsInv::INsInv()
{
    ;
}

INsInv::~INsInv()
{
    ;
}

void INsInv::Init()
{
    int nEqu = inscom.nEqu;
    prim.resize(nEqu);
    prim1.resize( nEqu );
    prim2.resize( nEqu );

    q.resize( nEqu );
    q1.resize( nEqu );
    q2.resize( nEqu );

    dq.resize( nEqu );

    flux.resize( nEqu );
    flux1.resize( nEqu );
    flux2.resize( nEqu );

}

INsInvterm::INsInvterm()
{
    ;
}

INsInvterm::~INsInvterm()
{
    ;
}

void INsInvterm::Solve()
{
}

void INsInvterm::CmpINsinvFlux()
{

		INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);

		INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);


		iinv.rf[ug.fId] = (iinv.rl + iinv.rr) * half;    //初始界面上的值（u、v、w ）

		iinv.uf[ug.fId] = (iinv.ul + iinv.ur) * half;

		iinv.vf[ug.fId] = (iinv.vl + iinv.vr) * half;

		iinv.wf[ug.fId] = (iinv.wl + iinv.wr) * half;

		iinv.vnflow[ug.fId] = (*ug.xfn)[ug.fId] * iinv.uf[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vf[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wf[ug.fId] - gcom.vfn;  //初始界面上 V*n

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow[ug.fId] * gcom.farea; //初始界面上的质量通量

}

void INsInvterm::CmpINsBcinvFlux()
{

	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);

	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	int bcType = ug.bcRecord->bcType[ug.fId];


	iinv.rf[ug.fId] = (iinv.rl + iinv.rr) * half;    //初始界面上的值（u、v、w ）

	iinv.uf[ug.fId] = (iinv.ul + iinv.ur) * half;

	iinv.vf[ug.fId] = (iinv.vl + iinv.vr) * half;

	iinv.wf[ug.fId] = (iinv.wl + iinv.wr) * half;

	iinv.vnflow[ug.fId] = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;  //初始界面上 V*n

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow[ug.fId] * gcom.farea; //初始界面上的质量通量

}

void INsInvterm::CmpINsinvTerm()
{ 
		Real clr = MAX(0, iinv.fq[ug.fId]);  //从界面左侧单元流入右侧单元的初始质量流量

		Real crl = clr - iinv.fq[ug.fId];   //从界面右侧单元流入左侧单元的初始质量流量

		iinv.aii1[ug.fId] = crl;   //该面流向左单元的流量

		iinv.aii2[ug.fId] = clr;   //该面流向右单元的流量

		iinv.ai1[ug.lc] += crl;   //流入单元的流量

		iinv.ai2[ug.rc] += clr;   //流出单元的流量
}

void INsInvterm::CmpINsBcinvTerm()
{
	
		Real clr = MAX(0, iinv.fq[ug.fId]);  //从界面左侧单元流入右侧单元的初始质量流量
		Real crl = clr - iinv.fq[ug.fId];   //从界面右侧单元流入左侧单元的初始质量流量


		iinv.aii1[ug.fId] = crl;   //该面流向左单元的流量
		iinv.aii2[ug.fId] = clr;   //该面流向右单元的流量

		iinv.ai1[ug.lc] += crl;   //流入单元的流量
		iinv.ai2[ug.rc] += clr;   //流出单元的流量
}

void INsInvterm::CmpINsFaceflux()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	iinv.Vau[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * (*ug.xfn)[ug.fId] / (iinv.spu[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * (*ug.xfn)[ug.fId] / (iinv.spu[ug.rc])); //Df*n
	iinv.Vav[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * (*ug.yfn)[ug.fId] / (iinv.spv[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * (*ug.yfn)[ug.fId] / (iinv.spv[ug.rc]));
	iinv.Vaw[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * (*ug.zfn)[ug.fId] / (iinv.spw[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * (*ug.zfn)[ug.fId] / (iinv.spw[ug.rc]));

	iinv.dist[ug.fId] = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	iinv.dsrl[ug.fId] = sqrt(((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc])*((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc])*((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc])*((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]));

	iinv.elrn[ug.fId] = iinv.dist[ug.fId] / iinv.dsrl[ug.fId];

	iinv.Deun[ug.fId] = iinv.Vau[ug.fId] / iinv.elrn[ug.fId];   //Df*n/e*n
	iinv.Devn[ug.fId] = iinv.Vav[ug.fId] / iinv.elrn[ug.fId];
	iinv.Dewn[ug.fId] = iinv.Vaw[ug.fId] / iinv.elrn[ug.fId];


	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dx2 = (*ug.xcc)[ug.rc] - (*ug.xfc)[ug.fId];
	Real dy2 = (*ug.ycc)[ug.rc] - (*ug.yfc)[ug.fId];
	Real dz2 = (*ug.zcc)[ug.rc] - (*ug.zfc)[ug.fId];

	iinv.dlf[ug.fId] = sqrt(dx1*dx1 + dy1 * dy1 + dz1 * dz1);
	iinv.dfr[ug.fId] = sqrt(dx2*dx2 + dy2 * dy2 + dz2 * dz2);

	iinv.Bpe[ug.fId] = (*uinsf.dqdx)[IIDX::IIP][ug.lc] * (dx1/ iinv.dlf[ug.fId]) + (*uinsf.dqdy)[IIDX::IIP][ug.lc] * (dy1/ iinv.dlf[ug.fId]) + (*uinsf.dqdz)[IIDX::IIP][ug.lc] * (dz1/ iinv.dlf[ug.fId])+
					   (*uinsf.dqdx)[IIDX::IIP][ug.rc] * (dx2/ iinv.dfr[ug.fId]) + (*uinsf.dqdy)[IIDX::IIP][ug.rc] * (dy2/ iinv.dfr[ug.fId]) + (*uinsf.dqdz)[IIDX::IIP][ug.rc] * (dz2/ iinv.dfr[ug.fId])-
					   (iinv.pr - iinv.pl);


	iinv.uf[ug.fId] = (iinv.f1[ug.fId] * iinv.ul + iinv.f2[ug.fId] * iinv.ur);//iinv.Deun[ug.fId]* iinv.Bpe[ug.fId];  //下一时刻的界面预测速度
		iinv.vf[ug.fId] = (iinv.f1[ug.fId] * iinv.vl + iinv.f2[ug.fId] * iinv.vr);// iinv.Devn[ug.fId] * iinv.Bpe[ug.fId];
		iinv.wf[ug.fId] = (iinv.f1[ug.fId] * iinv.wl + iinv.f2[ug.fId] * iinv.wr);//iinv.Dewn[ug.fId] * iinv.Bpe[ug.fId];

	iinv.vnflow[ug.fId] = (*ug.xfn)[ug.fId] * iinv.uf[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vf[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wf[ug.fId];

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow[ug.fId] * (*ug.farea)[ug.fId];  //下一时刻界面预测通量

}

void INsInvterm::CmpINsBcFaceflux()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	iinv.dist[ug.fId] = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);
	
	int bcType = ug.bcRecord->bcType[ug.fId];

	iinv.uf[ug.fId] = (iinv.ul + iinv.ur) / half;
	iinv.vf[ug.fId] = (iinv.vl + iinv.vr) / half;
	iinv.wf[ug.fId] = (iinv.wl + iinv.wr) / half;

	iinv.vnflow[ug.fId] = (*ug.xfn)[ug.fId] * iinv.uf[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vf[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wf[ug.fId];

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow[ug.fId] * (*ug.farea)[ug.fId];  //下一时刻界面预测通量



	//if (bcType == BC::SOLID_SURFACE)
	//{
	//	iinv.uf[ug.fId] = 0;

	//	iinv.vf[ug.fId] = 0;

	//	iinv.wf[ug.fId] = 0;

	//	iinv.fq[ug.fId] = 0;
	//}

}

void INsInvterm::CmpINsFaceCorrectPresscoef()
{

	iinv.Vdvu[ug.fId] = -iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc]*(*ug.xfn)[ug.fId] /((1+1)*iinv.spu[ug.lc] - iinv.sju[ug.lc])) - iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * (*ug.xfn)[ug.fId] / ((1+1)*iinv.spu[ug.rc] - iinv.sju[ug.rc]));  // -Mf*n，用于求面速度修正量
	iinv.Vdvv[ug.fId] = -iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] *(*ug.yfn)[ug.fId] / ((1+1)*iinv.spv[ug.lc] - iinv.sjv[ug.lc])) - iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc]* (*ug.yfn)[ug.fId] / ((1+1)*iinv.spv[ug.rc] - iinv.sjv[ug.rc]));
	iinv.Vdvw[ug.fId] = -iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * (*ug.zfn)[ug.fId] / ((1+1)*iinv.spw[ug.lc] - iinv.sjw[ug.lc])) - iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * (*ug.zfn)[ug.fId] / ((1+1)*iinv.spw[ug.rc] - iinv.sjw[ug.rc]));
	
	//iinv.aju[ug.fId] = iinv.rf[ug.fId] * iinv.Vdvu[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / (iinv.dist[ug.fId]); 
	//iinv.ajv[ug.fId] = iinv.rf[ug.fId] * iinv.Vdvv[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / (iinv.dist[ug.fId]);
	//iinv.ajw[ug.fId] = iinv.rf[ug.fId] * iinv.Vdvw[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / (iinv.dist[ug.fId]);


	//iinv.ajp[ug.fId] = iinv.aju[ug.fId] * gcom.xfn + iinv.ajv[ug.fId] * gcom.yfn + iinv.ajw[ug.fId] * gcom.zfn;  //压力修正方程中，矩阵中非零系数
	iinv.ajp[ug.fId] = iinv.rf[ug.fId] * (iinv.Vdvu[ug.fId]* (*ug.xfn)[ug.fId] + iinv.Vdvv[ug.fId]* (*ug.yfn)[ug.fId] + iinv.Vdvw[ug.fId]* (*ug.zfn)[ug.fId]) * (*ug.farea)[ug.fId] / iinv.elrn[ug.fId]-half* iinv.fq[ug.fId]; //压力修正方程中，矩阵中非零系数
}

void INsInvterm::CmpINsBcFaceCorrectPresscoef()
{

	iinv.Vdvu[ug.fId] = 0;  // (Vp/dv)j，用于求面速度修正量
	iinv.Vdvv[ug.fId] = 0;
	iinv.Vdvw[ug.fId] = 0;

	iinv.ajp[ug.fId] = 0;
}


  














EndNameSpace