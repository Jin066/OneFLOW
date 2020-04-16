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
#include "UINsCom.h"
#include "Zone.h"
#include "DataBase.h"
#include "UCom.h"
#include "Com.h"
#include "INsCom.h"
#include "INsIDX.h"
#include "HXMath.h"
#include "Ctrl.h"

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
    prim.resize( nEqu );
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

void INsInvterm::CmpINsinvTerm()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	Real v2l = ONEFLOW::SQR(iinv.ul, iinv.vl, iinv.wl);
	Real v2r = ONEFLOW::SQR(iinv.ur, iinv.vr, iinv.wr);

	Real vnl = gcom.xfn * iinv.ul + gcom.yfn * iinv.vl + gcom.zfn * iinv.wl - gcom.vfn;       // V * n
	Real vnr = gcom.xfn * iinv.ur + gcom.yfn * iinv.vr + gcom.zfn * iinv.wr - gcom.vfn;

	Real rvnl = iinv.rl * vnl;   //ρ * V * n
	Real rvnr = iinv.rr * vnr;


	iinv.rf[ug.fId] = (iinv.rl + iinv.rr) * half;    //初始界面上的值（u、v、w ）
	iinv.uf[ug.fId] = (iinv.ul + iinv.ur ) * half;
	iinv.vf[ug.fId] = (iinv.vl + iinv.vr ) * half;
	iinv.wf[ug.fId] = (iinv.wl + iinv.wr ) * half;

	
	iinv.vnflow[ug.fId] = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId];  //初始界面上 V*n

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow[ug.fId] * gcom.farea; //初始界面上的质量通量

     Real clr = MAX(0, iinv.fq[ug.fId]);  //从界面左侧单元流入右侧单元的初始质量流量
     Real crl = clr - iinv.fq[ug.fId];   //从界面右侧单元流入左侧单元的初始质量流量


	iinv.aii1[ug.lc] = crl;   //该面流向左单元的流量
	iinv.aii2[ug.rc] = clr;   //该面流向右单元的流量
	
    iinv.ai1[ug.lc]+= crl;   //流入单元的流量
	iinv.ai2[ug.rc]+=clr;   //流出单元的流量
	
}

void INsInvterm::CmpINsFaceflux()
{
	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dx2 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
	Real dy2 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
	Real dz2 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

	Real de1 = DIST(dx1, dy1, dz1);
	Real de2 = DIST(dx2, dy2, dz2);
	Real de = 1.0 / (de1 + de2);

	iinv.f1[ug.fId] = de2 * de;  //左单元权重
    iinv.f2[ug.fId] = de1 * de;  //右单元权重
 
	Real Vau = iinv.f1[ug.fId] * (gcom.cvol1 / iinv.spu[ug.lc]) + iinv.f2[ug.fId] * (gcom.cvol2 / iinv.spu[ug.rc]);  //（Vj/a）
	Real Vav = iinv.f1[ug.fId] * (gcom.cvol1 / iinv.spv[ug.lc]) + iinv.f2[ug.fId] * (gcom.cvol2 / iinv.spv[ug.rc]);
	Real Vaw = iinv.f1[ug.fId] * (gcom.cvol1 / iinv.spw[ug.lc]) + iinv.f2[ug.fId] * (gcom.cvol2 / iinv.spw[ug.rc]);

    iinv.dist[ug.fId] = gcom.xfn * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + gcom.yfn * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + gcom.zfn * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);
	Real Pd1 = visQ.dqdx1[IIDX::IIP] * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + visQ.dqdy1[IIDX::IIP] * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) + visQ.dqdz1[IIDX::IIP] * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);  //压力梯度项
	Real Pd2 = visQ.dqdx2[IIDX::IIP] * ((*ug.xcc)[ug.rc] - (*ug.xfc)[ug.fId]) + visQ.dqdy2[IIDX::IIP] * ((*ug.ycc)[ug.rc] - (*ug.yfc)[ug.fId]) + visQ.dqdz2[IIDX::IIP] * ((*ug.zcc)[ug.rc] - (*ug.zfc)[ug.fId]);
	Real Pd = Pd1 + Pd2;

	
	iinv.uf[ug.fId] = (iinv.f1[ug.fId] *iinv.ul + iinv.f2[ug.fId] *iinv.ur) + (iinv.Vau*gcom.xfn / iinv.dist[ug.fId])*( Pd - (iinv.pr - iinv.pl));  //下一时刻的界面预测速度
	iinv.vf[ug.fId] = (iinv.f1[ug.fId] *iinv.vl + iinv.f2[ug.fId] *iinv.vr) + (iinv.Vav*gcom.yfn / iinv.dist[ug.fId])*(Pd - (iinv.pr - iinv.pl));
	iinv.wf[ug.fId] = (iinv.f1[ug.fId] *iinv.wl + iinv.f2[ug.fId] *iinv.wr) + (iinv.Vaw*gcom.zfn / iinv.dist[ug.fId])*(Pd - (iinv.pr - iinv.pl));

	iinv.vnflow[ug.fId] = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId];
	

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow[ug.fId] * gcom.farea;  //下一时刻界面预测通量
}

void INsInvterm::CmpINsFaceCorrectPresscoef()
{

	iinv.Vdvu[ug.fId] = iinv.f1[ug.fId] * ( gcom.cvol1/((1+1)*iinv.spu[ug.lc] - iinv.sj[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 2)*iinv.spu[ug.rc] - iinv.sj[ug.rc]));  // (Vp/dv)j，用于求面速度修正量
	iinv.Vdvv[ug.fId] = iinv.f1[ug.fId] * (gcom.cvol1 / ((1+1)*iinv.spv[ug.lc] - iinv.sj[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 2)*iinv.spv[ug.rc] - iinv.sj[ug.rc]));
	iinv.Vdvw[ug.fId] = iinv.f1[ug.fId] * (gcom.cvol1 / ((1+1)*iinv.spw[ug.lc] - iinv.sj[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 2)*iinv.spw[ug.rc] - iinv.sj[ug.rc]));
	

	iinv.aju[ug.fId] = iinv.rm * iinv.Vdvu[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId]; //压力修正方程中，矩阵中非零系数
	iinv.ajv[ug.fId] = iinv.rm * iinv.Vdvv[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId];
	iinv.ajw[ug.fId] = iinv.rm * iinv.Vdvw[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId];

	iinv.ajp[ug.fId] = iinv.aju[ug.fId] * gcom.xfn + iinv.ajv[ug.fId] * gcom.yfn + iinv.ajw[ug.fId] * gcom.zfn;
}


  














EndNameSpace