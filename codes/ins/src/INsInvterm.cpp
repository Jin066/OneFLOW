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

	Real rvnl = iinv.rl * vnl;   //�� * V * n
	Real rvnr = iinv.rr * vnr;

	Real ratio = sqrt(iinv.rr / iinv.rl);
	Real coef = 1.0 / (1.0 + ratio);

	iinv.rm = (iinv.rl + iinv.rr ) * half;    //��ʼ�����ϵ�ֵ��u��v��w ��
	iinv.um = (iinv.ul + iinv.ur ) * half;
	iinv.vm = (iinv.vl + iinv.vr ) * half;
	iinv.wm = (iinv.wl + iinv.wr ) * half;
	iinv.hm = (iinv.hl + iinv.hr * ratio) * coef;
	iinv.pm = (iinv.pl + iinv.pr * ratio) * coef;
	iinv.gama = (iinv.gama1 + iinv.gama2 * ratio) * coef;

	

	//Real v2 = ONEFLOW::SQR(iinv.um, iinv.vm, iinv.wm);


	iinv.vnflow = gcom.xfn * iinv.um + gcom.yfn * iinv.vm + gcom.zfn * iinv.wm;  //��ʼ������ V*n(�Ķ�)

	//iinv.fq0[ug.fId] = iinv.rl * iinv.vnflow * gcom.farea; //��ʼ�����ϵ�����ͨ��

	iinv.fq0 = iinv.rl * iinv.vnflow * gcom.farea;

	//Real gamm1 = iinv.gama - one;

	//Real c2 = gamm1 * (iinv.hm - half * v2);
	//iinv.cm = sqrt(ABS(c2));

	//iinv.aeig1 = ABS(iinv.vnrel);
	//iinv.aeig2 = ABS(iinv.vnrel + iinv.cm);
	//iinv.aeig3 = ABS(iinv.vnrel - iinv.cm);


	iinv.clr = MAX(0, iinv.fq0);  //�ӽ�����൥Ԫ�����Ҳ൥Ԫ�ĳ�ʼ��������
	iinv.crl = iinv.clr - iinv.fq0;   //�ӽ����Ҳ൥Ԫ������൥Ԫ�ĳ�ʼ��������

	//iinv.ai1[ug.lc] = iinv.crl;   //������൥Ԫ��ϵ��
	//iinv.ai2[ug.rc] = iinv.clr;   //�����Ҳ൥Ԫ��ϵ��

	iinv.ai1= 1+iinv.crl;   //������൥Ԫ��ϵ��
	iinv.ai2= 0.5+iinv.clr;   //�����Ҳ൥Ԫ��ϵ��
	
     //iinv.ai1[ug.lc] = iinv.ai1[ug.lc]+crl;   //������൥Ԫ��ϵ��
	//iinv.ai2[ug.rc] = iinv.ai2[ug.rc]+clr;   //�����Ҳ൥Ԫ��ϵ��

	//iinv.flux[IIDX::IIRU] = iinv.rm * gcom.xfn * half* (iinv.ul + iinv.ur) * gcom.farea ;  
	//iinv.flux[IIDX::IIRV] = iinv.rm * gcom.yfn * half* (iinv.vl + iinv.vr)* gcom.farea ; 
	//iinv.flux[IIDX::IIRW] = iinv.rm * gcom.zfn * half* (iinv.vl + iinv.vr)* gcom.farea ; 
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

	iinv.f1[ug.fId] = de2 * de;  //��ԪȨ��
    iinv.f2[ug.fId] = de1 * de;  //�ҵ�ԪȨ��
 
	Real Vau = iinv.f1[ug.fId] * (gcom.cvol1 / iinv.spu1[ug.lc]) + iinv.f2[ug.fId] * (gcom.cvol2 / iinv.spu2[ug.rc]);  //��Vj/a��
	Real Vav = iinv.f1[ug.fId] * (gcom.cvol1 / iinv.spv1[ug.lc]) + iinv.f2[ug.fId] * (gcom.cvol2 / iinv.spv2[ug.lc]);
	Real Vaw = iinv.f1[ug.fId] * (gcom.cvol1 / iinv.spw1[ug.lc]) + iinv.f2[ug.fId] * (gcom.cvol2 / iinv.spw2[ug.lc]);

    iinv.dist[ug.fId] = gcom.xfn * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + gcom.yfn * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + gcom.zfn * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);
	Real Pd1 = visQ.dqdx1[IIDX::IIP] * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + visQ.dqdy1[IIDX::IIP] * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) + visQ.dqdz1[IIDX::IIP] * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);  //ѹ���ݶ���
	Real Pd2 = visQ.dqdx2[IIDX::IIP] * ((*ug.xcc)[ug.rc] - (*ug.xfc)[ug.fId]) + visQ.dqdy2[IIDX::IIP] * ((*ug.ycc)[ug.rc] - (*ug.yfc)[ug.fId]) + visQ.dqdz2[IIDX::IIP] * ((*ug.zcc)[ug.rc] - (*ug.zfc)[ug.fId]);
	Real Pd = Pd1 + Pd2;

	iinv.rm = (iinv.rl + iinv.rr) * half;  //�����ܶ�
	iinv.um = (iinv.f1[ug.fId] *iinv.ul + iinv.f2[ug.fId] *iinv.ur) + (iinv.Vau*gcom.xfn / iinv.dist[ug.fId])*( Pd - (iinv.pr - iinv.pl));  //�����ܶ�
	iinv.vm = (iinv.f1[ug.fId] *iinv.vl + iinv.f2[ug.fId] *iinv.vr) + (iinv.Vav*gcom.yfn / iinv.dist[ug.fId])*(Pd - (iinv.pr - iinv.pl));   //�����ٶ�
	iinv.wm = (iinv.f1[ug.fId] *iinv.wl + iinv.f2[ug.fId] *iinv.wr) + (iinv.Vaw*gcom.zfn / iinv.dist[ug.fId])*(Pd - (iinv.pr - iinv.pl));

	iinv.flux[IIDX::IIRU] = iinv.rm * gcom.xfn * iinv.um * gcom.farea ;
	iinv.flux[IIDX::IIRV] = iinv.rm * gcom.yfn * iinv.vm * gcom.farea ;
	iinv.flux[IIDX::IIRW] = iinv.rm * gcom.zfn * iinv.wm * gcom.farea ;
}

void INsInvterm::CmpINsFaceCorrectPresscoef()
{

     //iinv.Vdvj[ug.fId] = iinv.f1[ug.fId] * ( gcom.cvol1/((1+1)*iinv.sp[ug.lc]- iinv.spj[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 1)*iinv.sp[ug.rc] - iinv.spj[ug.rc]));  // (Vp/dv)j�����������ٶ�������
	 //iinv.aji[ug.fId] = iinv.rm * iinv.Vdvj[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId]; //ajp

	iinv.Vdvu[ug.fId] = iinv.f1[ug.fId] * ( gcom.cvol1/((1+1)*iinv.spu1[ug.lc] - iinv.sp1[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 1)*iinv.spu2[ug.lc] - iinv.sp2[ug.lc]));  // (Vp/dv)j�����������ٶ�������
	iinv.Vdvv[ug.fId] = iinv.f1[ug.fId] * (gcom.cvol1 / ((1 + 1)*iinv.spv1[ug.lc] - iinv.sp1[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 1)*iinv.spv2[ug.lc] - iinv.sp2[ug.lc]));
	iinv.Vdvw[ug.fId] = iinv.f1[ug.fId] * (gcom.cvol1 / ((1 + 1)*iinv.spw1[ug.lc] - iinv.sp1[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 1)*iinv.spw2[ug.lc] - iinv.sp2[ug.lc]));
	
	iinv.aju[ug.fId] = iinv.rm * iinv.Vdvu[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId]; //ajp
	iinv.ajv[ug.fId] = iinv.rm * iinv.Vdvv[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId];
	iinv.ajw[ug.fId] = iinv.rm * iinv.Vdvw[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId];

	iinv.VdU[ug.lc] = -gcom.cvol / ((1 + 1)*iinv.spu1[ug.lc] + iinv.sp1[ug.lc]); //������Ԫ�����ٶ���;
	iinv.VdV[ug.lc] = -gcom.cvol / ((1 + 1)*iinv.spv1[ug.lc] + iinv.sp1[ug.lc]);
	iinv.VdW[ug.lc] = -gcom.cvol / ((1 + 1)*iinv.spw1[ug.lc] + iinv.sp1[ug.lc]);
}


  














EndNameSpace