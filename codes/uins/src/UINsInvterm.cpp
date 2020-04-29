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

#include "UINsInvterm.h"
#include "INsInvterm.h"
#include "UINsVisterm.h"
#include "UINsGrad.h"
#include "Zone.h"
#include "Atmosphere.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "UCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Multigrid.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "UINsLimiter.h"
#include "FieldImp.h"
#include "Iteration.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "Ctrl.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

UINsInvterm::UINsInvterm()
{
    limiter = new INsLimiter();
    limf = limiter->limf;
}

UINsInvterm::~UINsInvterm()
{
    delete limiter;
}

void UINsInvterm::CmpLimiter()
{
    limiter->CmpLimiter();
}

void UINsInvterm::CmpInvFace()  //单元数据重构
{
    //uins_grad.Init();
    //uins_grad.CmpGrad();


    this->CmpLimiter();   //不改

    this->GetQlQrField();  //不改

    this->ReconstructFaceValueField();  //不改

    this->BoundaryQlQrFixField();  //不改
}

void UINsInvterm::GetQlQrField()
{
    limf->GetQlQr();
}

void UINsInvterm::ReconstructFaceValueField()
{
    limf->CmpFaceValue();
    //limf->CmpFaceValueWeighted();
}

void UINsInvterm::BoundaryQlQrFixField()
{
    limf->BcQlQrFix();
}

void UINsInvterm::CmpInvcoff()
{
    if ( inscom.icmpInv == 0 ) return;
    iinv.Init();
    ug.Init();
    uinsf.Init();
    //Alloc();

    this->CmpInvFace();
    this->CmpInvMassFlux();  //需要改动

   //DeAlloc();
}

void UINsInvterm::CmpINsTimestep()
{
	iinv.timestep = 0.01;
}

void UINsInvterm::CmpInvMassFlux()
{
	if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
	{
		iinv.ai1.resize(ug.nTCell);
		iinv.ai2.resize(ug.nTCell);
		iinv.aii1.resize(ug.nFace);
		iinv.aii2.resize(ug.nFace);
		iinv.akk1.resize(ug.nFace);
		iinv.akk2.resize(ug.nFace);
		iinv.ak1.resize(ug.nTCell);
		iinv.ak2.resize(ug.nTCell);
		iinv.vnflow.resize(ug.nFace);
		//iinv.aji.resize(ug.nFace);
		//iinv.fq0.resize(ug.nFace);
		iinv.bm1.resize(ug.nTCell);
		iinv.bm2.resize(ug.nTCell);
		iinv.bmu1.resize(ug.nTCell);
		iinv.bmu2.resize(ug.nTCell);
		iinv.bmv1.resize(ug.nTCell);
		iinv.bmv2.resize(ug.nTCell);
		iinv.bmw1.resize(ug.nTCell);
		iinv.bmw2.resize(ug.nTCell);
		iinv.dist.resize(ug.nFace);
		iinv.f1.resize(ug.nFace);
		iinv.f2.resize(ug.nFace);
		iinv.rf.resize(ug.nFace);
		iinv.uf.resize(ug.nFace);
		iinv.vf.resize(ug.nFace);
		iinv.wf.resize(ug.nFace);
		//iinv.Vdvj.resize(ug.nFace);
		//iinv.Vdv.resize(ug.nFace);
		iinv.Vdvu.resize(ug.nFace);
		iinv.Vdvv.resize(ug.nFace);
		iinv.Vdvw.resize(ug.nFace);
		iinv.aju.resize(ug.nFace);
		iinv.ajv.resize(ug.nFace);
		iinv.ajw.resize(ug.nFace);
		iinv.VdU.resize(ug.nTCell);
		iinv.VdV.resize(ug.nTCell);
		iinv.VdW.resize(ug.nTCell);
		iinv.buc.resize(ug.nTCell);
		iinv.bvc.resize(ug.nTCell);
		iinv.bwc.resize(ug.nTCell);
		//iinv.bc.resize(ug.nCell);
		//iinv.sp.resize(ug.nFace);
		iinv.sp1.resize(ug.nFace);
		iinv.sp2.resize(ug.nFace);
		//iinv.spj.resize(ug.nFace);
		iinv.spu1.resize(ug.nFace);
		iinv.spv1.resize(ug.nFace);
		iinv.spw1.resize(ug.nFace);
		iinv.spu2.resize(ug.nFace);
		iinv.spv2.resize(ug.nFace);
		iinv.spw2.resize(ug.nFace);
		iinv.bpu.resize(ug.nFace);
		iinv.bpv.resize(ug.nFace);
		iinv.bpw.resize(ug.nFace);
		iinv.bp.resize(ug.nFace);
		iinv.ajp.resize(ug.nFace);
		iinv.sppu.resize(ug.nFace);
		iinv.sppv.resize(ug.nFace);
		iinv.sppw.resize(ug.nFace);
		iinv.sj.resize(ug.nTCell);
		iinv.fq.resize(ug.nFace);
		iinv.spu.resize(ug.nTCell);
		iinv.spv.resize(ug.nTCell);
		iinv.spw.resize(ug.nTCell);
		iinv.spuj.resize(ug.nTCell, ug.nTCell);
		iinv.spvj.resize(ug.nTCell, ug.nTCell);
		iinv.spwj.resize(ug.nTCell, ug.nTCell);
		iinv.sjp.resize(ug.nTCell, ug.nTCell);
		//iinv.ajp.resize(ug.nFace);
		//iinv.app.resize(ug.nCell);
		iinv.spp.resize(ug.nTCell);
		iinv.pp.resize(ug.nTCell);
		iinv.pp0.resize(ug.nTCell);
		iinv.pc.resize(ug.nTCell);
		iinv.pp1.resize(ug.nFace);
		iinv.pp2.resize(ug.nFace);
		iinv.uu.resize(ug.nTCell);
		iinv.vv.resize(ug.nTCell);
		iinv.ww.resize(ug.nTCell);
		iinv.uuj.resize(ug.nFace);
		iinv.vvj.resize(ug.nFace);
		iinv.wwj.resize(ug.nFace);
		iinv.muc.resize(ug.nTCell);
		iinv.mvc.resize(ug.nTCell);
		iinv.mwc.resize(ug.nTCell);
		iinv.mp.resize(ug.nTCell);
		iinv.uc.resize(ug.nTCell);
		iinv.vc.resize(ug.nTCell);
		iinv.wc.resize(ug.nTCell);
		iinv.up.resize(ug.nTCell);
		iinv.vp.resize(ug.nTCell);
		iinv.wp.resize(ug.nTCell);
		iinv.spt.resize(ug.nTCell);
		iinv.but.resize(ug.nTCell);
		iinv.bvt.resize(ug.nTCell);
		iinv.bwt.resize(ug.nTCell);
		iinv.ppr.resize(ug.nFace);
		iinv.dqqdx.resize(ug.nFace);
		iinv.dqqdy.resize(ug.nFace);
		iinv.dqqdz.resize(ug.nFace);
		iinv.fux.resize(ug.nFace);

		iinv.ai1 = 0;
		iinv.ai2 = 0;
		iinv.spu1 = 1;
		iinv.spv1 = 1;
		iinv.spw1 = 1;
		iinv.spu2 = 1;
		iinv.spv2 = 1;
		iinv.spw2 = 1;
		
		//iinv.bm = 0;
		iinv.buc = 0;
		iinv.bvc = 0;
		iinv.bwc = 0;
		iinv.sp1 = 0;
		iinv.sp2 = 0;
		iinv.spj = 0;
		iinv.spp = 0;
		iinv.sppu = 0;
		iinv.sppv = 0;
		iinv.sppw = 0;

		iinv.bpu = 0;
		iinv.bpv = 0;
		iinv.bpw = 0;
		iinv.bp = 0;
		iinv.pp = 0;
		iinv.pp0 = 0;
		iinv.pp1 = 0;
		iinv.pp2 = 0;

		iinv.muc = 0;
		iinv.mvc = 0;
		iinv.mwc = 0;

		iinv.uc = 0;
		iinv.vc = 0;
		iinv.wc = 0;

		iinv.uu = 0;
		iinv.vv = 0;
		iinv.ww = 0;


		for (int fId = 0; fId < ug.nFace; ++fId)
		{
			ug.fId = fId;

			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			this->PrepareFaceValue();

			this->CmpINsinvTerm();

		}
	}

	else
	{
	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareFaceValue();

		this->CmpINsinvTerm();

	}
	}
}

void UINsInvterm::PrepareFaceValue()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    inscom.gama1 = ( * uinsf.gama )[ 0 ][ ug.lc ];
    inscom.gama2 = ( * uinsf.gama )[ 0 ][ ug.rc ];

    iinv.gama1 = inscom.gama1;
    iinv.gama2 = inscom.gama2;

    for ( int iEqu = 0; iEqu < limf->nEqu; ++ iEqu )
    {
        iinv.prim1[ iEqu ] = (*limf->q)[iEqu][ug.lc];
        iinv.prim2[ iEqu ] = (*limf->q)[iEqu][ug.rc];
    }
}

void UINsInvterm::PrepareProFaceValue()
{
	gcom.xfn = (*ug.xfn)[ug.fId];
	gcom.yfn = (*ug.yfn)[ug.fId];
	gcom.zfn = (*ug.zfn)[ug.fId];
	gcom.vfn = (*ug.vfn)[ug.fId];
	gcom.farea = (*ug.farea)[ug.fId];


	for (int iEqu = 0; iEqu < limf->nEqu; ++iEqu)
	{
		iinv.prim1[iEqu] = (*uinsf.q)[iEqu][ug.lc];
		iinv.prim2[iEqu] = (*uinsf.q)[iEqu][ug.rc];
	}
}

void UINsInvterm::MomPre()
{

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;
		int fn = (*ug.c2f)[ug.cId].size();
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];
			if (ug.cId == ug.lc)
			{
				iinv.muc[ug.cId] += iinv.spuj[ug.cId][ug.rc]*iinv.ur;   //使用高斯赛德尔迭代时，相邻单元对其的影响，矩阵法不需要
				iinv.mvc[ug.cId] += iinv.spvj[ug.cId][ug.rc]*iinv.vr;
				iinv.mwc[ug.cId] += iinv.spwj[ug.cId][ug.rc]*iinv.wr;
			}
			else if (ug.cId == ug.rc)
			{
				iinv.muc[ug.cId] += iinv.spuj[ug.cId][ug.lc]*iinv.ul; //使用高斯赛德尔迭代时，相邻单元对其的影响，矩阵法不需要
				iinv.mvc[ug.cId] += iinv.spvj[ug.cId][ug.lc]*iinv.vl;
				iinv.mwc[ug.cId] += iinv.spwj[ug.cId][ug.lc]*iinv.wl;
			}
		}


		if (!iinv.spu[ug.cId] == 0)
		{
			iinv.uc[ug.cId] = (iinv.muc[ug.cId] + iinv.buc[ug.cId]) / (iinv.spu[ug.cId]);  //下一时刻速度的预测值
			(*uinsf.q)[IIDX::IIU][ug.cId] = iinv.uc[ug.cId];
		}
		else
		{
			//iinv.uc[ug.cId] = (iinv.muc[ug.cId] + iinv.buc[ug.cId]) / (0.01+iinv.spu[ug.cId]);
			iinv.uc[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];
		}
		if (!iinv.spv[ug.cId] == 0)
		{
			iinv.vc[ug.cId] = (iinv.mvc[ug.cId] + iinv.bvc[ug.cId]) / (iinv.spv[ug.cId]);
			(*uinsf.q)[IIDX::IIV][ug.cId] = iinv.vc[ug.cId];
		}
		else
		{
			//iinv.vc[ug.cId] = (iinv.mvc[ug.cId] + iinv.bvc[ug.cId]) / (0.01+iinv.spv[ug.cId]);
			iinv.vc[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
		}
		if (!iinv.spv[ug.cId] == 0)
		{
			iinv.wc[ug.cId] = (iinv.mwc[ug.cId] + iinv.bwc[ug.cId]) / (iinv.spw[ug.cId]);
			(*uinsf.q)[IIDX::IIW][ug.cId] = iinv.wc[ug.cId];
		}
		else
		{
			//iinv.wc[ug.cId] = (iinv.mwc[ug.cId] + iinv.bwc[ug.cId]) / (0.01+iinv.spw[ug.cId]);
			 iinv.wc[ug.cId]=(*uinsf.q)[IIDX::IIW][ug.cId];
		}

		//(*uinsf.q)[IIDX::IIU][ug.cId] = iinv.uc[ug.cId];
		//(*uinsf.q)[IIDX::IIV][ug.cId] = iinv.vc[ug.cId];
		//(*uinsf.q)[IIDX::IIW][ug.cId] = iinv.wc[ug.cId];

    }

	    this->CmpINsMomRes();
}



void UINsInvterm::CmpFaceflux()
{

	iinv.Init();
	ug.Init();
	uinsf.Init();
	Alloc();
	//this->CmpInvFace();  //边界处理
	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareProFaceValue();

		this->CmpINsFaceflux();
	}

}

void UINsInvterm::CmpINsMomRes()
{
	iinv.res_u = 0;
	iinv.res_v = 0;
	iinv.res_w = 0;

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.res_u += (iinv.buc[ug.cId]+ iinv.muc[ug.cId]- iinv.uc[ug.cId]* (iinv.spu[ug.cId]))*(iinv.buc[ug.cId] + iinv.muc[ug.cId] - iinv.uc[ug.cId] * (iinv.spu[ug.cId]));
		iinv.res_v += (iinv.bvc[ug.cId] + iinv.mvc[ug.cId] - iinv.vc[ug.cId] * (iinv.spv[ug.cId]))*(iinv.bvc[ug.cId] + iinv.mvc[ug.cId] - iinv.vc[ug.cId] * (iinv.spv[ug.cId]));
		iinv.res_w += (iinv.bwc[ug.cId] + iinv.mwc[ug.cId] - iinv.wc[ug.cId] * (iinv.spw[ug.cId]))*(iinv.bwc[ug.cId] + iinv.mwc[ug.cId] - iinv.wc[ug.cId] * (iinv.spw[ug.cId]));

	}

	iinv.res_u = sqrt(iinv.res_u);
	iinv.res_v = sqrt(iinv.res_v);
	iinv.res_w = sqrt(iinv.res_w);
	
}

void UINsInvterm::AddFlux()
{
	UnsGrid * grid = Zone::GetUnsGrid();
	MRField * res = GetFieldPointer< MRField >(grid, "res");
	int nEqu = res->GetNEqu();
	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];
		//if ( ug.lc == 0 ) cout << fId << endl;

		for (int iEqu = 0; iEqu < nEqu; ++iEqu)
		{
			(*res)[iEqu][ug.lc] -= (*iinvflux)[iEqu][ug.fId];
		}
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//if ( ug.lc == 0 || ug.rc == 0 ) cout << fId << endl;

		for (int iEqu = 0; iEqu < nEqu; ++iEqu)
		{
			(*res)[iEqu][ug.lc] -= (*iinvflux)[iEqu][ug.fId];
			(*res)[iEqu][ug.rc] += (*iinvflux)[iEqu][ug.fId];
		}
	}

	//ONEFLOW::AddF2CField(res, iinvflux);

}

void UINsInvterm::CmpCorrectPresscoef()
{

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->CmpINsFaceCorrectPresscoef();
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.VdU[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spu[ug.cId] + iinv.sj[ug.cId]); //用于求单元修正速度量;
		iinv.VdV[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spv[ug.cId] + iinv.sj[ug.cId]);
		iinv.VdW[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spw[ug.cId] + iinv.sj[ug.cId]);

		int fn = (*ug.c2f)[ug.cId].size();
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			iinv.bp[ug.cId] -= iinv.fq[ug.fId];  //压力修正方程源项

			if (ug.cId == ug.lc)
			{
				iinv.sjp[ug.cId][ug.rc] = -iinv.ajp[ug.fId]; //求解压力修正方程的非零系数
			}
			else if(ug.cId==ug.rc)
			{
				iinv.sjp[ug.cId][ug.lc] = -iinv.ajp[ug.fId];
			}
			
			iinv.spp[ug.cId] += iinv.ajp[ug.fId];   //压力修正方程的矩阵主对角项
		}
	}
}

void UINsInvterm::CmpPressCorrectEqu()
{
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;
		int fn = (*ug.c2f)[ug.cId].size();
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];
			if (ug.cId == ug.lc)
			{
				iinv.mp[ug.cId] += iinv.ajp[ug.fId] * iinv.pp0[ug.rc]; //高斯赛戴尔迭代求解时的相邻单元的值，矩阵法不需要
			}
			else if (ug.cId == ug.rc)
			{
				iinv.mp[ug.cId] += iinv.ajp[ug.fId] * iinv.pp0[ug.lc];
			}
		}
		iinv.pp0[ug.cId] = (iinv.bp[ug.cId]+iinv.mp[ug.cId]) / (0.01+iinv.spp[ug.cId]); //压力修正值
		iinv.pp[ug.cId] = iinv.pp0[ug.cId]; //当前时刻的压力修正值

		iinv.pc[ug.cId] = (*uinsf.q)[IIDX::IIP][ug.cId] + iinv.pp[ug.cId]; //下一时刻的压力值

		//(*uinsf.q)[IIDX::IIP][ug.cId] = (*uinsf.q)[IIDX::IIP][ug.cId] + iinv.pp[ug.cId];
	}
	this->CmpINsPreRes();
}

void UINsInvterm::CmpINsPreRes()
{
	iinv.res_p = 0;


	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.res_p += (iinv.bp[ug.cId] + iinv.mp[ug.cId] - iinv.pp0[ug.cId]* (0.01+iinv.spp[ug.cId]))*(iinv.bp[ug.cId] + iinv.mp[ug.cId] - iinv.pp0[ug.cId]* (0.01+iinv.spp[ug.cId]));
	}

	iinv.res_p = sqrt(iinv.res_p);
	iinv.res_p = 0;
}


void UINsInvterm::UpdateFaceflux()
{
	iinv.Init();
	ug.Init();
	uinsf.Init();
	Alloc();
	//this->CmpInvFace();  //边界处理
	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//this->PrepareFaceValue();

		this->CmpUpdateINsFaceflux();

	}

}

void UINsInvterm::CmpUpdateINsFaceflux()
{
	iinv.uuj[ug.fId] =  iinv.Vdvu[ug.fId] *(iinv.pp[ug.lc] - iinv.pp[ug.rc]) *gcom.xfn / iinv.dist[ug.fId]; //面速度修正量
	iinv.vvj[ug.fId] =  iinv.Vdvv[ug.fId] *(iinv.pp[ug.lc] - iinv.pp[ug.rc]) *gcom.yfn / iinv.dist[ug.fId];
	iinv.wwj[ug.fId] =  iinv.Vdvw[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) *gcom.zfn / iinv.dist[ug.fId];

	iinv.uf[ug.fId] = iinv.uf[ug.fId] + iinv.uuj[ug.fId]; //下一时刻面速度
	iinv.vf[ug.fId] = iinv.vf[ug.fId] + iinv.vvj[ug.fId];
	iinv.wf[ug.fId] = iinv.wf[ug.fId] + iinv.wwj[ug.fId];

	iinv.fux[ug.fId] = iinv.rf[ug.fId] * (gcom.xfn * iinv.uuj[ug.fId] + gcom.yfn * iinv.vvj[ug.fId] + gcom.zfn * iinv.wwj[ug.fId])*gcom.farea;
	iinv.fq[ug.fId] = iinv.fq[ug.fId] + iinv.fux[ug.fId];

}

void UINsInvterm::UpdateSpeed()
{
	this->CmpPreGrad();

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.uu[ug.cId] = iinv.VdU[ug.cId] * iinv.dqqdx[ug.cId]; //速度修正量
		iinv.vv[ug.cId] = iinv.VdV[ug.cId] * iinv.dqqdy[ug.cId];
		iinv.ww[ug.cId] = iinv.VdW[ug.cId] * iinv.dqqdz[ug.cId];

		iinv.up[ug.cId]= iinv.uc[cId] + iinv.uu[ug.cId];  //下一时刻的速度值
		iinv.vp[ug.cId]= iinv.vc[cId] + iinv.vv[ug.cId];
		iinv.wp[ug.cId]= iinv.wc[cId] + iinv.ww[ug.cId];

		//(*uinsf.q)[IIDX::IIU][ug.cId] = iinv.up[ug.cId];
		//(*uinsf.q)[IIDX::IIV][ug.cId] = iinv.vp[ug.cId];
		//(*uinsf.q)[IIDX::IIW][ug.cId] = iinv.wp[ug.cId];

	}	
}

void UINsInvterm::CmpPreGrad()
{
	iinv.dqqdx = 0;
	iinv.dqqdy = 0;
	iinv.dqqdz = 0;

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		if (fId == 432)
		{
			int kkk = 1;
		}
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real dxl = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dyl = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dzl = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real dxr = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
		Real dyr = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
		Real dzr = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

		Real delt1 = DIST(dxl, dyl, dzl);
		Real delt2 = DIST(dxr, dyr, dzr);
		Real delta = 1.0 / (delt1 + delt2 + SMALL);

		Real cl = delt2 * delta;
		Real cr = delt1 * delta;

		Real value = cl * iinv.pp[ug.lc] + cr * iinv.pp[ug.rc];

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		iinv.dqqdx[ug.lc] += fnxa * value;
		iinv.dqqdy[ug.lc] += fnya * value;
		iinv.dqqdz[ug.lc] += fnza * value;

		if (ug.fId < ug.nBFace) continue;
		iinv.dqqdx[ug.rc] -= fnxa * value;
		iinv.dqqdy[ug.rc] -= fnya * value;
		iinv.dqqdz[ug.rc] -= fnza * value;
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		Real ovol = one / (*ug.cvol)[cId];
		iinv.dqqdx[cId] *= ovol;
		iinv.dqqdy[cId] *= ovol;
		iinv.dqqdz[cId] *= ovol;
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.dqqdx[ug.rc] = iinv.dqqdx[ug.lc];
		iinv.dqqdy[ug.rc] = iinv.dqqdy[ug.lc];
		iinv.dqqdz[ug.rc] = iinv.dqqdz[ug.lc];
	}
}


void UINsInvterm::Alloc()
{
    iinvflux = new MRField( inscom.nEqu, ug.nFace );
}

void UINsInvterm::DeAlloc()
{
    delete iinvflux;
}


void UINsInvterm::ReadTmp()
{
    static int iii = 0;
    if ( iii ) return;
    iii = 1;
    fstream file;
    file.open( "nsflow.dat", ios_base::in | ios_base::binary );
    if ( ! file ) exit( 0 );

    uinsf.Init();

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        for ( int iEqu = 0; iEqu < 5; ++ iEqu )
        {
            file.read( reinterpret_cast< char * >( & ( * uinsf.q )[ iEqu ][ cId ] ), sizeof( double ) );
        }
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * uinsf.visl )[ 0 ][ cId ] ), sizeof( double ) );
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * uinsf.vist )[ 0 ][ cId ] ), sizeof( double ) );
    }

    vector< Real > tmp1( ug.nTCell ), tmp2( ug.nTCell );

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        tmp1[ cId ] = ( * uinsf.timestep )[ 0 ][ cId ];
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * uinsf.timestep )[ 0 ][ cId ] ), sizeof( double ) );
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        tmp2[ cId ] = ( * uinsf.timestep )[ 0 ][ cId ];
    }

    turbcom.Init();
    uturbf.Init();
    for ( int iCell = 0; iCell < ug.nTCell; ++ iCell )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            file.read( reinterpret_cast< char * >( & ( * uturbf.q )[ iEqu ][ iCell ] ), sizeof( double ) );
        }
    }
    file.close();
    file.clear();
}



EndNameSpace