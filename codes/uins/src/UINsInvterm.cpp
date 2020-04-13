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
#include "Boundary.h"
#include "BcRecord.h"
#include "UINsLimiter.h"
#include "FieldImp.h"
#include "Iteration.h"
#include "TurbCom.h"
#include "UTurbCom.h"
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

void UINsInvterm::CmpInvFace()  //不改动
{
    //uins_grad.Init();
    //uins_grad.CmpGrad();

    this->CmpLimiter();   //不改

    this->GetQlQrField();  //不改

    //this->ReconstructFaceValueField();  //不改

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
    Alloc();

    //this->SetPointer( inscom.ischeme );

    //ReadTmp();
    this->CmpInvFace();
    this->CmpInvMassFlux();  //需要改动
    //this->AddInv();

    DeAlloc();
}

void UINsInvterm::CmpInvMassFlux()
{
    if ( Iteration::outerSteps == 2 )
    {
        int kkk = 1;
    }
	iinv.ai1.resize(ug.nFace);
	iinv.ai2.resize(ug.nFace);
	//iinv.aji.resize(ug.nFace);
	//iinv.fq0.resize(ug.nFace);
	iinv.bm.resize(ug.nFace);
	iinv.dist.resize(ug.nFace);
	iinv.f1.resize(ug.nFace);
	iinv.f2.resize(ug.nFace);
	//iinv.Vdvj.resize(ug.nFace);
	//iinv.Vdv.resize(ug.nFace);
	iinv.Vdvu.resize(ug.nFace);
	iinv.Vdvv.resize(ug.nFace);
	iinv.Vdvw.resize(ug.nFace);
	iinv.aju.resize(ug.nFace);
	iinv.ajv.resize(ug.nFace);
	iinv.ajw.resize(ug.nFace);
	iinv.VdU.resize(ug.nFace);
	iinv.VdV.resize(ug.nFace);
	iinv.VdW.resize(ug.nFace);
	iinv.buc.resize(ug.nCell);
	iinv.bvc.resize(ug.nCell);
	iinv.bwc.resize(ug.nCell);
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
	//iinv.ajp.resize(ug.nFace);
	//iinv.app.resize(ug.nCell);
	iinv.spp.resize(ug.nFace);
	iinv.pp.resize(ug.nFace);
	iinv.pp1.resize(ug.nFace);
	iinv.pp2.resize(ug.nFace);
	iinv.uu.resize(ug.nCell);
	iinv.vv.resize(ug.nCell);
	iinv.ww.resize(ug.nCell);
	iinv.uuj.resize(ug.nFace);
	iinv.vvj.resize(ug.nFace);
	iinv.wwj.resize(ug.nFace);
	iinv.muc.resize(ug.nCell);
	iinv.mvc.resize(ug.nCell);
	iinv.mwc.resize(ug.nCell);
	iinv.mp.resize(ug.nCell);
	iinv.uc.resize(ug.nCell);
	iinv.vc.resize(ug.nCell);
	iinv.wc.resize(ug.nCell);
	iinv.ppr.resize(ug.nFace);
	iinv.dqqdx.resize(ug.nFace);
	iinv.dqqdy.resize(ug.nFace);
	iinv.dqqdz.resize(ug.nFace);

	iinv.ai1 = 0;
	iinv.ai2 = 0;
	iinv.spu1 = 1;
	iinv.spv1 = 1;
	iinv.spw1 = 1;
	iinv.spu2 = 1;
	iinv.spv2 = 1;
	iinv.spw2 = 1;

	iinv.bm = 0;
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
	iinv.bp= 0;
	iinv.pp = 0;
	iinv.pp1 = 0;
	iinv.pp2 = 0;
	iinv.uu = 0;
	iinv.vv = 0;
	iinv.ww = 0;

	iinv.muc = 0;
	iinv.mvc = 0;
	iinv.mwc = 0;
	

    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;

        if ( fId == 10127 )
        {
            int kkk = 1;
        }

        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareFaceValue();

        this->CmpINsinvTerm();

        //this->UpdateFaceInvFlux();
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
        iinv.prim1[ iEqu ] = ( * limf->qf1 )[ iEqu ][ ug.fId ];
        iinv.prim2[ iEqu ] = ( * limf->qf2 )[ iEqu ][ ug.fId ];
    }
}

void UINsInvterm::MomPre()
{

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		int fn = (*ug.c2f)[ug.cId].size();
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
	        ug.rc = (*ug.rcf)[ug.fId];

 			iinv.muc[ug.lc] += iinv.ai2[ug.rc] * iinv.ur;  
			iinv.mvc[ug.lc] += iinv.ai2[ug.rc] * iinv.vr;
			iinv.mwc[ug.lc] += iinv.ai2[ug.rc] * iinv.wr;
		}
     	iinv.uc[ug.cId] = (iinv.muc[ug.lc] + iinv.buc[ug.lc]) / iinv.spu1[ug.lc];
		iinv.vc[ug.cId] = (iinv.mvc[ug.lc] + iinv.bvc[ug.lc]) / iinv.spv1[ug.lc];
		iinv.wc[ug.cId] = (iinv.mwc[ug.lc] + iinv.bwc[ug.lc]) / iinv.spw1[ug.lc];

	    inscom.prim[IIDX::IIU] = iinv.uc[ug.cId];
	    inscom.prim[IIDX::IIV] = iinv.vc[ug.cId];
	    inscom.prim[IIDX::IIW] = iinv.wc[ug.cId];
    }

	    this->CmpINsRes();
}



void UINsInvterm::CmpFaceflux()
{

	iinv.Init();
	ug.Init();
	uinsf.Init();
	Alloc();
	this->CmpInvFace();  //边界处理
	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareFaceValue();

		this->CmpINsFaceflux();
	}

	   
	//this->AddFlux();
}

void UINsInvterm::CmpINsRes()
{
	
	UnsGrid * grid = Zone::GetUnsGrid();
	MRField * res = GetFieldPointer< MRField >(grid, "res");
	int nEqu = res->GetNEqu();
	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];
		//if ( ug.lc == 0 ) cout << fId << endl;

		
			(*res)[0][ug.lc] = 0;
			(*res)[1][ug.lc] = iinv.spu1[ug.lc]* iinv.ul+iinv.muc[ug.lc]- iinv.buc[ug.lc];
			(*res)[2][ug.lc] = iinv.spv1[ug.lc] * iinv.vl + iinv.mvc[ug.lc] - iinv.bvc[ug.lc];
			(*res)[3][ug.lc] = iinv.spw1[ug.lc] * iinv.wl + iinv.mwc[ug.lc] - iinv.bwc[ug.lc];
			(*res)[4][ug.lc] = 0;
		
	}
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

		//this->PrepareFaceValue();

		this->CmpINsFaceCorrectPresscoef();
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		int fn = (*ug.c2f)[ug.cId].size();
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			//iinv.bpu[ug.lc] -= iinv.flux[IIDX::IIRU];  //源项
			//iinv.bpv[ug.lc] -= iinv.flux[IIDX::IIRV];
			//iinv.bpw[ug.lc] -= iinv.flux[IIDX::IIRW];

			iinv.bp[ug.lc] -= iinv.flux[IIDX::IIRU] + iinv.flux[IIDX::IIRV] + iinv.flux[IIDX::IIRW];


			//iinv.sppu[ug.lc] += iinv.aju[ug.fId];  //主参数
			//iinv.sppv[ug.lc] += iinv.ajv[ug.fId];
			//iinv.sppw[ug.lc] += iinv.ajw[ug.fId];

			iinv.spp[ug.lc] += iinv.ajp[ug.fId];
		}
		//iinv.Vdv[ug.cId] = -gcom.cvol / ((1 + 1)*iinv.sp[ug.cId] + iinv.spj[ug.cId]); //用于求单元修正速度量
	}
}

void UINsInvterm::CmpPressCorrectEqu()
{
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		int fn = (*ug.c2f)[ug.cId].size();
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];


			iinv.mp[ug.cId] += iinv.ajp[ug.fId] * iinv.pp[ug.rc];

		}
		iinv.pp[ug.lc] = (iinv.mp[ug.cId] + iinv.bp[ug.lc]) / iinv.spp[ug.lc];

		inscom.prim[IIDX::IIP] = inscom.prim[IIDX::IIP] + iinv.pp[ug.lc];
	}
}


void UINsInvterm::UpdateFaceflux()
{
	iinv.Init();
	ug.Init();
	uinsf.Init();
	Alloc();
	this->CmpInvFace();  //边界处理
	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareFaceValue();

		this->CmpUpdateINsFaceflux();

	}
	//for (int fId = 0; fId < ug.nFace; ++fId)
	//{
	//		ug.fId = fId;
	//		ug.lc = (*ug.lcf)[ug.fId];
	//		ug.rc = (*ug.rcf)[ug.fId];
	//		iinv.uuj = 0*iinv.Vdvu *(iinv.pp1-iinv.pp2) *gcom.xfn / iinv.dist;
	//		iinv.vvj = 0*iinv.Vdvv *(iinv.pp1 - iinv.pp2) *gcom.yfn / iinv.dist;
	//		iinv.wwj = 0*iinv.Vdvw * (iinv.pp1 - iinv.pp2) *gcom.zfn / iinv.dist;

	//		(*iinvflux)[0][ug.fId] = iinv.flux[IIDX::IIRU]+ 0 * iinv.rm * gcom.xfn * iinv.uuj * gcom.farea;
	//		(*iinvflux)[1][ug.fId] = iinv.flux[IIDX::IIRV] +0 * iinv.rm * gcom.xfn * iinv.vvj * gcom.farea;
	//		(*iinvflux)[2][ug.fId] = iinv.flux[IIDX::IIRW] +0 * iinv.rm * gcom.xfn * iinv.wwj * gcom.farea;
	//		(*iinvflux)[3][ug.fId] = 0;
	//		(*iinvflux)[4][ug.fId] = 0;
//	}

	 //this->AddFlux();
}

void UINsInvterm::CmpUpdateINsFaceflux()
{
	iinv.uuj[ug.fId] =  iinv.Vdvu[ug.fId] *(iinv.pp[ug.lc] - iinv.pp[ug.rc]) *gcom.xfn / iinv.dist[ug.fId];
	iinv.vvj[ug.fId] =  iinv.Vdvv[ug.fId] *(iinv.pp[ug.lc] - iinv.pp[ug.rc]) *gcom.yfn / iinv.dist[ug.fId];
	iinv.wwj[ug.fId] =  iinv.Vdvw[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) *gcom.zfn / iinv.dist[ug.fId];

	//PrimToQ(iinv.prim1, iinv.gama1, iinv.q1);
	//PrimToQ(iinv.prim2, iinv.gama2, iinv.q2);

	//for (int iEqu = 0; iEqu < inscom.nEqu; ++iEqu)
	//{
	//	iinv.dq[iEqu] = iinv.q2[iEqu] - iinv.q1[iEqu];
	//}

	//iinv.flux[IIDX::IIR] = iinv.flux[IIDX::IIR];

	iinv.um = iinv.um + iinv.uuj[ug.fId];
	iinv.vm = iinv.vm + iinv.vvj[ug.fId];
	iinv.wm = iinv.wm + iinv.wwj[ug.fId];

	iinv.prim[IIDX::IIR] = iinv.rm;
	iinv.prim[IIDX::IIU] = iinv.um;
	iinv.prim[IIDX::IIV] = iinv.vm;
	iinv.prim[IIDX::IIW] = iinv.wm;
	iinv.prim[IIDX::IIP] = iinv.pm;

    iinv.flux[IIDX::IIRU] += iinv.rm * gcom.xfn * iinv.uuj[ug.fId] * gcom.farea;
	iinv.flux[IIDX::IIRV] += iinv.rm * gcom.yfn * iinv.vvj[ug.fId] * gcom.farea;
	iinv.flux[IIDX::IIRW] += iinv.rm * gcom.zfn * iinv.wwj[ug.fId] * gcom.farea;
	//iinv.flux[IIDX::IIRE] = iinv.flux[IIDX::IIRE];



	for (int iEqu = 0; iEqu < inscom.nTEqu; ++iEqu)
	{
		(*iinvflux)[iEqu][ug.fId] = iinv.flux[iEqu];
	}

	//(*iinvflux)[0][ug.fId] = iinv.flux[IIDX::IIRU] + iinv.rm * gcom.xfn * iinv.uuj[ug.fId] * gcom.farea;
	//(*iinvflux)[1][ug.fId] = iinv.flux[IIDX::IIRV] + iinv.rm * gcom.xfn * iinv.vvj[ug.fId] * gcom.farea;
	//(*iinvflux)[2][ug.fId] = iinv.flux[IIDX::IIRW] + iinv.rm * gcom.xfn * iinv.wwj[ug.fId] * gcom.farea;
	//(*iinvflux)[3][ug.fId] = 0;
	//(*iinvflux)[4][ug.fId] = 0;
}

void UINsInvterm::UpdateSpeed()
{
	this->CmpPreGrad();

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.uu[ug.cId] = iinv.VdU[ug.lc] * iinv.dqqdx[ug.lc]; //应该乘压力修正的梯度
		iinv.vv[ug.cId] = iinv.VdV[ug.lc] * iinv.dqqdy[ug.lc];
		iinv.ww[ug.cId] = iinv.VdW[ug.lc] * iinv.dqqdz[ug.lc];

		inscom.prim[IIDX::IIU] = inscom.prim[IIDX::IIU] + iinv.uu[ug.cId];
		inscom.prim[IIDX::IIV] = inscom.prim[IIDX::IIV] + iinv.vv[ug.cId];
		inscom.prim[IIDX::IIW] = inscom.prim[IIDX::IIW] + iinv.ww[ug.cId];
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