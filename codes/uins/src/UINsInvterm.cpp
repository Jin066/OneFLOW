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
	iinv.sppu.resize(ug.nFace);
	iinv.sppv.resize(ug.nFace);
	iinv.sppw.resize(ug.nFace);
	//iinv.ajp.resize(ug.nFace);
	//iinv.app.resize(ug.nCell);
	//iinv.spp.resize(ug.nCell);
	iinv.pp.resize(ug.nCell);
	iinv.pp1.resize(ug.nFace);
	iinv.pp2.resize(ug.nFace);
	iinv.uu.resize(ug.nCell);
	iinv.vv.resize(ug.nCell);
	iinv.ww.resize(ug.nCell);
	iinv.uuj.resize(ug.nFace);
	iinv.vvj.resize(ug.nFace);
	iinv.wwj.resize(ug.nFace);
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
	iinv.pp = 0;
	iinv.pp1 = 0;
	iinv.pp2 = 0;
	iinv.uu = 0;
	iinv.vv = 0;
	iinv.ww = 0;

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

void UINsInvterm::MomPred()
{
	;
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

void UINsInvterm::AddFlux()
{
	UnsGrid * grid = Zone::GetUnsGrid();
	MRField * res = GetFieldPointer< MRField >(grid, "res");

	ONEFLOW::AddF2CField(res, iinvflux);
//	//ONEFLOW::AddF2CFieldDebug( res, invflux );
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

			iinv.bpu[ug.lc] -= iinv.flux[IIDX::IIRU];  //源项
			iinv.bpv[ug.lc] -= iinv.flux[IIDX::IIRV];
			iinv.bpw[ug.lc] -= iinv.flux[IIDX::IIRW];
			iinv.sppu[ug.lc] += iinv.aju[ug.fId];  //主参数
			iinv.sppv[ug.lc] += iinv.ajv[ug.fId];
			iinv.sppw[ug.lc] += iinv.ajw[ug.fId];
		}
		//iinv.Vdv[ug.cId] = -gcom.cvol / ((1 + 1)*iinv.sp[ug.cId] + iinv.spj[ug.cId]); //用于求单元修正速度量
	}
}

void UINsInvterm::CmpPressCorrectEqu()
{
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		inscom.prim[IIDX::IIP] = inscom.prim[IIDX::IIP]+iinv.pp[ug.cId];
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
	this->AddFlux();
}

void UINsInvterm::CmpUpdateINsFaceflux()
{
	iinv.uuj[ug.fId] =  iinv.Vdvu[ug.fId] *(iinv.pp1[ug.lc] - iinv.pp2[ug.rc]) *gcom.xfn / iinv.dist[ug.fId];
	iinv.vvj[ug.fId] =  iinv.Vdvv[ug.fId] *(iinv.pp1[ug.lc] - iinv.pp2[ug.lc]) *gcom.yfn / iinv.dist[ug.fId];
	iinv.wwj[ug.fId] =  iinv.Vdvw[ug.fId] * (iinv.pp1[ug.lc] - iinv.pp2[ug.lc]) *gcom.zfn / iinv.dist[ug.fId];

	//PrimToQ(iinv.prim1, iinv.gama1, iinv.q1);
	//PrimToQ(iinv.prim2, iinv.gama2, iinv.q2);

	//for (int iEqu = 0; iEqu < inscom.nEqu; ++iEqu)
	//{
	//	iinv.dq[iEqu] = iinv.q2[iEqu] - iinv.q1[iEqu];
	//}

	//iinv.flux[IIDX::IIR] = iinv.flux[IIDX::IIR];
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
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.uu[ug.cId] = iinv.VdU[ug.lc] * 1; //应该乘压力修正的梯度
		iinv.vv[ug.cId] = iinv.VdV[ug.lc] * 1;
		iinv.ww[ug.cId] = iinv.VdW[ug.lc] * 1;

		inscom.prim[IIDX::IIU] = inscom.prim[IIDX::IIU] + iinv.uu[ug.cId];
		inscom.prim[IIDX::IIV] = inscom.prim[IIDX::IIV] + iinv.vv[ug.cId];
		inscom.prim[IIDX::IIW] = inscom.prim[IIDX::IIW] + iinv.ww[ug.cId];
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