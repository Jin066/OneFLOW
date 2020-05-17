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

#include "UINsVisterm.h"
#include "INsInvterm.h"
#include "UINsInvterm.h"
#include "INsVisterm.h"
#include "Iteration.h"
#include "HeatFlux.h"
#include "Zone.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "INsCtrl.h"
#include "UCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "VisGrad.h"
#include "UINsGrad.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "ULimiter.h"
#include "FieldImp.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )


UINsVisterm::UINsVisterm()
{
    ;
}

UINsVisterm::~UINsVisterm()
{
    ;
}


void UINsVisterm::CmpViscoff()
{
    if ( vis_model.vismodel == 0 ) return;
    ug.Init();
    uinsf.Init();
    visQ.Init( inscom.nEqu );

    Alloc();

    this->PrepareField();
    this->CmpVisterm();

    DeAlloc();
}

void UINsVisterm::Alloc()
{
    visflux = new MRField( inscom.nEqu, ug.nFace );
}

void UINsVisterm::DeAlloc()
{
    delete visflux;
}

void UINsVisterm::PrepareField()
{
	//uins_grad.Init();
	//uins_grad.CmpGrad();  //计算梯度
    //ut_grad.CmpGradDebug();
	this->CmpPreandVisGrad();
}

void UINsVisterm::CmpPreandVisGrad()
{
	(*uinsf.dqdx)[IIDX::IIR] = 0;
	(*uinsf.dqdy)[IIDX::IIR] = 0;
	(*uinsf.dqdz)[IIDX::IIR] = 0;

	(*uinsf.dqdx)[IIDX::IIU] = 0;
	(*uinsf.dqdy)[IIDX::IIU] = 0;
	(*uinsf.dqdz)[IIDX::IIU] = 0;

	(*uinsf.dqdx)[IIDX::IIV] = 0;
	(*uinsf.dqdy)[IIDX::IIV] = 0;
	(*uinsf.dqdz)[IIDX::IIV] = 0;

	(*uinsf.dqdx)[IIDX::IIW] = 0;
	(*uinsf.dqdy)[IIDX::IIW] = 0;
	(*uinsf.dqdz)[IIDX::IIW] = 0;

	(*uinsf.dqdx)[IIDX::IIP] = 0;
	(*uinsf.dqdy)[IIDX::IIP] = 0;
	(*uinsf.dqdz)[IIDX::IIP] = 0;

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
		Real delta = 1.0 / (delt1 + delt2);

		Real cl = delt2 * delta;
		Real cr = delt1 * delta;

		Real value1 = cl * (*uinsf.q)[IIDX::IIU][ug.lc] + cr * (*uinsf.q)[IIDX::IIU][ug.rc];
		Real value2 = cl * (*uinsf.q)[IIDX::IIV][ug.lc] + cr * (*uinsf.q)[IIDX::IIV][ug.rc];
		Real value3 = cl * (*uinsf.q)[IIDX::IIW][ug.lc] + cr * (*uinsf.q)[IIDX::IIW][ug.rc];
		Real value4 = cl * (*uinsf.q)[IIDX::IIP][ug.lc] + cr * (*uinsf.q)[IIDX::IIP][ug.rc];

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		(*uinsf.dqdx)[IIDX::IIU][ug.lc] += fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.lc] += fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.lc] += fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.lc] += fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.lc] += fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.lc] += fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.lc] += fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.lc] += fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.lc] += fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.lc] += fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.lc] += fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.lc] += fnza * value4;

		if (ug.fId < ug.nBFace) continue;
		(*uinsf.dqdx)[IIDX::IIU][ug.rc] -= fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.rc] -= fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.rc] -= fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.rc] -= fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.rc] -= fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.rc] -= fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.rc] -= fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.rc] -= fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.rc] -= fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] -= fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.rc] -= fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.rc] -= fnza * value4;
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		Real ovol = one / (*ug.cvol)[cId];
		(*uinsf.dqdx)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIP][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIP][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIP][ug.cId] *= ovol;
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		(*uinsf.dqdx)[IIDX::IIU][ug.rc] = (*uinsf.dqdx)[IIDX::IIU][ug.lc];
		(*uinsf.dqdy)[IIDX::IIU][ug.rc] = (*uinsf.dqdy)[IIDX::IIU][ug.lc];
		(*uinsf.dqdz)[IIDX::IIU][ug.rc] = (*uinsf.dqdz)[IIDX::IIU][ug.lc];
		(*uinsf.dqdx)[IIDX::IIV][ug.rc] = (*uinsf.dqdx)[IIDX::IIV][ug.lc];
		(*uinsf.dqdy)[IIDX::IIV][ug.rc] = (*uinsf.dqdy)[IIDX::IIV][ug.lc];
		(*uinsf.dqdz)[IIDX::IIV][ug.rc] = (*uinsf.dqdz)[IIDX::IIV][ug.lc];
		(*uinsf.dqdx)[IIDX::IIW][ug.rc] = (*uinsf.dqdx)[IIDX::IIW][ug.lc];
		(*uinsf.dqdy)[IIDX::IIW][ug.rc] = (*uinsf.dqdy)[IIDX::IIW][ug.lc];
		(*uinsf.dqdz)[IIDX::IIW][ug.rc] = (*uinsf.dqdz)[IIDX::IIW][ug.lc];
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] = (*uinsf.dqdx)[IIDX::IIP][ug.lc];
		(*uinsf.dqdy)[IIDX::IIP][ug.rc] = (*uinsf.dqdy)[IIDX::IIP][ug.lc];
		(*uinsf.dqdz)[IIDX::IIP][ug.rc] = (*uinsf.dqdz)[IIDX::IIP][ug.lc];
	}
}

void UINsVisterm::CmpVisterm()
{
    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;

        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        if ( fId == 147489 )
        {
            int kkk = 1;
        }
		iinv.ukl[ug.fId] = (*uinsf.q)[IIDX::IIU][ug.lc];
		iinv.ukr[ug.fId] = (*uinsf.q)[IIDX::IIU][ug.rc];
		iinv.vkl[ug.fId] = (*uinsf.q)[IIDX::IIV][ug.lc];
		iinv.vkr[ug.fId] = (*uinsf.q)[IIDX::IIV][ug.rc];
		iinv.wkl[ug.fId] = (*uinsf.q)[IIDX::IIW][ug.lc];
		iinv.wkr[ug.fId] = (*uinsf.q)[IIDX::IIW][ug.rc];

        this->CmpFaceVisterm();  //要改动

    }
}

void UINsVisterm::CmpFaceVisterm()
{

	Real l2rdx = (*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc];  //界面左右单元中心距
	Real l2rdy = (*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc];
	Real l2rdz = (*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc];

	iinv.c2d = sqrt(l2rdx*l2rdx+ l2rdy * l2rdy+ l2rdz * l2rdz);

	iinv.Fnu[ug.fId] = (*ug.farea)[ug.fId] *((*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId]+ (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId]+ (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId])/ iinv.c2d;//μ / (n * d) 法向扩散项系数(改动)
	iinv.Fnv[ug.fId] = (*ug.farea)[ug.fId] * ((*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.c2d;
	iinv.Fnw[ug.fId] = (*ug.farea)[ug.fId] * ((*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.c2d;

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

	iinv.Puf[ug.fId] = (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId];  //▽q*n

	iinv.Pvf[ug.fId] = (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId];
	
	iinv.Pwf[ug.fId] = (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId];
	

	iinv.Pdx[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*l2rdx +
		               (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*l2rdy +
		               (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*l2rdz)/ iinv.c2d;

	iinv.Pdy[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*l2rdx +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*l2rdy +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*l2rdz)/ iinv.c2d;
	
	iinv.Pdz[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*l2rdx +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*l2rdy +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*l2rdz)/ iinv.c2d;


	iinv.Pufd[ug.fId] = iinv.Pdx[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + iinv.Pdx[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Pdx[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId];
	iinv.Pvfd[ug.fId] = iinv.Pdy[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + iinv.Pdy[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Pdy[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId];
	iinv.Pwfd[ug.fId] = iinv.Pdz[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + iinv.Pdz[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Pdz[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId];

	iinv.Ftu[ug.fId] = (iinv.Puf[ug.fId] + iinv.Pufd[ug.fId])*(*ug.farea)[ug.fId];   //扩散项中归入源项的部分
    iinv.Ftv[ug.fId] = (iinv.Pvf[ug.fId] + iinv.Pvfd[ug.fId])*(*ug.farea)[ug.fId];
	iinv.Ftw[ug.fId] = (iinv.Pwf[ug.fId] + iinv.Pwfd[ug.fId])*(*ug.farea)[ug.fId];


	iinv.FuT[ug.fId] = ((iinv.ukr[ug.fId]- iinv.ukl[ug.fId])*(*ug.xfn)[ug.fId]* (*ug.xfn)[ug.fId]+
		               (iinv.vkr[ug.fId] - iinv.vkl[ug.fId])*(*ug.xfn)[ug.fId] * (*ug.yfn)[ug.fId]+
		               (iinv.wkr[ug.fId] - iinv.wkl[ug.fId])*(*ug.xfn)[ug.fId] * (*ug.zfn)[ug.fId])/iinv.c2d*(*ug.farea)[ug.fId];

	iinv.FvT[ug.fId] = ((iinv.ukr[ug.fId] - iinv.ukl[ug.fId])*(*ug.xfn)[ug.fId] * (*ug.yfn)[ug.fId] +
		                (iinv.vkr[ug.fId] - iinv.vkl[ug.fId])*(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] +
		                (iinv.wkr[ug.fId] - iinv.wkl[ug.fId])*(*ug.yfn)[ug.fId] * (*ug.zfn)[ug.fId])/iinv.c2d*(*ug.farea)[ug.fId];

	iinv.FwT[ug.fId] = ((iinv.ukr[ug.fId] - iinv.ukl[ug.fId])*(*ug.xfn)[ug.fId] * (*ug.zfn)[ug.fId] +
		                (iinv.vkr[ug.fId] - iinv.vkl[ug.fId])*(*ug.yfn)[ug.fId] * (*ug.zfn)[ug.fId] +
		                (iinv.wkr[ug.fId] - iinv.wkl[ug.fId])*(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId])/iinv.c2d*(*ug.farea)[ug.fId];


	iinv.PufT[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*(*ug.zfn)[ug.fId]);

	iinv.PvfT[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*(*ug.zfn)[ug.fId]);

	iinv.PwfT[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId]);


	iinv.PufdT[ug.fId] = (iinv.Pdx[ug.fId]* (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId]+ iinv.Pdy[ug.fId]* (*ug.xfn)[ug.fId] * (*ug.yfn)[ug.fId] + 
		                 iinv.Pdz[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.zfn)[ug.fId]);

	iinv.PvfdT[ug.fId] = (iinv.Pdx[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Pdy[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] +
		                  iinv.Pdz[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.zfn)[ug.fId]);

	iinv.PwfdT[ug.fId] = (iinv.Pdx[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.zfn)[ug.fId] + iinv.Pdy[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.zfn)[ug.fId] +
		                  iinv.Pdz[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]);

	iinv.Pud[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.xfn)[ug.fId];

    iinv.Pvd[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.yfn)[ug.fId];

	iinv.Pwd[ug.fId] =((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.zfn)[ug.fId];

	iinv.FtuT[ug.fId] = (iinv.PufT[ug.fId] + iinv.PufdT[ug.fId]+ iinv.Pud[ug.fId])* (*ug.farea)[ug.fId];
	iinv.FtvT[ug.fId] = (iinv.PvfT[ug.fId] + iinv.PvfdT[ug.fId]+ iinv.Pvd[ug.fId])* (*ug.farea)[ug.fId];
	iinv.FtwT[ug.fId] = (iinv.PwfT[ug.fId] + iinv.PwfdT[ug.fId]+ iinv.Pwd[ug.fId])* (*ug.farea)[ug.fId];


	iinv.Fpu[ug.fId] = -half * ((*uinsf.q)[IIDX::IIP][ug.lc] + (*uinsf.q)[IIDX::IIP][ug.rc])*(*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];  //压力梯度项（归入源项）
	iinv.Fpv[ug.fId] = -half * ((*uinsf.q)[IIDX::IIP][ug.lc] + (*uinsf.q)[IIDX::IIP][ug.rc])*(*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
	iinv.Fpw[ug.fId] = -half * ((*uinsf.q)[IIDX::IIP][ug.lc] + (*uinsf.q)[IIDX::IIP][ug.rc])*(*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];



	
	iinv.akku1[ug.fId] = iinv.Fnu[ug.fId];    //该界面上的扩散流
	iinv.akku2[ug.fId] = iinv.Fnu[ug.fId];
	iinv.akkv1[ug.fId] = iinv.Fnv[ug.fId];    //该界面上的扩散流
	iinv.akkv2[ug.fId] = iinv.Fnv[ug.fId];
	iinv.akkw1[ug.fId] = iinv.Fnw[ug.fId];    //该界面上的扩散流
	iinv.akkw2[ug.fId] = iinv.Fnw[ug.fId];

	iinv.aku1[ug.lc] += iinv.Fnu[ug.fId];   //归于动量方程中主对角线系数
	iinv.aku2[ug.rc] += iinv.Fnu[ug.fId];
	iinv.akv1[ug.lc] += iinv.Fnv[ug.fId];   //归于动量方程中主对角线系数
	iinv.akv2[ug.rc] += iinv.Fnv[ug.fId];
	iinv.akw1[ug.lc] += iinv.Fnw[ug.fId];   //归于动量方程中主对角线系数
	iinv.akw2[ug.rc] += iinv.Fnw[ug.fId];

	iinv.bmu1[ug.lc] += (iinv.Ftu[ug.fId] + iinv.Fpu[ug.fId] * ((*ug.xfc)[ug.fId]) + iinv.FtuT[ug.fId] - iinv.FuT[ug.fId]); //界面上归入源项的扩散项
	iinv.bmu2[ug.rc] += (iinv.Ftu[ug.fId] + iinv.Fpu[ug.fId] * (-(*ug.xfc)[ug.fId]) + iinv.FtuT[ug.fId] + iinv.FuT[ug.fId]);

	iinv.bmv1[ug.lc] += (iinv.Ftv[ug.fId] + iinv.Fpv[ug.fId] * ((*ug.yfc)[ug.fId]) + iinv.FtvT[ug.fId] - iinv.FvT[ug.fId]); //界面上归入源项的扩散项
	iinv.bmv2[ug.rc] += (iinv.Ftv[ug.fId] + iinv.Fpv[ug.fId] * (-(*ug.yfc)[ug.fId]) + iinv.FtvT[ug.fId] + iinv.FvT[ug.fId]);

	iinv.bmw1[ug.lc] += (iinv.Ftw[ug.fId] + iinv.Fpw[ug.fId] + iinv.FtwT[ug.fId] - iinv.FwT[ug.fId]); //界面上归入源项的扩散项
	iinv.bmw2[ug.rc] += (iinv.Ftw[ug.fId] + iinv.Fpw[ug.fId] + iinv.FtwT[ug.fId] + iinv.FwT[ug.fId]);



}

void UINsVisterm::CmpUnsteadcoff()
{

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.spt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]/ iinv.timestep;  //矩阵对角线元素的非稳态项
		if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
		{
			iinv.up[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];
			iinv.vp[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
			iinv.wp[ug.cId] = (*uinsf.q)[IIDX::IIW][ug.cId];

			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]*iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]*iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]* iinv.wp[ug.cId]/ iinv.timestep;
		}
		else
		{
			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.wp[ug.cId] / iinv.timestep;
		}
	}
	for (int cId = ug.nCell; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.spt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] / iinv.timestep;  //矩阵对角线元素的非稳态项
		
		if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
		{
			iinv.up[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];
			iinv.vp[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
			iinv.wp[ug.cId] = (*uinsf.q)[IIDX::IIW][ug.cId];

			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.wp[ug.cId] / iinv.timestep;
		}
		else
		{
			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.wp[ug.cId] / iinv.timestep;
		}

		//iinv.spt[ug.cId] = 0;
		//iinv.but[ug.cId] = 0;
		//iinv.bvt[ug.cId] = 0;
		//iinv.bwt[ug.cId] = 0;
	}
}



void UINsVisterm::CmpINsSrc()
{
		for (int cId = 0; cId < ug.nCell; ++cId)
		{
			ug.cId = cId;

			iinv.spu[ug.cId] = iinv.ai1[ug.cId] + iinv.ai2[ug.cId] + iinv.aku1[ug.cId] + iinv.aku2[ug.cId] + iinv.spt[ug.cId]; //矩阵主对角线系数，动量方程单元主系数
			iinv.spv[ug.cId] = iinv.ai1[ug.cId] + iinv.ai2[ug.cId] + iinv.akv1[ug.cId] + iinv.akv2[ug.cId] + iinv.spt[ug.cId];
			iinv.spw[ug.cId] = iinv.ai1[ug.cId] + iinv.ai2[ug.cId] + iinv.akw1[ug.cId] + iinv.akw2[ug.cId] + iinv.spt[ug.cId];


			iinv.buc[ug.cId] = iinv.bmu1[ug.cId] + iinv.bmu2[ug.cId]+ iinv.but[ug.cId];   //动量方程源项
			iinv.bvc[ug.cId] = iinv.bmv1[ug.cId] + iinv.bmv2[ug.cId]+iinv.bvt[ug.cId];
			iinv.bwc[ug.cId] = iinv.bmw1[ug.cId] + iinv.bmw2[ug.cId]+ iinv.bwt[ug.cId];


			int fn = (*ug.c2f)[ug.cId].size();
			for (int iFace = 0; iFace < fn; ++iFace)
			{
				int fId = (*ug.c2f)[ug.cId][iFace];
				ug.fId = fId;
				ug.lc = (*ug.lcf)[ug.fId];
				ug.rc = (*ug.rcf)[ug.fId];
				if (ug.cId == ug.lc)
				{
					iinv.spuj[ug.cId][ug.rc] = (iinv.aii2[ug.fId] + iinv.akku2[ug.fId]);  //矩阵非零系数，动量方程中与主单元相邻的单元面通量
					iinv.spvj[ug.cId][ug.rc] = (iinv.aii2[ug.fId] + iinv.akkv2[ug.fId]);
					iinv.spwj[ug.cId][ug.rc] = (iinv.aii2[ug.fId] + iinv.akkw2[ug.fId]);


					iinv.sju[ug.cId] += iinv.aii2[ug.fId] + iinv.akku2[ug.fId]; //矩阵非零系数相加，用于求压力修正方程
					iinv.sjv[ug.cId] += iinv.aii2[ug.fId] + iinv.akkv2[ug.fId];
					iinv.sjw[ug.cId] += iinv.aii2[ug.fId] + iinv.akkw2[ug.fId];


				}
				else if (ug.cId == ug.rc)
				{
					iinv.spuj[ug.cId][ug.lc] = (iinv.aii1[ug.fId] + iinv.akku1[ug.fId]);
					iinv.spvj[ug.cId][ug.lc] = (iinv.aii1[ug.fId] + iinv.akkv1[ug.fId]);
					iinv.spwj[ug.cId][ug.lc] = (iinv.aii1[ug.fId] + iinv.akkw1[ug.fId]);

					iinv.sju[ug.cId] += iinv.aii1[ug.fId] + iinv.akku1[ug.fId];
					iinv.sjv[ug.cId] += iinv.aii1[ug.fId] + iinv.akkv1[ug.fId];
					iinv.sjw[ug.cId] += iinv.aii1[ug.fId] + iinv.akkw1[ug.fId];

				}
			}
		}

		for (int cId = ug.nCell; cId < ug.nTCell; ++cId)
		{
			ug.cId = cId;

			ug.cId = cId;

			iinv.spu[ug.cId] = iinv.ai1[ug.cId] + iinv.ai2[ug.cId] + iinv.aku1[ug.cId] + iinv.aku2[ug.cId] + iinv.spt[ug.cId]; //矩阵主对角线系数，动量方程单元主系数
			iinv.spv[ug.cId] = iinv.ai1[ug.cId] + iinv.ai2[ug.cId] + iinv.akv1[ug.cId] + iinv.akv2[ug.cId] + iinv.spt[ug.cId];
			iinv.spw[ug.cId] = iinv.ai1[ug.cId] + iinv.ai2[ug.cId] + iinv.akw1[ug.cId] + iinv.akw2[ug.cId] + iinv.spt[ug.cId];


			iinv.buc[ug.cId] = iinv.bmu1[ug.cId] + iinv.bmu2[ug.cId] + iinv.but[ug.cId];   //动量方程源项
			iinv.bvc[ug.cId] = iinv.bmv1[ug.cId] + iinv.bmv2[ug.cId] + iinv.bvt[ug.cId];
			iinv.bwc[ug.cId] = iinv.bmw1[ug.cId] + iinv.bmw2[ug.cId] + iinv.bwt[ug.cId];


			int fn = (*ug.c2f)[ug.cId].size();
			for (int iFace = 0; iFace < fn; ++iFace)
			{
				int fId = (*ug.c2f)[ug.cId][iFace];
				ug.fId = fId;
				ug.lc = (*ug.lcf)[ug.fId];
				ug.rc = (*ug.rcf)[ug.fId];
				if (ug.cId == ug.lc)
				{
					iinv.spuj[ug.cId][ug.rc] = (iinv.aii2[ug.fId] + iinv.akku2[ug.fId]);  //矩阵非零系数，动量方程中与主单元相邻的单元面通量
					iinv.spvj[ug.cId][ug.rc] = (iinv.aii2[ug.fId] + iinv.akkv2[ug.fId]);
					iinv.spwj[ug.cId][ug.rc] = (iinv.aii2[ug.fId] + iinv.akkw2[ug.fId]);


					iinv.sju[ug.cId] += iinv.aii2[ug.fId] + iinv.akku2[ug.fId]; //矩阵非零系数相加，用于求压力修正方程
					iinv.sjv[ug.cId] += iinv.aii2[ug.fId] + iinv.akkv2[ug.fId];
					iinv.sjw[ug.cId] += iinv.aii2[ug.fId] + iinv.akkw2[ug.fId];


				}
				else if (ug.cId == ug.rc)
				{
					iinv.spuj[ug.cId][ug.lc] = (iinv.aii1[ug.fId] + iinv.akku1[ug.fId]);
					iinv.spvj[ug.cId][ug.lc] = (iinv.aii1[ug.fId] + iinv.akkv1[ug.fId]);
					iinv.spwj[ug.cId][ug.lc] = (iinv.aii1[ug.fId] + iinv.akkw1[ug.fId]);

					iinv.sju[ug.cId] += iinv.aii1[ug.fId] + iinv.akku1[ug.fId];
					iinv.sjv[ug.cId] += iinv.aii1[ug.fId] + iinv.akkv1[ug.fId];
					iinv.sjw[ug.cId] += iinv.aii1[ug.fId] + iinv.akkw1[ug.fId];

				}
			}

			//iinv.spu[ug.cId] = 1; //矩阵主对角线系数，动量方程单元主系数
			//iinv.spv[ug.cId] = 1;
			//iinv.spw[ug.cId] = 1;

			//iinv.buc[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];   //动量方程源项
			//iinv.bvc[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
			//iinv.bwc[ug.cId] = (*uinsf.q)[IIDX::IIW][ug.cId];

			//iinv.spuj[ug.cId] = 0;
			//iinv.spvj[ug.cId] = 0;
			//iinv.spwj[ug.cId] = 0;
		}

}

//void UINsVisterm::Addcoff()
//{
    //UnsGrid * grid = Zone::GetUnsGrid();
    //MRField * res = GetFieldPointer< MRField >( grid, "res" );

   // ONEFLOW::AddF2CField( res, visflux );
//}

void UINsVisterm::PrepareFaceValue()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    gcom.CmpTangent();

    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        visQ.dqdx1[ iEqu ] = ( * uinsf.dqdx )[ iEqu ][ ug.lc ];
        visQ.dqdy1[ iEqu ] = ( * uinsf.dqdy )[ iEqu ][ ug.lc ];
        visQ.dqdz1[ iEqu ] = ( * uinsf.dqdz )[ iEqu ][ ug.lc ];

        visQ.dqdx2[ iEqu ] = ( * uinsf.dqdx )[ iEqu ][ ug.rc ];
        visQ.dqdy2[ iEqu ] = ( * uinsf.dqdy )[ iEqu ][ ug.rc ];
        visQ.dqdz2[ iEqu ] = ( * uinsf.dqdz )[ iEqu ][ ug.rc ];
    }

   // for ( int iEqu = 0; iEqu < inscom.nTModel; ++ iEqu )
    //{
   //     visT.dqdx1[ iEqu ] = ( * uinsf.dtdx )[ iEqu ][ ug.lc ];
   //     visT.dqdy1[ iEqu ] = ( * uinsf.dtdy )[ iEqu ][ ug.lc ];
   //     visT.dqdz1[ iEqu ] = ( * uinsf.dtdz )[ iEqu ][ ug.lc ];

//        visT.dqdx2[ iEqu ] = ( * uinsf.dtdx )[ iEqu ][ ug.rc ];
//        visT.dqdy2[ iEqu ] = ( * uinsf.dtdy )[ iEqu ][ ug.rc ];
//        visT.dqdz2[ iEqu ] = ( * uinsf.dtdz )[ iEqu ][ ug.rc ];
//    }

    inscom.visl1 = ( * uinsf.visl )[ 0 ][ ug.lc ];
    inscom.visl2 = ( * uinsf.visl )[ 0 ][ ug.rc ];

   inscom.vist1 = ( * uinsf.vist )[ 0 ][ ug.lc ];
    inscom.vist2 = ( * uinsf.vist )[ 0 ][ ug.rc ];

    inscom.visl = half * ( inscom.visl1 + inscom.visl2 );
    inscom.vist = half * ( inscom.vist1 + inscom.vist2 );
    inscom.vis  = inscom.visl + inscom.vist;

    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        visQ.q1[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.lc ];
        visQ.q2[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        visQ.q11[ iEqu ] = visQ.q1[ iEqu ];
        visQ.q22[ iEqu ] = visQ.q2[ iEqu ];
    }

   // for ( int iEqu = 0; iEqu < inscom.nTModel; ++ iEqu )
   // {
   //     visT.q1[ iEqu ] = ( * uinsf.tempr )[ iEqu ][ ug.lc ];
   //     visT.q2[ iEqu ] = ( * uinsf.tempr )[ iEqu ][ ug.rc ];
   // }

    for ( int iEqu = 0; iEqu < inscom.nTModel; ++ iEqu )
    {
        visT.q11[ iEqu ] = visT.q1[ iEqu ];
        visT.q22[ iEqu ] = visT.q2[ iEqu ];
    }

    this->AverGrad();
    this->CmpFaceWeight();
    this->SaveFacePara();
}

void UINsVisterm::SaveFacePara()
{
    Ivis.dudx  = visQ.dqdx[ IIDX::IIU ];
    Ivis.dudy  = visQ.dqdy[ IIDX::IIU ];
    Ivis.dudz  = visQ.dqdz[ IIDX::IIU ];

    Ivis.dvdx  = visQ.dqdx[ IIDX::IIV ];
    Ivis.dvdy  = visQ.dqdy[ IIDX::IIV ];
    Ivis.dvdz  = visQ.dqdz[ IIDX::IIV ];

    Ivis.dwdx  = visQ.dqdx[ IIDX::IIW ];
    Ivis.dwdy  = visQ.dqdy[ IIDX::IIW ];
    Ivis.dwdz  = visQ.dqdz[ IIDX::IIW ];

	Ivis.dpdx = visQ.dqdx[IIDX::IIP];
	Ivis.dpdy = visQ.dqdy[IIDX::IIP];
	Ivis.dpdz = visQ.dqdz[IIDX::IIP];

	Ivis.p1 = visQ.q1[IIDX::IIP];
	Ivis.p2 = visQ.q2[IIDX::IIP];

    Ivis.um  = visQ.q[ IIDX::IIU ];
    Ivis.vm  = visQ.q[ IIDX::IIV ];
    Ivis.wm  = visQ.q[ IIDX::IIW ];

    //Ivis.dtdn = visT.dqdn[ IIDX::IITT ];
    //Ivis.tmid = visT.q[ IIDX::IITT ];
}

void UINsVisterm::CmpFaceWeight()
{
    vgg.CmpFaceWeight();
}


void UINsVisterm::CmpGradCoef()
{
    vgg.CmpGradCoef();
}


void UINsVisterm::PrepareCellGeom()
{
    vgg.PrepareCellGeom();
}


void ICmpLaminarViscosity(int flag)
{
		ug.Init();
		uinsf.Init();
		ug.SetStEd(flag);

		Real minLimit = 0.0;

		for (int cId = ug.ist; cId < ug.ied; ++cId)
		{
			Real temperature = ( *uinsf.tempr )[ IIDX::IITT ][ cId ];
			Real visl = Iutherland::ICmpViscosity( temperature );
			//( *uinsf.visl )[ 0 ][ cId ] = MAX( minLimit, visl );
			(*uinsf.visl)[0][cId] = 0;
		}
}

EndNameSpace

