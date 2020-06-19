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
#include "UINsLimiter.h"
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

		//Real value1 = cl * (*uinsf.q)[IIDX::IIU][ug.lc] + cr * (*uinsf.q)[IIDX::IIU][ug.rc];
		//Real value2 = cl * (*uinsf.q)[IIDX::IIV][ug.lc] + cr * (*uinsf.q)[IIDX::IIV][ug.rc];
		//Real value3 = cl * (*uinsf.q)[IIDX::IIW][ug.lc] + cr * (*uinsf.q)[IIDX::IIW][ug.rc];
		//Real value4 = cl * (*uinsf.q)[IIDX::IIP][ug.lc] + cr * (*uinsf.q)[IIDX::IIP][ug.rc];

		Real value1 = 0.5 * (*uinsf.q)[IIDX::IIU][ug.lc] + 0.5 * (*uinsf.q)[IIDX::IIU][ug.rc];
		Real value2 = 0.5 * (*uinsf.q)[IIDX::IIV][ug.lc] + 0.5 * (*uinsf.q)[IIDX::IIV][ug.rc];
		Real value3 = 0.5 * (*uinsf.q)[IIDX::IIW][ug.lc] + 0.5 * (*uinsf.q)[IIDX::IIW][ug.rc];
		Real value4 = 0.5 * (*uinsf.q)[IIDX::IIP][ug.lc] + 0.5 * (*uinsf.q)[IIDX::IIP][ug.rc];

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
		(*uinsf.dqdx)[IIDX::IIU][ug.rc] += -fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.rc] += -fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.rc] += -fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.rc] += -fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.rc] += -fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.rc] += -fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.rc] += -fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.rc] += -fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.rc] += -fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] += -fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.rc] += -fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.rc] += -fnza * value4;
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		Real ovol = one / (*ug.cvol)[ug.cId];
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
		//if (ug.rc > ug.nCell)
		//{
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
		//}

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
		
		//iinv.ukl[ug.fId] = (*limf->qf1)[IIDX::IIU][ug.fId];
		//iinv.ukr[ug.fId] = (*limf->qf2)[IIDX::IIU][ug.fId];
		//iinv.vkl[ug.fId] = (*limf->qf1)[IIDX::IIV][ug.fId];
		//iinv.vkr[ug.fId] = (*limf->qf2)[IIDX::IIV][ug.fId];
		//iinv.wkl[ug.fId] = (*limf->qf1)[IIDX::IIW][ug.fId];
		//iinv.wkr[ug.fId] = (*limf->qf2)[IIDX::IIW][ug.fId];

        this->CmpFaceVisterm();  //要改动

    }
}

void UINsVisterm::CmpFaceVisterm()
{

	iinv.l2rdx[ug.fId] = (*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc];  //界面左右单元中心距
	iinv.l2rdy[ug.fId] = (*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc];
	iinv.l2rdz[ug.fId] = (*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc];

	iinv.c2d = sqrt(iinv.l2rdx[ug.fId] * iinv.l2rdx[ug.fId] + iinv.l2rdy[ug.fId] * iinv.l2rdy[ug.fId] + iinv.l2rdz[ug.fId] * iinv.l2rdz[ug.fId]);

	//iinv.visu[ug.fId] = 1 / inscom.reynolds;  //动力粘度
	//iinv.visv[ug.fId] = 1 / inscom.reynolds;
	//iinv.visw[ug.fId] = 1 / inscom.reynolds;

	iinv.vis = 1 / inscom.reynolds;  //动力粘度

	iinv.dist[ug.fId] = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	//iinv.Fnu[ug.fId] = iinv.visu[ug.fId] *(*ug.farea)[ug.fId] *((*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId]
	//	               + (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId]+ (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId])/ iinv.c2d;  //归入系数矩阵的部分
	//iinv.Fnv[ug.fId] = iinv.visv[ug.fId] * (*ug.farea)[ug.fId] * ((*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] 
	//	               + (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.c2d;
	//iinv.Fnw[ug.fId] = iinv.visw[ug.fId] * (*ug.farea)[ug.fId] * ((*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] 
	//	               + (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.c2d;


	iinv.Fn[ug.fId] = iinv.vis * (*ug.farea)[ug.fId] / iinv.dist[ug.fId];
	//iinv.Fnu[ug.fId] = iinv.visu[ug.fId] *(*ug.farea)[ug.fId] / iinv.dist[ug.fId];  //归入系数矩阵的部分
	//iinv.Fnv[ug.fId] = iinv.visv[ug.fId] * (*ug.farea)[ug.fId] / iinv.dist[ug.fId];
	//iinv.Fnw[ug.fId] = iinv.visw[ug.fId] * (*ug.farea)[ug.fId] / iinv.dist[ug.fId];

	//if ((*ug.xfn)[ug.fId] == 0)
	//{
	//	iinv.Fnu[ug.fId] = 0;
	//	iinv.Fnv[ug.fId] = 0;
	//	iinv.Fnw[ug.fId] = 0;
	//}
	//else
	//{
	//	iinv.Fnu[ug.fId] = iinv.visu[ug.fId] * (*ug.farea)[ug.fId] * (*ug.xfn)[ug.fId]/abs(iinv.l2rdx[ug.fId]);
	//	iinv.Fnv[ug.fId] = iinv.visv[ug.fId] * (*ug.farea)[ug.fId] * (*ug.xfn)[ug.fId]/abs(iinv.l2rdx[ug.fId]);
	//	iinv.Fnw[ug.fId] = iinv.visw[ug.fId] * (*ug.farea)[ug.fId] * (*ug.xfn)[ug.fId] / abs(iinv.l2rdx[ug.fId]);
	//}
	//if ((*ug.yfn)[ug.fId] == 0)
	//{
	//	iinv.Fnu[ug.fId] += 0;
	//	iinv.Fnv[ug.fId] += 0;
	//	iinv.Fnw[ug.fId] += 0;
	//}
	//else
	//{
	//	iinv.Fnu[ug.fId] += iinv.visu[ug.fId] * (*ug.farea)[ug.fId] * (*ug.yfn)[ug.fId] / abs(iinv.l2rdy[ug.fId]);
	//	iinv.Fnv[ug.fId] += iinv.visv[ug.fId] * (*ug.farea)[ug.fId] * (*ug.yfn)[ug.fId] / abs(iinv.l2rdy[ug.fId]);
	//	iinv.Fnw[ug.fId] += iinv.visw[ug.fId] * (*ug.farea)[ug.fId] * (*ug.yfn)[ug.fId] / abs(iinv.l2rdy[ug.fId]);
	//}
	//if ((*ug.zfn)[ug.fId] == 0)
	//{
	//	iinv.Fnu[ug.fId] += 0;
	//	iinv.Fnv[ug.fId] += 0;
	//	iinv.Fnw[ug.fId] += 0;
	//}
	//else
	//{
	//	iinv.Fnu[ug.fId] += iinv.visu[ug.fId] * (*ug.farea)[ug.fId] * (*ug.zfn)[ug.fId] / abs(iinv.l2rdz[ug.fId]);
	//	iinv.Fnv[ug.fId] += iinv.visv[ug.fId] * (*ug.farea)[ug.fId] * (*ug.zfn)[ug.fId] / abs(iinv.l2rdz[ug.fId]);
	//	iinv.Fnw[ug.fId] += iinv.visw[ug.fId] * (*ug.farea)[ug.fId] * (*ug.zfn)[ug.fId] / abs(iinv.l2rdz[ug.fId]);
	//}

	//iinv.Fnu[ug.fId] = iinv.visu[ug.fId] *(*ug.farea)[ug.fId] *((*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId]
		               //+ (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId]+ (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId])/ iinv.c2d;  //归入系数矩阵的部分
	//iinv.Fnv[ug.fId] = iinv.visv[ug.fId] * (*ug.farea)[ug.fId] * ((*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] 
		               //+ (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.c2d;
	//iinv.Fnw[ug.fId] = iinv.visw[ug.fId] * (*ug.farea)[ug.fId] * ((*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] 
		               //+ (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.c2d;

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
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*(*ug.zfn)[ug.fId];  //▽q*n

	iinv.Pvf[ug.fId] = (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*(*ug.zfn)[ug.fId];
	
	iinv.Pwf[ug.fId] = (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId];
	
	iinv.Pdu[ug.fId] = -((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*iinv.l2rdx[ug.fId]+
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*iinv.l2rdy[ug.fId]+
	 	                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*iinv.l2rdz[ug.fId])/ iinv.dist[ug.fId];

	iinv.Pdv[ug.fId] = -((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*iinv.l2rdx[ug.fId]  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*iinv.l2rdy[ug.fId]  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*iinv.l2rdz[ug.fId] ) / iinv.dist[ug.fId];

	iinv.Pdw[ug.fId] = -((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*iinv.l2rdx[ug.fId]  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*iinv.l2rdy[ug.fId]  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*iinv.l2rdz[ug.fId] ) / iinv.dist[ug.fId];

	//iinv.Pdu[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*iinv.l2rdx[ug.fId] *(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] +
	//	               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*iinv.l2rdy[ug.fId] *(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] +
	//	               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*iinv.l2rdz[ug.fId] *(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId])/ iinv.c2d;

	//iinv.Pdv[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*iinv.l2rdx[ug.fId] *(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] +
	//	               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*iinv.l2rdy[ug.fId] *(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] +
	//	               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*iinv.l2rdz[ug.fId] *(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId])/ iinv.c2d;
	
	//iinv.Pdw[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*iinv.l2rdx[ug.fId] *(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] +
	//	               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*iinv.l2rdy[ug.fId] *(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] +
	//	               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*iinv.l2rdz[ug.fId] *(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId])/ iinv.c2d;

	//iinv.Puf[ug.fId] = (half * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		//(half * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*(*ug.yfn)[ug.fId] +
		//(half * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*(*ug.zfn)[ug.fId];  //▽q*n

	//iinv.Pvf[ug.fId] = (half * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*(*ug.xfn)[ug.fId] +
		//(half * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		//(half * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*(*ug.zfn)[ug.fId];

	//iinv.Pwf[ug.fId] = (half * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*(*ug.xfn)[ug.fId] +
		//(half * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*(*ug.yfn)[ug.fId] +
		//(half * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId];


	//iinv.Pdu[ug.fId] = ((half * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*iinv.l2rdx[ug.fId]*(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] +
		//(half * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*iinv.l2rdy[ug.fId] *(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] +
		//(half * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*iinv.l2rdz[ug.fId] *(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.c2d;

	//iinv.Pdv[ug.fId] = ((half * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*iinv.l2rdx[ug.fId] *(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] +
		//(half * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*iinv.l2rdy[ug.fId] *(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] +
		//(half * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*iinv.l2rdz[ug.fId] *(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.c2d;

	//iinv.Pdw[ug.fId] = ((half * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*iinv.l2rdx[ug.fId] *(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] +
		//(half * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*iinv.l2rdy[ug.fId] *(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] +
		//(half * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*iinv.l2rdz[ug.fId] *(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.c2d;


	iinv.Ftu1[ug.fId] = iinv.Puf[ug.fId] *(*ug.farea)[ug.fId]*iinv.visu[ug.fId];   //扩散项中归入源项的部分1
    iinv.Ftv1[ug.fId] = iinv.Pvf[ug.fId] *(*ug.farea)[ug.fId]*iinv.visv[ug.fId];
	iinv.Ftw1[ug.fId] = iinv.Pwf[ug.fId] *(*ug.farea)[ug.fId]*iinv.visw[ug.fId];

	iinv.Ftu2[ug.fId] = iinv.Pdu[ug.fId]*(*ug.farea)[ug.fId] * iinv.visu[ug.fId];   //扩散项中归入源项的部分2
	iinv.Ftv2[ug.fId] = iinv.Pdv[ug.fId]*(*ug.farea)[ug.fId] * iinv.visv[ug.fId];
	iinv.Ftw2[ug.fId] = iinv.Pdw[ug.fId]*(*ug.farea)[ug.fId] * iinv.visw[ug.fId];



	iinv.PufT[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.zfn)[ug.fId]);

	iinv.PvfT[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId]);

	iinv.PwfT[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId]);



	iinv.Pud[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.xfn)[ug.fId];

    iinv.Pvd[ug.fId] = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.yfn)[ug.fId];

	iinv.Pwd[ug.fId] =((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.zfn)[ug.fId];


	//iinv.PufT[ug.fId] = ((half * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		//(half * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*(*ug.yfn)[ug.fId] +
		//(half * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*(*ug.zfn)[ug.fId]);

	//iinv.PvfT[ug.fId] = ((half * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*(*ug.xfn)[ug.fId] +
		//(half * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		//(half * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*(*ug.zfn)[ug.fId]);

	//iinv.PwfT[ug.fId] = ((half * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*(*ug.xfn)[ug.fId] +
		//(half * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*(*ug.yfn)[ug.fId] +
		//(half * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId]);



	//iinv.Pud[ug.fId] = ((half * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		//(half * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		//(half * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.xfn)[ug.fId];

	//iinv.Pvd[ug.fId] = ((half * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		//(half * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		//(half * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.yfn)[ug.fId];

	//iinv.Pwd[ug.fId] = ((half * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + half * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		//(half * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + half * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		//(half * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + half * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.zfn)[ug.fId];


	//iinv.Fu1[ug.fId] = (-2 / 3)*iinv.visu[ug.fId] * (iinv.Pud[ug.fId])*(*ug.farea)[ug.fId];  //λ(▽*V)，表面源项
	//iinv.Fv1[ug.fId] = (-2 / 3)*iinv.visv[ug.fId] * (iinv.Pvd[ug.fId])*(*ug.farea)[ug.fId];
	//iinv.Fw1[ug.fId] = (-2 / 3)*iinv.visw[ug.fId] * (iinv.Pwd[ug.fId])*(*ug.farea)[ug.fId];


	//if (iinv.l2rdx[ug.fId] == 0)
	//{
	//	iinv.FtuT[ug.fId] = 0;
	//	iinv.FtvT[ug.fId] = 0;
	//	iinv.FtwT[ug.fId] = 0;
	//}
	//else
	//{
	//	iinv.FtuT[ug.fId] = (iinv.ukr[ug.fId] - iinv.ukl[ug.fId])*(*ug.xfn)[ug.fId] / abs(iinv.l2rdx[ug.fId]);
	//	iinv.FtvT[ug.fId] = (iinv.vkr[ug.fId] - iinv.vkl[ug.fId])*(*ug.xfn)[ug.fId] / abs(iinv.l2rdx[ug.fId]);
	//	iinv.FtwT[ug.fId] = (iinv.wkr[ug.fId] - iinv.wkl[ug.fId])*(*ug.xfn)[ug.fId] / abs(iinv.l2rdx[ug.fId]);
	//}

	//if (iinv.l2rdy[ug.fId] == 0)
	//{
	//	iinv.FtuT[ug.fId] += 0;
	//	iinv.FtvT[ug.fId] += 0;
	//	iinv.FtwT[ug.fId] += 0;
	//}
	//else
	//{
	//	iinv.FtuT[ug.fId] += (iinv.ukr[ug.fId] - iinv.ukl[ug.fId])*(*ug.yfn)[ug.fId] / abs(iinv.l2rdy[ug.fId]);
	//	iinv.FtvT[ug.fId] += (iinv.vkr[ug.fId] - iinv.vkl[ug.fId])*(*ug.yfn)[ug.fId] / abs(iinv.l2rdy[ug.fId]);
	//	iinv.FtwT[ug.fId] += (iinv.wkr[ug.fId] - iinv.wkl[ug.fId])*(*ug.yfn)[ug.fId] / abs(iinv.l2rdy[ug.fId]);
	//}

	//if (iinv.l2rdz[ug.fId] == 0)
	//{
	//	iinv.FtuT[ug.fId] += 0;
	//	iinv.FtvT[ug.fId] += 0;
	//	iinv.FtwT[ug.fId] += 0;
	//}
	//else
	//{
	//	iinv.FtuT[ug.fId] += (iinv.ukr[ug.fId] - iinv.ukl[ug.fId])*(*ug.zfn)[ug.fId] / abs(iinv.l2rdz[ug.fId]);
	//	iinv.FtvT[ug.fId] += (iinv.vkr[ug.fId] - iinv.vkl[ug.fId])*(*ug.zfn)[ug.fId] / abs(iinv.l2rdz[ug.fId]);
	//	iinv.FtwT[ug.fId] += (iinv.wkr[ug.fId] - iinv.wkl[ug.fId])*(*ug.zfn)[ug.fId] / abs(iinv.l2rdz[ug.fId]);
	//}

	iinv.FtuT[ug.fId] = iinv.PufT[ug.fId] * (*ug.farea)[ug.fId] * iinv.visu[ug.fId];  //Г(▽V)T，表面源项
	iinv.FtvT[ug.fId] = iinv.PvfT[ug.fId] * (*ug.farea)[ug.fId] * iinv.visv[ug.fId];
	iinv.FtwT[ug.fId] = iinv.PwfT[ug.fId] * (*ug.farea)[ug.fId] * iinv.visw[ug.fId];



	//iinv.Fpu[ug.fId] = -half * ((*uinsf.q)[IIDX::IIP][ug.lc] + (*uinsf.q)[IIDX::IIP][ug.rc])*(*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];  //压力梯度项（归入源项）
	//iinv.Fpv[ug.fId] = -half * ((*uinsf.q)[IIDX::IIP][ug.lc] + (*uinsf.q)[IIDX::IIP][ug.rc])*(*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
	//iinv.Fpw[ug.fId] = -half * ((*uinsf.q)[IIDX::IIP][ug.lc] + (*uinsf.q)[IIDX::IIP][ug.rc])*(*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

	iinv.ai[0][ug.fId] += iinv.Fn[ug.fId];
	iinv.ai[1][ug.fId] += iinv.Fn[ug.fId];
	//iinv.akku1[ug.fId] = iinv.Fnu[ug.fId];    //该界面上的扩散流
	//iinv.akku2[ug.fId] = iinv.Fnu[ug.fId];
	//iinv.akkv1[ug.fId] = iinv.Fnv[ug.fId];    //该界面上的扩散流
	//iinv.akkv2[ug.fId] = iinv.Fnv[ug.fId];
	//iinv.akkw1[ug.fId] = iinv.Fnw[ug.fId];    //该界面上的扩散流
	//iinv.akkw2[ug.fId] = iinv.Fnw[ug.fId];

	//iinv.aku1[ug.lc] += iinv.Fnu[ug.fId];   //归于动量方程中主对角线系数
	//iinv.aku2[ug.rc] += iinv.Fnu[ug.fId];
	//iinv.akv1[ug.lc] += iinv.Fnv[ug.fId];   //归于动量方程中主对角线系数
	//iinv.akv2[ug.rc] += iinv.Fnv[ug.fId];
	//iinv.akw1[ug.lc] += iinv.Fnw[ug.fId];   //归于动量方程中主对角线系数
	//iinv.akw2[ug.rc] += iinv.Fnw[ug.fId];


	//iinv.bmu1[ug.lc] += iinv.Ftu1[ug.fId]+iinv.Ftu2[ug.fId] + iinv.FtuT[ug.fId];// +iinv.Fpu[ug.fId];//  -iinv.FuT[ug.fId]); //界面上归入源项的扩散项
	//iinv.bmu2[ug.rc] += -iinv.Ftu1[ug.fId] -iinv.Ftu2[ug.fId] - iinv.FtuT[ug.fId];// -iinv.Fpu[ug.fId]; // ;//  + iinv.FuT[ug.fId]);

	//iinv.bmv1[ug.lc] += iinv.Ftv1[ug.fId] +iinv.Ftv2[ug.fId] + iinv.FtvT[ug.fId];// +iinv.Fpv[ug.fId];// ;//  - iinv.FvT[ug.fId]); //界面上归入源项的扩散项
	//iinv.bmv2[ug.rc] += -iinv.Ftv1[ug.fId] -iinv.Ftv2[ug.fId] - iinv.FtvT[ug.fId];// -iinv.Fpv[ug.fId];// +iinv.Fpv[ug.fId]);// + iinv.FvT[ug.fId]);

	//iinv.bmw1[ug.lc] += iinv.Ftw1[ug.fId]+iinv.Ftw2[ug.fId] + iinv.FtwT[ug.fId]; //+iinv.Fpw[ug.fId];// +iinv.Fpw[ug.fId]);//  - iinv.FwT[ug.fId]); //界面上归入源项的扩散项
	//iinv.bmw2[ug.rc] += -iinv.Ftw1[ug.fId]-iinv.Ftw2[ug.fId] - iinv.FtwT[ug.fId]; //-iinv.Fpw[ug.fId];// );//  + iinv.FwT[ug.fId]);

	
	iinv.biu[0][ug.fId] = iinv.Ftu1[ug.fId] + iinv.Ftu2[ug.fId];
	iinv.biu[1][ug.fId] = -iinv.Ftu1[ug.fId] - iinv.Ftu2[ug.fId];

	iinv.biv[0][ug.fId] = iinv.Ftv1[ug.fId] + iinv.Ftv2[ug.fId];
	iinv.biv[1][ug.fId] = -iinv.Ftv1[ug.fId] - iinv.Ftv2[ug.fId];

	iinv.biw[0][ug.fId] = iinv.Ftw1[ug.fId] + iinv.Ftw2[ug.fId];
	iinv.biw[1][ug.fId] = -iinv.Ftw1[ug.fId] - iinv.Ftw2[ug.fId];

	//iinv.bmu1[ug.lc] += iinv.Ftu1[ug.fId]+ iinv.Ftu2[ug.fId];// +iinv.Fpu[ug.fId];//  -iinv.FuT[ug.fId]); //界面上归入源项的扩散项
	//iinv.bmu2[ug.rc] += -iinv.Ftu1[ug.fId] - iinv.Ftu2[ug.fId];// -iinv.Fpu[ug.fId]; // ;//  + iinv.FuT[ug.fId]);

	//iinv.bmv1[ug.lc] += iinv.Ftv1[ug.fId] + iinv.Ftv2[ug.fId];// +iinv.Fpv[ug.fId];// ;//  - iinv.FvT[ug.fId]); //界面上归入源项的扩散项
	//iinv.bmv2[ug.rc] += -iinv.Ftv1[ug.fId] - iinv.Ftv2[ug.fId];// -iinv.Fpv[ug.fId];// +iinv.Fpv[ug.fId]);// + iinv.FvT[ug.fId]);

	//iinv.bmw1[ug.lc] += iinv.Ftw1[ug.fId] + iinv.Ftw2[ug.fId]; //+iinv.Fpw[ug.fId];// +iinv.Fpw[ug.fId]);//  - iinv.FwT[ug.fId]); //界面上归入源项的扩散项
	//iinv.bmw2[ug.rc] += -iinv.Ftw1[ug.fId] - iinv.Ftw2[ug.fId]; //-iinv.Fpw[ug.fId];// );//  + iinv.FwT[ug.fId]);
}

void UINsVisterm::CmpUnsteadcoff()
{

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.spt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]/ iinv.timestep;  //矩阵对角线元素的非稳态项
		if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
		{
			//iinv.up[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];
			//iinv.vp[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
			//iinv.wp[ug.cId] = (*uinsf.q)[IIDX::IIW][ug.cId];

			iinv.up[ug.cId] = 0;
			iinv.vp[ug.cId] = 0;
			iinv.wp[ug.cId] = 0;

			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]*iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]*iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]* iinv.wp[ug.cId]/ iinv.timestep;

			//iinv.but[ug.cId] = 0; //源项的非稳态项
			//iinv.bvt[ug.cId] = 0;
			//iinv.bwt[ug.cId] = 0;
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
			//iinv.up[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];
			//iinv.vp[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
			//iinv.wp[ug.cId] = (*uinsf.q)[IIDX::IIW][ug.cId];

			iinv.up[ug.cId] = 0;
			iinv.vp[ug.cId] = 0;
			iinv.wp[ug.cId] = 0;

			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]*iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]*iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]* iinv.wp[ug.cId]/ iinv.timestep;

			//iinv.but[ug.cId] = 0; //源项的非稳态项
			//iinv.bvt[ug.cId] = 0;
			//iinv.bwt[ug.cId] = 0;
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
	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.spc[ug.lc] += iinv.ai[0][ug.fId];
		iinv.spc[ug.rc] += iinv.ai[1][ug.fId];

		iinv.buc[ug.lc] += iinv.biu[0][ug.fId];
		iinv.buc[ug.rc] += iinv.biu[1][ug.fId];

		iinv.bvc[ug.lc] += iinv.biv[0][ug.fId];
		iinv.bvc[ug.rc] += iinv.biv[1][ug.fId];

		iinv.bwc[ug.lc] += iinv.biw[0][ug.fId];
		iinv.bwc[ug.rc] += iinv.biw[1][ug.fId];

	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.spc[ug.cId] += iinv.spt[ug.cId];

		iinv.buc[ug.cId] += iinv.but[ug.cId];
		iinv.bvc[ug.cId] += iinv.bvt[ug.cId];
		iinv.bwc[ug.cId] += iinv.bwt[ug.cId];

		int fn = (*ug.c2f)[ug.cId].size();
		iinv.sj.resize(ug.nTCell, fn);
		iinv.sd.resize(ug.nTCell, fn);
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			if (ug.cId == ug.lc)
			{
				iinv.sj[ug.cId][iFace] = iinv.ai[0][ug.fId];  //矩阵非零系数，动量方程中与主单元相邻的单元面通量
				iinv.sd[ug.cId][iFace] = ug.rc;
			}
			else if (ug.cId == ug.rc)
			{
				iinv.sj[ug.cId][iFace] = iinv.ai[1][ug.fId];  //矩阵非零系数，动量方程中与主单元相邻的单元面通量
				iinv.sd[ug.cId][iFace] = ug.lc;
			}

		}
	}
}



EndNameSpace

