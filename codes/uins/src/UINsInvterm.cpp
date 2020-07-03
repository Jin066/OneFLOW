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
#include <iomanip>
using namespace std;

BeginNameSpace(ONEFLOW)

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
	if (inscom.icmpInv == 0) return;
	iinv.Init();
	ug.Init();
	uinsf.Init();
	//Alloc();

	//this->CmpInvFace();
	this->CmpInvMassFlux();  //需要改动

   //DeAlloc();
}

void UINsInvterm::CmpINsTimestep()
{
	iinv.timestep = GetDataValue< Real >("global_dt");
}
void UINsInvterm::CmpINsPreflux()
{
	if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
	{
		if (inscom.icmpInv == 0) return;
		iinv.Init();
		ug.Init();
		uinsf.Init();
		//Alloc();

		this->CmpInvFace();
		this->INsPreflux();
	}

}

void UINsInvterm::INsPreflux()
{
	this->Initflux();
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareFaceValue();

		this->CmpINsinvFlux();

	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareFaceValue();

		this->CmpINsBcinvFlux();
	}

}
void UINsInvterm::Initflux()
{
	iinv.ai1.resize(ug.nTCell);
	iinv.ai2.resize(ug.nTCell);
	iinv.aii1.resize(ug.nFace);
	iinv.aii2.resize(ug.nFace);
	iinv.akku1.resize(ug.nFace);
	iinv.akku2.resize(ug.nFace);
	iinv.akkv1.resize(ug.nFace);
	iinv.akkv2.resize(ug.nFace);
	iinv.akkw1.resize(ug.nFace);
	iinv.akkw2.resize(ug.nFace);
	iinv.aku1.resize(ug.nTCell);
	iinv.aku2.resize(ug.nTCell);
	iinv.akv1.resize(ug.nTCell);
	iinv.akv2.resize(ug.nTCell);
	iinv.akw1.resize(ug.nTCell);
	iinv.akw2.resize(ug.nTCell);
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
	iinv.bp.resize(ug.nTCell);
	iinv.ajp.resize(ug.nFace);
	iinv.sppu.resize(ug.nFace);
	iinv.sppv.resize(ug.nFace);
	iinv.sppw.resize(ug.nFace);
	iinv.sju.resize(ug.nTCell);
	iinv.sjv.resize(ug.nTCell);
	iinv.sjw.resize(ug.nTCell);
	iinv.fq.resize(ug.nFace);
	iinv.spc.resize(ug.nTCell);

	iinv.spu.resize(ug.nTCell);
	iinv.spv.resize(ug.nTCell);
	iinv.spw.resize(ug.nTCell);
	iinv.spuj.resize(ug.nTCell, ug.nTCell);
	iinv.spvj.resize(ug.nTCell, ug.nTCell);
	iinv.spwj.resize(ug.nTCell, ug.nTCell);
	/*iinv.sjp.resize(ug.nTCell, ug.nTCell);*/
	iinv.ai.resize(2, ug.nFace);
	iinv.biu.resize(2, ug.nFace);
	iinv.biv.resize(2, ug.nFace);
	iinv.biw.resize(2, ug.nFace);

	//iinv.ajp.resize(ug.nFace);
	//iinv.app.resize(ug.nCell);
	iinv.spp.resize(ug.nTCell);
	iinv.pp.resize(ug.nTCell);
	iinv.pp0.resize(ug.nTCell);
	iinv.pc.resize(ug.nTCell);
	iinv.pp1.resize(ug.nTCell);
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
	iinv.ump.resize(ug.nTCell);
	iinv.vmp.resize(ug.nTCell);
	iinv.wmp.resize(ug.nTCell);
	iinv.ppr.resize(ug.nFace);
	iinv.dqqdx.resize(ug.nTCell);
	iinv.dqqdy.resize(ug.nTCell);
	iinv.dqqdz.resize(ug.nTCell);
	iinv.fux.resize(ug.nFace);
	iinv.Fn.resize(ug.nFace);
	iinv.Fnu.resize(ug.nFace);
	iinv.Fnv.resize(ug.nFace);
	iinv.Fnw.resize(ug.nFace);
	iinv.Ftu1.resize(ug.nFace);
	iinv.Ftv1.resize(ug.nFace);
	iinv.Ftw1.resize(ug.nFace);
	iinv.Ftu2.resize(ug.nFace);
	iinv.Ftv2.resize(ug.nFace);
	iinv.Ftw2.resize(ug.nFace);
	iinv.ukl.resize(ug.nFace);
	iinv.ukr.resize(ug.nFace);
	iinv.vkl.resize(ug.nFace);
	iinv.vkr.resize(ug.nFace);
	iinv.wkl.resize(ug.nFace);
	iinv.wkr.resize(ug.nFace);
	iinv.uml.resize(ug.nFace);
	iinv.umr.resize(ug.nFace);
	iinv.vml.resize(ug.nFace);
	iinv.vmr.resize(ug.nFace);
	iinv.wml.resize(ug.nFace);
	iinv.wmr.resize(ug.nFace);
	iinv.Puf.resize(ug.nFace);
	iinv.Pvf.resize(ug.nFace);
	iinv.Pwf.resize(ug.nFace);
	iinv.Pufd.resize(ug.nFace);
	iinv.Pvfd.resize(ug.nFace);
	iinv.Pwfd.resize(ug.nFace);
	iinv.Fpu.resize(ug.nFace);
	iinv.Fpv.resize(ug.nFace);
	iinv.Fpw.resize(ug.nFace);
	iinv.tux.resize(ug.nFace);
	iinv.tvy.resize(ug.nFace);
	iinv.twz.resize(ug.nFace);
	iinv.Fqu.resize(ug.nFace);
	iinv.Fqv.resize(ug.nFace);
	iinv.Fqw.resize(ug.nFace);
	iinv.Pdu.resize(ug.nFace);
	iinv.Pdv.resize(ug.nFace);
	iinv.Pdw.resize(ug.nFace);
	iinv.FuT.resize(ug.nFace);
	iinv.FvT.resize(ug.nFace);
	iinv.FwT.resize(ug.nFace);
	iinv.PufT.resize(ug.nFace);
	iinv.PvfT.resize(ug.nFace);
	iinv.PwfT.resize(ug.nFace);
	iinv.PufdT.resize(ug.nFace);
	iinv.PvfdT.resize(ug.nFace);
	iinv.PwfdT.resize(ug.nFace);
	iinv.FtuT.resize(ug.nFace);
	iinv.FtvT.resize(ug.nFace);
	iinv.FtwT.resize(ug.nFace);
	iinv.Pud.resize(ug.nFace);
	iinv.Pvd.resize(ug.nFace);
	iinv.Pwd.resize(ug.nFace);
	iinv.dsrl.resize(ug.nFace);
	iinv.elrn.resize(ug.nFace);
	iinv.Deun.resize(ug.nFace);
	iinv.Devn.resize(ug.nFace);
	iinv.Dewn.resize(ug.nFace);
	iinv.dlf.resize(ug.nFace);
	iinv.dfr.resize(ug.nFace);
	iinv.Bpe.resize(ug.nFace);
	iinv.Vau.resize(ug.nFace);
	iinv.Vav.resize(ug.nFace);
	iinv.Vaw.resize(ug.nFace);
	iinv.value.resize(ug.nFace);
	iinv.visu.resize(ug.nFace);
	iinv.visv.resize(ug.nFace);
	iinv.visw.resize(ug.nFace);
	iinv.Fu1.resize(ug.nFace);
	iinv.Fv1.resize(ug.nFace);
	iinv.Fw1.resize(ug.nFace);
	iinv.bppu.resize(ug.nTCell);
	iinv.bppv.resize(ug.nTCell);
	iinv.bppw.resize(ug.nTCell);
	iinv.l2rdx.resize(ug.nFace);
	iinv.l2rdy.resize(ug.nFace);
	iinv.l2rdz.resize(ug.nFace);
	iinv.fqr.resize(ug.nFace);
	iinv.bi1.resize(ug.nTCell);
	iinv.bi2.resize(ug.nTCell);


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
	//iinv.pp0 = 0;

	iinv.muc = 0;
	iinv.mvc = 0;
	iinv.mwc = 0;

	iinv.uc = 0;
	iinv.vc = 0;
	iinv.wc = 0;

	iinv.uu = 0;
	iinv.vv = 0;
	iinv.ww = 0;
}

void UINsInvterm::CmpInvMassFlux()
{

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//this->PrepareFaceValue();

		this->CmpINsinvTerm();
	}
}

void UINsInvterm::PrepareFaceValue()
{
	gcom.xfn = (*ug.xfn)[ug.fId];
	gcom.yfn = (*ug.yfn)[ug.fId];
	gcom.zfn = (*ug.zfn)[ug.fId];
	gcom.vfn = (*ug.vfn)[ug.fId];
	gcom.farea = (*ug.farea)[ug.fId];

	inscom.gama1 = (*uinsf.gama)[0][ug.lc];
	inscom.gama2 = (*uinsf.gama)[0][ug.rc];

	iinv.gama1 = inscom.gama1;
	iinv.gama2 = inscom.gama2;

	for (int iEqu = 0; iEqu < limf->nEqu; ++iEqu)
	{
		iinv.prim1[iEqu] = (*limf->q)[iEqu][ug.lc];
		iinv.prim2[iEqu] = (*limf->q)[iEqu][ug.rc];
	}
}

void UINsInvterm::PrepareProFaceValue()
{
	gcom.xfn = (*ug.xfn)[ug.fId];
	gcom.yfn = (*ug.yfn)[ug.fId];
	gcom.zfn = (*ug.zfn)[ug.fId];
	gcom.vfn = (*ug.vfn)[ug.fId];
	gcom.farea = (*ug.farea)[ug.fId];


	//for (int iEqu = 0; iEqu < limf->nEqu; ++iEqu)
	//{
	iinv.prim1[IIDX::IIR] = (*uinsf.q)[IIDX::IIR][ug.lc];
	iinv.prim1[IIDX::IIU] = iinv.uc[ug.lc];
	iinv.prim1[IIDX::IIV] = iinv.vc[ug.lc];
	iinv.prim1[IIDX::IIW] = iinv.wc[ug.lc];
	iinv.prim1[IIDX::IIP] = (*uinsf.q)[IIDX::IIP][ug.lc];

	iinv.prim2[IIDX::IIR] = (*uinsf.q)[IIDX::IIR][ug.rc];
	iinv.prim2[IIDX::IIU] = iinv.uc[ug.rc];
	iinv.prim2[IIDX::IIV] = iinv.vc[ug.rc];
	iinv.prim2[IIDX::IIW] = iinv.vc[ug.rc];
	iinv.prim2[IIDX::IIP] = (*uinsf.q)[IIDX::IIP][ug.rc];
	//}
}

UINsInvterm NonZero;
void UINsInvterm::Init()
{
	int Number = 0;
}

void UINsInvterm::MomPre()
{
	this->CmpINsMomRes();

	iinv.muc = 0;
	iinv.mvc = 0;
	iinv.mwc = 0;

	/* for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.ump[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];
		iinv.vmp[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
		iinv.wmp[ug.cId] = (*uinsf.q)[IIDX::IIW][ug.cId];
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

			if (ug.cId == ug.lc)
			{
				iinv.muc[ug.cId] += iinv.sj[ug.cId][iFace] * iinv.ump[ug.rc];   //使用高斯赛德尔迭代时，相邻单元对其的影响，矩阵法不需要
				iinv.mvc[ug.cId] += iinv.sj[ug.cId][iFace] * iinv.vmp[ug.rc];
				iinv.mwc[ug.cId] += iinv.sj[ug.cId][iFace] * iinv.wmp[ug.rc];

			}
			else if (ug.cId == ug.rc)
			{

				iinv.muc[ug.cId] += iinv.sj[ug.cId][iFace] * iinv.ump[ug.lc]; //使用高斯赛德尔迭代时，相邻单元对其的影响，矩阵法不需要
				iinv.mvc[ug.cId] += iinv.sj[ug.cId][iFace] * iinv.vmp[ug.lc];
				iinv.mwc[ug.cId] += iinv.sj[ug.cId][iFace] * iinv.wmp[ug.lc];

			}
		}


		//	iinv.uc[ug.cId] = (iinv.muc[ug.cId] + iinv.buc[ug.cId]) / (iinv.spc[ug.cId]);  //下一时刻速度的预测值


		//	iinv.vc[ug.cId] = (iinv.mvc[ug.cId] + iinv.bvc[ug.cId]) / (iinv.spc[ug.cId]);

		//	iinv.wc[ug.cId] = (iinv.mwc[ug.cId] + iinv.bwc[ug.cId]) / (iinv.spc[ug.cId]);

	}

	//for (int cId = ug.nCell; cId < ug.nTCell; ++cId)
	//{
	//	ug.cId = cId;

	//	iinv.uc[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];
	//	iinv.vc[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
	//	iinv.wc[ug.cId] = (*uinsf.q)[IIDX::IIW][ug.cId];
	//} */

	//BGMRES求解
	NonZero.Number = 0;
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{                                          
		int fn = (*ug.c2f)[cId].size();                             //相邻单元的个数                                    
		NonZero.Number += fn;                                          //非对角线上非零元的个数
	}
	NonZero.Number = NonZero.Number + ug.nTCell;                     //非零元的总个数         
	Rank.RANKNUMBER = ug.nTCell;                                     // 矩阵的行列大小
	Rank.NUMBER = NonZero.Number;                                    // 矩阵非零元素个数传到计算程序中
	Rank.COLNUMBER = 1;                                              //右端项个数
	Rank.Init();                                                     //传入GMRES计算程序的中间变量
	double residual_u, residual_v, residual_w;
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		Rank.TempIA[0] = 0;
		int n = Rank.TempIA[cId];
		int fn = (*ug.c2f)[cId].size();
		Rank.TempIA[cId + 1] = Rank.TempIA[cId] + fn + 1;                                                  // 前n+1行非零元素的个数
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];                                                            // 相邻面的编号
			ug.lc = (*ug.lcf)[fId];                                                                     // 面左侧单元
			ug.rc = (*ug.rcf)[fId];                                                                     // 面右侧单元
			if (cId == ug.lc)
			{
				Rank.TempA[n + iFace] = -iinv.ai[0][fId];
				Rank.TempJA[n + iFace] = ug.rc;
			}
			else if (cId == ug.rc)
			{
				Rank.TempA[n + iFace] = -iinv.ai[1][fId];
				Rank.TempJA[n + iFace] = ug.lc;
			}
		}
		Rank.TempA[n + fn] = iinv.spc[cId];                          //主对角线元素值
		Rank.TempJA[n + fn] = cId;                                      //主对角线纵坐标

	}
	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		Rank.TempB[cId][0] = iinv.buc[cId];
	}
	bgx.BGMRES();
	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		iinv.uc[cId] = Rank.TempX[cId][0];                       // 解的输出
	}
	residual_u = Rank.residual;
	iinv.res_u = residual_u;
	cout << "residual_u:" << residual_u << endl;

	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		Rank.TempB[cId][0] = iinv.bvc[cId];
	}	
	bgx.BGMRES();
	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		iinv.vc[cId] = Rank.TempX[cId][0];
	}
	residual_v = Rank.residual;
	iinv.res_v = residual_v;
	cout << "residual_v:" << residual_v << endl;

	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		Rank.TempB[cId][0] = iinv.bwc[cId];
	}
	bgx.BGMRES();
	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		iinv.wc[cId] = Rank.TempX[cId][0];
	}
	residual_w = Rank.residual;
	iinv.res_w = residual_w;
	cout << "residual_w:" << residual_w << endl;

	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	inscom.prim[IIDX::IIU] = iinv.uc[cId];                      // 解的输出
	//	inscom.prim[IIDX::IIV] = iinv.vc[cId];                      // iinv.vc 为线性方程组求得的解
	//	inscom.prim[IIDX::IIW] = iinv.wc[cId];                      // 速度压力前面乘的系数，需修改
	//}
}

void UINsInvterm::CmpFaceflux()
{

	iinv.Init();
	ug.Init();
	uinsf.Init();
	Alloc();
	//this->CmpInvFace();  //边界处理
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
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

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareProFaceValue();

		this->CmpINsBcFaceflux();
	}

}

void UINsInvterm::CmpINsMomRes()
{
	iinv.res_u = 0;
	iinv.res_v = 0;
	iinv.res_w = 0;

	//判别迭代收敛的条件
	//double phiscale, temp;
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	phiscale = iinv.uc[0];
	//	if (phiscale < iinv.uc[cId])
	//	{
	//		phiscale = iinv.uc[cId];
	//	}
	//}
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	if (iinv.spc[cId] * phiscale - 0.0 > 1e-6)
	//	{
	//		temp = iinv.buc[cId]/(iinv.spc[cId]*phiscale);
	//		iinv.res_u += temp * temp;
	//	}

	//}
	//iinv.res_u = sqrt(iinv.res_u);

	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	phiscale = iinv.vc[0];
	//	if (phiscale < iinv.vc[cId])
	//	{
	//		phiscale = iinv.vc[cId];
	//	}
	//}
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	if (iinv.spc[cId] * phiscale - 0.0 > 1e-6)
	//	{
	//		temp = iinv.bvc[cId] / (iinv.spc[cId] * phiscale);
	//		iinv.res_v += temp * temp;
	//	}

	//}
	//iinv.res_v = sqrt(iinv.res_v);

	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	phiscale = iinv.wc[0];
	//	if (phiscale < iinv.wc[cId])
	//	{
	//		phiscale = iinv.wc[cId];
	//	}
	//}
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	if (iinv.spc[cId] * phiscale - 0.0 > 1e-6)
	//	{
	//		temp = iinv.bwc[cId] / (iinv.spc[cId] * phiscale);
	//		iinv.res_w += temp * temp;
	//	}

	//}
	//iinv.res_w = sqrt(iinv.res_w);


	/*for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.res_u += (iinv.buc[ug.cId]+iinv.muc[ug.cId] - iinv.ump[ug.cId]* (iinv.spu[ug.cId]))*(iinv.buc[ug.cId]+iinv.muc[ug.cId]  - iinv.ump[ug.cId] * (iinv.spu[ug.cId]));
		iinv.res_v += (iinv.bvc[ug.cId]+iinv.mvc[ug.cId] - iinv.vmp[ug.cId] * (iinv.spv[ug.cId]))*(iinv.bvc[ug.cId]+iinv.mvc[ug.cId] - iinv.vmp[ug.cId] * (iinv.spv[ug.cId]));
		iinv.res_w += (iinv.bwc[ug.cId]+iinv.mwc[ug.cId] - iinv.wmp[ug.cId] * (iinv.spw[ug.cId]))*(iinv.bwc[ug.cId]+iinv.mwc[ug.cId] - iinv.wmp[ug.cId] * (iinv.spw[ug.cId]));

	}

	iinv.res_u = sqrt(iinv.res_u);
	iinv.res_v = sqrt(iinv.res_v);
	iinv.res_w = sqrt(iinv.res_w);*/

}

void UINsInvterm::AddFlux()
{
	UnsGrid* grid = Zone::GetUnsGrid();
	MRField* res = GetFieldPointer< MRField >(grid, "res");
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
	this->CmpNewMomCoe();
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
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

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->CmpINsBcFaceCorrectPresscoef();
	}

	iinv.spp = 0;
	iinv.bp = 0;

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.spp[ug.lc] += iinv.ajp[ug.fId];
		iinv.spp[ug.rc] += iinv.ajp[ug.fId];

		iinv.bp[ug.lc] += -iinv.fq[ug.fId];
		iinv.bp[ug.rc] += iinv.fq[ug.fId];

		if (ug.fId < ug.nBFace)
		{
			//iinv.spp[ug.rc] = 0.001;
			iinv.spp[ug.rc] = 1;
		}
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		//iinv.VdU[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spu[ug.cId] - iinv.sju[ug.cId]); //用于求单元修正速度量;
		//iinv.VdV[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spv[ug.cId] - iinv.sjv[ug.cId]);
		//iinv.VdW[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spw[ug.cId] - iinv.sjw[ug.cId]);

		iinv.VdU[ug.cId] = -(*ug.cvol)[ug.cId] / (iinv.spc[ug.cId]); //用于求单元修正速度量;
		iinv.VdV[ug.cId] = -(*ug.cvol)[ug.cId] / (iinv.spc[ug.cId]);
		iinv.VdW[ug.cId] = -(*ug.cvol)[ug.cId] / (iinv.spc[ug.cId]);

		//iinv.spp[ug.cId] = (*ug.cvol)[ug.cId] / iinv.timestep;

		//iinv.bp[ug.cId] = iinv.bi1[ug.cId]+ iinv.bi2[ug.cId];

		int fn = (*ug.c2f)[ug.cId].size();
		iinv.sjp.resize(ug.nTCell, fn);
		iinv.sjd.resize(ug.nTCell, fn);
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			if (ug.cId == ug.lc)
			{
				iinv.sjp[ug.cId][iFace] = -iinv.ajp[ug.fId]; //求解压力修正方程的非零系数
				iinv.sjd[ug.cId][iFace] = ug.rc;

				//cout << "iinv.sjp=" << iinv.sjp[ug.cId][iFace] << "iinv.sjd=" << ug.rc << "\n";
			}
			else if (ug.cId == ug.rc)
			{
				iinv.sjp[ug.cId][iFace] = -iinv.ajp[ug.fId];
				iinv.sjd[ug.cId][iFace] = ug.lc;

				//cout << "iinv.sjp=" << iinv.sjp[ug.cId][iFace] << "iinv.sjd=" << ug.lc << "\n";
			}
		}
	}

	/*for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		cout << "iinv.bp=" << iinv.buc[ug.cId] << "\n";
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		cout << "iinv.spp=" << iinv.spp[ug.cId] << "\n";
	}*/

	//iinv.spp[0] = 3.996004185733362E-003;
	//iinv.spp[1] = 3.996004185733362E-003;
	//iinv.spp[2] = 3.996004185733362E-003;
	//iinv.spp[3] = 3.996004185733362E-003;
	//iinv.spp[4] = 3.996004185733362E-003;
	//iinv.spp[5] = 3.996004185733362E-003; 
	//iinv.spp[6] = 3.996004185733362E-003;
	//iinv.spp[7] = 3.996004185733362E-003;
	//iinv.spp[8] = 3.996004185733362E-003; 
	//iinv.spp[9] = 3.996004185733362E-003; 
	//iinv.spp[10] = 3.996004185733362E-003;
	//iinv.spp[11] = 3.996004185733362E-003;
	//iinv.spp[12] = 3.996004185733362E-003;
	//iinv.spp[13] = 3.996004185733362E-003; 
	//iinv.spp[14] = 3.996004185733362E-003;
	//iinv.spp[15] = 3.996003562604947E-003; 
	//	iinv.spp[16] = 3.996004185576957E-003; 
	//	iinv.spp[17] = 3.996004185733322E-003; 
	//	iinv.spp[18] = 3.996004185733362E-003; 
	//	iinv.spp[19] = 3.996004809018999E-003; 
	//	iinv.spp[20] = 3.993512901138565E-003; 
	//	iinv.spp[21] = 3.996003562135695E-003; 
	//	iinv.spp[22] = 3.996004185576840E-003;
	//	iinv.spp[23] = 3.996004809018960E-003;
	//	iinv.spp[24] = 3.998507933456209E-003;

		//for (int cId = 0; cId < 20; ++cId)
		//{
		//	ug.cId = cId;

		//	iinv.bp[ug.cId] = 0;
		//}

		//iinv.bp[20] = 9.990010315470651E-005;
		//iinv.bp[21] = 2.507471755350644E-008;
		//iinv.bp[22] = 6.309350054906251E-012;
		//iinv.bp[23] = 1.587575580783794E-015;
		//iinv.bp[24] = -9.992518418319765E-005;


}

void UINsInvterm::CmpNewMomCoe()
{
	iinv.spc = 0;

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.spc[ug.lc] += iinv.ai[0][ug.fId] + iinv.Fn[ug.fId];
		iinv.spc[ug.rc] += iinv.ai[1][ug.fId] + iinv.Fn[ug.fId];
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;
		iinv.spc[ug.cId] += iinv.spt[ug.cId];
	}

	//for (int cId = 0; cId < ug.nTCell; ++cId)
	//{
	//	ug.cId = cId;

	//	iinv.spu[ug.cId] = iinv.bi1[ug.cId] + iinv.bi2[ug.cId] + iinv.aku1[ug.cId] + iinv.aku2[ug.cId] + iinv.spt[ug.cId]; //矩阵主对角线系数，动量方程单元主系数
	//	iinv.spv[ug.cId] = iinv.bi1[ug.cId] + iinv.bi2[ug.cId] + iinv.akv1[ug.cId] + iinv.akv2[ug.cId] + iinv.spt[ug.cId];
	//	iinv.spw[ug.cId] = iinv.bi1[ug.cId] + iinv.bi2[ug.cId] + iinv.akw1[ug.cId] + iinv.akw2[ug.cId] + iinv.spt[ug.cId];
	//}

}

void UINsInvterm::CmpPressCorrectEqu()
{
	//double rhs_p = 1e-8;
	//iinv.res_p = 1;
	//iinv.mp = 0;
    //iinv.pp = 0.5;
	//while (iinv.res_p >= rhs_p)
	//{
	//	iinv.res_p = 0.0;

	//	for (int cId = 0; cId < ug.nCell; ++cId)
	//	{
	//		ug.cId = cId;

	//		iinv.ppd = iinv.pp[ug.cId];
	//		int fn = (*ug.c2f)[ug.cId].size();
	//		for (int iFace = 0; iFace < fn; ++iFace)
	//		{
	//			int fId = (*ug.c2f)[ug.cId][iFace];
	//			ug.fId = fId;
	//			if (ug.fId < ug.nBFace) continue;

	//			ug.lc = (*ug.lcf)[ug.fId];
	//			ug.rc = (*ug.rcf)[ug.fId];
	//			if (ug.cId == ug.lc)
	//			{
	//				iinv.mp[ug.cId] += iinv.sjp[ug.cId][iFace] * iinv.pp[ug.rc]; //高斯赛戴尔迭代求解时的相邻单元的值，矩阵法不需要
	//			}
	//			else if (ug.cId == ug.rc)
	//			{
	//				iinv.mp[ug.cId] += iinv.sjp[ug.cId][iFace] * iinv.pp[ug.lc];
	//			}
	//		}
	//		iinv.pp[ug.cId] = (iinv.bp[ug.cId] + iinv.mp[ug.cId]) / (iinv.spp[ug.cId]); //压力修正值

	//		iinv.res_p = MAX(iinv.res_p, abs(iinv.ppd - iinv.pp[ug.cId]));

	//	}

	//}

		//BGMRES求解
	NonZero.Number = 0;
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{   
		//ug.cId = cId;                                                                  // 主单元编号
		int fn = (*ug.c2f)[cId].size();                                                                 // 单元相邻面的个数
		NonZero.Number += fn;
	}
	NonZero.Number = NonZero.Number + ug.nTCell;                                                        // 非零元素的计数
	Rank.RANKNUMBER = ug.nTCell;                                                                        // 矩阵的行列
	Rank.COLNUMBER = 1;
	Rank.NUMBER = NonZero.Number;                                                                      // 矩阵非零元素个数
	Rank.Init();
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		iinv.ppd = iinv.pp[cId];
		Rank.TempIA[0] = 0;
		int n = Rank.TempIA[cId];
		int fn = (*ug.c2f)[cId].size();
		Rank.TempIA[cId + 1] = Rank.TempIA[cId] + fn + 1;                  // 前n+1行非零元素的个数
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];                           // 相邻面的编号
			ug.fId = fId;
			ug.lc = (*ug.lcf)[fId];                                    // 面左侧单元
			ug.rc = (*ug.rcf)[fId];                                    // 面右侧单元
			if (cId == ug.lc)
			{
				Rank.TempA[n + iFace] = iinv.sjp[cId][iFace];          //非对角线元素值
				Rank.TempJA[n + iFace] = ug.rc;                           //非对角线元素纵坐标
			}
			else if (cId == ug.rc)
			{
				Rank.TempA[n + iFace] = iinv.sjp[cId][iFace];          //非对角线元素值
				Rank.TempJA[n + iFace] = ug.lc;                           //非对角线元素纵坐标
			}
		}
		Rank.TempA[n + fn] = iinv.spp[cId];                            //主对角线元素
		Rank.TempJA[n + fn] = cId;                                        //主对角线纵坐标

		Rank.TempB[cId][0] = iinv.bp[cId];                             //右端项
	}
	bgx.BGMRES();
	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		//ug.cId = cId;
		iinv.pp[cId] = Rank.TempX[cId][0]; //当前时刻的压力修正值
	}

	//iinv.res_p = 0;
	//iinv.res_p = MAX(iinv.res_p, abs(iinv.ppd - iinv.pp[ug.cId]));

	//边界单元
	/*for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.pp[ug.rc] = iinv.pp[ug.lc];
	}*/

	/*if (Iteration::outerSteps != 1)
	{
		iinv.pp = 0;
	}*/

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;
		(*uinsf.q)[IIDX::IIP][ug.cId] = (*uinsf.q)[IIDX::IIP][ug.cId] +iinv.pp[ug.cId];
	}



	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	iinv.pc[ug.cId] = inscom.prim[IIDX::IIP] + iinv.pp[ug.cId]; //下一时刻的压力值
	//}
}


void UINsInvterm::CmpINsPreRes()
{
	//iinv.res_p = 0;


	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		if (ug.cId == 0)
		{
			iinv.res_p = abs(iinv.bp[ug.cId]);
		}
		else
		{
			iinv.res_p = MAX(abs(iinv.bp[ug.cId]), abs(iinv.bp[ug.cId - 1]));
		}

		//iinv.res_p += (iinv.bp[ug.cId]+iinv.mp[ug.cId] - iinv.pp1[ug.cId]* (0.01+iinv.spp[ug.cId]))*(iinv.bp[ug.cId]+iinv.mp[ug.cId] - iinv.pp1[ug.cId]* (0.01+iinv.spp[ug.cId]));
	}

	//iinv.res_p = sqrt(iinv.res_p);
	//iinv.res_p = 0;
}


void UINsInvterm::UpdateFaceflux()
{
	iinv.Init();
	ug.Init();
	uinsf.Init();
	Alloc();
	//this->CmpInvFace();  //边界处理
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
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

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//this->PrepareFaceValue();

		this->CmpUpdateINsBcFaceflux();
	}

}

void UINsInvterm::CmpUpdateINsBcFaceflux()
{
	iinv.uuj[ug.fId] = 0; //面速度修正量
	iinv.vvj[ug.fId] = 0;
	iinv.wwj[ug.fId] = 0;

	iinv.uf[ug.fId] = iinv.uf[ug.fId] + iinv.uuj[ug.fId]; //下一时刻面速度
	iinv.vf[ug.fId] = iinv.vf[ug.fId] + iinv.vvj[ug.fId];
	iinv.wf[ug.fId] = iinv.wf[ug.fId] + iinv.wwj[ug.fId];

	iinv.fux[ug.fId] = iinv.rf[ug.fId] * (gcom.xfn * iinv.uuj[ug.fId] + gcom.yfn * iinv.vvj[ug.fId] + gcom.zfn * iinv.wwj[ug.fId]) * gcom.farea;
	iinv.fq[ug.fId] = iinv.fq[ug.fId] + iinv.fux[ug.fId];
}


void UINsInvterm::CmpUpdateINsFaceflux()
{
	iinv.uuj[ug.fId] = iinv.Vdvu[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) * (*ug.xfn)[ug.fId] / iinv.dist[ug.fId]; //面速度修正量
	iinv.vvj[ug.fId] = iinv.Vdvv[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) * (*ug.yfn)[ug.fId] / iinv.dist[ug.fId];
	iinv.wwj[ug.fId] = iinv.Vdvw[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) * (*ug.zfn)[ug.fId] / iinv.dist[ug.fId];

	iinv.uf[ug.fId] = iinv.uf[ug.fId] + iinv.uuj[ug.fId]; //下一时刻面速度
	iinv.vf[ug.fId] = iinv.vf[ug.fId] + iinv.vvj[ug.fId];
	iinv.wf[ug.fId] = iinv.wf[ug.fId] + iinv.wwj[ug.fId];

	iinv.fux[ug.fId] = iinv.rf[ug.fId] * ((*ug.xfn)[ug.fId] * iinv.uuj[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vvj[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wwj[ug.fId]) * gcom.farea;
	iinv.fq[ug.fId] = iinv.fq[ug.fId] + iinv.fux[ug.fId];

}

void UINsInvterm::UpdateSpeed()
{
	this->CmpPreGrad();

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.uu[ug.cId] = iinv.VdU[ug.cId] * iinv.dqqdx[ug.cId]; //速度修正量
		iinv.vv[ug.cId] = iinv.VdV[ug.cId] * iinv.dqqdy[ug.cId];
		iinv.ww[ug.cId] = iinv.VdW[ug.cId] * iinv.dqqdz[ug.cId];

		iinv.up[ug.cId] = iinv.uc[cId] + iinv.uu[ug.cId];  //下一时刻的速度值
		iinv.vp[ug.cId] = iinv.vc[cId] + iinv.vv[ug.cId];
		iinv.wp[ug.cId] = iinv.wc[cId] + iinv.ww[ug.cId];

		(*uinsf.q)[IIDX::IIU][ug.cId] = iinv.up[ug.cId];
		(*uinsf.q)[IIDX::IIV][ug.cId] = iinv.vp[ug.cId];
		(*uinsf.q)[IIDX::IIW][ug.cId] = iinv.wp[ug.cId];

	}

	for (int cId = ug.nCell; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.uu[ug.cId] = 0; //速度修正量
		iinv.vv[ug.cId] = 0;
		iinv.ww[ug.cId] = 0;

		iinv.up[ug.cId] = iinv.uc[cId] + iinv.uu[ug.cId];  //下一时刻的速度值
		iinv.vp[ug.cId] = iinv.vc[cId] + iinv.vv[ug.cId];
		iinv.wp[ug.cId] = iinv.wc[cId] + iinv.ww[ug.cId];

		(*uinsf.q)[IIDX::IIU][ug.cId] = iinv.up[ug.cId];
		(*uinsf.q)[IIDX::IIV][ug.cId] = iinv.vp[ug.cId];
		(*uinsf.q)[IIDX::IIW][ug.cId] = iinv.wp[ug.cId];
	}
}

void UINsInvterm::CmpPreGrad()
{
	iinv.dqqdx = 0;
	iinv.dqqdy = 0;
	iinv.dqqdz = 0;

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
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
		//if (ug.fId < ug.nBFace)
		//{
		//	iinv.value[ug.fId] = iinv.pp[ug.lc] + iinv.pp[ug.rc];
		//}
		//else
		//{
		iinv.value[ug.fId] = cl * iinv.pp[ug.lc] + cr * iinv.pp[ug.rc];
		//}

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		iinv.dqqdx[ug.lc] += fnxa * iinv.value[ug.fId];
		iinv.dqqdy[ug.lc] += fnya * iinv.value[ug.fId];
		iinv.dqqdz[ug.lc] += fnza * iinv.value[ug.fId];

		if (ug.fId < ug.nBFace) continue;

		iinv.dqqdx[ug.rc] += -fnxa * iinv.value[ug.fId];
		iinv.dqqdy[ug.rc] += -fnya * iinv.value[ug.fId];
		iinv.dqqdz[ug.rc] += -fnza * iinv.value[ug.fId];
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		Real ovol = one / (*ug.cvol)[ug.cId];
		iinv.dqqdx[ug.cId] *= ovol;
		iinv.dqqdy[ug.cId] *= ovol;
		iinv.dqqdz[ug.cId] *= ovol;
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//if (ug.rc > ug.nCell)
		//{
		iinv.dqqdx[ug.rc] = iinv.dqqdx[ug.lc];
		iinv.dqqdy[ug.rc] = iinv.dqqdy[ug.lc];
		iinv.dqqdz[ug.rc] = iinv.dqqdz[ug.lc];
		//}

	}

}


void UINsInvterm::Alloc()
{
	iinvflux = new MRField(inscom.nEqu, ug.nFace);
}

void UINsInvterm::DeAlloc()
{
	delete iinvflux;
}


void UINsInvterm::ReadTmp()
{
	static int iii = 0;
	if (iii) return;
	iii = 1;
	fstream file;
	file.open("nsflow.dat", ios_base::in | ios_base::binary);
	if (!file) exit(0);

	uinsf.Init();

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		for (int iEqu = 0; iEqu < 5; ++iEqu)
		{
			file.read(reinterpret_cast<char*>(&(*uinsf.q)[iEqu][cId]), sizeof(double));
		}
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		file.read(reinterpret_cast<char*>(&(*uinsf.visl)[0][cId]), sizeof(double));
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		file.read(reinterpret_cast<char*>(&(*uinsf.vist)[0][cId]), sizeof(double));
	}

	vector< Real > tmp1(ug.nTCell), tmp2(ug.nTCell);

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		tmp1[cId] = (*uinsf.timestep)[0][cId];
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		file.read(reinterpret_cast<char*>(&(*uinsf.timestep)[0][cId]), sizeof(double));
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		tmp2[cId] = (*uinsf.timestep)[0][cId];
	}

	turbcom.Init();
	uturbf.Init();
	for (int iCell = 0; iCell < ug.nTCell; ++iCell)
	{
		for (int iEqu = 0; iEqu < turbcom.nEqu; ++iEqu)
		{
			file.read(reinterpret_cast<char*>(&(*uturbf.q)[iEqu][iCell]), sizeof(double));
		}
	}
	file.close();
	file.clear();
}



EndNameSpace