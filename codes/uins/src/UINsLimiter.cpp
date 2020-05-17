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

#include "UINsLimiter.h"
#include "UnsGrid.h"
#include "Ctrl.h"
#include "Zone.h"
#include "DataBase.h"
#include "UCom.h"
#include "INsCom.h"
#include "UINsCom.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "INsIdx.h"
#include "Iteration.h"
#include "INsInvterm.h"

BeginNameSpace( ONEFLOW )

INsLimField::INsLimField()
{
    qf1 = 0;
    qf2 = 0;
    this->nEqu = inscom.nEqu;
}

INsLimField::~INsLimField()
{
    delete qf1;
    delete qf2;
}

void INsLimField::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    q       = GetFieldPointer< MRField > ( grid, "q" );
    dqdx    = GetFieldPointer< MRField > ( grid, "dqdx" );
    dqdy    = GetFieldPointer< MRField > ( grid, "dqdy" );
    dqdz    = GetFieldPointer< MRField > ( grid, "dqdz" );
    limiter = GetFieldPointer< MRField > ( grid, "limiter" );

    this->nEqu = q->GetNEqu();

    qf1 = new MRField( this->nEqu, grid->nFace );
    qf2 = new MRField( this->nEqu, grid->nFace );

    this->ckfun = & INsCheckFunction;
}

void INsLimField::BcQlQrFix()
{
	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		int bcType = ug.bcRecord->bcType[fId];
		if (bcType == BC::INTERFACE) continue;
		if (bcType == BC::PERIODIC) continue;

		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real tmp1 = half * ((*this->q)[IIDX::IIR][ug.lc] + (*this->q)[IIDX::IIR][ug.rc]);
		Real tmp2 = half * ((*this->q)[IIDX::IIU][ug.lc] + (*this->q)[IIDX::IIU][ug.rc]);
		Real tmp3 = half * ((*this->q)[IIDX::IIV][ug.lc] + (*this->q)[IIDX::IIV][ug.rc]);
		Real tmp4 = half * ((*this->q)[IIDX::IIW][ug.lc] + (*this->q)[IIDX::IIW][ug.rc]);
		Real tmp5 = half * ((*this->q)[IIDX::IIP][ug.lc] + (*this->q)[IIDX::IIP][ug.rc]);

		(*this->qf1)[IIDX::IIR][ug.fId] = tmp1;
		(*this->qf2)[IIDX::IIR][ug.fId] = tmp1;
		(*this->qf1)[IIDX::IIP][ug.fId] = tmp5;
		(*this->qf2)[IIDX::IIP][ug.fId] = tmp5;

		if (ug.lc < ug.nCell)
		{
			(*this->q)[IIDX::IIU][ug.lc] += tmp2;
			(*this->q)[IIDX::IIV][ug.lc] += tmp3;
			(*this->q)[IIDX::IIW][ug.lc] += tmp4;

			(*this->q)[IIDX::IIU][ug.rc] = tmp2;
			(*this->q)[IIDX::IIV][ug.rc] = tmp3;
			(*this->q)[IIDX::IIW][ug.rc] = tmp4;
		}
		else if (ug.rc < ug.nCell)
		{
			(*this->q)[IIDX::IIU][ug.rc] += tmp2;
			(*this->q)[IIDX::IIV][ug.rc] += tmp3;
			(*this->q)[IIDX::IIW][ug.rc] += tmp4;

			(*this->q)[IIDX::IIU][ug.lc] = tmp2;
			(*this->q)[IIDX::IIV][ug.lc] = tmp3;
			(*this->q)[IIDX::IIW][ug.lc] = tmp4;
		}

	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{

		int bcType = ug.bcRecord->bcType[fId];
		if (bcType == BC::INTERFACE) continue;
		if (bcType == BC::PERIODIC) continue;

		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		if (bcType == BC::SOLID_SURFACE)
		{
			if (ug.lc < ug.nCell)
			{
				for (int iEqu = 0; iEqu < this->nEqu; ++iEqu)
				{
					(*this->qf2)[iEqu][ug.fId] = (*uinsf.bc_q)[iEqu][ug.fId];
					(*this->q)[iEqu][ug.rc] = (*this->qf2)[iEqu][ug.fId];
				}
			}
			else
			{
				for (int iEqu = 0; iEqu < this->nEqu; ++iEqu)
				{
					(*this->qf1)[iEqu][ug.fId] = (*uinsf.bc_q)[iEqu][ug.fId];
					(*this->q)[iEqu][ug.lc] = (*this->qf1)[iEqu][ug.fId];
				}
			}
		}

	}

}


INsLimiter::INsLimiter()
{
    limf = new INsLimField();
    limflag = ctrl.ilim;
}

INsLimiter::~INsLimiter()
{
    delete limf;
}

EndNameSpace