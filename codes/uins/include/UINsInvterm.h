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


#pragma once
#include "INsInvterm.h"

BeginNameSpace( ONEFLOW )

class UINsFField;
class Limiter;
class LimField;

class UINsInvterm : public INsInvterm
{
public:
    UINsInvterm();
    ~UINsInvterm();
public:
    void Alloc();
    void DeAlloc();
	void CmpINsTimestep();
    void CmpInvcoff();
    void CmpInvMassFlux();
    void CmpInvFace();
    void CmpLimiter();
	void CmpFaceflux();
	void CmpINsMomRes();
	void CmpINsPreRes();
	void CmpCorrectPresscoef();
	void CmpPressCorrectEqu();
	void UpdateFaceflux();
	void CmpUpdateINsFaceflux();
	void UpdateSpeed();
    void AddFlux();
    void PrepareFaceValue();
	void CmpPreGrad();
	//void CmpINsinvTerm();
    //void UpdateFaceInvFlux();
    void ReadTmp();
public:
    void GetQlQrField();
    void ReconstructFaceValueField();
    void BoundaryQlQrFixField();
	void MomPre();
public:
    Limiter * limiter;
    LimField * limf;
    MRField * iinvflux;
};
//void PrimToQ(RealField & prim, Real gama, RealField & q);

EndNameSpace