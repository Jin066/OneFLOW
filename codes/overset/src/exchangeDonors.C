//
// This file is part of the Tioga software library
//
// Tioga  is a tool for overset grid assembly on parallel distributed systems
// Copyright (C) 2015 Jay Sitaraman
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#include <algorithm>
#include "codetypes.h"
#include "tioga.h"
using namespace TIOGA;
void tioga::exchangeDonors(void)
{
  int nsend,nrecv;
  int *sndMap;
  int *rcvMap;
  PACKET *sndPack,*rcvPack;
  //
  // get the processor map for sending 获取用于发送和接收的处理器映射
  // and receiving
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend == 0) return;  
  //
  // create packets to send and receive创建数据包以发送和接收并将其初始化为零
  // and initialize them to zero
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  //
  pc->initPackets(sndPack,rcvPack);

  std::vector<int> nintsSend(nsend,0), nrealsSend(nsend,0);

  // First pass through all mesh blocks and determine data sizes for each
  // intersection pair 首先通过所有网格块并确定每个交叉口对的数据大小
  for (int ib=0;ib<nblocks;ib++) {
    auto &mb =mblocks[ib];
    mb->getMBDonorPktSizes(nintsSend, nrealsSend);
  }
  //
  // Allocate sndPack 分配发送包
  //
  for(int k=0;k<nsend;k++)
    {
      sndPack[k].nints=nintsSend[k];
      sndPack[k].nreals=nrealsSend[k];
      sndPack[k].intData = (int*)malloc(sizeof(int) * sndPack[k].nints);
      sndPack[k].realData = (double*)malloc(sizeof(double) * sndPack[k].nreals);    
    }
  //
  // Populate send packets with data from each mesh block in this partition
  //用此分区中每个网格块的数据填充发送数据包
  std::vector<int> ixOffset(nsend,0), rxOffset(nsend,0);
  for (int ib=0;ib<nblocks;ib++) {
    auto &mb = mblocks[ib];
    mb->getMBDonorPackets(ixOffset, rxOffset, sndPack); //得到网格块贡献单元的发送包
  }
  //
  // communicate donors (comm1)  联络贡献单元（通讯1）
  //
  pc->sendRecvPackets(sndPack,rcvPack); //接收和发送包
  // Initialize linked lists and populate donor data from rcvPack 从接收包中初始化链接列表并填充贡献者数据
  for (int ib=0;ib<nblocks;ib++) {
    auto& mb = mblocks[ib];
    mb->initializeDonorList(); //初始化贡献者列表
  }
  //
  for (int k=0; k<nrecv; k++) {
    int m=0;
    int l=0;
    for(int i=0;i< rcvPack[k].nints/4;i++)
      {
        int meshtag = rcvPack[k].intData[m++];
        int pointid = rcvPack[k].intData[m++];
        int remoteid = rcvPack[k].intData[m++];
	    int ib       = rcvPack[k].intData[m++];
        double donorRes = rcvPack[k].realData[l++];
        double receptorRes = rcvPack[k].realData[l++];
	auto& mb = mblocks[ib];
        mb->insertAndSort(pointid, k, meshtag, remoteid, donorRes,receptorRes); //插入和排序
      }
  }
  //
  // Figure out the state of each point (i.e., if it is a hole, fringe, or a
  // field point) 找出每个点的状态（即，如果它是一个洞、条纹或一个场点）
  //
  std::vector<int> nrecords(nblocks,0);
  int** donorRecords = (int**)malloc(sizeof(int*)*nblocks);
  double** receptorResolution = (double**)malloc(sizeof(double*)*nblocks);
  for (int ib=0; ib<nblocks; ib++) {
    auto& mb = mblocks[ib];
    mb->processDonors(holeMap, nmesh, &(donorRecords[ib]),
                      &(receptorResolution[ib]),&(nrecords[ib])); //进程贡献者
  }
  //
  // Reset all send/recv data structures 重置所有发送/接收数据结构
  //
  pc->clearPackets(sndPack, rcvPack);
  for (int i=0; i<nsend; i++) {
    sndPack[i].nints=0;
    sndPack[i].nreals=0;
    nintsSend[i] = 0;
    nrealsSend[i] = 0;
    ixOffset[i] = 0;
    rxOffset[i] = 0;
  }

  for (int n=0; n<nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int k = donorRecords[n][3*i];
      sndPack[k].nints+=2;
      sndPack[k].nreals++;
    }
  }

  for(int k=0;k<nsend;k++)
    {
      sndPack[k].intData = (int*)malloc(sizeof(int) * sndPack[k].nints);
      sndPack[k].realData = (double*)malloc(sizeof(double) * sndPack[k].nreals);
    }

  for (int n=0; n<nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int k = donorRecords[n][3*i];      
      sndPack[k].intData[ixOffset[k]++]=donorRecords[n][3*i+1];
      sndPack[k].intData[ixOffset[k]++]=donorRecords[n][3*i+2];
      sndPack[k].realData[rxOffset[k]++]=receptorResolution[n][i];
    }
  }
  //
  // comm 2 通讯2
  // notify donors that they are accepted (for now)通知捐助者他们已被接受（目前）
  //
  pc->sendRecvPackets(sndPack,rcvPack);

  std::vector<int> ninterp(nblocks,0);
  for(int k=0;k<nrecv;k++)
    {
      int m=0;
      for(int j=0; j < rcvPack[k].nints/2;j++)
	{
	  m++;   // skip over point id 跳过点id
	  int ib=tag_iblk_map[rcvPack[k].intData[m++]];
	  ninterp[ib]++;
	}
    }

  for (int i=0; i<nblocks; i++) {
    mblocks[i]->initializeInterpList(ninterp[i]); //初始化插值列表
  }

  std::fill(ninterp.begin(), ninterp.end(), 0);//填充

  for(int k=0; k < nrecv;k++)
    {
      int l=0;
      int m=0;
      for(int j=0;j< rcvPack[k].nints/2;j++)
	{
	  int recid=rcvPack[k].intData[m++];
	  int ib = tag_iblk_map[rcvPack[k].intData[m++]];
	  double receptorRes=rcvPack[k].realData[l++];
	  mblocks[ib]->findInterpData(&(ninterp[ib]),recid,receptorRes); //寻找插值数据
	}
    }
  
  for (int ib=0; ib<nblocks; ib++) {
    mblocks[ib]->set_ninterp(ninterp[ib]);
  }

  pc->clearPackets(sndPack, rcvPack);
  //
  // Find cancellation data (based on donor quality)查找注销数据（根据捐助者的质量）
  //
  for (int i=0; i<nblocks; i++) {
    if (donorRecords[i]) {
      TIOGA_FREE(donorRecords[i]);
      donorRecords[i] = NULL;
    }
    nrecords[i] = 0;
    mblocks[i]->getCancellationData(&(nrecords[i]),
                                    &(donorRecords[i]));
  }
  std::fill(nintsSend.begin(), nintsSend.end(), 0);
  std::fill(ixOffset.begin(), ixOffset.end(), 0);

  for (int n=0; n < nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int k = donorRecords[n][3*i];
      sndPack[k].nints+=2;
    }
  }
  for(int k=0;k<nsend;k++)
    sndPack[k].intData = (int*)malloc(sizeof(int) * sndPack[k].nints);

  for (int n=0; n < nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int k= donorRecords[n][3*i];
      sndPack[k].intData[ixOffset[k]++]=donorRecords[n][3*i+1];
      sndPack[k].intData[ixOffset[k]++]=donorRecords[n][3*i+2];
    }
  }
  //
  // communciate cancellation data comm 3 通讯 取消数据 （通讯3）
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  
  for (int k=0; k<nrecv; k++) {
    int m = 0;
    for (int j=0;j<rcvPack[k].nints/2;j++) {
      int recid=rcvPack[k].intData[m++];
      int ib = tag_iblk_map[rcvPack[k].intData[m++]];
      mblocks[ib]->cancelDonor(recid);
    }
  }

  for (int ib=0;ib<nblocks;ib++) {
    auto &mb =mblocks[ib];
    mb->resetCoincident();
  }

  //
  pc->clearPackets(sndPack, rcvPack);
  //
  // Find final interpolation data查找最终插值数据
  //
  for (int i=0; i<nblocks; i++) {
    if (donorRecords[i]) {
      TIOGA_FREE(donorRecords[i]);
      donorRecords[i] = NULL;
    }
    nrecords[i] = 0;    
    mblocks[i]->getInterpData(&(nrecords[i]),
                              &(donorRecords[i]));
  }  
  std::fill(nintsSend.begin(), nintsSend.end(), 0);
  std::fill(ixOffset.begin(), ixOffset.end(), 0);
  for (int n=0; n < nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int k = donorRecords[n][3*i];
      sndPack[k].nints+=2;
    }
  }
  for(int k=0;k<nsend;k++)
    sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
  for (int n=0; n < nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int k = donorRecords[n][3*i];
      sndPack[k].intData[ixOffset[k]++]=donorRecords[n][3*i+1];
      sndPack[k].intData[ixOffset[k]++]=donorRecords[n][3*i+2];
    }
  }
  //
  // comm 4 通讯4
  // final receptor data to set iblanks最终受体数据以设置iblank
  //     
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  for(int ib=0;ib<nblocks;ib++)
    mblocks[ib]->clearIblanks();
  
  for (int k=0; k<nrecv; k++) {
    int m = 0;
    for(int j=0;j< rcvPack[k].nints/2;j++)
      {
	int pointid=rcvPack[k].intData[m++];
	int ib=rcvPack[k].intData[m++];
	mblocks[ib]->setIblanks(pointid);
      }
  }
  pc->clearPackets(sndPack,rcvPack);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
  
  if (donorRecords) {
    for (int i=0; i<nblocks; i++) {
      if (donorRecords[i]) TIOGA_FREE(donorRecords[i]);
    }
    TIOGA_FREE(donorRecords);
  }
  if (receptorResolution) {
    for (int i=0; i<nblocks; i++) {
      if (receptorResolution[i]) TIOGA_FREE(receptorResolution[i]);
    }
    TIOGA_FREE(receptorResolution);
  }   
}
  
void tioga::outputStatistics(void)
{
#ifdef TIOGA_OUTPUT_STATS
  int mstats[2], mstats_sum[2], mstats_global[2];
  mstats_sum[0] = 0;
  mstats_sum[1] = 0;
  for (int ib=0;ib< nblocks;ib++) {
    auto& mb = mblocks[ib];
    mb->getStats(mstats);
    mstats_sum[0] += mstats[0];
    mstats_sum[1] += mstats[1];
  }
  MPI_Reduce(mstats_sum,mstats_global,2,MPI_INT,MPI_SUM,0,scomm);
  if (myid==0) {
    printf("#tioga -----------------------------------------\n");
    printf("#tioga : total receptors:\t%d\n",mstats_global[1]);
    printf("#tioga : total holes    :\t%d\n",mstats_global[0]);
    printf("#tioga -----------------------------------------\n");
  }
#endif
}
