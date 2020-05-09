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
#include <vector>
#include <cstring>
#include "codetypes.h"
#include "tioga.h"
using namespace TIOGA;
void tioga::exchangeSearchData(int at_points)
{
  int i;
  int nsend, nrecv;
  PACKET *sndPack, *rcvPack;
  int* sndMap;
  int* rcvMap;
  //
  // get the processor map for sending 获取用于发送/和接收的处理器映射
  // and receiving
  //
  pc->getMap(&nsend, &nrecv, &sndMap, &rcvMap); //得到map
  //
  // create packets to send and receive 创建发送和接收数据包并将其初始化为零
  // and initialize them to zero
  //
  sndPack = (PACKET*)malloc(sizeof(PACKET) * nsend);
  rcvPack = (PACKET*)malloc(sizeof(PACKET) * nrecv);
  //
  for (i = 0; i < nsend; i++) {
    sndPack[i].nints = sndPack[i].nreals = 0;
    sndPack[i].intData = NULL;
    sndPack[i].realData = NULL;
  }
  //
  for (i = 0; i < nrecv; i++) {
    rcvPack[i].nints = rcvPack[i].nreals = 0;
    rcvPack[i].intData = NULL;
    rcvPack[i].realData = NULL;
  }

  // Process each intersection pair and determine the total data that needs to
  // be sent 处理每个相交对（重叠区），并确定需要发送的总数据
  int nobb = obblist.size();
  std::vector<int> nintsSend(nobb);
  std::vector<int> nrealsSend(nobb);
  int** int_data = (int**)malloc(sizeof(int*) * nobb);
  double** real_data = (double**)malloc(sizeof(double*) * nobb);

  for (int ii=0; ii < nobb; ii++) {
    int ib = obblist[ii].iblk_local; //OBB框里的网格块索引（本地）
    auto& mb = mblocks[ib];
    if (at_points==0) 
    {
      mb->getQueryPoints2(
      &obblist[ii], &nintsSend[ii], &int_data[ii], &nrealsSend[ii],
      &real_data[ii]);  //得到查询点2
    } 
    else
    {
     mb->getExtraQueryPoints(&obblist[ii],
                            &nintsSend[ii],&int_data[ii],
                            &nrealsSend[ii],&real_data[ii]); //获取额外的查询点
    }
  }

  // Populate send packets and exchange data with other processors 填充发送数据包并与其他处理器交换数据
  for (int k=0; k<nsend; k++) {
    sndPack[k].nints = 3 * ibsPerProc[k];
    sndPack[k].nreals = 0;

    for (int i=0; i < ibsPerProc[k]; i++) {
      int ii = ibProcMap[k][i];

      sndPack[k].nints += nintsSend[ii];
      sndPack[k].nreals += nrealsSend[ii];
    }
    sndPack[k].intData = (int*)malloc(sizeof(int) * sndPack[k].nints);
    sndPack[k].realData = (double*)malloc(sizeof(double) * sndPack[k].nreals);

    int n = 0;
    int m = 0;

    for (int i=0; i < ibsPerProc[k]; i++) {
      int ii = ibProcMap[k][i];

      sndPack[k].intData[n++] = obblist[ii].send_tag;
      sndPack[k].intData[n++] = nintsSend[ii];
      sndPack[k].intData[n++] = nrealsSend[ii];

      for (int j=0; j < nintsSend[ii]; j++)
        sndPack[k].intData[n++] = int_data[ii][j];

      for (int j=0; j < nrealsSend[ii]; j++)
        sndPack[k].realData[m++] = real_data[ii][j];
    }
  }
  pc->sendRecvPackets(sndPack, rcvPack);  //发送和接收包

  // Reset MeshBlock data structures 重置网格块数据结构（置零）
  for (int ib=0;ib<nblocks;ib++) {
    auto &mb = mblocks[ib];
    mb->nsearch = 0;

    if (mb->xsearch) {
      TIOGA_FREE(mb->xsearch);
      mb->xsearch=NULL;
    }
    if (mb->isearch) {
      TIOGA_FREE(mb->isearch);
      mb->isearch=NULL;
    }
    if (mb->tagsearch) {
      TIOGA_FREE(mb->tagsearch);
      mb->tagsearch=NULL;
    }
    if (mb->donorId) {
      TIOGA_FREE(mb->donorId);
      mb->donorId=NULL;
    }
    if (mb->res_search) {
      TIOGA_FREE(mb->res_search);
      mb->res_search=NULL;
    }
   if (at_points==1) {
     if (mb->rst) {
       TIOGA_FREE(mb->rst);
       mb->rst=NULL;
     }
    }
#ifdef TIOGA_HAS_NODEGID
   mb->gid_search.clear();
#endif
  }

#ifdef TIOGA_HAS_NODEGID
  int nintsPerNode = 3;
#else
  int nintsPerNode = 1;
#endif

  // Loop through recv packets and estimate search array sizes in MeshBlock data
  // structures 通过recv数据包进行循环，并估计网格块数据结构中的搜索数组大小
  std::vector<int> nintsRecv(nobb);
  std::vector<int> nrealsRecv(nobb);
  for (int k=0; k<nrecv; k++) {
    int m = 0;

    for (int i=0; i < ibsPerProc[k]; i++) {
      int key = rcvPack[k].intData[m++];
      int ii = intBoxMap[key];
      int ib = obblist[ii].iblk_local;
      auto& mb = mblocks[ib];

      nintsRecv[ii] = rcvPack[k].intData[m++];
      nrealsRecv[ii] = rcvPack[k].intData[m++];

      // Skip the nodal data (indices and global IDs) 跳过节点数据（索引和全局ID）
      m += nintsRecv[ii];
      // Adjust for node GIDs if available 如果可用，请调整节点gid
      nintsRecv[ii] /= nintsPerNode;
      mb->nsearch += nintsRecv[ii];
    }
  }

  // Resize MeshBlock array sizes 调整网格块数组大小
  for (int ib=0;ib<nblocks;ib++) {
    auto &mb = mblocks[ib];
    if (mb->nsearch < 1) continue;
    mb->xsearch = (double*)malloc(sizeof(double) * 3 * mb->nsearch);
    mb->res_search = (double*)malloc(sizeof(double) * mb->nsearch);
    mb->isearch = (int*)malloc(3 * sizeof(int) * mb->nsearch);
    mb->tagsearch = (int*)malloc(sizeof(int) * mb->nsearch);
    mb->donorId = (int*)malloc(sizeof(int) * mb->nsearch);
#ifdef TIOGA_HAS_NODEGID
    mb->gid_search.resize(mb->nsearch);
#endif
    if (at_points==1) mb->rst = (double*)malloc(sizeof(double) * 3 * mb->nsearch);
  }

  //
  // Update search arrays in mesh blocks from recv packets从recv数据包更新网格块中的搜索数组
  //
  std::vector<int> icOffset(nblocks,0); // Index of isearch arrays where next fill happens下一个填充发生的isearch数组的索引
  std::vector<int> dcOffset(nblocks, 0); // Index of xsearch arrays where next fill happens下一次填充发生的xsearch数组索引
  std::vector<int> rcOffset(nblocks, 0); // Index of res_search arrays where next fill happens下一个填充发生的res_search数组的索引
  std::vector<int> igOffset(nblocks, 0); // Index of gid_search where next fill happens下一次填充发生的gid_search索引
  for (int k=0; k < nrecv; k++) {
    int l = 0;
    int m = 0;

    for (int i=0; i < ibsPerProc[k]; i++) {
      int key = rcvPack[k].intData[m++];
      int ii = intBoxMap[key];
      int ib = obblist[ii].iblk_local;
      auto& mb = mblocks[ib];

      int ioff = icOffset[ib];
      int doff = dcOffset[ib];
      int roff = rcOffset[ib];
      int goff = igOffset[ib];

      m += 2; // Skip nints and nreals information 跳过nints和nreals信息
      for (int j=0; j < nintsRecv[ii]; j++) {
        mb->isearch[ioff++] = k;
        mb->isearch[ioff++] = rcvPack[k].intData[m++];
        mb->isearch[ioff++] = obblist[ii].iblk_remote;
        mb->tagsearch[ioff/3-1]=obblist[ii].tag_remote;

#ifdef TIOGA_HAS_NODEGID
        std::memcpy(&(mb->gid_search[goff]), &(rcvPack[k].intData[m]), sizeof(uint64_t));
        m += 2;
        goff++;
#endif
      }

      for (int j = 0; j < nrealsRecv[ii] / 4; j++) {
        for (int mm = 0; mm < 3; mm++)
          mb->xsearch[doff++] = rcvPack[k].realData[l++];
        mb->res_search[roff++] = rcvPack[k].realData[l++];
      }

      icOffset[ib] = ioff;
      dcOffset[ib] = doff;
      rcOffset[ib] = roff;
      igOffset[ib] = goff;
    }
  }

  pc->clearPackets(sndPack, rcvPack);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
  // printf("%d %d\n",myid,mb->nsearch);

  if (int_data) {
    for (int i=0; i<nobb; i++) {
      if (int_data[i]) TIOGA_FREE(int_data[i]);
    }
    TIOGA_FREE(int_data);
  }
  if (real_data) {
    for (int i=0; i<nobb; i++) {
      if (real_data[i]) TIOGA_FREE(real_data[i]);
    }
    TIOGA_FREE(real_data);
  }
}
