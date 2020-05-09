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

#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <utility>
#include <cassert>
#include "codetypes.h"
#include "tioga.h"
using namespace TIOGA;
extern "C"{
	;
}

int obbIntersectCheck(double vA[3][3], double xA[3], double dxA[3],
	double vB[3][3], double xB[3], double dxB[3]);

void tioga::exchangeBoxes(void)
{
  int *sndMap; //发送图
  int *rcvMap; //接收图
  int nsend; //发送的数量
  int nrecv; //接受的数量
  PACKET *sndPack,*rcvPack; //发送和接收打包

  std::vector<int> nbPerProc(numprocs); // Number of chunks per processor 每个处理器的块数（进程个数）
  std::vector<int> obSizePerProc(numprocs); // Number of real data (OBBs) per processor 每个处理器的real数据(OBBs)个数（进程个数）

  MPI_Allgather(&nblocks, 1, MPI_INT, nbPerProc.data(), 1, MPI_INT, scomm); //并行

  // Total number mesh chunks across all procs 所有进程的总网格块数
  int ntotalblks = std::accumulate(nbPerProc.begin(),nbPerProc.end(),0); //从头到尾检查每个进程的网格块

  std::vector<int> alltags(ntotalblks); // Mesh tags for all blocks across all procs 所有处理器的所有块的网格标记
  std::vector<int> displs(numprocs+1);  // Offsets for tags per proc 每个处理器标签的偏移量
  std::vector<bool> sendFlag(numprocs,false); // Flag indicating send/recv from this proc 从这个处理器发送/接收的标志
  displs[0] = 0; //偏移量=0
  obSizePerProc[0]=nbPerProc[0]*16; //？ 每个处理器的实际数据数（obbs）
  for (int i=1; i <= numprocs; i++) {
    displs[i] = displs[i-1] + nbPerProc[i-1]; //前一个处理器标签的偏移量+前一个处理器的块数
    if (i < numprocs) obSizePerProc[i]=nbPerProc[i]*16;
  }

  MPI_Allgatherv(mytag.data(), nblocks, MPI_INT, alltags.data(),
                 nbPerProc.data(), displs.data(), MPI_INT, scomm); //并行

  int maxtag = -1;//网格块数
  //for (auto itag: alltags)
  for(int i=0;i<ntotalblks;i++) { //ntotalblks：所有进程的总网格块数目
    int itag=abs(alltags[i]); //网格块标签个数
    if (maxtag < itag) maxtag = itag;
  }
  int mxtgsqr = maxtag * maxtag;

  displs[0] = 0;
  for (int i=1; i <= numprocs; i++) {
    displs[i] = displs[i-1] + nbPerProc[i-1]*16; //i进程的偏移量
  }
  
  std::vector<double> myOBBdata(nblocks*16);//我的OBB数据（网格块数）
  std::vector<double> allOBBdata(ntotalblks*16); //全部的OBB数据（从头到尾检查每个进程的网格块）

  int m=0;
  for (int ib=0; ib < nblocks; ib++) {
    myOBBdata[m++] = (double)mytag[ib];
    for (int i=0; i < 3; i++)
      for (int j=0; j< 3; j++)
	myOBBdata[m++] = mblocks[ib]->obb->vec[i][j];
    for (int i=0; i< 3; i++)
      myOBBdata[m++] = mblocks[ib]->obb->xc[i];
    for (int i=0; i< 3; i++)
      myOBBdata[m++] = mblocks[ib]->obb->dxc[i];
  }

  MPI_Allgatherv(myOBBdata.data(), nblocks*16, MPI_DOUBLE, allOBBdata.data(),
                 obSizePerProc.data(), displs.data(), MPI_DOUBLE, scomm);

  // Determine total number of OBBs received 确定收到的obbs总数
  int nobb = ntotalblks;

  // Store all received OBBs in a temporary list将所有收到的obbs存储在临时列表中
  std::vector<OBB> obbRecv(nobb);
  std::vector<int> obbID(nobb); // Mesh tag corresponding to the OBB  对应于OBB的网格标记
  std::vector<int> obbProc(nobb); // Proc ID corresponding to OBB  与OBB相对应的处理器id
  m=0;
  for(int k=0,ix=0;k < numprocs;k++) {
    for (int n=0;n < nbPerProc[k];n++)
      {
	obbProc[ix]=k;
	obbID[ix]=(int)(allOBBdata[m++]+0.5);
	for (int i=0; i<3; i++)
	  for (int j=0; j<3; j++)
	    obbRecv[ix].vec[i][j] = allOBBdata[m++]; //顶点

      for (int i=0; i<3; i++)
        obbRecv[ix].xc[i] = allOBBdata[m++]; //

      for (int i=0; i<3; i++)
        obbRecv[ix].dxc[i] = allOBBdata[m++];//

      ix++;
    }
  }

  // Mapping of (local block_id, remote OBB block_id) for every intersected pair 每个相交对的（本地块_id，远程obb块_id）映射
  std::vector<std::pair<int, int>> intersectIDs;
  // Counter tracking number of intersected blocks per partition每个分区交叉块的计数器跟踪数
  std::vector<int> obPerProc(numprocs,0);

  // Reset sendflags 重置发送标志
  std::fill(sendFlag.begin(), sendFlag.end(), false);
  //
  // Check for intersection of OBBs检查obbs的相交
  //
  nsend=nrecv=numprocs; 
  for (int ob=0; ob < nobb; ob++) {
    for (int ib=0; ib < nblocks; ib++) {
      auto& mb = mblocks[ib];
      int meshtag = mb->getMeshTag();
        if (abs(obbID[ob]) == meshtag) continue;   //保证为不同网格块上的OBB

        if ( obbIntersectCheck(
               mb->obb->vec, mb->obb->xc, mb->obb->dxc,
               obbRecv[ob].vec, obbRecv[ob].xc, obbRecv[ob].dxc) ||
             obbIntersectCheck(
               obbRecv[ob].vec, obbRecv[ob].xc, obbRecv[ob].dxc,
               mb->obb->vec, mb->obb->xc, mb->obb->dxc)) {
          int overlap_present=1; //重叠
          if (obbID[ob] < 0 || mytag[ib] < 0) {
            mb->check_intersect_p4est(&obbProc[ob],&overlap_present);
          }
          // If there is an intersection, store the index pair, increment number 如果有交集，则存储该处理器的索引对、增加交叉口数量，并激活发送标志
          // of intersections for this processor, and activate send flag
          if (overlap_present) {
          intersectIDs.push_back(std::make_pair(ib, ob)); //相交的ID
          obPerProc[obbProc[ob]]++; //相交对个数
          sendFlag[obbProc[ob]] = true; //激活发送标志
          }
      }
    }
  }
  int new_send = std::count(sendFlag.begin(), sendFlag.end(), true);
  assert(new_send <= nsend);
  // Populate send and recv maps  填充发送和发送映射 
  std::map<int,int> invMap;
  nsend = nrecv = new_send;

  // allocate sndMap and rcvMap  分配sndmap和rcvmap
  sndMap = (int*)malloc(sizeof(int) * nsend);
  rcvMap = (int*)malloc(sizeof(int) * nrecv);
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);

  for (int p=0, ip=0; p < numprocs; p++)
    if (sendFlag[p]) {
      sndMap[ip] = p;
      rcvMap[ip] = p;
      invMap[p] = ip;  //逆映射
      ip++;
    }

  // clear packets before nsend and nrecv are modified in pc->setMap    在pc->setmap中修改nsend和nrecv之前的清除数据包
  pc->setMap(nsend,nrecv,sndMap,rcvMap); //设置发送和接收MAP
  pc->initPackets(sndPack,rcvPack); //发送包和接收包初始化


  // if (obblist) TIOGA_FREE(obblist);
  // obblist = (OBB*) malloc(sizeof(OBB) * intersectIDs.size());
  intBoxMap.clear();
  ibsPerProc.clear();
  ibProcMap.clear();
  obblist.clear();
  obblist.resize(intersectIDs.size());
  ibsPerProc.resize(nsend);
  ibProcMap.resize(nsend);

  // Determine packet sizes and reallocate arrays 确定分组大小和重新分配数组
  for(int k=0; k < nsend; k++) {
    sndPack[k].nints = 3 * obPerProc[sndMap[k]];
    sndPack[k].nreals = sndPack[k].nints * 6;
    sndPack[k].intData = (int*) malloc(sizeof(int)*sndPack[k].nints);
    sndPack[k].realData = (double*) malloc(sizeof(double)*sndPack[k].nreals);
    ibsPerProc[k] = obPerProc[sndMap[k]];
    ibProcMap[k].resize(ibsPerProc[k]);
  }

  // Array tracking indices for populating reduced OBBs  用于填充减少的obbs的阵列跟踪索引
  std::vector<int> idxOffset(nsend,0);
  for (size_t i=0; i<intersectIDs.size(); i++){
    auto ids = intersectIDs[i]; //相交对的ID(网格块ib,OBB块的ob)
    int ib = ids.first;           // Block ID of the local mesh block  局部网格块的块id
    int ob = ids.second;          // Index of the intersected block in OBB list  obb列表中相交块的索引
    int k = invMap[obbProc[ob]];  // Index in sndMap for this proc ID  此proc ID在sndMap的索引
    auto& mb = mblocks[ib];       // Mesh block data object 网格块数据对象
    int ip = obbProc[ob];  //ob所在进程的Id

    int ioff = idxOffset[k];      // Index to fill in sndPack  填写sndPack的索引
    int roff = idxOffset[k] * 6;

    int key_recv = mxtgsqr * ip + maxtag * (mtags[ib] - 1) + obbID[ob] - 1;
    int key_send = mxtgsqr * myid + maxtag * (obbID[ob]-1) + (mtags[ib]-1);
    intBoxMap[key_recv] = i; //i是相交对的id
    ibProcMap[k][ioff] = i;
    obblist[i].comm_idx = k; //k：此proc ID在sndMap的索引
    obblist[i].iblk_local = ib; //ib：局部网格块的块id
    // obblist[i].iblk_remote = obbID[ob];
    obblist[i].send_tag = key_send; //发送标志
    obblist[i].recv_tag = key_recv; //接收标志

    sndPack[k].intData[3*ioff] = key_send; // mb->getMeshTag(); 发送标志
    sndPack[k].intData[3*ioff+1] = ib; 
    sndPack[k].intData[3*ioff+2] = mb->getMeshTag();
    //mb->getReducedOBB2(&obbRecv[ob], &(sndPack[k].realData[roff]));
    mb->getReducedOBB(&obbRecv[ob], &(sndPack[k].realData[roff]));  //得到简化obb

    for(int ii=0; ii<3; ii++)
      for(int j=0; j<3; j++)
        obblist[i].vec[ii][j] = obbRecv[ob].vec[ii][j];

    // Increment index offset for next fill  下一次填充的增量索引偏移量
    idxOffset[k]++;
    // std::cout << "# " << myid << " " << obbProc[ob] << " "
    //           << mtags[ib] << " " << obbID[ob] << " "
    //           << std::setw(3) << key_recv << " "
    //           << std::setw(3) << key_send << std::endl;
  }
  pc->sendRecvPackets(sndPack,rcvPack);  //发送recv包

  for (int k=0; k<nrecv; k++) {
    int m=0;

    for (int n=0; n< rcvPack[k].nints; n+=3) {
      int key = rcvPack[k].intData[n]; //接收包中的接收标志
      int ii = intBoxMap[key]; //接收索引
      obblist[ii].iblk_remote = rcvPack[k].intData[n+1]; //接收包中OBB的局部网格块的块id
      obblist[ii].tag_remote = rcvPack[k].intData[n+2];

      for (int i=0; i<3; i++)
        obblist[ii].xc[i] = rcvPack[k].realData[m++]; //接收包中OBB的边框中心

      for (int i=0; i<3; i++)
        obblist[ii].dxc[i] = rcvPack[k].realData[m++]; //接收包中OBB的边框对角线（一半）

    }
  }

  pc->clearPackets(sndPack,rcvPack);
  //
  // Free local memory
  //
  TIOGA_FREE(sndMap);
  TIOGA_FREE(rcvMap);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
}
