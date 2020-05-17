
#include "mpif.h"


//initialize tioga 初始化Tioga

//ntypes = 1;
//nv1 = 6;
//nv2 = 8;
do ib = 1, 2
g = > gr(ib)
if (g%n6 > 0)  then  !g中的n6成员
call tioga_registergrid_data_mb(ib, g%bodytag(1), g%nv, g%x, g%iblank, g%nwbc, g%nobc, g%wbcnode, g%obcnode, &
	ntypes, nv1, g%n6, g%ndc6)
else if (g%n8 > 0) then
call tioga_registergrid_data_mb(ib, g%bodytag(1), g%nv, g%x, g%iblank, g%nwbc, g%nobc, g%wbcnode, g%obcnode, &
	ntypes, nv2, g%n8, g%ndc8)
	endif
	enddo
	!
	!example call to tioga_registergrid_data should be like this
	!
	!call tioga_registergrid_data(meshtag, !< mesh tag for this partition(scalar)此分区的网格标记（标量）
		!nnodes, !< number of nodes in this partition(scalar)此分区中的节点数（标量）
		!x, !< coordinates(1 - D real array with size = (3 * nnodes))坐标（1维的实数数组，大小 = 3 * 节点个数）
		!iblank, !< iblank(integer array with size = (nnodes)) 节点标志（整数数组大小 = 节点个数）
		!nwbc, !< number of wall boundary nodes   壁边界节点数
		!nobc, !< number of outer boundary nodes  外部边界节点数
		!wbcnode, !< index of wall boundary nodes size = (nwbc)壁边界节点索引
		!obcnode, !< index of overset boundary nodes size = (nobc)外边界节点索引
		!ntypes, !< number of type of cells(scalar) 单元类型个数（标量）
		!nv, !< number of vertices per cell for the first cell type 第一个单元类型中每个单元的顶点数
		!ncells, !< number cells of the first type 第一种类型的单元个数
		!connectivity, !< connectivity of the first type of cells 第一种类型单元的连通性
		!.., !< number of vertices per cell for the second cell type
		!.., !< number of cells of second type
		!..)            !< connectivity of the second type of cells
	!!< ..third, fourth etc
	call tioga_preprocess_grids                  !< preprocess the grids(call again if dynamic) 网格预处理
	call cpu_time(t1)
	call tioga_performconnectivity               !< determine iblanking and interpolation patterns确定消隐和插值模式
	call cpu_time(t2)

	call mpi_barrier(mpi_comm_world, ierr)
	if (myid == 0) write(6, *) 'connectivity time=', t2 - t1   !连接的时间

		call cpu_time(t1)
		do ib = 1, 2
			g = >gr(ib)
			m = 1
			do i = 1, g%nv
				xt(:) = g % x(3 * i - 2:3 * i)
				do n = 1, g%nvar
					g%q(m) = (xt(1) + xt(2) + xt(3))  !?
					m = m + 1
					enddo
					enddo
					enddo

					do ib = 1, 2
						call tioga_getdonorcount(ib, dcount, fcount)  !搜索贡献单元, tiogaInterface.C
						allocate(receptorInfo(4 * dcount))  !分配插值信息
						allocate(inode(fcount), frac(fcount))  !
						write(6, *) dcount, fcount
						!
						!use this API if you need to interpolate the field variables yourself 如果需要自行插入字段变量，请使用此API
						!
						call tioga_getdonorinfo(gr(ib) % bodytag(1), receptorInfo, inode, frac, dcount)  !tiogaInterface.C
						!
						!> dcount = number of donors 贡献者个数
						!> receptorInfo = { receptorProcess Id, receptor Index, receptor Block id, number of fractions }*dcount 受体信息 = { 受体处理ID，受体索引，受体块ID，分数个数 }*贡献个数
						!> inode = indices for each receptor one group after the other 每个受体的索引一组接一组
						!> frac = weights for each receptor one group after the other 每个受体的权重一组接一组
						!>
						m = 0
						do i = 1, dcount
							write(1000 + myid, "(20(1x,I6))") ib, (inode(m + j), j = 1, receptorInfo(4 * i) + 1), receptorInfo(4 * i - 3:4 * i)
							m = m + (receptorInfo(4 * i) + 1)
							enddo
							call flush(1000 + myid) !?
							deallocate(receptorInfo, inode, frac) !释放
							end do

							do ib = 1, 2
								g = > gr(ib)
								call tioga_registersolution(g%bodytag(1), g%q)  !tiogaInterface.C
								enddo
								call tioga_dataupdate_mb(gr(1) % nvar, 'row')    !< update the q - variables(can be called anywhere)更新q变量（可以在任何地方调用）
								!< nvar = number of field variables per node每个节点的场变量个数
								!< if fields are different arrays, you can also 如果场是不同的数组，可以为每个场调用多次
								!< call this multiple times for each field
								call mpi_barrier(mpi_comm_world, ierr)
								call cpu_time(t2)
								if (myid == 0) write(6, *) 'data update time=', t2 - t1

									!
									!compute error in interpolation 计算插值误差
									!should be machine - zero for the linear
									!problem here
									!
									rnorm = 0d0
									do ib = 1, 2
										g = >gr(ib)
										m = 1
										do i = 1, g%nv
											xt(:) = g % x(3 * i - 2:3 * i)
											do n = 1, g%nvar
												rnorm = rnorm + (g%q(m) - (xt(1) + xt(2) + xt(3)))**2
												m = m + 1
												enddo
												enddo
												enddo
												if (myid == 0) write(6, "(A36)") '-- Interpolation error statistics --' !插值误差统计
													if (myid == 0) write(6, "(A15)") 'ProcId    Error'
														call flush(6)
														call mpi_barrier(mpi_comm_world, ierr)
														write(6, "(I4,3x,E15.7)") myid, sqrt(rnorm / g % nv / g % nvar)
														call mpi_barrier(mpi_comm_world, ierr)

														call tioga_writeoutputfiles(g%nvar, 'row') !< write output files, if need be 如果需要，请编写输出文件。

														200 continue
														call mpi_barrier(mpi_comm_world, ierr)
														call tioga_delete
														call mpi_finalize(ierr)
														!
														end program testTioga
