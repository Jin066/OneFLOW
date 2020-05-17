
#include "mpif.h"


//initialize tioga ��ʼ��Tioga

//ntypes = 1;
//nv1 = 6;
//nv2 = 8;
do ib = 1, 2
g = > gr(ib)
if (g%n6 > 0)  then  !g�е�n6��Ա
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
	!call tioga_registergrid_data(meshtag, !< mesh tag for this partition(scalar)�˷����������ǣ�������
		!nnodes, !< number of nodes in this partition(scalar)�˷����еĽڵ�����������
		!x, !< coordinates(1 - D real array with size = (3 * nnodes))���꣨1ά��ʵ�����飬��С = 3 * �ڵ������
		!iblank, !< iblank(integer array with size = (nnodes)) �ڵ��־�����������С = �ڵ������
		!nwbc, !< number of wall boundary nodes   �ڱ߽�ڵ���
		!nobc, !< number of outer boundary nodes  �ⲿ�߽�ڵ���
		!wbcnode, !< index of wall boundary nodes size = (nwbc)�ڱ߽�ڵ�����
		!obcnode, !< index of overset boundary nodes size = (nobc)��߽�ڵ�����
		!ntypes, !< number of type of cells(scalar) ��Ԫ���͸�����������
		!nv, !< number of vertices per cell for the first cell type ��һ����Ԫ������ÿ����Ԫ�Ķ�����
		!ncells, !< number cells of the first type ��һ�����͵ĵ�Ԫ����
		!connectivity, !< connectivity of the first type of cells ��һ�����͵�Ԫ����ͨ��
		!.., !< number of vertices per cell for the second cell type
		!.., !< number of cells of second type
		!..)            !< connectivity of the second type of cells
	!!< ..third, fourth etc
	call tioga_preprocess_grids                  !< preprocess the grids(call again if dynamic) ����Ԥ����
	call cpu_time(t1)
	call tioga_performconnectivity               !< determine iblanking and interpolation patternsȷ�������Ͳ�ֵģʽ
	call cpu_time(t2)

	call mpi_barrier(mpi_comm_world, ierr)
	if (myid == 0) write(6, *) 'connectivity time=', t2 - t1   !���ӵ�ʱ��

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
						call tioga_getdonorcount(ib, dcount, fcount)  !�������׵�Ԫ, tiogaInterface.C
						allocate(receptorInfo(4 * dcount))  !�����ֵ��Ϣ
						allocate(inode(fcount), frac(fcount))  !
						write(6, *) dcount, fcount
						!
						!use this API if you need to interpolate the field variables yourself �����Ҫ���в����ֶα�������ʹ�ô�API
						!
						call tioga_getdonorinfo(gr(ib) % bodytag(1), receptorInfo, inode, frac, dcount)  !tiogaInterface.C
						!
						!> dcount = number of donors �����߸���
						!> receptorInfo = { receptorProcess Id, receptor Index, receptor Block id, number of fractions }*dcount ������Ϣ = { ���崦��ID�����������������ID���������� }*���׸���
						!> inode = indices for each receptor one group after the other ÿ�����������һ���һ��
						!> frac = weights for each receptor one group after the other ÿ�������Ȩ��һ���һ��
						!>
						m = 0
						do i = 1, dcount
							write(1000 + myid, "(20(1x,I6))") ib, (inode(m + j), j = 1, receptorInfo(4 * i) + 1), receptorInfo(4 * i - 3:4 * i)
							m = m + (receptorInfo(4 * i) + 1)
							enddo
							call flush(1000 + myid) !?
							deallocate(receptorInfo, inode, frac) !�ͷ�
							end do

							do ib = 1, 2
								g = > gr(ib)
								call tioga_registersolution(g%bodytag(1), g%q)  !tiogaInterface.C
								enddo
								call tioga_dataupdate_mb(gr(1) % nvar, 'row')    !< update the q - variables(can be called anywhere)����q�������������κεط����ã�
								!< nvar = number of field variables per nodeÿ���ڵ�ĳ���������
								!< if fields are different arrays, you can also ������ǲ�ͬ�����飬����Ϊÿ�������ö��
								!< call this multiple times for each field
								call mpi_barrier(mpi_comm_world, ierr)
								call cpu_time(t2)
								if (myid == 0) write(6, *) 'data update time=', t2 - t1

									!
									!compute error in interpolation �����ֵ���
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
												if (myid == 0) write(6, "(A36)") '-- Interpolation error statistics --' !��ֵ���ͳ��
													if (myid == 0) write(6, "(A15)") 'ProcId    Error'
														call flush(6)
														call mpi_barrier(mpi_comm_world, ierr)
														write(6, "(I4,3x,E15.7)") myid, sqrt(rnorm / g % nv / g % nvar)
														call mpi_barrier(mpi_comm_world, ierr)

														call tioga_writeoutputfiles(g%nvar, 'row') !< write output files, if need be �����Ҫ�����д����ļ���

														200 continue
														call mpi_barrier(mpi_comm_world, ierr)
														call tioga_delete
														call mpi_finalize(ierr)
														!
														end program testTioga
