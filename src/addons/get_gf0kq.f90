subroutine get_gf0kq(indx,gf0_tau,gf0kq_tau)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
!
integer, intent(in) :: indx
real(8),intent(in) :: gf0_tau(nstsv,lr_nw,nkptnrloc)
!
real(8),intent(out) :: gf0kq_tau(nstsv,lr_nw,nkptnrloc)
!
integer,allocatable :: jkmap(:,:)
logical :: need_to_receive
integer :: ik,jk,ikstep,nkstep,jkloc,i,j,tag
!
!indx: +1   -- k+q
!      +3   -- k-q
nkstep=nkptnrloc
call mpi_grid_bcast(nkstep,dims=(/dim_k/))
allocate(jkmap(2,0:mpi_grid_dim_size(dim_k)-1))
do ikstep=1,nkstep
  jkmap=-1
  need_to_receive=.false.
  ! if this processor has a k-point for this step
  if (ikstep.le.nkptnrloc) then
    ! k-point global index
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikstep)
    ! k+q point global index
    jk=idxkq(indx,ik)
!    write(*,*) "jk,ik,cpu_i:",jk,ik,mpi_grid_dim_pos(dim_k)
    ! find index of processor and a local jk index at CPU j
    jkloc=mpi_grid_map(nkptnr,dim_k,x=j,glob=jk)
    ! store index of processor j, and local index of k+q point jkloc
    jkmap(1,mpi_grid_dim_pos(dim_k))=j
    jkmap(2,mpi_grid_dim_pos(dim_k))=jkloc
    ! make a local copy if jk is on the same processor
    if (j.eq.mpi_grid_dim_pos(dim_k)) then
      gf0kq_tau(:,:,ikstep)=gf0_tau(:,:,jkloc)
    else
      need_to_receive=.true.
    endif
  endif
  call mpi_grid_reduce(jkmap(1,0),2*mpi_grid_dim_size(dim_k),dims=(/dim_k/),all=.true.,op=op_max)
  ! check who needs k-point which is stored on this processor
  do i=0,mpi_grid_dim_size(dim_k)-1
    if (jkmap(1,i).eq.mpi_grid_dim_pos(dim_k).and.mpi_grid_dim_pos(dim_k).ne.i) then
     jkloc=jkmap(2,i)
     ! send to proc i
     tag=(ikstep*mpi_grid_dim_size(dim_k)+i)*10
     call mpi_grid_send(gf0_tau(1,1,jkloc),lr_nw*nstsv,(/dim_k/),(/i/),tag)
    endif
  enddo

  if (need_to_receive) then
    j=jkmap(1,mpi_grid_dim_pos(dim_k))
    tag=(ikstep*mpi_grid_dim_size(dim_k)+mpi_grid_dim_pos(dim_k))*10
    call mpi_grid_recieve(gf0kq_tau(1,1,ikstep),lr_nw*nstsv,(/dim_k/),(/j/),tag)
  endif
enddo !ikstep
deallocate(jkmap)
!call mpi_grid_barrier((/dim_k/))
call mpi_grid_barrier()

return
end subroutine
