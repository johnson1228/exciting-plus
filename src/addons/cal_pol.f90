subroutine cal_pol(iq,nwloc,gf0_tau,pol)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
!
integer,intent(in) :: iq
integer,intent(in) :: nwloc
real(8),intent(in) :: gf0_tau(nstsv,lr_nw,nkptnrloc)
complex(8),intent(inout) :: pol(ngq(iq),ngq(iq),lr_nw)
!
integer :: ig1,ig2,iw,i,ikloc,ist1,ist2,ik
integer :: ierr,n,iw1,iwloc
real(8),allocatable :: gf0kq_tau(:,:,:)
real(8),allocatable :: wt(:)
character*100 :: fname
!
if (allocated(megqblh2)) deallocate(megqblh2)
allocate(megqblh2(nstsv*nstsv,ngq(iq)))
allocate(gf0kq_tau(nstsv,lr_nw,nkptnrloc))
allocate(wt(nstsv*nstsv))
megqblh2=zzero
gf0kq_tau=0.d0
wt=0.d0
pol=zzero
!
! get G(k+q) from other processors
call get_gf0kq(1,gf0_tau,gf0kq_tau)
!
! calculate polarization function
do ikloc=1,nkptnrloc
 do iwloc=1,nwloc
  iw=mpi_grid_map(lr_nw,dim_q,loc=iwloc)
  iw1=lr_nw-iw+1
  do i=1,nmegqblh(ikloc)
   ist1=bmegqblh(1,i,ikloc)
   ist2=bmegqblh(2,i,ikloc)
   wt(i)=-gf0_tau(ist1,iw1,ikloc)*gf0kq_tau(ist2,iw,ikloc)
  enddo    
  do ig1=1,ngq(iq)
   megqblh2(:,ig1)=dconjg(megqblh(:,ig1,ikloc))*wt(:)
  enddo
  call zgemm('T','N',ngq(iq),ngq(iq),nstsv*nstsv,zone,&
           &megqblh(1,1,ikloc),nstsv*nstsv,megqblh2(1,1),nstsv*nstsv,&
           &zone,pol(1,1,iw),ngq(iq))
 enddo !iwloc
enddo !ikloc

call mpi_grid_reduce(pol(1,1,1),ngq(iq)*ngq(iq)*lr_nw,dims=(/dim_k,dim_q/), &
                     & all=.true.)

!when nspinor=1, occmax=2
pol(:,:,:)=pol(:,:,:)*occmax/(nkptnr*omega)

!print out pol(G,G',2) for test
!if (mpi_grid_root()) then
!  write(fname,'("pol_tau2_iq",I3.3)') iq
!  open(166,file=trim(adjustl(fname)),form="FORMATTED",status="REPLACE")
!  do ig2=1,ngq(iq)
!   do ig1=1,ngq(iq)
!    write(166,'(4(G16.6))') ig1,ig2,dreal(pol(ig1,ig2,2)),dimag(pol(ig1,ig2,2))
!   enddo
!  enddo
!  close(166)
!endif

deallocate(megqblh2,wt)
deallocate(gf0kq_tau)

return
end subroutine
