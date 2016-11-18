subroutine init_gf_tau(bndrg,gf0)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
!
implicit none
integer, intent(in) :: bndrg(2,nspinor)
real(8), intent(out) :: gf0(nstsv,lr_nw,nkptnrloc)
integer :: iw,ist,isp1,ikloc,ik
real(8) :: nf,arg1,arg2
character*100 :: fname,fname_tot,fspn
!test
integer :: i
real(8) :: tt(lr_nw),y(lr_nw),z(lr_nw),x,d
real(8),external :: spline3_eval
!end test
!
gf0=0.d0
! suppose tau > 0
! assume the green function is diagonal at this moment
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do iw=1,lr_nw
   do ist=1,nstsv
     arg1=evalsvnr(ist,ik)*dreal(lr_w(iw)) 
     arg2=evalsvnr(ist,ik)*(dreal(lr_w(iw))-bhbar) 
     gf0(ist,iw,ikloc)=-1.d0/(exp(arg1)+exp(arg2))
     if (abs(gf0(ist,iw,ikloc)).lt.1.d-16) gf0(ist,iw,ikloc)=0.d0
   enddo !ist
  enddo !iw
 enddo !ik


if (mpi_grid_root((/dim_q/))) then
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do isp1=1,nspinor
   write(fname,'("gf_tau_I0_k",I4.4,"_isp",I1.1)') ik,isp1
   do ist=bndrg(1,isp1),bndrg(2,isp1)
    write(fspn,'("_ist",I3.3)') ist
    fname_tot="Green_func/"//trim(adjustl(fname))//trim(adjustl(fspn))
    open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
    do iw=1,lr_nw
      write(168,'(2G20.8)') dreal(lr_w(iw))/ha2ev,gf0(ist,iw,ikloc)
    enddo
    close(168)
   enddo !ist
  enddo !isp1
 enddo !ikloc
endif

! test
!if (mpi_grid_root()) then
! tt(:)=dreal(lr_w(:))
! y(:)=gf0(4,:,1)
! call spline3_coef(lr_nw-1,tt,y,z)
! write(*,*) "cubic spline:"

! do i=1,4*lr_nw-3
!  x=dreal(lr_w(lr_nw))/(4*lr_nw-4)*(i-1)
!  d=spline3_eval(lr_nw-1,tt,y,z,x)
!  write(*,*) x,d
! enddo
!endif

return
end subroutine
