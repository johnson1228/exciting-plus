subroutine find_degenency(bndrg,evalmap,neval)
use modmain
use mod_nrkp
use mod_linresp
!
implicit none
!
integer,intent(in) :: bndrg(2,nspinor)
integer,intent(out) :: evalmap(nbnd,nspinor,nkptnrloc)
integer,intent(out) :: neval(nspinor,nkptnrloc)
integer :: ikloc,ik,ist,isp1,i,ibnd,fbnd
real(8) :: t1
!

do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 do isp1=1,nspinor
  ibnd=bndrg(1,isp1)
  fbnd=bndrg(2,isp1)
! check lower bound
  if (ibnd.gt.1) then
   t1=abs(evalsvnr(ibnd,ik)-evalsvnr(ibnd-1,ik))
   if (t1.lt.epslat) then
    if (mpi_grid_root()) then
     write(*,'("warning! Degenerate states are not completely included!")')
     write(*,'("ik",I4,1X,"isp1:",I2,"ib",I4)') ik,isp1,ibnd
    endif
   endif
  endif
! check upper bound
  t1=abs(evalsvnr(fbnd,ik)-evalsvnr(fbnd+1,ik))
  if (t1.lt.epslat) then
   if (mpi_grid_root()) then
    write(*,'("warning! Degenerate states are not completely included!")')
    write(*,'("ik",I4,1X,"isp1:",I2,"ib",I4)') ik,isp1,fbnd
   endif
  endif

 enddo !isp1
enddo !ikloc

! if there is only one band, then no need to check
if (bndrg(2,1).eq.bndrg(1,1)) then
 neval(:,:)=1
 do isp1=1,nspinor
  evalmap(:,isp1,:)=bndrg(1,isp1)
 enddo
 return
endif

! initial value
evalmap(:,:,:)=0
neval(:,:)=1

do isp1=1,nspinor
 evalmap(1,isp1,:)=bndrg(1,isp1)
enddo

do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 do isp1=1,nspinor
  ibnd=bndrg(1,isp1)
  fbnd=bndrg(2,isp1)
  do ist=ibnd+1,fbnd
   i=ist-ibnd+1
   t1=abs(evalsvnr(ist,ik)-evalsvnr(ist-1,ik))
   if (t1.gt.epslat) then
    neval(isp1,ikloc)=neval(isp1,ikloc)+1
    evalmap(neval(isp1,ikloc),isp1,ikloc)=ist
   endif
  enddo
  if (evalmap(neval(isp1,ikloc),isp1,ikloc).lt.bndrg(2,isp1)) then
   neval(isp1,ikloc)=neval(isp1,ikloc)+1
   evalmap(neval(isp1,ikloc),isp1,ikloc)=bndrg(2,isp1)+1
  endif
 enddo !isp1
enddo !ikloc

if (mpi_grid_root()) then
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do isp1=1,nspinor
   do i=1,neval(isp1,ikloc)
    write(*,*) "i,isp1,ik,evalmap:",i,isp1,ik,evalmap(i,isp1,ikloc)
   enddo
  enddo
 enddo
endif

return
end subroutine
