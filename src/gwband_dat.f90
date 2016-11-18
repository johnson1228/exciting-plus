subroutine gwband_dat(fnum)
use modmain
use mod_linresp
!
implicit none
!
integer,intent(in) :: fnum
real(8), allocatable :: e_qp(:,:,:),e_ks(:,:,:),dk(:,:),vk(:,:)
integer :: ierr1,ik,ikloc,isp1,i,ib,tmp,exst,ios
real(8) :: dt,vc(3),d1k,d2k,dc12
character*100 :: fspn,fname,fname_tot
integer :: ndk

nbnd=qpnb(2)-qpnb(1)+1
allocate(e_qp(nbnd,nspinor,nkptnr))
allocate(e_ks(nbnd,nspinor,nkptnr))
e_qp=0.d0
e_ks=0.d0
ierr1=0

do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 write(fname,'("Eqp_k",I4.4)')ik

 do isp1=1,nspinor
  if (nspinor.eq.1) then
   fspn=''
  else
   write(fspn,'("_spn",I1.1)') isp1
  endif
  fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
  inquire(file=fname_tot, exist=exst)
  if (exst) then
   open(163,file=fname_tot,action='read',form="FORMATTED",status="old")

   do i=1,nbnd
    read(163,*,iostat=ios) tmp,e_qp(i,isp1,ik),tmp,e_ks(i,isp1,ik)
   enddo !i
   close(163)
  else
   write(*,'("Eqp files are not found!",I4)') ik
   ierr1=1
  endif
 enddo !isp1
enddo !ikloc

call mpi_grid_reduce(ierr1,1,dims=(/dim_k/),all=.true.,op=op_max)
if (ierr1.ne.0) call pstop

call mpi_grid_reduce(e_qp(1,1,1),nbnd*nspinor*nkptnr,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(e_ks(1,1,1),nbnd*nspinor*nkptnr,dims=(/dim_k/),all=.true.)

! generate the gwband_plot using the root cpu
if (mpi_grid_root()) then
 allocate(dk(2,npp1d))
 allocate(vk(3,nkptnr))
 write(fname,'("gwband.dat")')
 do isp1=1,nspinor
  if (nspinor.eq.1) then
   fspn=''
  else
   write(fspn,'("_spn",I1.1)') isp1
  endif
  fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
  open(163,file=fname_tot,form='formatted',status='replace')
  dt=0.d0
  ndk=0
  do i=1,nvp1d-1
   call r3mv(bvec,vvlp1d(:,i+1)-vvlp1d(:,i),vc)
   dc12=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
   do ik=1,nkptnr
    call r3mv(bvec,vklnr(:,ik)-vvlp1d(:,i),vc)
    d1k=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
    call r3mv(bvec,vvlp1d(:,i+1)-vklnr(:,ik),vc)
    d2k=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
    if (abs(d1k+d2k-dc12).lt.1.d-8) then
     ndk=ndk+1
     dk(1,ndk)=dt+d1k
     dk(2,ndk)=dble(ik)
    endif
   enddo !ik
   dt=dt+dc12
   write(fnum,'("dt:",f16.8)') dt
  enddo !i

  if (ndk.eq.0) then
   write(fnum,'("ndk = 0!")')
   call flushifc(fnum)
  else
   write(fnum,'("dt:",f16.8)'),dt
   do i=1,nbnd
    do ik=1,ndk
     write(163,'(3G18.10)') dk(1,ik),e_qp(i,isp1,int(dk(2,ik))),&
          &                   e_ks(i,isp1,int(dk(2,ik)))
    enddo
    write(163,'(" ")')
   enddo
  endif
  close(163)
 enddo !isp1
 write(fnum,'("done with GW band!")')
 deallocate(dk,vk)
endif !root_cpu

deallocate(e_qp,e_ks)
return
end subroutine
