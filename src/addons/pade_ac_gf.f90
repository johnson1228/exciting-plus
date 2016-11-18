subroutine pade_ac_gf(bndrg,gfniw,gf_tau)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
!
integer,intent(in) :: bndrg(2,nspinor)
integer,intent(in) :: gfniw
real(8),intent(in) :: gf_tau(nstsv,lr_nw,nkptnrloc)
!
integer :: iw1,iw2,ist1,isp1,ikloc,n,ik
complex(8),allocatable :: gf_iw(:,:,:)
complex(8),allocatable :: gf_w(:,:,:)
complex(8),allocatable :: gn(:,:)
complex(8),allocatable :: sigcoe(:,:,:)
complex(8) :: gf_tauc(nstsv,lr_nw,nkptnrloc)
complex(8) :: enk,enk1,dSe,wn(gfniw)
character*100 :: fname,fspn,fname_tot
complex(8),external :: ac_func
!
allocate(gn(gfniw,gfniw))
allocate(sigcoe(gfniw,nstsv,nkptnrloc))
allocate(gf_iw(nstsv,gfniw,nkptnrloc))
allocate(gf_w(nstsv,gfniw,nkptnrloc))
!
sigcoe=zzero
gf_w=zzero
gf_iw=zzero
gn=zzero
wn=zzero
gf_tauc=zzero

gf_tauc(:,:,:)=dcmplx(gf_tau(:,:,:))

!gf(w)
call ft_tw3(1,-1,nstsv,lr_nw,nkptnrloc,gfniw,gf_tauc,gf_iw)

! iw mesh
n=0
do iw1=-int((gfniw-1)/2),int(gfniw/2)
 n=n+1
 wn(n)=zi*(2.d0*iw1+1)*pi/bhbar
enddo

do ikloc=1,nkptnrloc
  do ist1=qpnb(1),qpnb(2)
   do iw1=1,gfniw
    do iw2=iw1,gfniw
     if (iw1.eq.1) then
      gn(1,iw2)=gf_iw(ist1,iw2,ikloc)
     else
      gn(iw1,iw2)=(gn(iw1-1,iw1-1)-gn(iw1-1,iw2))/&
               &(wn(iw2)-wn(iw1-1))/gn(iw1-1,iw2)
     endif
    enddo !iw2
    ! gn(i,i) are the coefficients of the fitted self-energy
    sigcoe(iw1,ist1,ikloc)=gn(iw1,iw1)
   enddo !iw1
  enddo !ist1
enddo !ikloc

do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do iw1=1,gfniw
   enk=dcmplx(dimag(wn(iw1)),lr_eta)
   do ist1=qpnb(1),qpnb(2)
    gf_w(ist1,iw1,ikloc)=ac_func(gfniw,sigcoe(:,ist1,ikloc),wn,enk)
   enddo !ist1
  enddo !iw1
enddo !ikloc

if (mpi_grid_root((/dim_q/))) then
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do isp1=1,nspinor
   write(fname,'("gf_w1_k",I4.4,"_isp",I1.1)') ik,isp1
   do ist1=bndrg(1,isp1),bndrg(2,isp1)
     write(fspn,'("_ist",I3.3)') ist1
     fname_tot="Green_func/"//trim(adjustl(fname))//trim(adjustl(fspn))
     open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
     do iw1=1,gfniw
!    enk=dcmplx(evalsvnr(ist1,ik)+dimag(wn(iw1)),lr_eta)
      enk=dcmplx(dimag(wn(iw1)),lr_eta)
      write(168,'(3G18.8)') dreal(enk)*ha2ev,&
                  &dreal(gf_w(ist1,iw1,ikloc)),&
                  &dimag(gf_w(ist1,iw1,ikloc))
     enddo
     close(168)
   enddo !ist1
  enddo !isp1
 enddo !ikloc
endif

deallocate(gn,sigcoe)
deallocate(gf_iw,gf_w)
return
end subroutine
