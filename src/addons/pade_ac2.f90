subroutine pade_ac2(gf0_tau)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
!
complex(8),intent(in) :: gf0_tau(lr_nw,nstsv,nkptnrloc)
!
integer :: iw1,iw2,ist1,isp1,ikloc,n,ik
complex(8),allocatable :: gf0_iw(:,:,:)
complex(8),allocatable :: gf0_w(:,:,:)
complex(8),allocatable :: gn(:,:),ctmp1(:,:,:),ctmp2(:,:,:)
complex(8),allocatable :: sigcoe(:,:,:)
complex(8) :: enk,enk1,dSe,wn(niw)
character*100 :: fname,fspn,fname_tot
complex(8),external :: ac_func
!
allocate(gn(niw,niw))
allocate(sigcoe(niw,nstsv,nkptnrloc))
allocate(ctmp1(niw,nstsv,nkptnrloc))
allocate(ctmp2(nstsv,nkptnrloc,lr_nw))
allocate(gf0_iw(nstsv,nkptnrloc,niw))
allocate(gf0_w(nstsv,nkptnrloc,niw))
!
sigcoe=zzero
gf0_w=zzero
gf0_iw=zzero
gn=zzero
ctmp1=zzero
ctmp2=zzero
wn=zzero

!change order of indices
do iw1=1,lr_nw
 ctmp2(:,:,iw1)=gf0_tau(iw1,:,:)
enddo

!do ikloc=1,nkptnrloc
! ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! write(fname,'("gf_tau_k",I4.4)') ik
! do ist1=qpnb(1),qpnb(2)
!  write(fspn,'("_ist",I3.3)') ist1
!  fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
!  open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
!  do iw1=1,lr_nw
!    write(168,'(3G18.8)') dreal(lr_w(iw1))/ha2ev,&
!                &dreal(gf0_tau(iw1,ist1,ikloc)),&
!                &dimag(gf0_tau(iw1,ist1,ikloc))
!  enddo
!  close(168)
! enddo !ist1
!enddo !ikloc

call ft_tw3(1,-1,nstsv,nkptnrloc,lr_nw,niw,ctmp2,gf0_iw)


! iw mesh
n=0
do iw1=-int((niw-1)/2),int(niw/2)
 n=n+1
 wn(n)=zi*(2.d0*iw1+1)*pi/bhbar
enddo

do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 write(fname,'("gf_iw_k",I4.4)') ik
 do ist1=qpnb(1),qpnb(2)
  write(fspn,'("_ist",I3.3)') ist1
  fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
  open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
  do iw1=1,niw
    enk1=zone/(wn(iw1)-evalsvnr(ist1,ik))
    write(168,'(5G18.8)') dimag(wn(iw1))*ha2ev,&
                &dreal(gf0_iw(ist1,ikloc,iw1))/ha2ev,&
                &dimag(gf0_iw(ist1,ikloc,iw1))/ha2ev,&
                &dreal(enk1)/ha2ev,dimag(enk1)/ha2ev
  enddo
  close(168)
 enddo !ist1
enddo !ikloc



do ikloc=1,nkptnrloc
  do ist1=1,nstsv
   do iw1=1,niw
    do iw2=iw1,niw
     if (iw1.eq.1) then
      gn(1,iw2)=gf0_iw(ist1,ikloc,iw2)
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
  do ist1=1,nstsv
   do iw1=1,niw
    enk=dcmplx(evalsvnr(ist1,ik)+dimag(wn(iw1)),lr_eta)
    gf0_w(ist1,ikloc,iw1)=ac_func(niw,sigcoe(:,ist1,ikloc),wn,enk)
!    if (ik.eq.11.and.ist1.eq.qpnb(1)) write(*,*) "gf0_w:",gf0_w(ist1,ikloc,iw1)
     if (ik.eq.11.and.ist1.eq.qpnb(1)) write(*,*) "sigcoe:",sigcoe(iw1,ist1,ikloc)
   enddo
  enddo
enddo

do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 write(fname,'("gf_w_k",I4.4)') ik
 do ist1=qpnb(1),qpnb(2)
  write(fspn,'("_ist",I3.3)') ist1
  fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
  open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
  do iw1=1,niw
    enk=dcmplx(evalsvnr(ist1,ik)+dimag(wn(iw1)),lr_eta)
    enk1=zone/(enk-evalsvnr(ist1,ik)+lr_eta*zi)
    write(168,'(5G18.8)') dreal(enk)*ha2ev,&
                &dreal(gf0_w(ist1,ikloc,iw1))/ha2ev,&
                &dimag(gf0_w(ist1,ikloc,iw1))/ha2ev,&
                &dreal(enk1)/ha2ev,dimag(enk1)/ha2ev
  enddo
  close(168)
 enddo !ist1
enddo !ikloc

deallocate(gn,ctmp1,ctmp2,sigcoe)
deallocate(gf0_iw,gf0_w)
return
end subroutine
