subroutine genvxcnk(nbnd,vxcnk)
use modmain
use mod_nrkp
use mod_linresp,  only : qpnb
implicit none
integer, intent(in) :: nbnd
real(8), intent(out) :: vxcnk(nbnd,nspinor,nkptnrloc) 
!
integer jst,ik,ikloc,ias,is,ic,l1,l2,l3,io1,io2,ispn,lm1,lm2,lm3,ig,ir
complex(8) zt
real(8), allocatable :: rf(:)
real(8), allocatable :: rhonkmt(:,:,:)
real(8), allocatable :: rhonkir(:)
complex(8), allocatable :: gntmp(:,:)
complex(8), allocatable :: zfft(:)
!
real(8), external :: rfinp
integer :: isp1, ibnd, fbnd
!
vxcnk=0.d0
allocate(gntmp(lmmaxapw,lmmaxapw)) 
allocate(zfft(ngrtot))
allocate(rhonkmt(lmmaxvr,nrmtmax,natmtot))
allocate(rhonkir(ngrtot))
allocate(rf(nrmtmax))

do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 do isp1 = 1, nspinor
  if (isp1.eq.1) then
   ibnd=qpnb(1)
   fbnd=qpnb(2)
  elseif (isp1.eq.2) then   !spin-polarized
   ibnd=qpnb(1)+int(nstsv/2)
   fbnd=qpnb(2)+int(nstsv/2)
  endif
  
  do jst=ibnd,fbnd
    rhonkmt=0.d0
    rhonkir=0.d0
    do lm3=1,lmmaxvr
      gntmp(:,:)=gntyry(lm3,:,:)   !Gaunt coefficient
      l3=lm2l(lm3)
      do ias=1,natmtot
        is=ias2is(ias)
        ic=ias2ic(ias)
        rf=0.d0
        do l1=0,lmaxapw; do io1=1,nufr(l1,is)
          do l2=0,lmaxapw; do io2=1,nufr(l2,is)
            if (mod(l1+l2+l3,2).eq.0) then
              zt=zzero
              do ispn=1,nspinor
                do lm2=l2**2+1,(l2+1)**2
                  do lm1=l1**2+1,(l1+1)**2
                    zt=zt+dconjg(wfsvmtnrloc(lm1,io1,ias,ispn,jst,ikloc))*&
                      &wfsvmtnrloc(lm2,io2,ias,ispn,jst,ikloc)*gntmp(lm1,lm2)
                  enddo !lm1
                enddo !lm2
              enddo !ispn
              rf(:)=rf(:)+dreal(zt)*ufr(:,l1,io1,ic)*ufr(:,l2,io2,ic)
            endif
          enddo; enddo !l2, io2
        enddo; enddo !l1, io1
        rhonkmt(lm3,:,ias)=rf(:)
      enddo !ias
    enddo !lm3
    do ispn=1,nspinor
      zfft=zzero
      do ig=1,ngknr(ikloc) 
        zfft(igfft(igkignr(ig,ikloc)))=wfsvitnrloc(ig,ispn,jst,ikloc)
      enddo
      call zfftifc(3,ngrid,1,zfft) 
      do ir=1,ngrtot
        rhonkir(ir)=rhonkir(ir)+(abs(zfft(ir))**2)/omega
      enddo
    enddo !ispn
    if (rho_val.or.pt_core) then
     vxcnk(jst-ibnd+1,isp1,ikloc)=rfinp(1,rhonkmt,vxcmt_val,rhonkir,vxcir_val)
    else
     vxcnk(jst-ibnd+1,isp1,ikloc)=rfinp(1,rhonkmt,vxcmt,rhonkir,vxcir)
    endif
  enddo !jst
 enddo ! isp1
enddo !ikloc
deallocate(gntmp,zfft,rhonkmt,rhonkir,rf)

return
end subroutine
