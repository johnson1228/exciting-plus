module mod_linresp
implicit none

! dimension for q-vectors
integer, parameter :: dim_q=2

! type of linear response calculation
!   0 : charge response
!   1 : magnetic response (not yet implemented)
integer lrtype
data lrtype/0/

! temporary array used in genchi0blh
complex(8), allocatable :: megqblh2(:,:)

! interval of bands to include in chi0
real(8) chi0_include_bands(2)
data chi0_include_bands/-100.1d0,100.1d0/

! interval of bands to exclude from chi0
real(8) chi0_exclude_bands(2)
data chi0_exclude_bands/100.1d0,-100.1d0/

complex(8), allocatable :: wann_cc(:,:,:)
complex(8), allocatable :: wann_cc2(:,:)

complex(8), allocatable :: chi0loc(:,:,:)
complex(8), allocatable :: chi0wanloc(:,:,:)

! number of energy-mesh points or Matsubara time grid points
integer lr_nw
data lr_nw/201/
! number of energy-mesh points
integer niw
data niw/201/
! number of energy-mesh points for dyson_gf
integer gf_niw
data gf_niw/10001/
! first energy point (Ha)
real(8) lr_w0
data lr_w0/0.d0/
! last energy point (Ha)
real(8) lr_w1
data lr_w1/1.d0/
! energy step
real(8) lr_dw
! energy mesh
complex(8), allocatable :: lr_w(:)
! energy mesh for GW integral
complex(8), allocatable :: lr_wc(:)
! broadening parameter (Ha)
real(8) lr_eta
data lr_eta/0.002d0/
! QP correction: 0 for GW0, 1 for GW for tau-domain GW
integer gw_mode
data gw_mode/0/
! number of steps in GW0 calculation; 1 if G0W0 calculation
integer gw0step
data gw0step/1/
! scGW number of iterations
integer scgwni
data scgwni/1/
! QP correction to bands between qpnb(1) and qpnb(2)
integer qpnb(2)
! number of bands in self-energy matrix
integer nbnd
! number of empty bands used in chi
integer nebd_chi
! number of empty bands used in self-energy
integer nebd_se
! GW calculation type: 
! 0: PPA  1:real-axis integration 2:full self-energy 3:new type of rai
integer caltype
data caltype/0/
! Type of calculating the valence contribution to exxnk:
! 0: in real-space; 1: in reciprocal space
integer exxtype
data exxtype/0/
! beta*hbar in \tau GW method
real(8) bhbar
data bhbar/1046.5846d0/    ! kT = 26 meV
! uniform power mesh for GW method, upm(1)=p, upm(2)=u
! defined in Stan et al. JCP 130,114105 (2009)
integer upm(2)
! energy step used in G0W0 calculation
real(8) del_e
data del_e/0.01d0/
! low and high frequency cutoffs in real-axis integration
! i.e. caltype = 3
real(8) raicut(2)
!restart flag in GW calculations
logical gw_restart
data gw_restart/.false./
!flag for GW band
logical gwband
data gwband/.false./
! spin index in GW band structure
integer gwband_sp
data gwband_sp/1/
! do not update GW energy in writewann.f90
logical eval_udt
data eval_udt/.true./
! core contribution to exx
logical core_exx
data core_exx/.true./
! algorithm of analytic continuation for self-energy
! 1: Pade approximation; 2: continuous pole approximation with Monte-Carlo (3) CPE original
integer ac_sigma
data ac_sigma/1/
! initial guess for alpha in CPE
real(8) cpe_alp0
data cpe_alp0/0.001d0/
! initial guess for sig0 in CPE
real(8) cpe_sig0
data cpe_sig0/0.0d0/
! tolerance for gradient square in CPE
real(8) cpe_gtor
data cpe_gtor/5.d-6/
! range of real-axis integration (-cpe_delta,cpe_delta) in CPE
real(8) cpe_delta 
data cpe_delta/20.d0/
! number of basis functions in (0,cpe_delta)
integer cpe_N
data cpe_N/4096/
! chemical potential used in dyson_gf when updating G
real(8) dyson_mu
data dyson_mu/0.d0/
! tolerance for occupation number in GW-tau calculation
real(8) occ_tor
data occ_tor/0.02d0/
! inverse temperature for the matsubara frequency in eV^-1
real(8) lr_beta
data lr_beta/30.d0/
! .true. if imaginary frequency mesh is required
logical timgw
data timgw/.false./
! first imaginary frequency
real(8) lr_iw0
data lr_iw0/0.d0/
! last imaginary frequency
real(8) lr_iw1
data lr_iw1/80.d0/

real(8) fxca0
data fxca0/0.d0/
real(8) fxca1
data fxca1/0.d0/
integer nfxca
data nfxca/1/
integer fxctype
data fxctype/0/

! high-level switch: compute chi0 and chi in Wannier functions basis
logical wannier_chi0_chi 
data wannier_chi0_chi/.false./
! high-level switch: .true. if chi0 should be multiplied by 2
logical wannier_chi0_afm
data wannier_chi0_afm/.false./

! indices of response functions in global array f_response(:,:,:)
integer, parameter :: f_chi0                 = 1
integer, parameter :: f_chi                  = 2
integer, parameter :: f_chi_scalar           = 3
integer, parameter :: f_chi_pseudo_scalar    = 4
integer, parameter :: f_epsilon_matrix_GqGq  = 5
integer, parameter :: f_epsilon_scalar_GqGq  = 6
integer, parameter :: f_inv_epsilon_inv_GqGq = 7
integer, parameter :: f_epsilon_eff          = 8
integer, parameter :: f_epsilon_eff_scalar   = 9
integer, parameter :: f_sigma                = 10
integer, parameter :: f_sigma_scalar         = 11
integer, parameter :: f_loss                 = 12
integer, parameter :: f_loss_scalar          = 13
integer, parameter :: f_chi0_wann            = 14
integer, parameter :: f_chi_wann             = 15
integer, parameter :: f_epsilon_eff_wann     = 16
integer, parameter :: f_sigma_wann           = 17
integer, parameter :: f_loss_wann            = 18
integer, parameter :: f_epsilon_inv_GqGq     = 19

integer, parameter :: nf_response            = 19
complex(8), allocatable :: f_response(:,:,:)

complex(8), allocatable :: u4(:,:,:,:)
logical screenu4
data screenu4/.true./

complex(8), allocatable :: gw_self_energy(:,:,:,:)
complex(8), allocatable :: self_energy_x(:,:)
contains

subroutine genchi0blh(ikloc,ngq,w,chi0w)
use modmain
use mod_nrkp
use mod_expigqr
implicit none
! arguments
integer, intent(in) :: ikloc
integer, intent(in) :: ngq
complex(8), intent(in) :: w
complex(8), intent(out) :: chi0w(ngq,ngq)
! local variables
logical l1
logical gw_flag,flag
integer i,ist1,ist2,ik,jk,ig
integer nst,itask
real(8) t1,t2
complex(8), allocatable :: wt(:)
logical, external :: bndint
! 
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
allocate(wt(nmegqblh(ikloc)))
if (allocated(megqblh2)) deallocate(megqblh2)
allocate(megqblh2(nstsv*nstsv,ngq))
wt(:)=zzero

! check if G0W0 calculation is performed
do itask=1,ntasks
 if (tasks(itask).eq.802) then
  gw_flag=.true.
  goto 222
 endif
enddo

222 continue
nst=(int(chgval/2.0)+nebd_chi)*nspinor
nst=min(nst,nstsv)

do i=1,nmegqblh(ikloc)
  ist1=bmegqblh(1,i,ikloc)
  ist2=bmegqblh(2,i,ikloc)

  flag=(((ist1.gt.nst/nspinor).and.(ist1.le.nstsv/nspinor)).or. &
        ((ist2.gt.nst/nspinor).and.(ist2.le.nstsv/nspinor)))
  if (gw_flag.and.flag) cycle
  flag=((ist1.gt.(nst/nspinor+nstsv/nspinor)).or. &
        (ist2.gt.(nst/nspinor+nstsv/nspinor)))
  if (gw_flag.and.flag) cycle

! default : include all interband transitions         
  l1=.true.
! cRPA case : don't include bands in energy window [crpa_e1,crpa_e2]
  if (bndint(ist1,evalsvnr(ist1,ik),chi0_exclude_bands(1),&
      chi0_exclude_bands(2)).and.bndint(ist2,evalsvnr(ist2,jk),&
      chi0_exclude_bands(1),chi0_exclude_bands(2))) l1=.false.
  if (l1) then
    t1=occsvnr(ist1,ik)-occsvnr(ist2,jk)
    if (abs(t1).gt.1d-6) then
      t2=sign(scissor,t1)
      wt(i)=t1/(evalsvnr(ist1,ik)-evalsvnr(ist2,jk)-t2+w)
    endif
  endif
enddo !i
call papi_timer_start(pt_megqblh2)
do ig=1,ngq
  megqblh2(1:nmegqblh(ikloc),ig)=dconjg(megqblh(1:nmegqblh(ikloc),ig,ikloc))*wt(1:nmegqblh(ikloc))
enddo
call papi_timer_stop(pt_megqblh2)
call papi_timer_start(pt_chi0_zgemm)
call zgemm('T','N',ngq,ngq,nmegqblh(ikloc),zone,megqblh(1,1,ikloc),nstsv*nstsv,&
  &megqblh2(1,1),nstsv*nstsv,zone,chi0w(1,1),ngq)
call papi_timer_stop(pt_chi0_zgemm)
deallocate(wt,megqblh2)
return
end subroutine

subroutine read_exxvxc(nst,exxnk,exxvc,vxcnk,vclnk,exind,ierr1)
use modmain
use mod_addons_q
!
implicit none
!
integer,intent(in) :: nst
real(8),intent(out) :: exxnk(nst,nspinor,nkptnrloc)
real(8),intent(out) :: exxvc(nst,nspinor,nkptnrloc)
real(8),intent(out) :: vxcnk(nst,nspinor,nkptnrloc)
real(8),intent(out) :: vclnk(nst,nspinor,nkptnrloc)
integer,intent(out) :: exind
integer,intent(out) :: ierr1
!
integer :: ikloc,isp1,ik,i,tmp,ios
logical :: exst
character*100 :: fname,fspn,fname_tot
!
!
ierr1=0
exind=0

do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 write(fname,'("Exx_Vxc_k",I4.4)')ik

 do isp1=1,nspinor
  if (nspinor.eq.1) then
   fspn=''
  else
   write(fspn,'("_spn",I1.1)') isp1
  endif
  fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
  inquire(file=fname_tot, exist=exst)
  if (exst) then
   open(164,file=fname_tot,action='read',form="FORMATTED",status="old")
   read(164,*,iostat=ios) exind
   do i=1,nst
    read(164,*,iostat=ios) tmp,exxnk(i,isp1,ikloc),exxvc(i,isp1,ikloc),vxcnk(i,isp1,ikloc),&
                           &vclnk(i,isp1,ikloc)
   enddo !i
   exxnk(:,isp1,ikloc)=exxnk(:,isp1,ikloc)/ha2ev  ! convert to Ha
   exxvc(:,isp1,ikloc)=exxvc(:,isp1,ikloc)/ha2ev  ! convert to Ha
   vxcnk(:,isp1,ikloc)=vxcnk(:,isp1,ikloc)/ha2ev
   vclnk(:,isp1,ikloc)=vclnk(:,isp1,ikloc)/ha2ev
   close(164)
  else
   write(*,*) "Exx files are not found!",ik
   ierr1=1
   exind=0
  endif
 enddo !isp1
enddo !ikloc

write(*,*) mpi_grid_dim_pos(dim_k),mpi_grid_dim_pos(dim_q),exind

call mpi_grid_reduce(exind,1,dims=(/dim_k,dim_q/),all=.true.,op=op_min)
call mpi_grid_reduce(ierr1,1,dims=(/dim_k,dim_q/),all=.true.,op=op_max)

if (ierr1.ne.0) then
 ierr1=1
 write(*,*) "ierr .ne. 0:",ierr1
 exind=0
 exxnk=0.d0
 exxvc=0.d0
 vxcnk=0.d0
 vclnk=0.d0
endif

return
end subroutine

subroutine write_exxvxc(iter,nst,bndrg,exxnk,exxvc,vxcnk,vclnk)
use modmain
use mod_addons_q
!
implicit none
!
integer,intent(in) :: iter
integer,intent(in) :: nst
integer,intent(in) :: bndrg(2,nspinor)
real(8),intent(in) :: exxnk(nst,nspinor,nkptnrloc)
real(8),intent(in) :: exxvc(nst,nspinor,nkptnrloc)
real(8),intent(in) :: vxcnk(nst,nspinor,nkptnrloc)
real(8),intent(in) :: vclnk(nst,nspinor,nkptnrloc)
!
integer :: ikloc,isp1,ik,i,ibnd,fbnd
character*100 :: fname,fspn,fname_tot
!
!outputing Exx_Vxc files
do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 do isp1=1,nspinor
  ibnd=bndrg(1,isp1)
  fbnd=bndrg(2,isp1)
  if (nspinor.eq.1) then
   fspn=''
  else
   write(fspn,'("_spn",I1.1)') isp1
  endif

  write(fname,'("Exx_Vxc_k",I4.4)')ik
  fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
  open(164,file=fname_tot,form="FORMATTED",status="REPLACE")
  write(164,'(I4)') iter
  do i=ibnd,fbnd
    write(164,'(I4,1X,4(f22.16,1X))') i,exxnk(i-ibnd+1,isp1,ikloc)*ha2ev,&
                                        exxvc(i-ibnd+1,isp1,ikloc)*ha2ev,&
                                        vxcnk(i-ibnd+1,isp1,ikloc)*ha2ev,&
                                        vclnk(i-ibnd+1,isp1,ikloc)*ha2ev
  enddo  !i
  close(164)
 enddo !isp1
enddo !ikloc

!delete temp files for exxnk
if (mpi_grid_root()) call system('rm -f Temp_files/tmp_Sx_ip*')

return
end subroutine

subroutine read_gwse(nw_se,nst,corr_se_aux,prev_Eqp,ierr2,step0,q0)
use modmain
use mod_addons_q
!
implicit none
!
integer,intent(in) :: nw_se
integer,intent(in) :: nst
complex(8),intent(inout) :: corr_se_aux(nw_se,nst,nspinor,nkptnr,nkpt+nvq0-1)
real(8),intent(inout) :: prev_Eqp(nst,nspinor,nkptnr)
integer,intent(inout) :: ierr2
integer,intent(inout) :: step0
integer,intent(inout) :: q0
!
integer :: ikloc,isp1,ik,i,tmp,ios,iz
logical :: exst
character*100 :: fname,fspn,fname_tot
real(8) :: prev_Eqp_(nst,nspinor,nkptnr)
!
!
ierr2=0
write(fname,'("Temp_files/tmp_Sc_ik",I4.4,"_iq",I4.4)') &
      &mpi_grid_dim_pos(dim_k),mpi_grid_dim_pos(dim_q)
inquire(file=trim(adjustl(fname)),exist=exst)

if (exst) then
 open(165,file=trim(adjustl(fname)),action='read',form='unformatted',&
    & status='old')
 read(165) step0,q0
 read(165) corr_se_aux
 read(165) prev_Eqp_
 close(165)
else
 write(*,'("tmp_Sc file not found!",I4)'),mpi_grid_dim_pos(dim_k)
 ierr2=1
endif

call mpi_grid_reduce(ierr2,1,dims=(/dim_k,dim_q/),all=.true.,op=op_max)
if (ierr2.ne.0) then
 q0=0
 step0=0
 corr_se_aux=zzero
 if (mpi_grid_root()) then
  write(151,'("Exx_Vxc arrays regenerated!")')
  write(151,'("But tmp_Sc files corrupted, need to recalculate!")')
  call timestamp(151)
  call flushifc(151)
 endif
else ! ierr2=0
 prev_Eqp(:,:,:)=prev_Eqp_(:,:,:)
endif

return
end subroutine

subroutine read_gwsetau(fnum,nst,sig_c_aux,sig_x_aux,gf_tau,sig_c,exxnk,&
                       ierr2,iter0,q0,iter1,maxdel_prev)
use modmain
use mod_addons_q
!
implicit none
!
integer,intent(in) :: fnum
integer,intent(in) :: nst
real(8),intent(inout) :: sig_c_aux(nst,nspinor,lr_nw,nkptnr,nkpt+nvq0-1)
real(8),intent(inout) :: sig_x_aux(nst,nspinor,nkptnr,nkpt+nvq0-1)
real(8),intent(inout) :: gf_tau(nstsv,lr_nw,nkptnrloc)
real(8),intent(inout) :: sig_c(nst,nspinor,lr_nw,nkptnrloc)
real(8),intent(inout) :: exxnk(nst,nspinor,nkptnrloc)
integer,intent(inout) :: ierr2
integer,intent(inout) :: iter0
integer,intent(inout) :: q0
integer,intent(inout) :: iter1
real(8),intent(inout) :: maxdel_prev
!
integer :: ikloc,isp1,ik,i,tmp,ios,iz
logical :: exst
character*100 :: fname,fspn,fname_tot
real(8) :: sig_c_aux_(nst,nspinor,lr_nw,nkptnr,nkpt+nvq0-1)
real(8) :: sig_x_aux_(nst,nspinor,nkptnr,nkpt+nvq0-1)
real(8) :: sig_c_(nst,nspinor,lr_nw,nkptnrloc)
real(8) :: exxnk_(nst,nspinor,nkptnrloc)
real(8) :: gf_tau_(nstsv,lr_nw,nkptnrloc)
!
!
ierr2=0
write(fname,'("Temp_files/tmp_Sc_ik",I4.4,"_iq",I4.4)') &
      & mpi_grid_dim_pos(dim_k),mpi_grid_dim_pos(dim_q)
inquire(file=trim(adjustl(fname)),exist=exst)

if (exst) then
  open(165,file=trim(adjustl(fname)),action='read',form='unformatted',&
       & status='old')
  read(165) iter0,q0,iter1
  read(165) sig_c_aux_
  read(165) gf_tau_
  read(165) sig_c_
  read(165) dyson_mu,maxdel_prev
  if (exxtype.eq.1) read(165) sig_x_aux_,exxnk_
  close(165)

  write(*,*) "iter1:",iter1,mpi_grid_dim_pos(dim_k),mpi_grid_dim_pos(dim_q)
else
  write(*,'("tmp_Sc file not found!",I4)') mpi_grid_dim_pos(dim_k)
  ierr2=1
endif

call mpi_grid_reduce(ierr2,1,dims=(/dim_k,dim_q/),all=.true.,&
                     &op=op_max)
if (ierr2.ne.0) then
  q0=0
  iter0=0
  iter1=0
  dyson_mu=0.d0
  maxdel_prev=0.d0
  if (mpi_grid_root()) then
    write(fnum,'("Exx_Vxc arrays regenerated!")')
    write(fnum,'("But tmp_Sc files corrupted, need to recalculate!")')
    call timestamp(fnum)
    call flushifc(fnum)
  endif
else ! ierr2=0
  sig_c_aux(:,:,:,:,:)=sig_c_aux_(:,:,:,:,:)
  gf_tau(:,:,:)=gf_tau_(:,:,:)
  sig_c(:,:,:,:)=sig_c_(:,:,:,:)
  if (exxtype.eq.1) then
    sig_x_aux(:,:,:,:)=sig_x_aux_(:,:,:,:)
    exxnk(:,:,:)=exxnk_(:,:,:)
  endif
  ! make sure all the processors receive identical values
  call mpi_grid_bcast(iter0,dims=(/dim_k/))
  call mpi_grid_bcast(q0,dims=(/dim_k/))
  call mpi_grid_bcast(iter1,dims=(/dim_k/))
  call mpi_grid_bcast(dyson_mu,dims=(/dim_k/))
  call mpi_grid_bcast(maxdel_prev,dims=(/dim_k/))
endif

return
end subroutine

subroutine ft_tw3(indx,bf,n1,n3,nk,n4,s0,s1)
use modmain
use mod_addons_q
!
implicit none
!
integer, intent(in) :: indx
integer, intent(in) :: bf
integer, intent(in) :: n1
integer, intent(in) :: n3
integer, intent(in) :: nk
integer, intent(in) :: n4
complex(8),intent(in) :: s0(n1,n3,nk)
complex(8),intent(out) :: s1(n1,n4,nk)
!
integer :: i1,i2,i3,ntmp,n,ikloc,nw,iw,w1
complex(8) :: ff
real(8) :: wn,dw0,z1(n3),yb,ya
real(8),allocatable :: wd(:)
real(8),external :: spline3_eval
complex(8),allocatable :: a1(:),a2(:)
!
! Fourier transformation between \tau domain and Matsubara frequency domain
!indx = +1    \tau -> iwn
!       -1     iwn -> \tau
! bf  =  1    bosonic
!     = -1    fermion
s1=zzero
n=0
!
dw0=dreal(lr_w(2)-lr_w(1))
nw=int(dreal(lr_w(lr_nw))/dw0+0.1d0)+1
allocate(wd(nw))
allocate(a1(nw-1))
allocate(a2(nw-1))

! define a dense uniform grid
do iw=1,nw
 wd(iw)=dw0*(iw-1)
enddo

if (indx.eq.1) then
 do ikloc=1,nk
  do i1=1,n1
   !cubic spline
   call spline3_coef(n3-1,dreal(lr_w(:)),dreal(s0(i1,:,ikloc)),z1(:))

   do iw=2,nw
     ya=spline3_eval(n3-1,dreal(lr_w(:)),dreal(s0(i1,:,ikloc)),&
        & z1(:),wd(iw-1))
     yb=spline3_eval(n3-1,dreal(lr_w(:)),dreal(s0(i1,:,ikloc)),&
        & z1(:),wd(iw))
     a1(iw-1)=(yb-ya)/dw0
     a2(iw-1)=yb-a1(iw-1)*wd(iw)
   enddo

   n=0
   do i2=-int((n4-1)/2),int(n4/2) !n4=niw
     n=n+1
     if (bf.eq.-1) then
       wn=(2.d0*i2+1)*pi/bhbar
     elseif (bf.eq.1) then
       wn=2.d0*i2*pi/bhbar
     endif
   
     if (abs(wn).gt.1.d-10) then
       do i3=2,nw
        ff=exp(zi*wn*wd(i3))*wd(i3)-exp(zi*wn*wd(i3-1))*wd(i3-1)
        ff=ff/zi/wn
        s1(i1,n,ikloc)=s1(i1,n,ikloc)+ff*a1(i3-1)

        ff=exp(zi*wn*wd(i3))-exp(zi*wn*wd(i3-1))
        ff=ff/zi/wn
        s1(i1,n,ikloc)=s1(i1,n,ikloc)+ff*(a2(i3-1)-a1(i3-1)/zi/wn)
       enddo
     else
       do i3=2,nw
        ff=0.5d0*(wd(i3)**2-wd(i3-1)**2)
        s1(i1,n,ikloc)=s1(i1,n,ikloc)+ff*a1(i3-1)
        ff=wd(i3)-wd(i3-1)
        s1(i1,n,ikloc)=s1(i1,n,ikloc)+ff*a2(i3-1)
       enddo
     endif

   enddo !i2
  enddo !i1
 enddo !ikloc

elseif (indx.eq.-1) then
 do ikloc=1,nk
  do i2=1,n4 !tau
   n=0
   do i1=-int((n3-1)/2),int(n3/2) !iwn
    if (bf.eq.-1) then
     wn=(2.d0*i1+1)*pi/bhbar
    elseif (bf.eq.1) then
     wn=2.d0*i1*pi/bhbar
    endif

    ff=exp(-zi*wn*lr_w(i2))
    n=n+1
    s1(:,i2,ikloc)=s1(:,i2,ikloc)+ff*s0(:,n,ikloc)
   enddo
  enddo
 enddo  !ikloc

 s1(:,:,:)=s1(:,:,:)/bhbar
endif

deallocate(wd,a1,a2)
return
end subroutine

subroutine ft_tw4(indx,bf,n1,n2,n3,nk,n4,s0,s1)
use modmain
use mod_addons_q
!
implicit none
!
integer, intent(in) :: indx
integer, intent(in) :: bf
integer, intent(in) :: n1
integer, intent(in) :: n2
integer, intent(in) :: n3
integer, intent(in) :: nk
integer, intent(in) :: n4
complex(8),intent(in) :: s0(n1,n2,n3,nk)
complex(8),intent(out) :: s1(n1,n2,n4,nk)
!
integer :: i1,i2,i3,i4,ntmp,n,ikloc,nw,iw,w1
complex(8) :: ff
real(8) :: wn,dw0,z1(n3),yb,ya
real(8),allocatable :: wd(:)
real(8),external :: spline3_eval
complex(8),allocatable :: a1(:),a2(:)
!
! Fourier transformation between \tau domain and Matsubara frequency domain
!indx = +1    \tau -> iwn
!       -1     iwn -> \tau
! bf  =  1    bosonic
!     = -1    fermion
s1=zzero
n=0
!
dw0=dreal(lr_w(2)-lr_w(1))
nw=int(dreal(lr_w(lr_nw))/dw0+0.1d0)+1
allocate(wd(nw))
allocate(a1(nw-1))
allocate(a2(nw-1))

! define a dense uniform grid
do iw=1,nw
 wd(iw)=dw0*(iw-1)
enddo

if (indx.eq.1) then
 do ikloc=1,nk
  do i2=1,n2
   do i1=1,n1
    !cubic spline
    call spline3_coef(n3-1,dreal(lr_w(:)),dreal(s0(i1,i2,:,ikloc)),z1(:))
    
    do iw=2,nw
     ya=spline3_eval(n3-1,dreal(lr_w(:)),dreal(s0(i1,i2,:,ikloc)),&
        & z1(:),wd(iw-1))
     yb=spline3_eval(n3-1,dreal(lr_w(:)),dreal(s0(i1,i2,:,ikloc)),&
        & z1(:),wd(iw))
     a1(iw-1)=(yb-ya)/dw0
     a2(iw-1)=yb-a1(iw-1)*wd(iw)
    enddo

    n=0
    do i4=-int((n4-1)/2),int(n4/2) !n4=niw
     n=n+1
     if (bf.eq.-1) then
       wn=(2.d0*i4+1)*pi/bhbar
     elseif (bf.eq.1) then
       wn=2.d0*i4*pi/bhbar
     endif

     if (abs(wn).gt.1.d-10) then
       do i3=2,nw
        ff=exp(zi*wn*wd(i3))*wd(i3)-exp(zi*wn*wd(i3-1))*wd(i3-1)
        ff=ff/zi/wn    
        s1(i1,i2,n,ikloc)=s1(i1,i2,n,ikloc)+ff*a1(i3-1)

        ff=exp(zi*wn*wd(i3))-exp(zi*wn*wd(i3-1))
        ff=ff/zi/wn
        s1(i1,i2,n,ikloc)=s1(i1,i2,n,ikloc)+ff*(a2(i3-1)-a1(i3-1)/zi/wn)
       enddo
     else
       do i3=2,nw
        ff=0.5d0*(wd(i3)**2-wd(i3-1)**2)
        s1(i1,i2,n,ikloc)=s1(i1,i2,n,ikloc)+ff*a1(i3-1)
        ff=wd(i3)-wd(i3-1)
        s1(i1,i2,n,ikloc)=s1(i1,i2,n,ikloc)+ff*a2(i3-1)
       enddo
     endif  

    enddo !i4
   enddo !i1
  enddo !i2
 enddo !ikloc

elseif (indx.eq.-1) then
 do ikloc=1,nk
  do i2=1,n4 !tau
   n=0
   do i1=-int((n3-1)/2),int(n3/2) !iwn
    if (bf.eq.-1) then
     wn=(2.d0*i1+1)*pi/bhbar
    elseif (bf.eq.1) then
     wn=2.d0*i1*pi/bhbar
    endif

    ff=exp(-zi*wn*lr_w(i2))
    n=n+1
    s1(:,:,i2,ikloc)=s1(:,:,i2,ikloc)+ff*s0(:,:,n,ikloc)
   enddo
  enddo
 enddo  !ikloc

 s1(:,:,:,:)=s1(:,:,:,:)/bhbar
endif

deallocate(wd,a1,a2)
return
end subroutine

subroutine update_occ(gf)
use modmain
use mod_nrkp
!
implicit none
!
real(8),intent(in) :: gf(nstsv,lr_nw,nkptnrloc)
!
!local variables
integer :: ist, ik, ikloc

! reset occupations
occsvnr(:,:)=0.d0

do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 do ist=1,nstsv
  occsvnr(ist,ik)=-gf(ist,lr_nw,ikloc)
 enddo
enddo

call mpi_grid_reduce(occsvnr(1,1),nkptnr*nstsv,dims=(/dim_k/),all=.true.)

occsvnr(:,:)=occsvnr(:,:)*occmax

return
end subroutine

subroutine get_delf(bndrg,gf0,gf,delf)
use modmain
use mod_nrkp
!
implicit none
!
integer,intent(in) :: bndrg(2,nspinor)
real(8),intent(in) :: gf0(nstsv,lr_nw,nkptnrloc)
real(8),intent(in) :: gf(nstsv,lr_nw,nkptnrloc)
real(8),intent(inout) :: delf(nstsv,nkptnr)
!
!local variables
integer :: ist, ik, ikloc, isp1

! calculate the occupation change in GW calculations
! Iek-Heng Chu, 2014
!
delf(:,:)=0.d0

do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 do isp1=1,nspinor  
  do ist=bndrg(1,isp1),bndrg(2,isp1)
    delf(ist,ik)=-gf(ist,lr_nw,ikloc)+gf0(ist,lr_nw,ikloc)
  enddo
 enddo
enddo

call mpi_grid_reduce(delf(1,1),nkptnr*nstsv,dims=(/dim_k/),all=.true.)

delf(:,:)=delf(:,:)*occmax

if (mpi_grid_root()) then
 do ik=1,nkptnr
  do ist=1,nstsv
    write(*,*) "delf,ist,ik:",delf(ist,ik),ist,ik
  enddo
 enddo
endif

return
end subroutine

end module
