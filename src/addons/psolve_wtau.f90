subroutine psolve_wtau(iq,nwloc,gf0,svq0)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
!
integer,intent(in) :: iq
integer,intent(in) :: nwloc
real(8),intent(in) :: gf0(nstsv,lr_nw,nkptnrloc)
complex(8),intent(out) :: svq0(ngq(iq),ngq(iq),nwloc)
!
! local variable
integer :: nxprocs,nyprocs,nxprocs_0,nyprocs_0
integer :: ictxt, my_row, my_col
integer :: mb, nb, m, n, mpdim, npdim, nbdim, mblock,nblock
integer :: mp, np, mc, mr, igloc, jgloc,indx, iloc,jloc
integer :: ig1,ig2,ierr,info,tmp1,iw,iw1,iw2,i,j
integer :: iwloc
integer :: iwmap(lr_nw,lr_nw)
! scalapack dimensions
integer,parameter :: IDESCLEN_=9
! descriptors a and b
integer :: idesca(IDESCLEN_), idescb(IDESCLEN_)
complex(8),allocatable :: krnl(:,:)
complex(8),allocatable :: aloc(:),bloc(:)
complex(8),allocatable :: pw_aux(:,:,:)
integer,allocatable :: ipiv(:),tmp2(:)
real(8) :: dw(lr_nw),dwr(lr_nw,lr_nw)
complex(8) :: lint,aux1,aux2
integer,external :: NUMROC
character*100 :: fname
!
iwmap=0
svq0=zzero

! amat: m x m ;   bmat: m x n
m=ngq(iq)*lr_nw
n=ngq(iq)
! desired block size of block-cyclic
mb=64
nb=64
! process row and column
if (nproc.gt.1) then
 do i=1,int(sqrt(dble(nproc)))
  if (mod(nproc,i).eq.0) tmp1=i
 enddo
 nxprocs=tmp1
 nyprocs=nproc/tmp1
else
 nxprocs=1
 nyprocs=1
endif
! 
ictxt=0
call sl_init(ictxt,nxprocs,nyprocs)
call blacs_gridinfo(ictxt,nxprocs_0,nyprocs_0,my_row,my_col)

! create the scalapack array descriptors for 2d-block-cyclic
! (1) for m x m amat
tmp1=int(m/nxprocs)
if (mod(m,nxprocs).ne.0) tmp1=tmp1+1
mblock=max(1,min(tmp1,mb))

tmp1=int(m/nyprocs)
if (mod(m,nyprocs).ne.0) tmp1=tmp1+1
nblock=max(1,min(tmp1,mb))

mpdim=NUMROC(m,mblock,my_row,0,nxprocs)
npdim=NUMROC(m,nblock,my_col,0,nyprocs)

call DESCINIT(idesca,m,m,mblock,nblock,0,0,ictxt,mpdim,info)

write(*,*) "mpdim,npdim,info:",mpdim,npdim,info,my_row,my_col

! (2) for nrhs of bmat
tmp1=int(n/nyprocs)
if (mod(n,nyprocs).ne.0) tmp1=tmp1+1
nblock=max(1,min(tmp1,nb))
nbdim=NUMROC(n,nblock,my_col,0,nyprocs)

call DESCINIT(idescb,m,n,mblock,nblock,0,0,ictxt,mpdim,info)

! now allocate local matrices
allocate(aloc(mpdim*npdim))
allocate(bloc(mpdim*nbdim))
allocate(ipiv(mpdim+mblock)) ! swap-info from pivoting
allocate(tmp2(ngq(iq)*lr_nw))
allocate(pw_aux(ngq(iq),ngq(iq),lr_nw))
allocate(krnl(ngq(iq),2))
krnl=zzero
pw_aux=zzero
aloc=zzero
bloc=zzero
ipiv=0
tmp2=0

! v(q)=4\pi/|q|^2
do ig1=1,ngq(iq)
  krnl(ig1,1)=dcmplx(vhgq(ig1,iq))
! different weight at Gamma
  krnl(ig1,2)=dcmplx(wtvhgq(ig1,iq))
enddo

! calculate the polarization function
call cal_pol(iq,nwloc,gf0,pw_aux)
!
! P -> vP
!
do iw=1,lr_nw
 do ig2=1,ngq(iq)
  do ig1=1,ngq(iq)
    pw_aux(ig1,ig2,iw)=pw_aux(ig1,ig2,iw)*krnl(ig1,1)
  enddo
 enddo
enddo

if (mpi_grid_root()) then
 write(151,'(" ")')
 write(151,'("Solving the screened Coulomb potential...")')
 write(151,'("m,n:",2I6)') m,n
 call flushifc(151)
endif

! now assign aloc!  igloc=(ig,iw)
mblock=idesca(5)
nblock=idesca(6)

! first calculate dw^{r}
dw(1)=dreal(lr_w(2)-lr_w(1))
dw(lr_nw)=dreal(lr_w(lr_nw)-lr_w(lr_nw-1))

do iw1=2,lr_nw-1
 dw(iw1)=dreal(lr_w(iw1+1)-lr_w(iw1-1))
enddo
dw(:)=dw(:)*0.5d0

! calculate tmp2
do igloc=1,m
 tmp2(igloc)=(igloc-1)/(mblock*nxprocs)
enddo

! mapping between iw and iw1-iw2, linear interpolation is applied
call findiw_new(lr_nw,iwmap,dwr)

do jgloc=1,m
 mc=mod((jgloc-1)/nblock,nyprocs)
 if (mc.ne.my_col) cycle

 ig2=mod(jgloc-1,ngq(iq))+1
 iw2=(jgloc-1)/ngq(iq)+1
 tmp1=(jgloc-1)/(nblock*nyprocs)
 jloc=tmp1*nblock+mod((jgloc-1),nblock)+1
 do igloc=1,m
  mr=mod((igloc-1)/mblock,nxprocs)
  if (mr.ne.my_row) cycle

  ig1=mod(igloc-1,ngq(iq))+1
  iw1=(igloc-1)/ngq(iq)+1
  iloc=tmp2(igloc)*mblock+mod(igloc-1,mblock)+1
  indx=iloc+(jloc-1)*mpdim 
  iw=iwmap(iw1,iw2)
  aux1=pw_aux(ig1,ig2,iw)
  aux2=pw_aux(ig1,ig2,iw+1)
  lint=aux1*(zone-dwr(iw1,iw2))+aux2*dwr(iw1,iw2)
  aloc(indx)=-lint*dw(iw2)
  if (ig1.eq.ig2.and.iw1.eq.iw2) aloc(indx)=aloc(indx)+zone
 enddo
enddo

!bloc
mblock=idescb(5)
nblock=idescb(6)

do jgloc=1,n
 mc=mod((jgloc-1)/nblock,nyprocs)
 if (mc.ne.my_col) cycle

 tmp1=(jgloc-1)/(nblock*nyprocs)
 jloc=tmp1*nblock+mod((jgloc-1),nblock)+1
 do igloc=1,m
  mr=mod((igloc-1)/mblock,nxprocs)
  if (mr.ne.my_row) cycle

  ig1=mod(igloc-1,ngq(iq))+1
  iw1=(igloc-1)/ngq(iq)+1
  tmp1=(igloc-1)/(mblock*nxprocs)
  iloc=tmp1*mblock+mod(igloc-1,mblock)+1
  indx=iloc+(jloc-1)*mpdim
  bloc(indx)=pw_aux(ig1,jgloc,iw1)*krnl(jgloc,2)
 enddo
enddo

!reset zero
pw_aux(:,:,:)=zzero

! wait for all the processes
call mpi_grid_barrier()

! solve the linear equation using pdgesv
call pzgesv(m,n,aloc,1,1,idesca,ipiv,bloc,1,1,idescb,info)
if (info.ne.0) write(*,*) "info.ne.0:",info

! store solution in svq0
mblock=idescb(5)
nblock=idescb(6)

do jloc=1,nbdim
  tmp1=(jloc-1)/nblock
  jgloc=tmp1*nblock*nyprocs+my_col*nblock+mod(jloc-1,nblock)+1
  if (jgloc.gt.n) cycle
  do iloc=1,mpdim
   tmp1=(iloc-1)/mblock
   igloc=tmp1*mblock*nxprocs+my_row*mblock+mod(iloc-1,mblock)+1
   if (igloc.gt.m) cycle
   indx=iloc+(jloc-1)*mpdim
   ig1=mod(igloc-1,ngq(iq))+1
   iw=(igloc-1)/ngq(iq)+1
   pw_aux(ig1,jgloc,iw)=bloc(indx)
  enddo
enddo

write(*,*) "mpi_reduce for svq matrix:",mpi_grid_dim_pos(dim_k),&
           &mpi_grid_dim_pos(dim_q)

call mpi_grid_reduce(pw_aux(1,1,1),ngq(iq)*ngq(iq)*lr_nw,all=.true.)

do iwloc=1,nwloc
 iw=mpi_grid_map(lr_nw,dim_q,loc=iwloc)
 svq0(:,:,iwloc)=pw_aux(:,:,iw)
enddo

!---test
!if (mpi_grid_root()) then
! write(fname,'("Pv_iq",I4.4)') iq
! open(888,file=fname,status='replace')
! do iw1=1,lr_nw-1
!  write(888,'(3(f18.6,1X))') dreal(lr_w(iw1)),dreal(aux(2,2,iw1)),dimag(aux(2,2,iw1))
!  iw=iwmap(iw1,2)
!  lint=aux(2,2,iw)*(zone-dwr(iw1,2))+aux(2,2,iw+1)*dwr(iw1,2)
!  write(888,'(3(f18.6,1X))') dreal(lr_w(iw)+dwr(iw1,2)*(lr_w(iw+1)-lr_w(iw))),&
!                           &dreal(lint),dimag(lint)
! enddo
!endif
!---end test

deallocate(krnl,ipiv,aloc,bloc,pw_aux)
deallocate(tmp2)
return
end subroutine
