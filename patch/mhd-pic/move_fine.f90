subroutine move_fine(ilevel)
  use amr_commons
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  use pm_commons
  use mpi_mod
  implicit none
  integer::ilevel,xtondim
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,icpu,ind,iskip,ivar,i
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  if(trajectories(1)>0 .or. ntrajectories>0)then
   filename='trajectory.dat'
   call title(myid,nchar)
   fileloc=TRIM(filename)//TRIM(nchar)
   open(25+myid, file = fileloc, status = 'unknown', access = 'append')
  endif

#ifdef TSC
    xtondim=threetondim
#else
    xtondim=twotondim
#endif

  ! Set unew = uold in the active region
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=2,4
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar)=&
           &uold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
  end do
  ! Set unew reception cells to zero
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=2,4
           do i=1,reception(icpu,ilevel)%ngrid
              unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0D0
           end do
        end do
     end do
  end do

  ! Update particles position and velocity
  ig=0
  ip=0
  ! Loop over grids
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     npart1=numbp(igrid)  ! Number of particles in the grid
     if(npart1>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle  <---- Very important !!!
           next_part=nextp(ipart)
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ip=ip+1
           ind_part(ip)=ipart
           ind_grid_part(ip)=ig
           if(ip==nvector)then
              call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)
              ip=0
              ig=0
           end if
           ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0)call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)

  ! Update MPI boundary conditions for unew for dust mass and momentum densities
  do ivar=2,4 ! Gas momentum indices
     call make_virtual_reverse_dp(unew(1,ivar),ilevel)
  end do

  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=2,4
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar)=&
           &unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
  end do

  do ivar=2,4 ! Gas momentum indices
     call make_virtual_fine_dp   (uold(1,ivar),ilevel)
  end do

  close(25+myid)
!!!!!!!!!!!!!!!!!!!!!!! NEW !!!!!!!!!!!!!!!!!!!!!!!!!!!
if(simple_boundary)call make_boundary_hydro(ilevel)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

111 format('   Entering move_fine for level ',I2)

end subroutine move_fine
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move_fine_static(ilevel)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  integer::ilevel,xtondim
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,npart2
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
#ifdef TSC
    xtondim=threetondim
#else
    xtondim=twotondim
#endif

  ! Update particles position and velocity
  ig=0
  ip=0
  ! Loop over grids
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     npart1=numbp(igrid)  ! Number of particles in the grid
     npart2=0

     ! Count particles
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           if(star) then
              if ( (.not. static_DM .and. is_DM(typep(ipart))) .or. &
                   & (.not. static_stars .and. is_not_DM(typep(ipart)) )  ) then
                 ! FIXME: there should be a static_sink as well
                 ! FIXME: what about debris?
                 npart2=npart2+1
              endif
           else
              if(.not.static_DM) then
                 npart2=npart2+1
              endif
           endif
           ipart=next_part  ! Go to next particle
        end do
     endif

     ! Gather star particles
     if(npart2>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           ! Select particles
           if(star) then
              if ( (.not. static_DM .and. is_DM(typep(ipart))) .or. &
                   & (.not. static_stars .and. is_not_DM(typep(ipart)) )  ) then
                 ! FIXME: there should be a static_sink as well
                 ! FIXME: what about debris?
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
           else
              if(.not.static_dm) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
           endif
           if(ip==nvector)then
              call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)
              ip=0
              ig=0
           end if
           ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0)call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)

111 format('   Entering move_fine for level ',I2)

end subroutine move_fine_static
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move1(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,xtondim)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  use amr_parameters, ONLY: cr_c_fraction
  implicit none
  integer::ng,np,ilevel,xtondim
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !------------------------------------------------------------
  ! This routine computes the force on each particle B2
  ! inverse CIC and computes new positions for all particles.
  ! If particle sits entirely in fine level, then CIC is performed
  ! at level ilevel. Otherwise, it is performed at level ilevel-1.
  ! This routine is called B2 move_fine.
  !------------------------------------------------------------
  logical::error
  integer::i,j,ind,idim,nx_loc,isink,index_part,ivar_mhd_pic,iskip,icpu
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::ctm,tempc,maxciso

  ! Grid-based arrays
  integer ,dimension(1:nvector),save::father_cell
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  real(dp),dimension(1:nvector),save:: sizes,charges! ERM: grain sizes and pic charges
  ! Particle-based arrays
#ifndef TSC
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_xp,new_vp,dd,dg
  real(dp),dimension(1:nvector,1:ndim),save::vv
  real(dp),dimension(1:nvector,1:ndim),save::bb,uu
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::big_vv,big_ww
  real(dp),dimension(1:nvector),save:: nu_stop,mov,dgr,ddgr,ciso
  real(dp),dimension(1:nvector),save:: lorentzf,new_lorentzf,new_mp
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
#else
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_xp,new_vp
  real(dp),dimension(1:nvector,1:ndim),save::vv,cl,cr,cc,wl,wr,wc
  real(dp),dimension(1:nvector,1:ndim),save::bb,uu
  real(dp),dimension(1:nvector,1:threetondim,1:ndim),save::big_vv,big_ww
  real(dp),dimension(1:nvector),save:: nu_stop,mov,dgr,ddgr,ciso ! ERM: fluid variables and stopping times
  real(dp),dimension(1:nvector),save:: lorentzf,new_lorentzf,new_mp
  integer ,dimension(1:nvector,1:ndim),save::igl,igr,igc,icl,icr,icc
  real(dp),dimension(1:nvector,1:threetondim),save::vol
  integer ,dimension(1:nvector,1:threetondim),save::igrid,icell,indp,kg
#endif
  real(dp),dimension(1:3)::skip_loc
  real(dp)::den_dust,den_gas,mom_dust,mom_gas,velocity_com, crsol
  ! Family
  logical,dimension(1:nvector),save :: classical_tracer, dust, cosr, mhd_pic
  ! ERM: w is the cell dust-gas drift, B the mag field.
  ctm = du_charge_to_mass
  crsol=cr_c_fraction*2.9979246d+10*units_time/units_length ! Reduced speed of light.
  !ts = t_stop!  ERM: Not used if constant_t_stop==.false.
  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather neighboring father cells (should be present anytime !)
  do i=1,ng
     father_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,&
       & ng,ilevel)

  ! Rescale particle position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

#ifndef TSC
  ! Check for illegal moves. Is this different for CIC vs TSC?
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in move'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=int(dd(j,idim))
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do ind=1,twotondim
     do j=1,np
        ok(j)=ok(j).and.igrid(j,ind)>0
     end do
  end do

  ! If not, rescale position at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           x(j,idim)=x(j,idim)/2.0D0
        end if
     end do
  end do
  ! If not, redo CIC at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=int(dd(j,idim))
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end if
     end do
  end do

 ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        else
           icg(j,idim)=ig(j,idim)
           icd(j,idim)=id(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,5)=1+icg(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,6)=1+icd(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,7)=1+icg(j,1)+3*icd(j,2)+9*icd(j,3)
        icell(j,8)=1+icd(j,1)+3*icd(j,2)+9*icd(j,3)
     end if
  end do
#endif

  ! Compute parent cell adresses
  do ind=1,twotondim
     do j=1,np
        if(ok(j))then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        else
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
        end if
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif

#else
!#include "tsc_fine.F90"
if (ndim .ne. 3)then
   write(*,*)'TSC not supported for ndim neq 3'
   call clean_stop
end if

xtondim=threetondim

! Check for illegal moves
error=.false.
do idim=1,ndim
   do j=1,np
      if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0)error=.true.
   end do
end do
if(error)then
   write(*,*)'problem in tsc_fine'
   do idim=1,ndim
      do j=1,np
         if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0)then
            write(*,*)x(j,1:ndim)
         endif
      end do
   end do
   stop
end if

! TSC at level ilevel; a particle contributes
!     to three cells in each dimension
! cl: position of leftmost cell centre
! cc: position of central cell centre
! cr: position of rightmost cell centre
! wl: weighting function for leftmost cell
! wc: weighting function for central cell
! wr: weighting function for rightmost cell
do idim=1,ndim
   do j=1,np
      cl(j,idim)=dble(int(x(j,idim)))-0.5D0
      cc(j,idim)=dble(int(x(j,idim)))+0.5D0
      cr(j,idim)=dble(int(x(j,idim)))+1.5D0
      wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
      wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
      wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
   end do
end do

! Compute parent grids
do idim=1,ndim
   do j=1,np
      igl(j,idim)=(int(cl(j,idim)))/2
      igc(j,idim)=(int(cc(j,idim)))/2
      igr(j,idim)=(int(cr(j,idim)))/2
   end do
end do
! #if NDIM==3
do j=1,np
   kg(j,1 )=1+igl(j,1)+3*igl(j,2)+9*igl(j,3)
   kg(j,2 )=1+igc(j,1)+3*igl(j,2)+9*igl(j,3)
   kg(j,3 )=1+igr(j,1)+3*igl(j,2)+9*igl(j,3)
   kg(j,4 )=1+igl(j,1)+3*igc(j,2)+9*igl(j,3)
   kg(j,5 )=1+igc(j,1)+3*igc(j,2)+9*igl(j,3)
   kg(j,6 )=1+igr(j,1)+3*igc(j,2)+9*igl(j,3)
   kg(j,7 )=1+igl(j,1)+3*igr(j,2)+9*igl(j,3)
   kg(j,8 )=1+igc(j,1)+3*igr(j,2)+9*igl(j,3)
   kg(j,9 )=1+igr(j,1)+3*igr(j,2)+9*igl(j,3)
   kg(j,10)=1+igl(j,1)+3*igl(j,2)+9*igc(j,3)
   kg(j,11)=1+igc(j,1)+3*igl(j,2)+9*igc(j,3)
   kg(j,12)=1+igr(j,1)+3*igl(j,2)+9*igc(j,3)
   kg(j,13)=1+igl(j,1)+3*igc(j,2)+9*igc(j,3)
   kg(j,14)=1+igc(j,1)+3*igc(j,2)+9*igc(j,3)
   kg(j,15)=1+igr(j,1)+3*igc(j,2)+9*igc(j,3)
   kg(j,16)=1+igl(j,1)+3*igr(j,2)+9*igc(j,3)
   kg(j,17)=1+igc(j,1)+3*igr(j,2)+9*igc(j,3)
   kg(j,18)=1+igr(j,1)+3*igr(j,2)+9*igc(j,3)
   kg(j,19)=1+igl(j,1)+3*igl(j,2)+9*igr(j,3)
   kg(j,20)=1+igc(j,1)+3*igl(j,2)+9*igr(j,3)
   kg(j,21)=1+igr(j,1)+3*igl(j,2)+9*igr(j,3)
   kg(j,22)=1+igl(j,1)+3*igc(j,2)+9*igr(j,3)
   kg(j,23)=1+igc(j,1)+3*igc(j,2)+9*igr(j,3)
   kg(j,24)=1+igr(j,1)+3*igc(j,2)+9*igr(j,3)
   kg(j,25)=1+igl(j,1)+3*igr(j,2)+9*igr(j,3)
   kg(j,26)=1+igc(j,1)+3*igr(j,2)+9*igr(j,3)
   kg(j,27)=1+igr(j,1)+3*igr(j,2)+9*igr(j,3)
end do

do ind=1,threetondim
   do j=1,np
      igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
   end do
end do

! Check if particles are entirely in level ilevel
ok(1:np)=.true.
do ind=1,threetondim
   do j=1,np
      ok(j)=ok(j).and.igrid(j,ind)>0
   end do
end do

! If not, rescale position at level ilevel-1
do idim=1,ndim
   do j=1,np
      if(.not.ok(j))then
         x(j,idim)=x(j,idim)/2.0D0
      end if
   end do
end do
! If not, redo TSC at level ilevel-1
do idim=1,ndim
   do j=1,np
      if(.not.ok(j))then
        cl(j,idim)=dble(int(x(j,idim)))-0.5D0
        cc(j,idim)=dble(int(x(j,idim)))+0.5D0
        cr(j,idim)=dble(int(x(j,idim)))+1.5D0
        wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
        wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
        wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
      end if
   end do
end do

! Compute parent cell position
do idim=1,ndim
   do j=1,np
     if(ok(j))then
      icl(j,idim)=int(cl(j,idim))-2*igl(j,idim)
      icc(j,idim)=int(cc(j,idim))-2*igc(j,idim)
      icr(j,idim)=int(cr(j,idim))-2*igr(j,idim)
     else ! ERM: this else may or may not be correct? But I believe it is.
      icl(j,idim)=int(cl(j,idim))
      icc(j,idim)=int(cc(j,idim))
      icr(j,idim)=int(cr(j,idim))
     endif
   end do
end do

! #if NDIM==3
do j=1,np
  if(ok(j))then
   icell(j,1 )=1+icl(j,1)+2*icl(j,2)+4*icl(j,3)
   icell(j,2 )=1+icc(j,1)+2*icl(j,2)+4*icl(j,3)
   icell(j,3 )=1+icr(j,1)+2*icl(j,2)+4*icl(j,3)
   icell(j,4 )=1+icl(j,1)+2*icc(j,2)+4*icl(j,3)
   icell(j,5 )=1+icc(j,1)+2*icc(j,2)+4*icl(j,3)
   icell(j,6 )=1+icr(j,1)+2*icc(j,2)+4*icl(j,3)
   icell(j,7 )=1+icl(j,1)+2*icr(j,2)+4*icl(j,3)
   icell(j,8 )=1+icc(j,1)+2*icr(j,2)+4*icl(j,3)
   icell(j,9 )=1+icr(j,1)+2*icr(j,2)+4*icl(j,3)
   icell(j,10)=1+icl(j,1)+2*icl(j,2)+4*icc(j,3)
   icell(j,11)=1+icc(j,1)+2*icl(j,2)+4*icc(j,3)
   icell(j,12)=1+icr(j,1)+2*icl(j,2)+4*icc(j,3)
   icell(j,13)=1+icl(j,1)+2*icc(j,2)+4*icc(j,3)
   icell(j,14)=1+icc(j,1)+2*icc(j,2)+4*icc(j,3)
   icell(j,15)=1+icr(j,1)+2*icc(j,2)+4*icc(j,3)
   icell(j,16)=1+icl(j,1)+2*icr(j,2)+4*icc(j,3)
   icell(j,17)=1+icc(j,1)+2*icr(j,2)+4*icc(j,3)
   icell(j,18)=1+icr(j,1)+2*icr(j,2)+4*icc(j,3)
   icell(j,19)=1+icl(j,1)+2*icl(j,2)+4*icr(j,3)
   icell(j,20)=1+icc(j,1)+2*icl(j,2)+4*icr(j,3)
   icell(j,21)=1+icr(j,1)+2*icl(j,2)+4*icr(j,3)
   icell(j,22)=1+icl(j,1)+2*icc(j,2)+4*icr(j,3)
   icell(j,23)=1+icc(j,1)+2*icc(j,2)+4*icr(j,3)
   icell(j,24)=1+icr(j,1)+2*icc(j,2)+4*icr(j,3)
   icell(j,25)=1+icl(j,1)+2*icr(j,2)+4*icr(j,3)
   icell(j,26)=1+icc(j,1)+2*icr(j,2)+4*icr(j,3)
   icell(j,27)=1+icr(j,1)+2*icr(j,2)+4*icr(j,3)
 else
   icell(j,1 )=1+icl(j,1)+3*icl(j,2)+9*icl(j,3)
   icell(j,2 )=1+icc(j,1)+3*icl(j,2)+9*icl(j,3)
   icell(j,3 )=1+icr(j,1)+3*icl(j,2)+9*icl(j,3)
   icell(j,4 )=1+icl(j,1)+3*icc(j,2)+9*icl(j,3)
   icell(j,5 )=1+icc(j,1)+3*icc(j,2)+9*icl(j,3)
   icell(j,6 )=1+icr(j,1)+3*icc(j,2)+9*icl(j,3)
   icell(j,7 )=1+icl(j,1)+3*icr(j,2)+9*icl(j,3)
   icell(j,8 )=1+icc(j,1)+3*icr(j,2)+9*icl(j,3)
   icell(j,9 )=1+icr(j,1)+3*icr(j,2)+9*icl(j,3)
   icell(j,10)=1+icl(j,1)+3*icl(j,2)+9*icc(j,3)
   icell(j,11)=1+icc(j,1)+3*icl(j,2)+9*icc(j,3)
   icell(j,12)=1+icr(j,1)+3*icl(j,2)+9*icc(j,3)
   icell(j,13)=1+icl(j,1)+3*icc(j,2)+9*icc(j,3)
   icell(j,14)=1+icc(j,1)+3*icc(j,2)+9*icc(j,3)
   icell(j,15)=1+icr(j,1)+3*icc(j,2)+9*icc(j,3)
   icell(j,16)=1+icl(j,1)+3*icr(j,2)+9*icc(j,3)
   icell(j,17)=1+icc(j,1)+3*icr(j,2)+9*icc(j,3)
   icell(j,18)=1+icr(j,1)+3*icr(j,2)+9*icc(j,3)
   icell(j,19)=1+icl(j,1)+3*icl(j,2)+9*icr(j,3)
   icell(j,20)=1+icc(j,1)+3*icl(j,2)+9*icr(j,3)
   icell(j,21)=1+icr(j,1)+3*icl(j,2)+9*icr(j,3)
   icell(j,22)=1+icl(j,1)+3*icc(j,2)+9*icr(j,3)
   icell(j,23)=1+icc(j,1)+3*icc(j,2)+9*icr(j,3)
   icell(j,24)=1+icr(j,1)+3*icc(j,2)+9*icr(j,3)
   icell(j,25)=1+icl(j,1)+3*icr(j,2)+9*icr(j,3)
   icell(j,26)=1+icc(j,1)+3*icr(j,2)+9*icr(j,3)
   icell(j,27)=1+icr(j,1)+3*icr(j,2)+9*icr(j,3)
 endif
end do

! Compute parent cell adress
do ind=1,threetondim
   do j=1,np
     if(ok(j))then
      indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
     else ! ERM: for AMR(?) there may be an issue with ind_grid_part(j) being used here.
       indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
     endif
   end do
end do

! Compute cloud volumes (NDIM==3)
do j=1,np
   vol(j,1 )=wl(j,1)*wl(j,2)*wl(j,3)
   vol(j,2 )=wc(j,1)*wl(j,2)*wl(j,3)
   vol(j,3 )=wr(j,1)*wl(j,2)*wl(j,3)
   vol(j,4 )=wl(j,1)*wc(j,2)*wl(j,3)
   vol(j,5 )=wc(j,1)*wc(j,2)*wl(j,3)
   vol(j,6 )=wr(j,1)*wc(j,2)*wl(j,3)
   vol(j,7 )=wl(j,1)*wr(j,2)*wl(j,3)
   vol(j,8 )=wc(j,1)*wr(j,2)*wl(j,3)
   vol(j,9 )=wr(j,1)*wr(j,2)*wl(j,3)
   vol(j,10)=wl(j,1)*wl(j,2)*wc(j,3)
   vol(j,11)=wc(j,1)*wl(j,2)*wc(j,3)
   vol(j,12)=wr(j,1)*wl(j,2)*wc(j,3)
   vol(j,13)=wl(j,1)*wc(j,2)*wc(j,3)
   vol(j,14)=wc(j,1)*wc(j,2)*wc(j,3)
   vol(j,15)=wr(j,1)*wc(j,2)*wc(j,3)
   vol(j,16)=wl(j,1)*wr(j,2)*wc(j,3)
   vol(j,17)=wc(j,1)*wr(j,2)*wc(j,3)
   vol(j,18)=wr(j,1)*wr(j,2)*wc(j,3)
   vol(j,19)=wl(j,1)*wl(j,2)*wr(j,3)
   vol(j,20)=wc(j,1)*wl(j,2)*wr(j,3)
   vol(j,21)=wr(j,1)*wl(j,2)*wr(j,3)
   vol(j,22)=wl(j,1)*wc(j,2)*wr(j,3)
   vol(j,23)=wc(j,1)*wc(j,2)*wr(j,3)
   vol(j,24)=wr(j,1)*wc(j,2)*wr(j,3)
   vol(j,25)=wl(j,1)*wr(j,2)*wr(j,3)
   vol(j,26)=wc(j,1)*wr(j,2)*wr(j,3)
   vol(j,27)=wr(j,1)*wr(j,2)*wr(j,3)
end do

#endif
  ! Gather center of mass 3-velocity
  ivar_mhd_pic=9
  if((pic_dust .or. pic_cr) .and. nvar<ivar_mhd_pic+ndim)then
     write(*,*)'You need to compile ramses with nvar=',ivar_mhd_pic+ndim
     stop
  endif

  ! ERM: need to put appropriate flags in place,
  !test to see if code still works the same, then try to add tracers,
  !then MC tracers, then CRs. Also enable alternate mass spectra for dust.

  ! Boolean flag for classical tracers (velocity advected ones)
  ! Note: this is enough as the subroutine is only called if `MC_tracer` is false (is this true??? Then MC_tracer won't work with this.)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do j = 1, np
      classical_tracer(j) = is_gas_tracer(typep(ind_part(j)))
      dust(j) = is_dust(typep(ind_part(j)))
      cosr(j) = is_cr(typep(ind_part(j)))
      mhd_pic(j) = is_cr_or_dust(typep(ind_part(j)))
   end do

  ! Various fields interpolated to particle positions
  ! Gather 3-velocity and 3-magnetic field
  uu(1:np,1:ndim)=0.0D0
  bb(1:np,1:ndim)=0.0D0
  if((pic_dust .or. pic_cr).and.hydro)then
     do ind=1,xtondim
        do idim=1,ndim
           do j=1,np
              uu(j,idim)=uu(j,idim)+uold(indp(j,ind),idim+1)/max(uold(indp(j,ind),1),smallr)*vol(j,ind)
              bb(j,idim)=bb(j,idim)+0.5D0*(uold(indp(j,ind),idim+5)+uold(indp(j,ind),idim+nvar))*vol(j,ind)
           end do
        end do
     end do
  endif


  ! ciso(1:np) = 0.0D0 ! Isothermal sound speed (squared).
  ! dgr(1:np) = 0.0D0 ! Gas density. Only used for trajectory file.
  ! if(pic_dust)then
  !    do ind=1,xtondim
  !        do j=1,np
  !          if(dust(j))then
  !             dgr(j)=dgr(j)+uold(indp(j,ind),1)*vol(j,ind)
  !             ciso(j)=ciso(j)+vol(j,ind)*&
  !             & sqrt((uold(indp(j,ind),5) - &
  !             &0.125D0*(uold(indp(j,ind),1+5)+uold(indp(j,ind),1+nvar))*(uold(indp(j,ind),1+5)+uold(indp(j,ind),1+nvar)) &
  !             &-0.125D0*(uold(indp(j,ind),2+5)+uold(indp(j,ind),2+nvar))*(uold(indp(j,ind),2+5)+uold(indp(j,ind),2+nvar)) &
  !             &-0.125D0*(uold(indp(j,ind),3+5)+uold(indp(j,ind),3+nvar))*(uold(indp(j,ind),3+5)+uold(indp(j,ind),3+nvar)) &
  !             &- 0.5D0*(&
  !             & uold(indp(j,ind),1+1)*uold(indp(j,ind),1+1)+ &
  !             & uold(indp(j,ind),2+1)*uold(indp(j,ind),2+1)+ &
  !             & uold(indp(j,ind),3+1)*uold(indp(j,ind),3+1) &
  !             &)/max(uold(indp(j,ind),1),smallr))&! Subtract from total energy agnetic and kinetic energies
  !             &*(gamma-1.0D0)/max(uold(indp(j,ind),1),smallr))! P/rho = ciso**2, to be interpolated.
  !           endif
  !       end do
  !    end do
  ! endif
  ciso(1:np) = 0.0D0 ! Isothermal sound speed (squared).
  dgr(1:np) = 0.0D0 ! Gas density. Only used for trajectory file.
  if(pic_dust)then
     do ind=1,xtondim
         do j=1,np
           if(dust(j))then
              dgr(j)=dgr(j)+uold(indp(j,ind),1)*vol(j,ind)
              tempc = uold(indp(j,ind),5) - &
              &0.125D0*(uold(indp(j,ind),1+5)+uold(indp(j,ind),1+nvar))**2 &
              &-0.125D0*(uold(indp(j,ind),2+5)+uold(indp(j,ind),2+nvar))**2&
              &-0.125D0*(uold(indp(j,ind),3+5)+uold(indp(j,ind),3+nvar))**2 &
              &- 0.5D0*(&
              & uold(indp(j,ind),1+1)**2+ &
              & uold(indp(j,ind),2+1)**2+ &
              & uold(indp(j,ind),3+1)**2)/max(uold(indp(j,ind),1),smallr) ! energy minus magnetic and kinetic.

              tempc=tempc*(gamma-1.0d0)/max(uold(indp(j,ind),1),smallr)
              tempc=min(1/smallr,max(tempc,smallr))
              ciso(j)=ciso(j)+vol(j,ind)*sqrt(tempc) ! P/rho = ciso**2, to be interpolated.
            endif
        end do
     end do
  endif
  ! tempc=1d100
  ! maxciso=-1d100
  ! if(pic_dust)then
  !   do j=1,np
  !     tempc=min(ciso(j),tempc)
  !     maxciso=max(ciso(j),maxciso)
  !   end do
  ! endif

  !write(*,*)'cisomin/max=',tempc,maxciso

  ddgr(1:np) = 0.0D0 ! Dust(+cr) density.
  if((pic_dust.or.pic_cr) .and. back_reaction .and. (trajectories(1) > 0 .or. ntrajectories > 0))then
     do ind=1,xtondim
         do j=1,np
           if(mhd_pic(j))then
              ddgr(j)=ddgr(j)+uold(indp(j,ind),ivar_mhd_pic)*vol(j,ind)
            endif
        end do
     end do
  endif

 if (trajectories(1)>0)then
    i=1
    do while(trajectories(i) .ne. 0)!Various fields interpolated to particle positions
         do j=1,np
             if( trajectories(i) .eq. idp(ind_part(j)) )then
                write(25+myid,*)t-dtnew(ilevel),idp(ind_part(j)),dgr(j),ddgr(j), & ! Old time
                     & xp(ind_part(j),1),xp(ind_part(j),2),xp(ind_part(j),3),& ! Old particle position
                     & vp(ind_part(j),1),vp(ind_part(j),2),vp(ind_part(j),3),& ! Old particle velocity
                     &  uu(j,1),uu(j,2),uu(j,3),& ! Old fluid velocity
                     &  bb(j,1),bb(j,2),bb(j,3)! Old magnetic field.
                     ! & new_vp(j,1),new_vp(j,2),new_vp(j,3) ! NEW particle velocity (for comparison)
             endif
        end do
       i=i+1
    end do
 endif

 if (ntrajectories>0)then
   do i=1,ntrajectories!Various fields interpolated to particle positions
        do j=1,np
            if(1+(i-1)*int(((ndust+ncr)*2**(3*levelmin))/ntrajectories) .eq. idp(ind_part(j)) )then
               write(25+myid,*)t-dtnew(ilevel),idp(ind_part(j)),dgr(j),ddgr(j), & ! Old time
                    & xp(ind_part(j),1),xp(ind_part(j),2),xp(ind_part(j),3),& ! Old particle position
                    & vp(ind_part(j),1),vp(ind_part(j),2),vp(ind_part(j),3),& ! Old particle velocity
                    &  uu(j,1),uu(j,2),uu(j,3),& ! Old fluid velocity
                    &  bb(j,1),bb(j,2),bb(j,3)! Old magnetic field.
                    ! & new_vp(j,1),new_vp(j,2),new_vp(j,3) ! NEW particle velocity (for comparison)
            endif
       end do
   end do
endif

  ! Gather 3-velocity
  ff(1:np,1:ndim)=0.0D0
  if(tracer.and.hydro)then
     do ind=1,xtondim
        do idim=1,ndim
          do j=1,np
             if (classical_tracer(j)) then
                ff(j,idim)=ff(j,idim) + &
                     uold(indp(j,ind),idim+1)/max(uold(indp(j,ind),1),smallr)*vol(j,ind)
             end if
          end do
        end do
     end do
  endif
  ! Gather 3-force
  if(poisson)then
     do ind=1,xtondim
        do idim=1,ndim
          do j=1,np
             if (.not. classical_tracer(j)) then
             ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*vol(j,ind)
             end if
          end do
        end do
#ifdef OUTPUT_PARTICLE_POTENTIAL
        do j=1,np
           ptcl_phi(ind_part(j)) = phi(indp(j,ind))
        end do
#endif
     end do
  endif

  if(pic_dust .and. ((accel_gr(1).ne.0).or.(accel_gr(2).ne.0).or.(accel_gr(3).ne.0)))then
     do idim=1,ndim
        do j=1,np
          if(dust(j))then
            ff(j,idim)=ff(j,idim)+accel_gr(idim)
          endif
        end do
     end do
  endif

  ! effective lorentz factors. vp is the 4-velocity for crs.
  lorentzf(1:np)=1.0d0
  if(pic_cr)then
    do j=1,np
      if(cosr(j))then
        lorentzf(j)=cr_c_fraction*sqrt(1.0d0+&
        &(vp(ind_part(j),1)*vp(ind_part(j),1)&
        &+vp(ind_part(j),2)*vp(ind_part(j),2)&
        &+vp(ind_part(j),3)*vp(ind_part(j),3))/(crsol*crsol))
      endif
    end do
  endif

  ! Update velocity
  do idim=1,ndim
     if(static)then
        do j=1,np
           new_vp(j,idim)=ff(j,idim)
        end do
     else
        do j=1,np
           if (classical_tracer(j)) then
              new_vp(j,idim)=ff(j,idim)
           else
              new_vp(j,idim)=vp(ind_part(j),idim)+ff(j,idim)*0.5D0*dtnew(ilevel)
           end if
        end do
     endif
  end do

  vv(1:np,1:ndim)=new_vp(1:np,1:ndim) ! post-external forces

  if((pic_dust.or.pic_cr).and.hydro)then
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! GRAIN SIZES AND CHARGES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sizes(1:np)=0.0d0
    charges(1:np)=0.0d0
    call PicParams(sizes,charges,dust,cosr,lorentzf,ind_part,np)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STOPPING RATE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We compute this before the EM kick for the sole reason that we had to do
  ! that in init_mhd_pic_fine. For second order accuracy, things will be more
  ! complicated.
  ! See about turning this off when no dust.
  call StoppingRate(np,dtnew(ilevel),indp,dust,vol,vv,nu_stop,ciso,dgr,sizes,xtondim)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! LORENTZ KICK
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ERM: Determining the evolved versions of the drift, then interpolating those
  ! For now, do the 8N calculations, but in the future, would be good to loop
  ! over CELLS rather than particles (esp. when doing TSC).
  ! Set unew's dust momentum slots to be the gas velocity.

    !call ResetUnewToFluidVel(ilevel)
    !call reset_unew(ilevel)
    big_vv(1:np,1:xtondim,1:ndim)=0.0D0 ! contains actual sub-cloud velocities.
    big_ww(1:np,1:xtondim,1:ndim)=0.0D0 ! Contains net mean drift velocity
    ! might want a "big_ww"? I think that's how I'll approach it.
    ! We want to evolve each of the subclouds. Knowing the new w will
    ! allow us to compute the evolution of the sub-clouds with the drag too.
if(pic_dust)then
    call EMKick(np,dtnew(ilevel),indp,mhd_pic,cosr,charges,lorentzf,&
    &ok,vol,vv,big_vv,big_ww,crsol,xtondim)
endif
if(pic_cr)then
#ifdef BORIS
    call BorisKick(np,dtnew(ilevel),indp,mhd_pic,cosr,charges,lorentzf,&
    &ok,vol,vv,big_vv,big_ww,crsol,xtondim)
#endif
#ifdef GORIS
    call GorisKick(np,dtnew(ilevel),indp,mhd_pic,cosr,charges,lorentzf,&
    &ok,vol,vv,big_vv,big_ww,crsol,xtondim)
! #else
!     call EMKick(np,dtnew(ilevel),indp,mhd_pic,cosr,charges,lorentzf,&
!     &ok,vol,vv,big_vv,big_ww,crsol,xtondim)
#endif
endif

    ! big_vv now contains changes to sub-cloud velocities. vv is still the old
    ! velocity. As well, unew's dust slot contains u**n+du**EM
    !write(*,*)'big_vv=',big_vv(1,1,1),big_vv(1,1,2),big_vv(1,1,3)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DRAG KICK
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call DragKick(np,dtnew(ilevel),indp,dust,ok,vol,nu_stop,big_vv,big_ww,vv,xtondim)
    ! Was just "DragKick"
    !write(*,*)'big_vv+=',big_vv(1,1,1),big_vv(1,1,2),big_vv(1,1,3)
    ! DragKick will modify big_ww as well as big_vv, but not vv.
    ! Now kick the dust given these quantities.
    vv(1:np,1:ndim)=new_vp(1:np,1:ndim)
    do ind=1,xtondim
      do idim=1,ndim
        do j=1,np ! DOES THIS NEED TO BE DIVIDED BY VOL_LOC? No, only if mp is involved.
          if(mhd_pic(j))then ! add up the cahnges to get the velocity post gas forces
            vv(j,idim)=vv(j,idim)+vol(j,ind)*(big_vv(j,ind,idim)-new_vp(j,idim))
          endif
        end do
      end do
    end do

    !call DragKick(np,dtnew(ilevel),indp,ok,vol,mov,nu_stop,big_vv,vv)
    vv(1:np,1:ndim)=vv(1:np,1:ndim)-new_vp(1:np,1:ndim) ! change in velocity
    new_vp(1:np,1:ndim)=new_vp(1:np,1:ndim)+vv(1:np,1:ndim)
     ! finally let new_vp be the post-gas-force velocity
  endif

  ! For sink cloud particle only
  if(sink)then
     ! Overwrite cloud particle velocity with sink velocity
     do idim=1,ndim
        do j=1,np
           if( is_cloud(typep(ind_part(j))) ) then
              isink=-idp(ind_part(j))
              new_vp(j,idim)=vsnew(isink,idim,ilevel)
           end if
        end do
     end do
  end if

  ! effective new lorentz factors. vp is the 4-velocity for crs.
  new_lorentzf(1:np)=1.0d0
  if(pic_cr)then
    do j=1,np ! NOTICE THE CR_C_FRACTION PRESENT HERE BUT NOT IN OTHER PLACES.
      if(cosr(j))then
        new_lorentzf(j)=cr_c_fraction*sqrt(1.0d0+&
        &(new_vp(j,1)*new_vp(j,1)&
        &+new_vp(j,2)*new_vp(j,2)&
        &+new_vp(j,3)*new_vp(j,3))/(crsol*crsol))
      endif
    end do
  endif


  ! Update position BEFORE setting new velocity using trapezoidal rule.
  do idim=1,ndim
     if(static)then
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)
        end do
     elseif(pic_cr)then
       do j=1,np
         if(cosr(j))then
           new_xp(j,idim)=xp(ind_part(j),idim)+0.5d0*cr_c_fraction*&
           &(new_vp(j,idim)/new_lorentzf(j)&
           &+vp(ind_part(j),idim)/lorentzf(j))*dtnew(ilevel)
         else
           new_xp(j,idim)=xp(ind_part(j),idim)+0.5d0*(new_vp(j,idim)+vp(ind_part(j),idim))*dtnew(ilevel)
         endif
       end do
     else
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)+0.5d0*(new_vp(j,idim)+vp(ind_part(j),idim))*dtnew(ilevel)
        end do
     endif
  end do

  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=new_xp(j,idim)
     end do
  end do



 ! Update effective cr masses
  ! do j=1,np
  !    new_mp(j)=mp(ind_part(j))*new_lorentzf(j)/lorentzf(j)
  ! end do

  ! Deposit minus final dust momentum to new gas momentum
if(back_reaction)then
  do ind=1,xtondim
     do idim=1,ndim
        do j=1,np
           if(ok(j).and.mhd_pic(j))then ! May need to change slightly for crs.
              unew(indp(j,ind),1+idim)=unew(indp(j,ind),1+idim)&
              &-mp(ind_part(j))*vv(j,idim)*vol(j,ind)/vol_loc
           ! elseif(ok(j).and.cosr(j))then
           !   unew(indp(j,ind),1+idim)=unew(indp(j,ind),1+idim)&
           !   &-(mp(ind_part(j))/lorentzf(j))*vv(j,idim)*vol(j,ind)/vol_loc
           !   !&-(cr_c_fraction*mp(ind_part(j))/lorentzf(j))*vv(j,idim)*vol(j,ind)/vol_loc ! Previous code.
           end if
        end do
     end do
  end do

  ! Deposit minus final dust and cr energy onto gas energy.
  do ind=1,xtondim
        do j=1,np
           if(ok(j).and.dust(j))then
              unew(indp(j,ind),5)=unew(indp(j,ind),5)&
              &-0.5*mp(ind_part(j))*&
              &(new_vp(j,idim)*vv(j,idim)-0.5d0*vv(j,idim)*vv(j,idim))&
              &*vol(j,ind)/vol_loc
           elseif(ok(j).and.cosr(j))then
             unew(indp(j,ind),5)=unew(indp(j,ind),5)&
             &-(new_lorentzf(j)-lorentzf(j))*mp(ind_part(j))*(crsol*crsol/(cr_c_fraction*cr_c_fraction))&
             &*vol(j,ind)/vol_loc
           end if
        end do
  end do
endif ! only need to do feedback steps if you have back_reaction


  ! ! Update effective cr masses
  ! do j=1,np
  !    if(cosr(j))then
  !       mp(ind_part(j))=new_mp(j)
  !    end if
  ! end do

  ! Store new velocity
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim)=new_vp(j,idim)
     end do
  end do



end subroutine move1
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine CRKick(nn,dt,indp,mhd_pic,cosm,ctms,lf,ok,vol,v,big_v,big_w,crsol,xtondim)
  ! This subroutine will compute changes to sub-cloud velocity in big_v,
  ! as well as set unew's dust momentum slot to being u+du**EM.
  ! need to add "lorentzf" array for lorentz factors
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer ::ivar_mhd_pic,xtondim ! cell-centered dust variables start.
  integer ::nn ! number of cells
  real(dp) ::dt,ctmav! timestep, average charge to mass ratio
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector)::ctms,lf! charge to mass ratios, lorentz factors, half-step lorentz factors
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector,1:xtondim,1:ndim)::big_v,big_w
  real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
  real(dp) ::den_dust,den_gas,den_tot,mu,den_cr,den_true,ctm,lfh,crsol,sol,piclf,dtp
  real(dp) ::den_dustp,den_gasp,den_totp,mup,den_crp,den_truep,ctmp,boost_mag,boost_lf,ulf,vcom_mag
  real(dp),dimension(1:3) ::vtemp,w,B,vcom,vh,utemp,vpic,Bp,u,vi,boost
  logical,dimension(1:nvector)::mhd_pic,cosm
  integer ::i,j,ind,idim! Just an -index
  sol=crsol/cr_c_fraction
  ivar_mhd_pic=9

if(back_reaction)then
  do ind=1,xtondim
     do i=1,nn
       if(mhd_pic(i))then
          den_gas=uold(indp(i,ind),1)
          den_dust=uold(indp(i,ind),ivar_mhd_pic) !dust mass only
          den_cr=unew(indp(i,ind),ivar_mhd_pic+2) ! mass of CRs /cr_c_fraction
          den_true=den_dust+cr_c_fraction*den_cr
          den_tot=den_dust+den_cr
          if(den_tot.eq.0)then
            write(*,*)'den_dust: ',den_dust
          endif
          !ctmav=unew(indp(i,ind),ivar_mhd_pic+1)/max(den_tot,smallr)! average with crs too.

                      ! if(.not.cosm(i))then ! Dust case... need to add this.
                      !
                      !   big_v(i,ind,1)=v(i,1)+& ! subcloud velocity update
                      !   &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
                      !   & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
                      !   &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)
                      !
                      !   big_v(i,ind,2)=v(i,2)+&
                      !   &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
                      !   & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
                      !   &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)
                      !
                      !   big_v(i,ind,3)=v(i,3)+&
                      !   &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
                      !   & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
                      !   &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)
                      ! endif
          ! you'll have to
            !mu=den_dust/max(den_gas,smallr)
          do idim=1,3
            B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
            utemp(idim)=uold(indp(i,ind),idim+1)/max(den_gas,smallr)
            vpic(idim)=uold(indp(i,ind),idim+ivar_mhd_pic)/max(den_true,smallr)
            !w(idim)=uold(indp(i,ind),idim+ivar_mhd_pic)/max(trum,smallr)&
            !&-uold(indp(i,ind),idim+1)/max(den_gas,smallr) ! True drift velocity
            vcom(idim)= (uold(indp(i,ind),idim+ivar_mhd_pic)&
            &+uold(indp(i,ind),idim+1))&
            &/max(den_gas+den_true,smallr)! effective CoM velocity
          end do
          u(1:3)=utemp(1:3)
          vi(1:3)=vpic(1:3)
          ! treat utemp as being very small with respect to the speed of light.
          ! Boost into CoM frame. Transform gas vel, CR vel, masses
          boost(1:3)=-1.0d0*vcom(1:3)/sol ! Definition of boost is opposite what it shoudl be below.
          boost_mag=sqrt(boost(1)**2+boost(2)**2+boost(3)**2)
          vcom_mag=sqrt(vcom(1)**2+vcom(2)**2+vcom(3)**2)
          boost_lf=1.0d0/sqrt(1.0d0-boost_mag**2)

          den_gasp=den_gasp*boost_lf&
          &-(boost_lf*boost(1)*uold(indp(i,ind),1+1)&
          &+boost_lf*boost(2)*uold(indp(i,ind),2+1)&
          &+boost_lf*boost(3)*uold(indp(i,ind),3+1))/crsol

          den_totp=den_tot*boost_lf&
          &-(boost_lf*boost(1)*uold(indp(i,ind),1+ivar_mhd_pic)&
          &+boost_lf*boost(2)*uold(indp(i,ind),2+ivar_mhd_pic)&
          &+boost_lf*boost(3)*uold(indp(i,ind),3+ivar_mhd_pic))/crsol

          mup=den_totp/den_gasp
          ctmav=unew(indp(i,ind),ivar_mhd_pic+1)/max(den_totp,smallr) !Works because we've effectively deposited a charge.
          if((vpic(1)**2+vpic(2)**2+vpic(3)**2).ge.sol**2)then
            write(*,*)'gammas'
            write(*,*)vpic(1),vpic(2),vpic(3)
            write(*,*)utemp(1),utemp(2),utemp(3)
            write(*,*)den_true,den_tot
            write(*,*)uold(indp(i,ind),1+ivar_mhd_pic),uold(indp(i,ind),2+ivar_mhd_pic),uold(indp(i,ind),3+ivar_mhd_pic)
          endif
          piclf=1./sqrt(1.-(vpic(1)**2+vpic(2)**2+vpic(3)**2)/sol**2)
          ulf=1./sqrt(1.-(utemp(1)**2+utemp(2)**2+utemp(3)**2)/sol**2)

          vpic(1:3)=piclf*vpic(1:3)
          utemp(1:3)=ulf*utemp(1:3)


          if(vcom_mag>0.0)then

            utemp(1)=-boost_lf*ulf*vcom(1)+&
            &(1.0d0+(boost_lf-1.0d0)*vcom(1)*vcom(1)/(vcom_mag**2))*utemp(1)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(1)*vcom(2)/(vcom_mag**2))*utemp(2)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(1)*vcom(3)/(vcom_mag**2))*utemp(3)


            utemp(2)=-boost_lf*ulf*vcom(2)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(2)*vcom(1)/(vcom_mag**2))*utemp(1)+&
            &(1.0d0+(boost_lf-1.0d0)*vcom(2)*vcom(2)/(vcom_mag**2))*utemp(2)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(2)*vcom(3)/(vcom_mag**2))*utemp(3)

            utemp(3)=-boost_lf*ulf*vcom(3)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(3)*vcom(1)/(vcom_mag**2))*utemp(1)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(3)*vcom(2)/(vcom_mag**2))*utemp(2)+&
            &(1.0d0+(boost_lf-1.0d0)*vcom(3)*vcom(3)/(vcom_mag**2))*utemp(3)

            vpic(1)=-boost_lf*piclf*vcom(1)+&
            &(1.0d0+(boost_lf-1.0d0)*vcom(1)*vcom(1)/(vcom_mag**2))*vpic(1)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(1)*vcom(2)/(vcom_mag**2))*vpic(2)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(1)*vcom(3)/(vcom_mag**2))*vpic(3)

            vpic(2)=-boost_lf*piclf*vcom(2)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(2)*vcom(1)/(vcom_mag**2))*vpic(1)+&
            &(1.0d0+(boost_lf-1.0d0)*vcom(2)*vcom(2)/(vcom_mag**2))*vpic(2)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(2)*vcom(3)/(vcom_mag**2))*vpic(3)

            vpic(3)=-boost_lf*piclf*vcom(3)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(3)*vcom(1)/(vcom_mag**2))*vpic(1)+&
            &(0.0d0+(boost_lf-1.0d0)*vcom(3)*vcom(2)/(vcom_mag**2))*vpic(2)+&
            &(1.0d0+(boost_lf-1.0d0)*vcom(3)*vcom(3)/(vcom_mag**2))*vpic(3)


            ! Transform the magnetic field too.
            Bp(1)=(-(vcom(2)*B(2)+vcom(3)*B(3))*&
            &(boost_lf*boost_mag**2*u(1)-(1.-boost_lf)*vcom(1))+&
            &B(1)*(vcom(1)**2+boost_lf*(vcom(2)**2+vcom(3)**2)&
            &+boost_lf*boost_mag**2*(vcom(2)*u(2)+vcom(3)*u(3))))/vcom_mag**2

            Bp(2)=(-(vcom(3)*B(3)+vcom(1)*B(1))*&
            &(boost_lf*boost_mag**2*u(2)-(1.-boost_lf)*vcom(2))+&
            &B(2)*(vcom(2)**2+boost_lf*(vcom(3)**2+vcom(1)**2)&
            &+boost_lf*boost_mag**2*(vcom(3)*u(3)+vcom(1)*u(1))))/vcom_mag**2

            Bp(3)=(-(vcom(1)*B(1)+vcom(2)*B(2))*&
            &(boost_lf*boost_mag**2*u(3)-(1.-boost_lf)*vcom(3))+&
            &B(3)*(vcom(3)**2+boost_lf*(vcom(1)**2+vcom(2)**2)&
            &+boost_lf*boost_mag**2*(vcom(1)*u(1)+vcom(2)*u(2))))/vcom_mag**2

          else
            Bp(1:3)=B(1:3)
          endif

          piclf=sqrt(1.0d0+(vpic(1)**2+vpic(2)**2+vpic(3)**2)/sol**2)
          ulf=sqrt(1.0d0+(utemp(1)**2+utemp(2)**2+utemp(3)**2)/sol**2)
          w(1:3)=vpic(1:3)/piclf-utemp(1:3)/ulf ! drift velocity in the COM frame.
          dtp=dt/boost_lf ! time-dilation matters for the timestep.




          ! Perform kick on gas velocity. I believe dt will have to be Lorentz transformed too.
          ! dteff = dt/gamma. Can try to figure this out in the simple case of a uniform force.
          ! We know how velocities transform and their relationship to one another.

          ! These velocity changes are drift velocity norm-preserving
          big_w(i,ind,1)=-1.*& !velocity changes to drift
          &(ctmav*dtp*(2.*Bp(2)**2*ctmav*dtp*(1.+mu)**2*w(1)+&
          &Bp(2)*(-2.*Bp(1)*ctmav*dtp*(1.+mu)**2*w(2) + 4.*(1.+mu)*w(3)) +&
          & Bp(3)*(2.*Bp(3)*ctmav*dtp*(1.+mu)**2*w(1)-4.*(1.+mu)*w(2) - 2*Bp(1)*ctmav*dtp*(1.+mu)**2*w(3))))/&
          &((4.+(Bp(1)**2+Bp(2)**2+Bp(3)**2)*ctmav**2*dtp**2*(1.+mu)**2))

          big_w(i,ind,2)=-1.*& !velocity changes to drift
          &(ctmav*dtp*(2.*Bp(3)**2*ctmav*dtp*(1.+mu)**2*w(2)+&
          &Bp(3)*(-2.*Bp(2)*ctmav*dtp*(1.+mu)**2*w(3) + 4.*(1.+mu)*w(1)) +&
          & Bp(1)*(2.*Bp(1)*ctmav*dtp*(1.+mu)**2*w(2)-4.*(1.+mu)*w(3) - 2*Bp(2)*ctmav*dtp*(1.+mu)**2*w(1))))/&
          &((4.+(Bp(1)**2+Bp(2)**2+Bp(3)**2)*ctmav**2*dtp**2*(1.+mu)**2))

          big_w(i,ind,3)=-1.*& !velocity changes to drift
          &(ctmav*dtp*(2.*Bp(1)**2*ctmav*dtp*(1.+mu)**2*w(3)+&
          &Bp(1)*(-2.*Bp(3)*ctmav*dtp*(1.+mu)**2*w(1) + 4.*(1.+mu)*w(2)) +&
          & Bp(2)*(2.*Bp(2)*ctmav*dtp*(1.+mu)**2*w(3)-4.*(1.+mu)*w(1) - 2*Bp(3)*ctmav*dtp*(1.+mu)**2*w(2))))/&
          &((4.+(Bp(1)**2+Bp(2)**2+Bp(3)**2)*ctmav**2*dtp**2*(1.+mu)**2))

          utemp(1:3)=utemp(1:3)-mu*ulf*big_w(i,ind,1:3)/(1+mu)
          vpic(1:3)=vpic(1:3)+piclf*big_w(i,ind,1:3)/(1+mu)

          ! Transform back into the Lab frame
          boost(1:3)=vcom(1:3)/sol ! Reverse the transform.

          utemp(1)=boost_lf*ulf*boost(1)*sol+&
          &(1.0d0+(boost_lf-1.0d0)*boost(1)*boost(1)/(boost_mag**2))*utemp(1)+&
          &(0.0d0+(boost_lf-1.0d0)*boost(1)*boost(2)/(boost_mag**2))*utemp(2)+&
          &(0.0d0+(boost_lf-1.0d0)*boost(1)*boost(3)/(boost_mag**2))*utemp(3)

          utemp(2)=boost_lf*ulf*boost(2)*sol+&
          &(0.0d0+(boost_lf-1.0d0)*boost(2)*boost(1)/(boost_mag**2))*utemp(1)+&
          &(1.0d0+(boost_lf-1.0d0)*boost(2)*boost(2)/(boost_mag**2))*utemp(2)+&
          &(0.0d0+(boost_lf-1.0d0)*boost(2)*boost(3)/(boost_mag**2))*utemp(3)

          utemp(3)=boost_lf*ulf*boost(3)*sol+&
          &(0.0d0+(boost_lf-1.0d0)*boost(3)*boost(1)/(boost_mag**2))*utemp(1)+&
          &(0.0d0+(boost_lf-1.0d0)*boost(3)*boost(2)/(boost_mag**2))*utemp(2)+&
          &(1.0d0+(boost_lf-1.0d0)*boost(3)*boost(3)/(boost_mag**2))*utemp(3)

          vpic(1)=boost_lf*piclf*boost(1)*sol+&
          &(1.0d0+(boost_lf-1.0d0)*boost(1)*boost(1)/(boost_mag**2))*vpic(1)+&
          &(0.0d0+(boost_lf-1.0d0)*boost(1)*boost(2)/(boost_mag**2))*vpic(2)+&
          &(0.0d0+(boost_lf-1.0d0)*boost(1)*boost(3)/(boost_mag**2))*vpic(3)

          vpic(2)=boost_lf*piclf*boost(2)*sol+&
          &(0.0d0+(boost_lf-1.0d0)*boost(2)*boost(1)/(boost_mag**2))*vpic(1)+&
          &(1.0d0+(boost_lf-1.0d0)*boost(2)*boost(2)/(boost_mag**2))*vpic(2)+&
          &(0.0d0+(boost_lf-1.0d0)*boost(2)*boost(3)/(boost_mag**2))*vpic(3)

          vpic(3)=boost_lf*piclf*boost(3)*sol+&
          &(0.0d0+(boost_lf-1.0d0)*boost(3)*boost(1)/(boost_mag**2))*vpic(1)+&
          &(0.0d0+(boost_lf-1.0d0)*boost(3)*boost(2)/(boost_mag**2))*vpic(2)+&
          &(1.0d0+(boost_lf-1.0d0)*boost(3)*boost(3)/(boost_mag**2))*vpic(3)

          piclf=sqrt(1.0d0+(vpic(1)**2+vpic(2)**2+vpic(3)**2)/sol**2)
          ulf=sqrt(1.0d0+(utemp(1)**2+utemp(2)**2+utemp(3)**2)/sol**2)

          utemp(1:3)=utemp(1:3)/ulf
          do idim=1,3
            big_w(i,ind,idim)=vpic(idim)/piclf-utemp(idim)
          end do
          ! Do 1/2 electric kick to compute lfh
          ! Half u is...
          u(1:3)=0.5d0*(utemp(1:3)+u(1:3))
          vh(1)=v(i,1)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(u(3)*B(2)-u(2)*B(3))
          vh(2)=v(i,2)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(u(1)*B(3)-u(3)*B(1))
          vh(3)=v(i,3)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(u(2)*B(1)-u(1)*B(2))

          lfh=sqrt(1.0d0+(vh(1)**2+vh(2)**2+vh(3)**2)/crsol**2)
          ctm=cr_c_fraction*cr_charge_to_mass/lfh
          ! Magnetic rotation

          vh(1)=vh(1)+&
          &(ctm*dt*(-2.*B(2)**2*ctm*dt*vh(1) + B(2)*(2.*B(1)*ctm*dt*vh(2) - 4.*vh(3)) +&
          & B(3)*(-2.*B(3)*ctm*dt*vh(1) + 4.*vh(2) + 2.*B(1)*ctm*dt*vh(3))))/&
          &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

          vh(2)=vh(2)+&
          &(ctm*dt*(-2.*B(3)**2*ctm*dt*vh(2) + B(3)*(2.*B(2)*ctm*dt*vh(3) - 4.*vh(1)) +&
          & B(1)*(-2.*B(1)*ctm*dt*vh(2) + 4.*vh(3) + 2.*B(2)*ctm*dt*vh(1))))/&
          &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

          vh(3)=vh(3)+&
          &(ctm*dt*(-2.*B(1)**2*ctm*dt*vh(3) + B(1)*(2.*B(3)*ctm*dt*vh(1) - 4.*vh(2)) +&
          & B(2)*(-2.*B(2)*ctm*dt*vh(3) + 4.*vh(1) + 2.*B(3)*ctm*dt*vh(2))))/&
          &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

          ! Final 1/2 electric kick.

          big_v(i,ind,1)=vh(1)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(u(3)*B(2)-u(2)*B(3))
          big_v(i,ind,2)=vh(2)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(u(1)*B(3)-u(3)*B(1))
          big_v(i,ind,3)=vh(3)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(u(2)*B(1)-u(1)*B(2))


        endif
     end do
  end do
else ! Need to update to what you choose to do in the case of BR present.
  do ind=1,xtondim
     do i=1,nn
       if(mhd_pic(i))then
          den_gas=uold(indp(i,ind),1)
          if (.not. cosm(i))then
            do idim=1,3
              B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
              vtemp(idim) = v(i,idim)-lf(i)*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            big_v(i,ind,1)=v(i,1)+& ! subcloud velocity update
            &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            big_v(i,ind,2)=v(i,2)+&
            &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            big_v(i,ind,3)=v(i,3)+&
            &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)
          else ! Need to get intermediate steps for lorentz factor and ctm
            do idim=1,3
              B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
              vtemp(idim) = v(i,idim)-lf(i)*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            vh(1)=v(i,1)+& ! subcloud velocity update
            &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            vh(2)=v(i,2)+&
            &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            vh(3)=v(i,3)+&
            &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            lfh = 0.5*(lf(i)+cr_c_fraction*sqrt(1.+(vh(1)*vh(1)+vh(2)*vh(2)+vh(3)*vh(3))/(crsol*crsol)))
            ctm = ctms(i)*lf(i)/lfh

            do idim=1,3
              vtemp(idim) = v(i,idim)-lfh*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            big_v(i,ind,1)=v(i,1)+& ! true subcloud velocity update
            &(ctm*dt*(-2.*B(2)**2*ctm*dt*vtemp(1) + B(2)*(2.*B(1)*ctm*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctm*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctm*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

            big_v(i,ind,2)=v(i,2)+&
            &(ctm*dt*(-2.*B(3)**2*ctm*dt*vtemp(2) + B(3)*(2.*B(2)*ctm*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctm*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctm*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

            big_v(i,ind,3)=v(i,3)+&
            &(ctm*dt*(-2.*B(1)**2*ctm*dt*vtemp(3) + B(1)*(2.*B(3)*ctm*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctm*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctm*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

          endif
        endif
     end do
  end do
endif
end subroutine CRKick
!#########################################################################
!#########################################################################

!#########################################################################
!#########################################################################
subroutine GorisKick(nn,dt,indp,mhd_pic,cosm,ctms,lf,ok,vol,v,big_v,big_w,crsol,xtondim)
  ! Generalized boris algorithm...
  ! This subroutine will compute changes to sub-cloud velocity in big_v,
  ! as well as set unew's dust momentum slot to being u+du**EM.
  ! need to add "lorentzf" array for lorentz factors
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer ::ivar_mhd_pic,xtondim ! cell-centered dust variables start.
  integer ::nn ! number of cells
  real(dp) ::dt,ctmav! timestep, average charge to mass ratio
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector)::ctms,lf! charge to mass ratios, lorentz factors, half-step lorentz factors
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector,1:xtondim,1:ndim)::big_v,big_w
  real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
  real(dp) ::den_dust,den_gas,den_true,den_tot,den_cr,mu,ctm,lfh,crsol,ctmeff,ctmavm
  real(dp):: gmo,gmn,sol
  real(dp),dimension(1:3) ::vtemp,w,B,vcom,vh,u,eeff,vt,wm, wp,uh,un
  logical,dimension(1:nvector)::mhd_pic,cosm
  integer ::i,j,ind,idim! Just an -index

  ivar_mhd_pic=9
  sol=crsol/cr_c_fraction
if(back_reaction)then
  do ind=1,xtondim
     do i=1,nn
       if(cosm(i))then
         den_gas=uold(indp(i,ind),1)
         den_dust=uold(indp(i,ind),ivar_mhd_pic) ! Cr mass (over cr_c_fraction)
         ctmav=unew(indp(i,ind),ivar_mhd_pic+1)/max(den_dust,smallr)! average with crs too.
         mu=den_dust/max(den_gas,smallr)

          do idim=1,3
            B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
            w(idim)=uold(indp(i,ind),idim+ivar_mhd_pic)/max(den_dust*cr_c_fraction,smallr)&
            &-uold(indp(i,ind),idim+1)/max(den_gas,smallr) ! True drift velocity
            vt(idim) = mu*&
            &uold(indp(i,ind),idim+ivar_mhd_pic)/max(den_dust*cr_c_fraction,smallr)&
            &+uold(indp(i,ind),idim+1)/max(den_gas,smallr) ! eps*gamma*v + u
          end do
          ! "effective electric field"
          eeff(1)=(cr_charge_to_mass*cr_c_fraction-ctmav)*(vt(2)*B(3)-vt(3)*B(2))
          eeff(2)=(cr_charge_to_mass*cr_c_fraction-ctmav)*(vt(3)*B(1)-vt(1)*B(3))
          eeff(3)=(cr_charge_to_mass*cr_c_fraction-ctmav)*(vt(1)*B(2)-vt(2)*B(1))

          do idim=1,3
            wm(idim)=w(idim)-0.5*dt*eeff(idim) ! half "electric" kick.
            vtemp(idim)=(vt(idim)+wm(idim))/(1+mu) !new 4-vel
            vh(idim)=(vt(idim)+w(idim))/(1+mu) !old 4-vel
          end do
          ! Now need to apply gamma_old/gamma_new to the ctmav
          gmn=sqrt(1.+(vtemp(1)**2+vtemp(2)**2+vtemp(3)**2)/sol**2)
          gmo=sqrt(1.+(vh(1)**2+vh(2)**2+vh(3)**2)/sol**2)

          ctmavm=ctmav*gmo/gmn
          ctmeff=(ctmavm+mu*cr_charge_to_mass*cr_c_fraction)/(1.+mu)

          wp(1)=wm(1)-1.*& !velocity changes to drift
          &(ctmeff*dt*(2.*B(2)**2*ctmeff*dt*(1.+mu)**2*wm(1)+&
          &B(2)*(-2.*B(1)*ctmeff*dt*(1.+mu)**2*wm(2) + 4.*(1.+mu)*wm(3)) +&
          & B(3)*(2.*B(3)*ctmeff*dt*(1.+mu)**2*wm(1)-4.*(1.+mu)*wm(2) - 2*B(1)*ctmeff*dt*(1.+mu)**2*wm(3))))/&
          &((4.+(B(1)**2+B(2)**2+B(3)**2)*ctmeff**2*dt**2*(1.+mu)**2))

          wp(2)=wm(2)-1.*& !velocity changes to drift
          &(ctmeff*dt*(2.*B(3)**2*ctmeff*dt*(1.+mu)**2*wm(2)+&
          &B(3)*(-2.*B(2)*ctmeff*dt*(1.+mu)**2*wm(3) + 4.*(1.+mu)*wm(1)) +&
          & B(1)*(2.*B(1)*ctmeff*dt*(1.+mu)**2*wm(2)-4.*(1.+mu)*wm(3) - 2*B(2)*ctmeff*dt*(1.+mu)**2*wm(1))))/&
          &((4.+(B(1)**2+B(2)**2+B(3)**2)*ctmeff**2*dt**2*(1.+mu)**2))

          wp(3)=wm(3)-1.*& !velocity changes to drift
          &(ctmeff*dt*(2.*B(1)**2*ctmeff*dt*(1.+mu)**2*wm(3)+&
          &B(1)*(-2.*B(3)*ctmeff*dt*(1.+mu)**2*wm(1) + 4.*(1.+mu)*wm(2)) +&
          & B(2)*(2.*B(2)*ctmeff*dt*(1.+mu)**2*wm(3)-4.*(1.+mu)*wm(1) - 2*B(3)*ctmeff*dt*(1.+mu)**2*wm(2))))/&
          &((4.+(B(1)**2+B(2)**2+B(3)**2)*ctmeff**2*dt**2*(1.+mu)**2))

          write(*,*)wp(1)**2+wp(2)**2+wp(3)**2-wm(1)**2-wm(2)**2-wm(3)**2

          eeff(1)=(cr_charge_to_mass*cr_c_fraction-ctmavm)*(vt(2)*B(3)-vt(3)*B(2))
          eeff(2)=(cr_charge_to_mass*cr_c_fraction-ctmavm)*(vt(3)*B(1)-vt(1)*B(3))
          eeff(3)=(cr_charge_to_mass*cr_c_fraction-ctmavm)*(vt(1)*B(2)-vt(2)*B(1))

          do idim=1,3
            big_w(i,ind,idim)=wp(idim)-0.5*dt*eeff(idim) ! half "electric" kick, now w^(n+1)
          end do

          do idim=1,ndim
            uh(idim)=0.5*((vt(idim)-mu*big_w(i,ind,idim))/(1+mu)+&
            &uold(indp(i,ind),idim+1)/max(den_gas,smallr))
            un(idim)=uold(indp(i,ind),idim+1)/max(den_gas,smallr)
            !uh(idim)=(vt(idim)-mu*0.5*(big_w(i,ind,idim)+w(idim)))/(1+mu) ! u(n+1/2)
            !uh(idim)=(vt(idim)-mu*(big_w(i,ind,idim)))/(1+mu) ! u(n+1/2)
            vtemp(idim) = v(i,idim)/cr_c_fraction
          end do

              vh(1)=vtemp(1)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(3)*B(2)-uh(2)*B(3))
              vh(2)=vtemp(2)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(1)*B(3)-uh(3)*B(1))
              vh(3)=vtemp(3)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(2)*B(1)-uh(1)*B(2))

              lfh=sqrt(1.0d0+(vh(1)**2+vh(2)**2+vh(3)**2)/sol**2)
              ctm=cr_c_fraction*cr_charge_to_mass/lfh

              vtemp(1)=vh(1)+&
              &(ctm*dt*(-2.*B(2)**2*ctm*dt*vh(1) + B(2)*(2.*B(1)*ctm*dt*vh(2) - 4.*vh(3)) +&
              & B(3)*(-2.*B(3)*ctm*dt*vh(1) + 4.*vh(2) + 2.*B(1)*ctm*dt*vh(3))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

              vtemp(2)=vh(2)+&
              &(ctm*dt*(-2.*B(3)**2*ctm*dt*vh(2) + B(3)*(2.*B(2)*ctm*dt*vh(3) - 4.*vh(1)) +&
              & B(1)*(-2.*B(1)*ctm*dt*vh(2) + 4.*vh(3) + 2.*B(2)*ctm*dt*vh(1))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

              vtemp(3)=vh(3)+&
              &(ctm*dt*(-2.*B(1)**2*ctm*dt*vh(3) + B(1)*(2.*B(3)*ctm*dt*vh(1) - 4.*vh(2)) +&
              & B(2)*(-2.*B(2)*ctm*dt*vh(3) + 4.*vh(1) + 2.*B(3)*ctm*dt*vh(2))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

              ! Final 1/2 electric kick.

              big_v(i,ind,1)=(vtemp(1)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(3)*B(2)-uh(2)*B(3)))*cr_c_fraction
              big_v(i,ind,2)=(vtemp(2)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(1)*B(3)-uh(3)*B(1)))*cr_c_fraction
              big_v(i,ind,3)=(vtemp(3)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(2)*B(1)-uh(1)*B(2)))*cr_c_fraction
        endif
     end do
  end do
else ! Need to update to what you choose to do in the case of BR present.
  do ind=1,xtondim
     do i=1,nn
       if(mhd_pic(i))then
          den_gas=uold(indp(i,ind),1)
          if (.not. cosm(i))then
            do idim=1,3
              B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
              vtemp(idim) = v(i,idim)-lf(i)*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            big_v(i,ind,1)=v(i,1)+& ! subcloud velocity update
            &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            big_v(i,ind,2)=v(i,2)+&
            &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            big_v(i,ind,3)=v(i,3)+&
            &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)
          else ! Need to get intermediate steps for lorentz factor and ctm
            do idim=1,3
              B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
              vtemp(idim) = v(i,idim)-lf(i)*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            vh(1)=v(i,1)+& ! subcloud velocity update
            &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            vh(2)=v(i,2)+&
            &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            vh(3)=v(i,3)+&
            &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            lfh = 0.5*(lf(i)+cr_c_fraction*sqrt(1.+(vh(1)*vh(1)+vh(2)*vh(2)+vh(3)*vh(3))/(crsol*crsol)))
            ctm = ctms(i)*lf(i)/lfh

            do idim=1,3
              vtemp(idim) = v(i,idim)-lfh*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            big_v(i,ind,1)=v(i,1)+& ! true subcloud velocity update
            &(ctm*dt*(-2.*B(2)**2*ctm*dt*vtemp(1) + B(2)*(2.*B(1)*ctm*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctm*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctm*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

            big_v(i,ind,2)=v(i,2)+&
            &(ctm*dt*(-2.*B(3)**2*ctm*dt*vtemp(2) + B(3)*(2.*B(2)*ctm*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctm*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctm*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

            big_v(i,ind,3)=v(i,3)+&
            &(ctm*dt*(-2.*B(1)**2*ctm*dt*vtemp(3) + B(1)*(2.*B(3)*ctm*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctm*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctm*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

          endif
        endif
     end do
  end do
endif
end subroutine GorisKick
!#########################################################################
!#########################################################################

!#########################################################################
!#########################################################################
subroutine BorisKick(nn,dt,indp,mhd_pic,cosm,ctms,lf,ok,vol,v,big_v,big_w,crsol,xtondim)
  ! Generalized boris algorithm...
  ! This subroutine will compute changes to sub-cloud velocity in big_v,
  ! as well as set unew's dust momentum slot to being u+du**EM.
  ! need to add "lorentzf" array for lorentz factors
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer ::ivar_mhd_pic,xtondim ! cell-centered dust variables start.
  integer ::nn ! number of cells
  real(dp) ::dt,ctmav! timestep, average charge to mass ratio
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector)::ctms,lf! charge to mass ratios, lorentz factors, half-step lorentz factors
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector,1:xtondim,1:ndim)::big_v,big_w
  real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
  real(dp) ::den_dust,den_gas,den_true,den_tot,den_cr,mu,ctm,lfh,crsol,ctmeff,ctmavm
  real(dp):: gmo,gmn,sol
  real(dp),dimension(1:3) ::vtemp,w,B,vcom,vh,u,eeff,vt,wm, wp,uh
  logical,dimension(1:nvector)::mhd_pic,cosm
  integer ::i,j,ind,idim! Just an -index

  ivar_mhd_pic=9
  sol=crsol/cr_c_fraction
if(back_reaction)then
  do ind=1,xtondim
     do i=1,nn
       if(cosm(i))then
         den_gas=uold(indp(i,ind),1)
         den_dust=uold(indp(i,ind),ivar_mhd_pic) ! Cr mass (over cr_c_fraction)
         ctmav=unew(indp(i,ind),ivar_mhd_pic+1)/max(den_dust,smallr)! average with crs too.
         mu=den_dust/max(den_gas,smallr)

          do idim=1,3
            B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
            w(idim)=uold(indp(i,ind),idim+ivar_mhd_pic)/max(den_dust*cr_c_fraction,smallr)&
            &-uold(indp(i,ind),idim+1)/max(den_gas,smallr) ! True drift velocity
            uh(idim)=uold(indp(i,ind),idim+1)/max(den_gas,smallr)
            vt(idim) = (den_dust/max(den_gas,smallr))*&
            &uold(indp(i,ind),idim+ivar_mhd_pic)/max(den_dust*cr_c_fraction,smallr)&
            &+uold(indp(i,ind),idim+1)/max(den_gas,smallr) ! eps*gamma*v + u
            vtemp(idim) = v(i,idim)/cr_c_fraction ! "true" 4-vel (not reduced)
            big_w(i,ind,idim)=w(idim)
          end do

              vh(1)=vtemp(1)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(3)*B(2)-uh(2)*B(3))
              vh(2)=vtemp(2)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(1)*B(3)-uh(3)*B(1))
              vh(3)=vtemp(3)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(2)*B(1)-uh(1)*B(2))

              lfh=sqrt(1.0d0+(vh(1)**2+vh(2)**2+vh(3)**2)/sol**2)
              ctm=cr_c_fraction*cr_charge_to_mass/lfh

              vtemp(1)=vh(1)+&
              &(ctm*dt*(-2.*B(2)**2*ctm*dt*vh(1) + B(2)*(2.*B(1)*ctm*dt*vh(2) - 4.*vh(3)) +&
              & B(3)*(-2.*B(3)*ctm*dt*vh(1) + 4.*vh(2) + 2.*B(1)*ctm*dt*vh(3))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

              vtemp(2)=vh(2)+&
              &(ctm*dt*(-2.*B(3)**2*ctm*dt*vh(2) + B(3)*(2.*B(2)*ctm*dt*vh(3) - 4.*vh(1)) +&
              & B(1)*(-2.*B(1)*ctm*dt*vh(2) + 4.*vh(3) + 2.*B(2)*ctm*dt*vh(1))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

              vtemp(3)=vh(3)+&
              &(ctm*dt*(-2.*B(1)**2*ctm*dt*vh(3) + B(1)*(2.*B(3)*ctm*dt*vh(1) - 4.*vh(2)) +&
              & B(2)*(-2.*B(2)*ctm*dt*vh(3) + 4.*vh(1) + 2.*B(3)*ctm*dt*vh(2))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

              ! Final 1/2 electric kick. Restore original magnitude.

              big_v(i,ind,1)=(vtemp(1)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(3)*B(2)-uh(2)*B(3)))*cr_c_fraction
              big_v(i,ind,2)=(vtemp(2)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(1)*B(3)-uh(3)*B(1)))*cr_c_fraction
              big_v(i,ind,3)=(vtemp(3)+cr_c_fraction*cr_charge_to_mass*0.5*dt*(uh(2)*B(1)-uh(1)*B(2)))*cr_c_fraction
        endif
     end do
  end do
else ! Need to update to what you choose to do in the case of BR present.
  do ind=1,xtondim
     do i=1,nn
       if(mhd_pic(i))then
          den_gas=uold(indp(i,ind),1)
          if (.not. cosm(i))then
            do idim=1,3
              B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
              vtemp(idim) = v(i,idim)-lf(i)*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            big_v(i,ind,1)=v(i,1)+& ! subcloud velocity update
            &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            big_v(i,ind,2)=v(i,2)+&
            &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            big_v(i,ind,3)=v(i,3)+&
            &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)
          else ! Need to get intermediate steps for lorentz factor and ctm
            do idim=1,3
              B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
              vtemp(idim) = v(i,idim)-lf(i)*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            vh(1)=v(i,1)+& ! subcloud velocity update
            &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            vh(2)=v(i,2)+&
            &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            vh(3)=v(i,3)+&
            &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            lfh = 0.5*(lf(i)+cr_c_fraction*sqrt(1.+(vh(1)*vh(1)+vh(2)*vh(2)+vh(3)*vh(3))/(crsol*crsol)))
            ctm = ctms(i)*lf(i)/lfh

            do idim=1,3
              vtemp(idim) = v(i,idim)-lfh*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            big_v(i,ind,1)=v(i,1)+& ! true subcloud velocity update
            &(ctm*dt*(-2.*B(2)**2*ctm*dt*vtemp(1) + B(2)*(2.*B(1)*ctm*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctm*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctm*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

            big_v(i,ind,2)=v(i,2)+&
            &(ctm*dt*(-2.*B(3)**2*ctm*dt*vtemp(2) + B(3)*(2.*B(2)*ctm*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctm*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctm*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

            big_v(i,ind,3)=v(i,3)+&
            &(ctm*dt*(-2.*B(1)**2*ctm*dt*vtemp(3) + B(1)*(2.*B(3)*ctm*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctm*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctm*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

          endif
        endif
     end do
  end do
endif
end subroutine BorisKick
!#########################################################################
!#########################################################################

!#########################################################################
!#########################################################################
subroutine EMKick(nn,dt,indp,mhd_pic,cosm,ctms,lf,ok,vol,v,big_v,big_w,crsol,xtondim)
  ! This subroutine will compute changes to sub-cloud velocity in big_v,
  ! as well as set unew's dust momentum slot to being u+du**EM.
  ! need to add "lorentzf" array for lorentz factors
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer ::ivar_mhd_pic,xtondim ! cell-centered dust variables start.
  integer ::nn ! number of cells
  real(dp) ::dt,ctmav! timestep, average charge to mass ratio
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector)::ctms,lf! charge to mass ratios, lorentz factors, half-step lorentz factors
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector,1:xtondim,1:ndim)::big_v,big_w
  real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
  real(dp) ::den_dust,den_gas,mu,ctm,lfh,crsol
  real(dp),dimension(1:3) ::vtemp,w,B,vcom,vh
  logical,dimension(1:nvector)::mhd_pic,cosm
  integer ::i,j,ind,idim! Just an -index

  ivar_mhd_pic=9

if(back_reaction)then
  do ind=1,xtondim
     do i=1,nn
       if(mhd_pic(i))then
         den_gas=uold(indp(i,ind),1)
         den_dust=uold(indp(i,ind),ivar_mhd_pic) !dust mass only
         ctmav=unew(indp(i,ind),ivar_mhd_pic+1)/max(den_dust,smallr)! average with crs too.
         mu=den_dust/max(den_gas,smallr)

          do idim=1,3
            B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
            w(idim)=uold(indp(i,ind),idim+ivar_mhd_pic)/max(den_dust,smallr)&
            &-uold(indp(i,ind),idim+1)/max(den_gas,smallr) ! True drift velocity
          end do

          big_w(i,ind,1)=-1.*& !velocity changes to drift
          &(ctmav*dt*(2.*B(2)**2*ctmav*dt*(1.+mu)**2*w(1)+&
          &B(2)*(-2.*B(1)*ctmav*dt*(1.+mu)**2*w(2) + 4.*(1.+mu)*w(3)) +&
          & B(3)*(2.*B(3)*ctmav*dt*(1.+mu)**2*w(1)-4.*(1.+mu)*w(2) - 2*B(1)*ctmav*dt*(1.+mu)**2*w(3))))/&
          &((4.+(B(1)**2+B(2)**2+B(3)**2)*ctmav**2*dt**2*(1.+mu)**2))

          big_w(i,ind,2)=-1.*& !velocity changes to drift
          &(ctmav*dt*(2.*B(3)**2*ctmav*dt*(1.+mu)**2*w(2)+&
          &B(3)*(-2.*B(2)*ctmav*dt*(1.+mu)**2*w(3) + 4.*(1.+mu)*w(1)) +&
          & B(1)*(2.*B(1)*ctmav*dt*(1.+mu)**2*w(2)-4.*(1.+mu)*w(3) - 2*B(2)*ctmav*dt*(1.+mu)**2*w(1))))/&
          &((4.+(B(1)**2+B(2)**2+B(3)**2)*ctmav**2*dt**2*(1.+mu)**2))

          big_w(i,ind,3)=-1.*& !velocity changes to drift
          &(ctmav*dt*(2.*B(1)**2*ctmav*dt*(1.+mu)**2*w(3)+&
          &B(1)*(-2.*B(3)*ctmav*dt*(1.+mu)**2*w(1) + 4.*(1.+mu)*w(2)) +&
          & B(2)*(2.*B(2)*ctmav*dt*(1.+mu)**2*w(3)-4.*(1.+mu)*w(1) - 2*B(3)*ctmav*dt*(1.+mu)**2*w(2))))/&
          &((4.+(B(1)**2+B(2)**2+B(3)**2)*ctmav**2*dt**2*(1.+mu)**2))

          do idim=1,ndim
            vtemp(idim) = v(i,idim)-lf(i)*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)&
            &+0.5*mu*big_w(i,ind,idim)/(1.+mu)) ! v^n - u^(n+1)
            big_w(i,ind,idim)=w(idim)+big_w(i,ind,idim) !w^(n+1)
          end do
            if(.not.cosm(i))then

              big_v(i,ind,1)=v(i,1)+& ! subcloud velocity update
              &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
              & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

              big_v(i,ind,2)=v(i,2)+&
              &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
              & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

              big_v(i,ind,3)=v(i,3)+&
              &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
              & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)
            else
              ! vcom will be the COM velocity. vh will be the scratch space.
              ! Boost into the COM, then do the operation, then boost out of it.
              ! Make sure to account for ctm changing in the boost, as well as
              ! the lorentz factor!

              ! vp(ipart,1)=cr_boost_lf*cr_lorentzf*cr_boost(1)+&
              ! &(1.0d0+(cr_boost_lf-1.0d0)*cr_boost(1)*cr_boost(1)/(cr_boost_mag**2))*vp(ipart,1)+&
              ! &(0.0d0+(cr_boost_lf-1.0d0)*cr_boost(1)*cr_boost(2)/(cr_boost_mag**2))*vp(ipart,2)+&
              ! &(0.0d0+(cr_boost_lf-1.0d0)*cr_boost(1)*cr_boost(3)/(cr_boost_mag**2))*vp(ipart,3)
              !
              ! vp(ipart,2)=cr_boost_lf*cr_lorentzf*cr_boost(2)+&
              ! &(0.0d0+(cr_boost_lf-1.0d0)*cr_boost(2)*cr_boost(1)/(cr_boost_mag**2))*vp(ipart,1)+&
              ! &(1.0d0+(cr_boost_lf-1.0d0)*cr_boost(2)*cr_boost(2)/(cr_boost_mag**2))*vp(ipart,2)+&
              ! &(0.0d0+(cr_boost_lf-1.0d0)*cr_boost(2)*cr_boost(3)/(cr_boost_mag**2))*vp(ipart,3)
              !
              ! vp(ipart,3)=cr_boost_lf*cr_lorentzf*cr_boost(3)+&
              ! &(0.0d0+(cr_boost_lf-1.0d0)*cr_boost(3)*cr_boost(1)/(cr_boost_mag**2))*vp(ipart,1)+&
              ! &(0.0d0+(cr_boost_lf-1.0d0)*cr_boost(3)*cr_boost(2)/(cr_boost_mag**2))*vp(ipart,2)+&
              ! &(1.0d0+(boost_lf-1.0d0)*cr_boost(3)*cr_boost(3)/(cr_boost_mag**2))*vp(ipart,3)


              ! Compute half lorentz factor (g^n+g^(n+1))/2 and charge-to-mass ratio
              !lfh = 0.5*(lf(i)+cr_c_fraction*sqrt(1.+(vh(i,1)*vh(i,1)+vh(i,2)*vh(i,2)+vh(i,3)*vh(i,3))/(crsol*crsol)))
              vh(1:3)=0.5*(v(i,1:3)+vh(1:3)) ! Define a half-velocity to get a true predictor-corrector scheme
              ! Try out a boris-solver but using E = your average E from the gas-vel update. Then try a boris solver without that entirely.
              ! Then also with the lfh from that paper on your phone. Basically, do an electric kick.
              lfh = cr_c_fraction*sqrt(1.+(vh(1)*vh(1)+vh(2)*vh(2)+vh(3)*vh(3))/(crsol*crsol)) ! forward looking
              ctm = ctms(i)*lf(i)/lfh

              do idim=1,ndim ! Try not modifying the lf in front of fluid vel
                vtemp(idim) = v(i,idim)-lf(i)*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)&
                &+0.5*mu*big_w(i,ind,idim)/(1.+mu)) ! v^n - u^(n+1)
                big_w(i,ind,idim)=w(idim)+big_w(i,ind,idim) !w^(n+1)
              end do

              big_v(i,ind,1)=v(i,1)+& ! true subcloud velocity update
              &(ctm*dt*(-2.*B(2)**2*ctm*dt*vtemp(1) + B(2)*(2.*B(1)*ctm*dt*vtemp(2) - 4.*vtemp(3)) +&
              & B(3)*(-2.*B(3)*ctm*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctm*dt*vtemp(3))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

              big_v(i,ind,2)=v(i,2)+&
              &(ctm*dt*(-2.*B(3)**2*ctm*dt*vtemp(2) + B(3)*(2.*B(2)*ctm*dt*vtemp(3) - 4.*vtemp(1)) +&
              & B(1)*(-2.*B(1)*ctm*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctm*dt*vtemp(1))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

              big_v(i,ind,3)=v(i,3)+&
              &(ctm*dt*(-2.*B(1)**2*ctm*dt*vtemp(3) + B(1)*(2.*B(3)*ctm*dt*vtemp(1) - 4.*vtemp(2)) +&
              & B(2)*(-2.*B(2)*ctm*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctm*dt*vtemp(2))))/&
              &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)
            endif
        endif
     end do
  end do
else ! Need to update to what you choose to do in the case of BR present.
  do ind=1,xtondim
     do i=1,nn
       if(mhd_pic(i))then
          den_gas=uold(indp(i,ind),1)
          if (.not. cosm(i))then
            do idim=1,3
              B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
              vtemp(idim) = v(i,idim)-lf(i)*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            big_v(i,ind,1)=v(i,1)+& ! subcloud velocity update
            &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            big_v(i,ind,2)=v(i,2)+&
            &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            big_v(i,ind,3)=v(i,3)+&
            &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)
          else ! Need to get intermediate steps for lorentz factor and ctm
            do idim=1,3
              B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
              vtemp(idim) = v(i,idim)-lf(i)*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            vh(1)=v(i,1)+& ! subcloud velocity update
            &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            vh(2)=v(i,2)+&
            &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            vh(3)=v(i,3)+&
            &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

            lfh = 0.5*(lf(i)+cr_c_fraction*sqrt(1.+(vh(1)*vh(1)+vh(2)*vh(2)+vh(3)*vh(3))/(crsol*crsol)))
            ctm = ctms(i)*lf(i)/lfh

            do idim=1,3
              vtemp(idim) = v(i,idim)-lfh*(uold(indp(i,ind),1+idim)/max(den_gas,smallr)) ! v^n - u^n
            end do

            big_v(i,ind,1)=v(i,1)+& ! true subcloud velocity update
            &(ctm*dt*(-2.*B(2)**2*ctm*dt*vtemp(1) + B(2)*(2.*B(1)*ctm*dt*vtemp(2) - 4.*vtemp(3)) +&
            & B(3)*(-2.*B(3)*ctm*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctm*dt*vtemp(3))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

            big_v(i,ind,2)=v(i,2)+&
            &(ctm*dt*(-2.*B(3)**2*ctm*dt*vtemp(2) + B(3)*(2.*B(2)*ctm*dt*vtemp(3) - 4.*vtemp(1)) +&
            & B(1)*(-2.*B(1)*ctm*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctm*dt*vtemp(1))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

            big_v(i,ind,3)=v(i,3)+&
            &(ctm*dt*(-2.*B(1)**2*ctm*dt*vtemp(3) + B(1)*(2.*B(3)*ctm*dt*vtemp(1) - 4.*vtemp(2)) +&
            & B(2)*(-2.*B(2)*ctm*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctm*dt*vtemp(2))))/&
            &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctm**2*dt**2)

          endif
        endif
     end do
  end do
endif
end subroutine EMKick
!#########################################################################
!#########################################################################


!#########################################################################
!#########################################################################
subroutine StoppingRate(nn,dt,indp,dust,vol,v,nu,c,dgr,gs,xtondim)
  ! The following subroutine will alter its last argument, nu
  ! to be a half-step advanced. Because we are operator splitting,
  ! one must use the updated dust and gas velocities.
  ! "Large dust fractions can prevent the propagation of soundwaves"
  ! Above is a paper that we should use to test our code at high mu
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer ::nn ! number of cells
  integer ::ivar_mhd_pic,xtondim  ! cell-centered dust variables start.
  real(dp) ::dt ! timestep.
  real(dp)::rd,cs! ERM: Grain size parameter
  real(dp),dimension(1:nvector) ::nu,c
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector) ::dgr,gs! gas density at grain, grain size array
  real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
  real(dp),dimension(1:nvector,1:xtondim,1:ndim)::big_v
  real(dp),dimension(1:nvector,1:ndim),save ::wh! drift at half step.
  logical,dimension(1:nvector)::dust
  integer ::i,j,idim,ind
  ivar_mhd_pic=9
  rd = 0.62665706865775 !*sqrt(gamma) constant for epstein drag law.

  if ((constant_t_stop).and.(stopping_rate .lt. 0.0))then ! add a "constant_nu_stop" option so you can turn drag totally off.
    nu(1:nvector)=(1./t_stop)*grain_size/max(gs(1:nvector),smallr) ! Or better yet, add pre-processor directives to turn drag off.
  else if ((constant_t_stop) .and. (stopping_rate .ge. 0.0))then
    nu(1:nvector)=stopping_rate*grain_size/max(gs(1:nvector),smallr)
  else

     wh(1:nn,1:ndim) = 0.0D0 ! Set to the drift velocity post-Lorentz force
     if(pic_dust.and.supersonic_drag)then
        do ind=1,xtondim
          do idim=1,ndim
            do j=1,nn
              if(dust(j))then
                 wh(j,idim)=wh(j,idim)+vol(j,ind)*&
                 &(v(j,idim)-uold(indp(j,ind),1+idim)/&
                 &max(uold(indp(j,ind),1),smallr))
              endif
           end do
         end do
        end do
     endif
     do i=1,nn
       if(dust(i))then
         nu(i)=(dgr(i)*c(i)/(rd*gs(i)))*sqrt(1.+&
         &0.22089323345553233*&
         &(wh(i,1)**2+wh(i,2)**2+wh(i,3)**2)&
         &/(c(i)*c(i)))
      endif
     end do
  endif

  do i=1,nn
    if(.not.dust(i))then
      nu(i)=0.0d0
    endif
  end do

end subroutine StoppingRate
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine DragKick(nn,dt,indp,dust,ok,vol,nu,big_v,big_w,v,xtondim) ! mp is actually mov
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer::nn
  integer ::ivar_mhd_pic, xtondim ! cell-centered dust variables start.  integer ::nn ! number of cells
  real(dp) ::dt ! timestep
  real(dp) ::vol_loc ! cloud volume
  real(dp),dimension(1:nvector) ::nu,mp
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector,1:ndim) ::v ! grain velocity
  real(dp),dimension(1:nvector,1:xtondim,1:ndim) ::big_v,big_w
  real(dp) ::den_dust,den_gas,mu,nuj,up,den_cr
  logical,dimension(1:nvector):: dust
  integer ::i,j,ind,idim! Just an index
  ivar_mhd_pic=9

  if (second_order .and. back_reaction)then
    do ind=1,xtondim
       do i=1,nn
         if(dust(i))then
            den_gas=uold(indp(i,ind),1)
            den_dust=uold(indp(i,ind),ivar_mhd_pic)
            den_cr=unew(indp(i,ind),ivar_mhd_pic+2)
            mu=(den_dust+den_cr)/max(den_gas,smallr)
            nuj=(1.+mu)*unew(indp(i,ind),ivar_mhd_pic)/max(den_cr+uold(indp(i,ind),ivar_mhd_pic),smallr)
            do idim=1,ndim
              ! w = &
              ! &uold(indp(i,ind),ivar_mhd_pic+idim)/max(uold(indp(i,ind),ivar_mhd_pic),smallr)&
              ! &-uold(indp(i,ind),1+idim)/max(uold(indp(i,ind),1),smallr)

              big_w(i,ind,idim)=big_w(i,ind,idim)&
              &/(1.+nuj*dt+0.5*nuj*nuj*dt*dt)

              up = -mu*big_w(i,ind,idim)*(1.+0.5*dt*nuj)/(1.+mu)+(uold(indp(i,ind),1+idim)&
              &+uold(indp(i,ind),ivar_mhd_pic+idim))/&
              &max(den_gas+den_dust,smallr)

              big_v(i,ind,idim)=(big_v(i,ind,idim)+dt*nu(i)*(1.+0.5*dt*nu(i))*up)&
              &/(1.+dt*nu(i)*(1.+0.5*dt*nu(i)))
              ! up = -mu*big_w(i,ind,idim)*(1.+0.5*dt*nuj)/(1.+mu)+(uold(indp(i,ind),1+idim)&
              ! &+uold(indp(i,ind),ivar_mhd_pic+idim))/&
              ! &max(uold(indp(i,ind),1)+uold(indp(i,ind),ivar_mhd_pic),smallr)
              !
              ! big_v(i,ind,idim)=(big_v(i,ind,idim)+dt*nu(i)*(1.+0.5*dt*nu(i))*up)&
              ! &/(1.+dt*nu(i)*(1.+0.5*dt*nu(i)))
            end do
            ! big_w corresponds directly to a change in the gas velocity.
          endif
       end do
    end do
  elseif(back_reaction)then
     do ind=1,xtondim
        do i=1,nn
          if(dust(i))then
             den_gas=uold(indp(i,ind),1)
             den_dust=uold(indp(i,ind),ivar_mhd_pic)
             mu=den_dust/max(den_gas,smallr)
             nuj=(1.+mu)*unew(indp(i,ind),ivar_mhd_pic)/max(uold(indp(i,ind),ivar_mhd_pic),smallr)
             do idim=1,ndim
               ! w = &
               ! &uold(indp(i,ind),ivar_mhd_pic+idim)/max(uold(indp(i,ind),ivar_mhd_pic),smallr)&
               ! &-uold(indp(i,ind),1+idim)/max(uold(indp(i,ind),1),smallr)

               big_w(i,ind,idim)=big_w(i,ind,idim)&
               &/(1.+nuj*dt)

               up = -mu*big_w(i,ind,idim)/(1.+mu)+(uold(indp(i,ind),1+idim)&
               &+uold(indp(i,ind),ivar_mhd_pic+idim))/&
               &max(uold(indp(i,ind),1)+uold(indp(i,ind),ivar_mhd_pic),smallr)

               big_v(i,ind,idim)=(big_v(i,ind,idim)+dt*nu(i)*up)/(1.+dt*nu(i))
             end do
             ! big_w corresponds directly to a change in the gas velocity.
           endif
        end do
     end do
   else ! want to write down the second order version too. Strang split.
     do idim=1,ndim
       do ind=1,xtondim
          do i=1,nn
            if(dust(i))then
               den_gas=uold(indp(i,ind),1)
               up = uold(indp(i,ind),1+idim)/max(den_gas,smallr)
               big_v(i,ind,idim)=(big_v(i,ind,idim)+dt*nu(i)*up)/(1.+dt*nu(i))
             endif
          end do
       end do
    end do
   endif
end subroutine DragKick

!#########################################################################
!#########################################################################
subroutine PicParams(sizes,charges,dust,cosr,lorentzf,ind_part,np)
  ! This subroutine will compute changes to sub-cloud velocity in big_v,
  ! as well as set unew's dust momentum slot to being u+du**EM.
  ! need to add "lorentzf" array for lorentz factors
  use amr_commons
  use pm_commons
  use hydro_parameters
  use poisson_commons
  use amr_parameters, ONLY: cr_c_fraction
  implicit none
  real(dp),dimension(1:nvector)::charges,sizes,lorentzf ! charge to mass ratios, lorentz factors
  integer,dimension(1:nvector)::ind_part
  logical,dimension(1:nvector) :: dust, cosr
  integer ::i,j,ind,idim,np,ivar_mhd_pic,chint,drint,gint,chbins,drbins! Just an -index


  ivar_mhd_pic=9

  if (ddex.ne.0.0 .and. .not. (universal_drag.or.universal_charge.or.size_bins).and. .not.lognormal .and. .not. (astrodust2.or.astrodust4)) then ! if there is a spectrum of grain sizes...
    do j=1,np ! construct charges and grain sizes
      if(dust(j))then
        sizes(j)=grain_size*1.0d1**((ddex*(idp(ind_part(j))-1.0d0))/(ndust*2.0d0**(3.*levelmin)-1.0d0))
        charges(j)=du_charge_to_mass*(sizes(j)/grain_size)**charge_slope
      elseif(cosr(j))then
        charges(j)=cr_c_fraction*cr_c_fraction*cr_charge_to_mass/lorentzf(j)
        ! Also need to divide by cr_c_fraction.
      endif
    end do
    elseif (astrodust2) then ! if we use the Hensley & Draine size distribution...
      do j=1,np ! construct charges and grain sizes
        if(dust(j))then
          sizes(j)=grain_size*2.0d0**(2.0d0*((idp(ind_part(j))-1.0d0)/(ndust*2.0d0**(3.*levelmin)-1.0d0)-0.5d0))
          charges(j)=du_charge_to_mass*(sizes(j)/grain_size)**charge_slope
        elseif(cosr(j))then
          charges(j)=cr_c_fraction*cr_c_fraction*cr_charge_to_mass/lorentzf(j)
          ! Also need to divide by cr_c_fraction.
        endif
      end do
      elseif (astrodust4) then ! if we use the Hensley & Draine size distribution...
         do j=1,np ! construct charges and grain sizes
           if(dust(j))then
             sizes(j)=grain_size*2.0d0**(4.0d0*((idp(ind_part(j))-1.0d0)/(ndust*2.0d0**(3.*levelmin)-1.0d0)-0.5d0))
             charges(j)=du_charge_to_mass*(sizes(j)/grain_size)**charge_slope
           elseif(cosr(j))then
             charges(j)=cr_c_fraction*cr_c_fraction*cr_charge_to_mass/lorentzf(j)
             ! Also need to divide by cr_c_fraction.
           endif
         end do
   elseif (ddex.ne.0.0 .and. lognormal)then ! if there is a spectrum of grain sizes...
      do j=1,np ! construct charges and grain sizes
        if(dust(j))then
          sizes(j)=grain_size*1.0d1**(ddex*((idp(ind_part(j))-1.0d0)/(ndust*2.0d0**(3.*levelmin)-1.0d0)-0.5d0))
          charges(j)=du_charge_to_mass*(sizes(j)/grain_size)**charge_slope
        elseif(cosr(j))then
          charges(j)=cr_c_fraction*cr_c_fraction*cr_charge_to_mass/lorentzf(j)
          ! Also need to divide by cr_c_fraction.
        endif
      end do
  elseif (ddex.ne.0.0 .and. universal_drag.and.size_bins)then ! if there is a spectrum of grain sizes...
    do j=1,np ! construct charges and grain sizes
      if(dust(j))then
        sizes(j)=grain_size*1.0d1**(ddex*int((idp(ind_part(j))-1.0d0)/(2.0d0**(3.*levelmin)))/ndust)
        charges(j)=du_charge_to_mass*(sizes(j)/grain_size)**charge_slope
        sizes(j)=grain_size
      elseif(cosr(j))then
        charges(j)=cr_c_fraction*cr_c_fraction*cr_charge_to_mass/lorentzf(j)
        ! Also need to divide by cr_c_fraction.
      endif
    end do
  elseif (ddex.ne.0.0 .and. universal_charge.and.size_bins)then ! if there is a spectrum of grain sizes...
    do j=1,np ! construct charges and grain sizes
      if(dust(j))then
        sizes(j)=grain_size*1.0d1**(ddex*int((idp(ind_part(j))-1.0d0)/(2.0d0**(3.*levelmin)))/ndust)
        charges(j)=du_charge_to_mass
      elseif(cosr(j))then
        charges(j)=cr_c_fraction*cr_c_fraction*cr_charge_to_mass/lorentzf(j)
        ! Also need to divide by cr_c_fraction.
      endif
    end do
  elseif (ddex.ne.0.0 .and. size_bins .and. grain_sampling_rate > 1)then ! if there is a spectrum of grain sizes...
    chbins = int(ndust/grain_sampling_rate)
    drbins=grain_sampling_rate
    do j=1,np ! construct charges and grain sizes
      if(dust(j))then
        gint=int((idp(ind_part(j))-1.0d0)/(2.0d0**(3.*levelmin)))
        ! gint goes from zero to ndust-1
        chint = int(gint/grain_sampling_rate)
        drint = mod(gint,grain_sampling_rate)
        sizes(j)=grain_size*1.0d1**(chdex*chint/chbins)
        charges(j)=du_charge_to_mass*(sizes(j)/grain_size)**charge_slope
        sizes(j)=grain_size*1.0d1**(ddex*drint/drbins)
      elseif(cosr(j))then
        charges(j)=cr_c_fraction*cr_c_fraction*cr_charge_to_mass/lorentzf(j)
        ! Also need to divide by cr_c_fraction.
      endif
    end do
  elseif (ddex.ne.0.0 .and. size_bins)then ! if there is a spectrum of grain sizes...
    do j=1,np ! construct charges and grain sizes
      if(dust(j))then
        sizes(j)=grain_size*1.0d1**(ddex*int((idp(ind_part(j))-1.0d0)/(2.0d0**(3.*levelmin)))/ndust)
        charges(j)=du_charge_to_mass*(sizes(j)/grain_size)**charge_slope
      elseif(cosr(j))then
        charges(j)=cr_c_fraction*cr_c_fraction*cr_charge_to_mass/lorentzf(j)
        ! Also need to divide by cr_c_fraction.
      endif
    end do
  else
    do j=1,np ! construct charges and grain sizes
      if(dust(j))then
        sizes(j)=grain_size
        charges(j)=du_charge_to_mass
      elseif(cosr(j))then
        charges(j)=cr_c_fraction*cr_c_fraction*cr_charge_to_mass/lorentzf(j)
      endif
    end do
  endif

end subroutine PicParams
!#########################################################################
!#########################################################################
