subroutine synchro_fine(ilevel)
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,xtondim
  !--------------------------------------------------------------------
  ! This routine synchronizes particle velocity with particle
  ! position for ilevel particle only. If particle sits entirely
  ! in level ilevel, then use inverse CIC at fine level to compute
  ! the force. Otherwise, use coarse level force and coarse level CIC.
  !--------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart
  integer::ig,ip,npart1,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  if(sink)then
     fsink_new=0
  endif

#ifdef TSC
    xtondim=threetondim
#else
    xtondim=twotondim
#endif

  ! Synchronize velocity using CIC
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
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ip=ip+1
           ind_part(ip)=ipart
           ind_grid_part(ip)=ig
           if(ip==nvector)then
              call sync(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)
              ip=0
              ig=0
           end if
           ipart=nextp(ipart)  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0)call sync(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)

  !sink cloud particles are used to average the grav. acceleration
  if(sink)then
     if(nsink>0)then
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(fsink_new,fsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
        fsink_all=fsink_new
#endif
     endif
     do isink=1,nsink
        if (.not. direct_force_sink(isink))then
           fsink_partial(isink,1:ndim,ilevel)=fsink_all(isink,1:ndim)
        end if
     end do
  endif

111 format('   Entering synchro_fine for level ',I2)

end subroutine synchro_fine
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine synchro_fine_static(ilevel)
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,xtondim
  !--------------------------------------------------------------------
  ! This routine synchronizes particle velocity with particle
  ! position for ilevel particle only. If particle sits entirely
  ! in level ilevel, then use inverse CIC at fine level to compute
  ! the force. Otherwise, use coarse level force and coarse level CIC.
  !--------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart
  integer::ig,ip,next_part,npart1,npart2,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  if(sink)then
     fsink_new=0
     fsink_all=0
  endif

#ifdef TSC
    xtondim=threetondim
#else
    xtondim=twotondim
#endif

  ! Synchronize velocity using CIC
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
                   & (.not. static_stars .and. (is_star(typep(ipart)) .or. is_debris(typep(ipart))) )  ) then
                 ! FIXME: there should be a static_sink as well
                 npart2=npart2+1
              endif
           else
              if(.not.static_dm) then
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
                   & (.not. static_stars .and. (is_star(typep(ipart)) .or. is_debris(typep(ipart))) )  ) then
                 ! FIXME: what about sinks?
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
              call sync(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)
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
  if(ip>0)call sync(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)

  !sink cloud particles are used to average the grav. acceleration
  if(sink)then
     if(nsink>0)then
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(fsink_new,fsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
        fsink_all=fsink_new
#endif
     endif
     do isink=1,nsink
        if (.not. direct_force_sink(isink))then
           fsink_partial(isink,1:ndim,ilevel)=fsink_all(isink,1:ndim)
        end if
     end do
  endif

111 format('   Entering synchro_fine for level ',I2)

end subroutine synchro_fine_static
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine sync(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,xtondim)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ng,np,ilevel,xtondim
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !
  !
  !
  logical::error
  integer::i,j,ind,idim,nx_loc,isink
  real(dp)::dx,scale
  ! Grid-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::dteff
  real(dp),dimension(1:3)::skip_loc
  ! Particle-based arrays
#ifndef TSC
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_vp,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
#else
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_vp
  real(dp),dimension(1:nvector,1:ndim),save::cl,cr,cc,wl,wr,wc
  integer ,dimension(1:nvector,1:ndim),save::igl,igr,igc,icl,icr,icc
  real(dp),dimension(1:nvector,1:threetondim),save::vol
  integer ,dimension(1:nvector,1:threetondim),save::igrid,icell,indp,kg
#endif


  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
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
#ifdef TSC
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

#else

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in sync'
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

!#else
!#include "tsc_fine.F90"
#endif
  ! Gather 3-force
  ff(1:np,1:ndim)=0.0D0
  if(poisson)then
     do ind=1,xtondim
        do idim=1,ndim
           do j=1,np
              ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*vol(j,ind)
           end do
        end do
     end do
  endif


    ! Extra acceleration forces on dust added here
    if(pic_dust.and.((accel_gr(1).ne.0).or.(accel_gr(2).ne.0).or.(accel_gr(3).ne.0)))then
       do idim=1,ndim
          do j=1,np
             if(is_dust(typep(ind_part(j))))then
             ff(j,idim)=ff(j,idim)+accel_gr(idim)
             endif
          end do
       end do
    endif


  ! For sink particle only, store contribution to the sink force
  if(sink)then
     do idim=1,ndim
        do j=1,np
           if ( is_cloud(typep(ind_part(j))) ) then
              isink=-idp(ind_part(j))
              if(.not. direct_force_sink(isink))then
                 fsink_new(isink,idim)=fsink_new(isink,idim)+ff(j,idim)
              endif
           endif
        end do
     end do
  end if

  ! Compute individual time steps
  do j=1,np
     if(levelp(ind_part(j))>=ilevel)then
        dteff(j)=dtnew(levelp(ind_part(j)))
     else
        dteff(j)=dtold(levelp(ind_part(j)))
     endif
  end do

  ! Update particles level
  do j=1,np
     levelp(ind_part(j))=ilevel
  end do

  ! Update 3-velocity
  do idim=1,ndim
     if(static)then
        do j=1,np
           new_vp(j,idim)=ff(j,idim)
        end do
     else
        do j=1,np
           new_vp(j,idim)=vp(ind_part(j),idim)+ff(j,idim)*0.5D0*dteff(j)
        end do
     endif
  end do
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim)=new_vp(j,idim)
     end do
  end do

  ! For sink particle only, overwrite cloud particle velocity with sink velocity
  if(sink)then
     do idim=1,ndim
        do j=1,np
           if ( is_cloud(typep(ind_part(j))) ) then
              isink=-idp(ind_part(j))
              ! Remember that vsink is half time step older than other particles
              vp(ind_part(j),idim)=vsink(isink,idim)
           endif
        end do
     end do
  end if

end subroutine sync
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_mhd_pic_fine(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,xtondim
  ! First, reset uold to zero.
  ! Can remove gravity and sink particle related things.
  ! Can remove synchro_fine_static as well.
  ! In "sync", want to remove the gravity...
  ! Syncing up the velocity, get rid of that too.
  integer::igrid,jgrid,ipart,jpart
  integer::ig,ip,npart1,isink
  integer::i,iskip,icpu,ind,ibound,ivar,ivar_mhd_pic
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ivar_mhd_pic=9

#ifdef TSC
    xtondim=threetondim
#else
    xtondim=twotondim
#endif

  ! Reset unew to zero for dust stopping rate and charge
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=ivar_mhd_pic,ivar_mhd_pic+3
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar)=0.0D0
        end do
     end do
  end do
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=ivar_mhd_pic,ivar_mhd_pic+3
           do i=1,reception(icpu,ilevel)%ngrid
              unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0D0
           end do
        end do
     end do
  end do

  ! Reset uold to zero for dust mass and momentum densities
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=ivar_mhd_pic,ivar_mhd_pic+3
           do i=1,reception(icpu,ilevel)%ngrid
              uold(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0D0
           end do
        end do
     end do
  end do
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=ivar_mhd_pic,ivar_mhd_pic+3
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar)=0.0D0
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!!!!!! NEW !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Reset dust variables in physical boundaries
   do ibound=1,nboundary
      do ind=1,twotondim
         iskip=ncoarse+(ind-1)*ngridmax
         do ivar=ivar_mhd_pic,ivar_mhd_pic+ndim
            do i=1,boundary(ibound,ilevel)%ngrid
               uold(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)=0.0D0
               ! unew(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)=&
               ! &uold(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)
            end do
         end do
      end do
   end do

   do ibound=1,nboundary
      do ind=1,twotondim
         iskip=ncoarse+(ind-1)*ngridmax
         do ivar=ivar_mhd_pic,ivar_mhd_pic+ndim
            do i=1,boundary(ibound,ilevel)%ngrid
               unew(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)=0.0D0
               ! unew(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)=&
               ! &uold(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)
            end do
         end do
      end do
   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Synchronize velocity using CIC (No longer need velocity to be synced.)
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
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ip=ip+1
           ind_part(ip)=ipart
           ind_grid_part(ip)=ig
           if(ip==nvector)then
              call init_mhd_pic(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)
              ip=0
              ig=0
           end if
           ipart=nextp(ipart)  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0)call init_mhd_pic(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)

  ! Update MPI boundary conditions for uold for dust mass and momentum densities
  ! Update MPI boundary conditions for unew for dust "mass" and "momentum" densities

  do ivar=ivar_mhd_pic,ivar_mhd_pic+ndim
     call make_virtual_reverse_dp(uold(1,ivar),ilevel)
     call make_virtual_fine_dp   (uold(1,ivar),ilevel)
     call make_virtual_reverse_dp(unew(1,ivar),ilevel)
     call make_virtual_fine_dp   (unew(1,ivar),ilevel)
  end do

  ! ! Update MPI boundary conditions for unew for dust "mass" and "momentum" densities
  ! call make_virtual_reverse_dp(unew(1,ivar_mhd_pic),ilevel)
  ! call make_virtual_fine_dp   (unew(1,ivar_mhd_pic),ilevel)



111 format('   Entering init_mhd_pic_fine for level ',I2)

end subroutine init_mhd_pic_fine
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine init_mhd_pic(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,xtondim)
  use amr_commons
  !use amr_parameters ERM
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  use amr_parameters, ONLY: cr_c_fraction
  implicit none
  integer::ng,np,ilevel,xtondim
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !
  !
  !
  logical::error
  integer::i,j,ind,idim,nx_loc,isink,ivar_mhd_pic,err
  real(dp)::dx,scale,dx_loc,vol_loc
  real(dp)::ctm, ts, rd, crsol

  ! Grid-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mmm,dteff,nu_stop,lorentzf
  real(dp),dimension(1:nvector),save::dgr,tss,mm,sizes,charges ! ERM: density, (non-constant) stopping times
  real(dp),dimension(1:nvector,1:ndim),save::uu,bb,vv ! ERM: Added these arrays
#ifndef TSC
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_vp,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
#else
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_vp,cl,cr,cc,wl,wr,wc
  integer ,dimension(1:nvector,1:ndim),save::igl,igr,igc,icl,icr,icc
  real(dp),dimension(1:nvector,1:threetondim),save::vol
  integer ,dimension(1:nvector,1:threetondim),save::igrid,icell,indp,kg
#endif
  real(dp),dimension(1:3)::skip_loc
  logical,dimension(1:nvector),save:: dust,cosr,mhd_pic

  ctm = du_charge_to_mass
  ts = t_stop
  rd = 0.62665706865775*grain_size ! constant for epstein drag law. used to havesqrt(gamma)*
  crsol=cr_c_fraction*2.9979246d+10*units_time/units_length ! Reduced speed of light.

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

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
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
  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in init_mhd_pic'
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

  do idim=1,ndim ! set vv equal to the velocity.
     do j=1,np
        vv(j,idim)=vp(ind_part(j),idim)
     end do
  end do

  ! Update old dust mass and momentum density variables
  ivar_mhd_pic=9
  if(nvar<ivar_mhd_pic+ndim .and. back_reaction)then
     write(*,*)'You need to compile ramses with nvar=',ivar_mhd_pic+ndim
     stop
  endif


  do j = 1, np
     dust(j) = is_dust(typep(ind_part(j)))
     cosr(j) = is_cr(typep(ind_part(j)))
     mhd_pic(j) = (cosr(j).or.dust(j))
  end do

  do j=1,np
    if(cosr(j))then
      lorentzf(j)=sqrt(1.0d0+&
      &(vv(j,1)*vv(j,1)+vv(j,2)*vv(j,2)+vv(j,3)*vv(j,3))/(crsol*crsol))
    else
      lorentzf(j)=1.0d0
    endif
  end do

  do ind=1,xtondim
     do j=1,np ! deposit the dust mass density.
        if(ok(j).and. mhd_pic(j))then !changed to total density
           uold(indp(j,ind),ivar_mhd_pic)=uold(indp(j,ind),ivar_mhd_pic)&
           &+mp(ind_part(j))*vol(j,ind)/vol_loc
        end if ! ERM: else if cosmic ray... (need to do things a bit differently to account for gammas etc.)
     end do
     do idim=1,ndim
        do j=1,np ! deposit the dust/cr momentum density
           if(ok(j) .and. mhd_pic(j))then
              uold(indp(j,ind),ivar_mhd_pic+idim)=uold(indp(j,ind),ivar_mhd_pic+idim)&
              &+mp(ind_part(j))*vv(j,idim)&
              &*vol(j,ind)/vol_loc
           end if
        end do
     end do
  end do

  sizes(1:np)=0.0d0
  charges(1:np)=0.0d0
  call InitPicParams(sizes,charges,dust,cosr,lorentzf,ind_part,np)
  ! if (ddex.ne.0.0)then ! if there is a spectrum of grain sizes...
  !   do j=1,np ! construct charges and grain sizes
  !     if(dust(j))then
  !       sizes(j)=grain_size*1.0d1**((ddex*(idp(ind_part(j))-1.0d0))/(ndust*2.0d0**(3.*levelmin)-1.0d0))
  !       charges(j)=du_charge_to_mass*(sizes(j)/grain_size)**charge_slope
  !     elseif(cosr(j))then
  !       charges(j)=cr_c_fraction*cr_charge_to_mass/lorentzf(j)
  !     endif
  !   end do
  ! else
  !   do j=1,np ! construct charges and grain sizes
  !     if(dust(j))then
  !       sizes(j)=grain_size
  !       charges(j)=du_charge_to_mass
  !     elseif(cosr(j))then
  !       charges(j)=cr_c_fraction*cr_charge_to_mass/lorentzf(j) ! Not accurate fully... need to divide by gamma.
  !     endif
  !   end do
  ! endif

  ! I don't think we actually want to do this until after the Lorentz kick
  call InitStoppingRate(np,dtnew(ilevel),indp,dust,ind_part,vol,vv,nu_stop,sizes,xtondim)
  !call InitCharge(np,dtnew(ilevel),indp,vol,vv,charges,xtondim)
  ! Is it ok for me to use other "dust" values of unew?
  ! Will have to deposit the grain charges too
  do ind=1,xtondim
     do j=1,np ! deposit the pic mass weighted stopping rate, charge, and CR mass
        if(ok(j) .and. mhd_pic(j))then
           unew(indp(j,ind),ivar_mhd_pic)=unew(indp(j,ind),ivar_mhd_pic)+&
           &(mp(ind_part(j))*vol(j,ind)/vol_loc)*&!rho^d_ij
           &nu_stop(j)
           unew(indp(j,ind),ivar_mhd_pic+1)=unew(indp(j,ind),ivar_mhd_pic+1)+&
           &(mp(ind_part(j))*vol(j,ind)/vol_loc)*&!rho^d_ij
           &charges(j)
        end if
        ! Need another statement to compute the CR mass (it's got a factor of 1/cr_c_fraction)
        ! if (ok(j).and.cosr(j))then
        !   unew(indp(j,ind),ivar_mhd_pic+2)=unew(indp(j,ind),ivar_mhd_pic+2)+&
        !   &(mp(ind_part(j))*vol(j,ind)/vol_loc)
        ! endif
     end do
  end do

end subroutine init_mhd_pic
!#########################################################################
!#########################################################################
subroutine InitStoppingRate(nn,dt,indp,dust,ind_part,vol,v,nu,gs,xtondim)
  ! The following subroutine will alter its last argument, nu
  ! to be a half-step advanced. Because we are operator splitting,
  ! one must use the updated dust and gas velocities.
  ! "Large dust fractions can prevent the propagation of soundwaves"
  ! Above is a paper that we should use to test our code at high mu
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer ::nn,xtondim ! number of cells
  integer ::ivar_mhd_pic ! cell-centered dust variables start.
  real(dp) ::dt ! timestep.
  real(dp)::rd,cs! ERM: Grain size parameter
  real(dp),dimension(1:nvector) ::nu,c,gs
  integer,dimension(1:nvector) :: ind_part
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector),save ::dgr! gas density at grain.
  real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
  real(dp),dimension(1:nvector,1:xtondim,1:ndim)::big_v
  real(dp),dimension(1:nvector,1:ndim),save ::w! drift at half step.
  logical,dimension(1:nvector)::dust
  integer ::i,j,idim,ind
  ivar_mhd_pic=9
  rd = 0.62665706865775 !constant for epstein drag law. #used to have *sqrt(gamma)
   ! isothermal sound speed... Need to get this right. This works for now,
         ! but only if you have scaled things so that the sound speed is 1.

     if ((constant_t_stop).and.(stopping_rate .lt. 0.0))then ! add a "constant_nu_stop" option so you can turn drag totally off.
       nu(1:nvector)=(1./t_stop)*grain_size/max(gs(1:nvector),smallr) ! Or better yet, add pre-processor directives to turn drag off.
     else if ((constant_t_stop) .and. (stopping_rate .ge. 0.0))then
       nu(1:nvector)=stopping_rate*grain_size/max(gs(1:nvector),smallr)
     else
     dgr(1:nn) = 0.0D0
     c(1:nn) = 0.0D0
     if(pic_dust)then
        do ind=1,xtondim
            do j=1,nn
              if(dust(j))then
               dgr(j)=dgr(j)+uold(indp(j,ind),1)*vol(j,ind)
               cs= uold(indp(j,ind),5) - &
               &0.125D0*(uold(indp(j,ind),1+5)+uold(indp(j,ind),1+nvar))**2 &
               &-0.125D0*(uold(indp(j,ind),2+5)+uold(indp(j,ind),2+nvar))**2&
               &-0.125D0*(uold(indp(j,ind),3+5)+uold(indp(j,ind),3+nvar))**2 &
               &- 0.5D0*(&
               & uold(indp(j,ind),1+1)**2+ &
               & uold(indp(j,ind),2+1)**2+ &
               & uold(indp(j,ind),3+1)**2)/max(uold(indp(j,ind),1),smallr) ! energy minus magnetic and kinetic.

               cs = cs*(gamma-1.0D0)/max(uold(indp(j,ind),1),smallr)
               cs = min(1/smallr,max(cs,smallr))
               c(j)=c(j)+vol(j,ind)*sqrt(cs)
              endif
           end do
        end do
     endif

     w(1:nn,1:ndim) = 0.0D0 ! Set to the drift velocity post-Lorentz force
     if(pic_dust .and. supersonic_drag)then
        do ind=1,xtondim
          do idim=1,ndim
            do j=1,nn
              if(dust(j))then
                 w(j,idim)=w(j,idim)+vol(j,ind)*&
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
         &(w(i,1)**2+w(i,2)**2+w(i,3)**2)&
         &/(c(i)*c(i)))
      endif
     end do
  endif

  do i=1,nn
    if(.not.dust(i))then
      nu(i)=0.0d0
    endif
  end do
end subroutine InitStoppingRate
!#########################################################################
!#########################################################################
subroutine InitPicParams(sizes,charges,dust,cosr,lorentzf,ind_part,np)
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
  integer ::i,j,ind,idim,np,ivar_mhd_pic,drint,chint,gint,chbins,drbins! Just an -index


  ivar_mhd_pic=9

  if (ddex.ne.0.0 .and. .not. (universal_drag.or.universal_charge.or.size_bins).and. .not. lognormal .and. .not. (astrodust2.or.astrodust4))then ! if there is a spectrum of grain sizes...
    do j=1,np ! construct charges and grain sizes
      if(dust(j))then
        sizes(j)=grain_size*1.0d1**((ddex*(idp(ind_part(j))-1.0d0))/(ndust*2.0d0**(3.*levelmin)-1.0d0))
        charges(j)=du_charge_to_mass*(sizes(j)/grain_size)**charge_slope
      elseif(cosr(j))then
        charges(j)=cr_c_fraction*cr_charge_to_mass/lorentzf(j)
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
             sizes(j)=grain_size*2.0d0**(2.0d0*((idp(ind_part(j))-1.0d0)/(ndust*2.0d0**(3.*levelmin)-1.0d0)-0.5d0))
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
        charges(j)=cr_c_fraction*cr_charge_to_mass/lorentzf(j)
      endif
    end do
  elseif (ddex.ne.0.0 .and. universal_charge.and.size_bins)then ! if there is a spectrum of grain sizes...
    do j=1,np ! construct charges and grain sizes
      if(dust(j))then
        sizes(j)=grain_size*1.0d1**(ddex*int((idp(ind_part(j))-1.0d0)/(2.0d0**(3.*levelmin)))/ndust)
        charges(j)=du_charge_to_mass
      elseif(cosr(j))then
        charges(j)=cr_c_fraction*cr_charge_to_mass/lorentzf(j)
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
        charges(j)=cr_c_fraction*cr_charge_to_mass/lorentzf(j)
      endif
    end do
  else
    do j=1,np ! construct charges and grain sizes
      if(dust(j))then
        sizes(j)=grain_size
        charges(j)=du_charge_to_mass
      elseif(cosr(j))then
        charges(j)=cr_c_fraction*cr_charge_to_mass/lorentzf(j)
      endif
    end do
  endif

end subroutine InitPicParams
!#########################################################################
!#########################################################################
