      SUBROUTINE ana_grid (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets model grid using an analytical expressions.       !
!                                                                      !
!  On Output:  stored in common blocks:                                !
!                                                                      !
!                           "grid"    (file grid.h)                    !
!                           "scalars" (file scalar.h)                  !
!                                                                      !
!     el       Length (m) of domain box in the ETA-direction.          !
!     f        Coriolis parameter (1/seconds) at RHO-points.           !
!     h        Bathymetry (meters; positive) at RHO-points.            !
!     hmin     Minimum depth of bathymetry (m).                        !
!     hmax     Maximum depth of bathymetry (m).                        !
!     pm       Coordinate transformation metric "m" (1/meters)         !
!              associated with the differential distances in XI        !
!              at RHO-points.                                          !
!     pn       Coordinate transformation metric "n" (1/meters)         !
!              associated with the differential distances in ETA.      !
!              at RHO-points.                                          !
!     xl       Length (m) of domain box in the XI-direction.           !
!     xp       XI-coordinates (m) at PSI-points.                       !
!     xr       XI-coordinates (m) at RHO-points.                       !
!     yp       ETA-coordinates (m) at PSI-points.                      !
!     yr       ETA-coordinates (m) at RHO-points.                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_grid_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    GRID(ng) % angler,                            &
#if defined CURVGRID && defined UV_ADV
     &                    GRID(ng) % dmde,                              &
     &                    GRID(ng) % dndx,                              &
#endif
#ifdef ICESHELF
     &                    GRID(ng) % zice,                              &
#endif
#ifdef SPHERICAL
     &                    GRID(ng) % lonp,                              &
     &                    GRID(ng) % lonr,                              &
     &                    GRID(ng) % lonu,                              &
     &                    GRID(ng) % lonv,                              &
     &                    GRID(ng) % latp,                              &
     &                    GRID(ng) % latr,                              &
     &                    GRID(ng) % latu,                              &
     &                    GRID(ng) % latv,                              &
#else
     &                    GRID(ng) % xp,                                &
     &                    GRID(ng) % xr,                                &
     &                    GRID(ng) % xu,                                &
     &                    GRID(ng) % xv,                                &
     &                    GRID(ng) % yp,                                &
     &                    GRID(ng) % yr,                                &
     &                    GRID(ng) % yu,                                &
     &                    GRID(ng) % yv,                                &
#endif
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % f,                                 &
     &                    GRID(ng) % h)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 7)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_grid
!
!***********************************************************************
      SUBROUTINE ana_grid_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          angler,                                 &
#if defined CURVGRID && defined UV_ADV
     &                          dmde, dndx,                             &
#endif
#ifdef ICESHELF
     &                          zice,                                   &
#endif
#ifdef SPHERICAL
     &                          lonp, lonr, lonu, lonv,                 &
     &                          latp, latr, latu, latv,                 &
#else
     &                          xp, xr, xu, xv,                         &
     &                          yp, yr, yu, yv,                         &
#endif
     &                          pn, pm, f, h)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
#endif
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: angler(LBi:,LBj:)
# if defined CURVGRID && defined UV_ADV
      real(r8), intent(out) :: dmde(LBi:,LBj:)
      real(r8), intent(out) :: dndx(LBi:,LBj:)
# endif
# ifdef ICESHELF
      real(r8), intent(out) :: zice(LBi:,LBj:)
# endif
# ifdef SPHERICAL
      real(r8), intent(out) :: lonp(LBi:,LBj:)
      real(r8), intent(out) :: lonr(LBi:,LBj:)
      real(r8), intent(out) :: lonu(LBi:,LBj:)
      real(r8), intent(out) :: lonv(LBi:,LBj:)
      real(r8), intent(out) :: latp(LBi:,LBj:)
      real(r8), intent(out) :: latr(LBi:,LBj:)
      real(r8), intent(out) :: latu(LBi:,LBj:)
      real(r8), intent(out) :: latv(LBi:,LBj:)
# else
      real(r8), intent(out) :: xp(LBi:,LBj:)
      real(r8), intent(out) :: xr(LBi:,LBj:)
      real(r8), intent(out) :: xu(LBi:,LBj:)
      real(r8), intent(out) :: xv(LBi:,LBj:)
      real(r8), intent(out) :: yp(LBi:,LBj:)
      real(r8), intent(out) :: yr(LBi:,LBj:)
      real(r8), intent(out) :: yu(LBi:,LBj:)
      real(r8), intent(out) :: yv(LBi:,LBj:)
# endif
      real(r8), intent(out) :: pn(LBi:,LBj:)
      real(r8), intent(out) :: pm(LBi:,LBj:)
      real(r8), intent(out) :: f(LBi:,LBj:)
      real(r8), intent(out) :: h(LBi:,LBj:)
#else
      real(r8), intent(out) :: angler(LBi:UBi,LBj:UBj)
# if defined CURVGRID && defined UV_ADV
      real(r8), intent(out) :: dmde(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: dndx(LBi:UBi,LBj:UBj)
# endif
# ifdef ICESHELF
      real(r8), intent(out) :: zice(LBi:UBi,LBj:UBj)
# endif
# ifdef SPHERICAL
      real(r8), intent(out) :: lonp(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: lonu(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: lonv(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: latp(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: latr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: latu(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: latv(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(out) :: xp(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: xu(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: xv(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: yp(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: yr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: yu(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: yv(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: h(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: NSUB, i, ival, j, k

      real(r8), parameter :: twopi = 2.0_r8*pi

      real(r8) :: Esize, Xsize, beta, cff, depth, dth
      real(r8) :: dx, dy, f0, my_min, my_max, r, theta, val1, val2

#ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len=3), dimension(2) :: op_handle
#endif
#ifdef WEDDELL
      real(r8) :: hwrk(-1:235), xwrk(-1:235), zwrk
#endif
      real(r8) :: wrkX(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: wrkY(IminS:ImaxS,JminS:JmaxS)

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set grid parameters:
!
!     Xsize    Length (m) of domain box in the XI-direction.
!     Esize    Length (m) of domain box in the ETA-direction.
!     depth    Maximum depth of bathymetry (m).
!     f0       Coriolis parameter, f-plane constant (1/s).
!     beta     Coriolis parameter, beta-plane constant (1/s/m).
!-----------------------------------------------------------------------
!
#if defined BASIN
      Xsize=3600.0E+03_r8
      Esize=2800.0E+03_r8
      depth=5000.0_r8
      f0=1.0E-04_r8
      beta=2.0E-11_r8
#elif defined BENCHMARK
      Xsize=360.0_r8              ! degrees of longitude
      Esize=20.0_r8               ! degrees of latitude
      depth=4000.0_r8
      f0=-1.0E-04_r8
      beta=2.0E-11_r8
#elif defined BL_TEST
      Xsize=100.0E+03_r8
      Esize=5.0E+03_r8
      depth=47.5_r8
      f0=9.25E-04_r8
      beta=0.0_r8
#elif defined CHANNEL
      Xsize=600.0E+03_r8
      Esize=360.0E+03_r8
      depth=500.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
#elif defined CANYON
      Xsize=128.0E+03_r8
      Esize=96.0E+03_r8
      depth=4000.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
#elif defined COUPLING_TEST
      Xsize=6000.0_r8*REAL(Lm(ng),r8)
      Esize=6000.0_r8*REAL(Mm(ng),r8)
      depth=1500.0_r8
      f0=5.0E-05_r8
      beta=0.0_r8
#elif defined DOUBLE_GYRE
      Xsize=1000.0E+03_r8
      Esize=2000.0E+03_r8
      depth=500.0_r8
!!    depth=5000.0_r8
      f0=7.3E-05_r8
      beta=2.0E-11_r8
#elif defined ESTUARY_TEST
      Xsize=100000.0_r8
      Esize=300.0_r8
      depth=10.0_r8
      f0=0.0_r8
      beta=0.0_r8
#elif defined KELVIN
      Xsize=20000.0_r8*REAL(Lm(ng),r8)
      Esize=20000.0_r8*REAL(Mm(ng),r8)
      depth=100.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
#elif defined FLT_TEST
      Xsize=1.0E+03_r8*REAL(Lm(ng),r8)
      Esize=1.0E+03_r8*REAL(Mm(ng),r8)
      depth=10.0_r8
      f0=0.0_r8
      beta=0.0_r8
#elif defined GRAV_ADJ
      Xsize=64.0E+03_r8
      Esize=2.0E+03_r8
      depth=20.0_r8
      f0=0.0_r8
      beta=0.0_r8
#elif defined LAB_CANYON
      Xsize=0.55_r8                  ! width of annulus
      Esize=2.0_r8*pi                ! azimuthal length (radians)
      f0=4.0_r8*pi/25.0_r8
      beta=0.0_r8
#elif defined LAKE_SIGNELL
      Xsize=50.0e3_r8
      Esize=10.0e3_r8
      depth=18.0_r8
      f0=0.0E-04_r8
      beta=0.0_r8
#elif defined LMD_TEST
      Xsize=100.0E+03_r8
      Esize=100.0E+03_r8
      depth=50.0_r8
      f0=1.09E-04_r8
      beta=0.0_r8
# elif defined MIXED_LAYER
      Xsize=500.0_r8
      Esize=400.0_r8
      depth=50.0_r8
      f0=0.0_r8
      beta=0.0_r8
#elif defined OVERFLOW
      Xsize=4.0E+03_r8
      Esize=200.0E+03_r8
      depth=4000.0_r8
      f0=0.0_r8
      beta=0.0_r8
#elif defined RIVERPLUME1
      Xsize=58.5E+03_r8
      Esize=201.0E+03_r8
      depth=150.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
#elif defined RIVERPLUME2
      Xsize=100.0E+03_r8
      Esize=210.0E+03_r8
      depth=190.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
#elif defined SEAMOUNT
      Xsize=320.0E+03_r8
      Esize=320.0E+03_r8
      depth=5000.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
#elif defined SOLITON
!!    Xsize=0.5_r8*REAL(Lm(ng),r8)
!!    Esize=0.5_r8*REAL(Mm(ng),r8)
      Xsize=48.0_r8
      Esize=16.0_r8
      depth=1.0_r8
      f0=0.0_r8
      beta=1.0_r8
      g=1.0_r8
#elif defined SED_TEST1
      Xsize=300.0_r8
      Esize=36.0_r8
      depth=10.0_r8
      f0=0.0_r8
      beta=0.0_r8
#elif defined SED_TOY
      Xsize=40.0_r8
      Esize=30.0_r8
      depth=0.5_r8
      f0=0.0_r8
      beta=0.0_r8
# elif defined SHOREFACE
      Xsize=1180.0_r8
      Esize=140.0_r8
      depth=15.0_r8
      f0=0.0E-04_r8
      beta=0.0_r8
#elif defined TEST_CHAN
      Xsize=10000.0_r8
      Esize=1000.0_r8
      depth=10.0_r8
      f0=0.0_r8
      beta=0.0_r8
#elif defined UPWELLING
      Xsize=1000.0_r8*REAL(Lm(ng),r8)
      Esize=1000.0_r8*REAL(Mm(ng),r8)
      depth=150.0_r8
      f0=-8.26E-05_r8
      beta=0.0_r8
#elif defined WEDDELL
      Xsize=4000.0_r8*REAL(Lm(ng),r8)
      Esize=4000.0_r8*REAL(Mm(ng),r8)
      depth=4500.0_r8
      f0=0.0_r8
      beta=0.0_r8
#elif defined WINDBASIN
      Xsize=2000.0_r8*REAL(Lm(ng),r8)
      Esize=1000.0_r8*REAL(Mm(ng),r8)
      depth=50.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
#elif defined ICETEST
      Xsize=570.466E+03_r8
      Esize=1120.0E+03_r8
      depth=900.0_r8
      f0=0.0_r8
      beta=0.0_r8
#elif defined ICESHELF2D
      Xsize=20.0E+03_r8
      Esize=500.0E+03_r8
      depth=980.0_r8
      f0=0.0_r8
      beta=0.0_r8
#elif defined MERTZ_TEST
      Xsize=350.0E+03_r8
      Esize=390.0E+03_r8
      depth=500.0_r8
      f0=0.0_r8
      beta=0.0_r8
#else
      ana_grid.h: no values provided for Xsize, Esize, depth, f0, beta.
#endif
!
!  Load grid parameters to global storage.
!
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        xl(ng)=Xsize
        el(ng)=Esize
      END IF
!
!-----------------------------------------------------------------------
!  Compute the (XI,ETA) coordinates at PSI- and RHO-points.
!  Set grid spacing (m).
!-----------------------------------------------------------------------
!
!  Determine I- and J-ranges for computing grid data.  These ranges
!  are special in periodic boundary conditons since periodicity cannot
!  be imposed in the grid coordinates.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=Istr-1
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=Iend+1
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=Jstr-1
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=Jend+1
      ELSE
        Jmax=Jend
      END IF

#if defined BENCHMARK
!
!  Spherical coordinates set-up.
!
      dx=Xsize/REAL(Lm(ng),r8)
      dy=Esize/REAL(Mm(ng),r8)
      spherical=.TRUE.
      DO j=Jmin,Jmax
        val1=-70.0_r8+dy*(REAL(j,r8)-0.5_r8)
        val2=-70.0_r8+dy*REAL(j,r8)
        DO i=Imin,Imax
          lonr(i,j)=dx*(REAL(i,r8)-0.5_r8)
          latr(i,j)=val1
          lonu(i,j)=dx*REAL(i,r8)
          lonp(i,j)=lonu(i,j)
          latu(i,j)=latr(i,j)
          lonv(i,j)=lonr(i,j)
          latv(i,j)=val2
          latp(i,j)=latv(i,j)
        END DO
      END DO
#elif defined LAB_CANYON
!
!  Polar coordinates set-up.
!
      dx=Xsize/REAL(Lm(ng),r8)
      dy=Esize/REAL(Mm(ng),r8)
!!    dth=twopi/REAL(Mm(ng),r8)               ! equal azimultal spacing
      dth=0.01_r8                             ! azimultal spacing
      cff=(4.0_r8*pi/(dth*REAL(Mm(ng),r8)))-1.0_r8   ! F
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          r=0.35_r8+dx*REAL(i-1,r8)
          theta=-pi+                                                    &
     &          0.5_r8*dth*((cff+1.0_r8)*REAL(j-1,r8)+                  &
     &                      (cff-1.0_r8)*(REAL(Mm(ng),r8)/twopi)*       &
     &                      SIN(twopi*REAL(j-1,r8)/REAL(Mm(ng),r8)))
          xp(i,j)=r*COS(theta)
          yp(i,j)=r*SIN(theta)
          r=0.35_r8+dx*(REAL(i-1,r8)+0.5_r8)
          theta=-pi+                                                    &
     &          0.5_r8*dth*((cff+1.0_r8)*(REAL(j-1,r8)+0.5_r8)+         &
     &                      (cff-1.0_r8)*(REAL(Mm(ng),r8)/twopi)*       &
     &                      SIN(twopi*(REAL(j-1,r8)+0.5_r8)/            &
     &                          REAL(Mm(ng),r8)))
          xr(i,j)=r*COS(theta)
          yr(i,j)=r*SIN(theta)
          xu(i,j)=xp(i,j)
          yu(i,j)=yr(i,j)
          xv(i,j)=xr(i,j)
          yv(i,j)=yp(i,j)
        END DO
      END DO
#elif defined ICETEST
!
!  Spherical coordinates set-up.
!
      spherical=.TRUE.
      DO j=JstrR,JendR
        cff=-80.0_r8+0.1_r8*REAL(j-1,r8)
        DO i=IstrR,IendR
          lonr(i,j)=0.3_r8*REAL(i-1,r8)
          latr(i,j)=cff
          lonu(i,j)=0.3_r8*REAL(i-1,r8)+0.15_r8
          lonp(i,j)=lonu(i,j)
          latu(i,j)=latr(i,j)
          lonv(i,j)=lonr(i,j)
          latv(i,j)=cff+0.05_r8
          latp(i,j)=latv(i,j)
        END DO
      END DO
#elif defined ICESHELF2D
      dx=Xsize/REAL(Lm(ng),r8)
      dy=Esize/REAL(Mm(ng),r8)
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          xp(i,j)=dx*REAL(i-1,r8)
          xr(i,j)=dx*(REAL(i-1,r8)+0.5_r8)
          xu(i,j)=xp(i,j)
          xv(i,j)=xr(i,j)
          yp(i,j)=dy*REAL(j-1,r8)
          yr(i,j)=dy*(REAL(j-1,r8)+0.5_r8)
          yu(i,j)=yr(i,j)
          yv(i,j)=yp(i,j)
        END DO
      END DO
#elif defined MERTZ_TEST
!
!  Spherical coordinates set-up.
!
      spherical=.TRUE.
      DO j=JstrR,JendR
        cff=-68.0_r8+0.1_r8*REAL(j-1,r8)
        DO i=IstrR,IendR
          lonr(i,j)=0.3_r8*REAL(i-1,r8)
          latr(i,j)=cff
          lonu(i,j)=0.3_r8*REAL(i-1,r8)+0.15_r8
          lonp(i,j)=lonu(i,j)
          latu(i,j)=latr(i,j)
          lonv(i,j)=lonr(i,j)
          latv(i,j)=cff+0.05_r8
          latp(i,j)=latv(i,j)
        END DO
      END DO
#else
      dx=Xsize/REAL(Lm(ng),r8)
      dy=Esize/REAL(Mm(ng),r8)
      DO j=Jmin,Jmax
        DO i=Imin,Imax
# ifdef BL_TEST
          dx=0.5_r8*(4000.0_r8/REAL(Lm(ng)+1,r8))*REAL(i,r8)+675.0_r8
# endif
          xp(i,j)=dx*REAL(i-1,r8)
          xr(i,j)=dx*(REAL(i-1,r8)+0.5_r8)
          xu(i,j)=xp(i,j)
          xv(i,j)=xr(i,j)
          yp(i,j)=dy*REAL(j-1,r8)
          yr(i,j)=dy*(REAL(j-1,r8)+0.5_r8)
          yu(i,j)=yr(i,j)
          yv(i,j)=yp(i,j)
        END DO
      END DO
#endif

#ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
# ifdef SPHERICAL
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    lonp, lonr, lonu, lonv)
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    latp, latr, latu, latv)
# else
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    xp, xr, xu, xv)
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    yp, yr, yu, yv)
# endif
#endif
!
!-----------------------------------------------------------------------
! Compute coordinate transformation metrics at RHO-points "pm" and
! "pn"  (1/m) associated with the differential distances in XI and
! ETA, respectively.
!-----------------------------------------------------------------------
!
#define J_RANGE MIN(JstrT,Jstr-1),MAX(Jend+1,JendT)
#define I_RANGE MIN(IstrT,Istr-1),MAX(Iend+1,IendT)

#if defined BENCHMARK
!
!  Spherical coordinates set-up.
!
      val1=REAL(Lm(ng),r8)/(2.0_r8*pi*Eradius)
      val2=REAL(Mm(ng),r8)*360.0_r8/(2.0_r8*pi*Eradius*Esize)
      DO j=J_RANGE
         cff=1.0_r8/COS((-70.0_r8+dy*(REAL(j,r8)-0.5_r8))*deg2rad)
        DO i=I_RANGE
          wrkX(i,j)=val1*cff
          wrkY(i,j)=val2
        END DO
      END DO
#elif defined LAB_CANYON
!
!  Polar coordinates set-up.
!
      DO j=J_RANGE
        DO i=I_RANGE
          r=0.35_r8+dx*(REAL(i-1,r8)+0.5_r8)
          theta=0.5_r8*dth*((cff+1.0_r8)+                               &
     &                      (cff-1.0_r8)*                               &
     &                      COS(twopi*REAL(j-1,r8)/REAL(Mm(ng),r8)))
          wrkX(i,j)=1.0_r8/dx
          wrkY(i,j)=1.0_r8/(r*theta)
        END DO
      END DO
# elif defined ICETEST
!
!  Spherical coordinates set-up.
!
      Eradius=6371020.0_r8
      val1=360.0_r8/(0.3_r8*(2.0_r8*pi*Eradius))
      val2=360.0_r8/(0.1_r8*(2.0_r8*pi*Eradius))
      DO j=J_RANGE
        DO i=I_RANGE
          wrkX(i,j)=val1/COS(latr(i,j)*deg2rad)
          wrkY(i,j)=val2
        END DO
      END DO
# elif defined ICESHELF2D
     DO j=J_RANGE
        DO i=I_RANGE
          wrkX(i,j)=1.0_r8/dx
          wrkY(i,j)=1.0_r8/dy
        END DO
      END DO
# elif defined MERTZ_TEST
!
!  Spherical coordinates set-up.
!
      Eradius=6371020.0_r8
      val1=360.0_r8/(0.3_r8*(2.0_r8*pi*Eradius))
      val2=360.0_r8/(0.1_r8*(2.0_r8*pi*Eradius))
      DO j=J_RANGE
        DO i=I_RANGE
          wrkX(i,j)=val1/COS(latr(i,j)*deg2rad)
          wrkY(i,j)=val2
        END DO
      END DO
#else
      DO j=J_RANGE
        DO i=I_RANGE
# ifdef BL_TEST
          dx=0.5_r8*(4000.0_r8/REAL(Lm(ng)+1,r8))*REAL(i,r8)+675.0_r8
# endif
          wrkX(i,j)=1.0_r8/dx
          wrkY(i,j)=1.0_r8/dy
        END DO
      END DO
#endif
#undef J_RANGE
#undef I_RANGE
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          pm(i,j)=wrkX(i,j)
          pn(i,j)=wrkY(i,j)
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pm)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pn)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pm, pn)
#endif

#if (defined CURVGRID && defined UV_ADV)
!
!-----------------------------------------------------------------------
!  Compute d(1/n)/d(xi) and d(1/m)/d(eta) at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          dndx(i,j)=0.5_r8*((1.0_r8/wrkY(i+1,j  ))-                     &
     &                      (1.0_r8/wrkY(i-1,j  )))
          dmde(i,j)=0.5_r8*((1.0_r8/wrkX(i  ,j+1))-                     &
     &                      (1.0_r8/wrkX(i  ,j-1)))
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          dndx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          dmde)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    dndx, dmde)
# endif
#endif
!
!-----------------------------------------------------------------------
! Angle (radians) between XI-axis and true EAST at RHO-points.
!-----------------------------------------------------------------------
!
#if defined LAB_CANYON
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          theta=-pi+                                                    &
     &          0.5_r8*dth*((cff+1.0_r8)*(REAL(j-1,r8)+0.5_r8)+         &
     &                      (cff-1.0_r8)*(REAL(Mm(ng),r8)/twopi)*       &
     &                      SIN(twopi*(REAL(j-1,r8)+0.5_r8)/            &
     &                          REAL(Mm(ng),r8)))
          angler(i,j)=theta
        END DO
      END DO
#elif defined WEDDELL
      val1=90.0_r8*deg2rad
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          angler(i,j)=val1
        END DO
      END DO
# elif defined ICETEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f(i,j)=(4.0_r8*pi/86164.1_r8)*                                &
     &           SIN(latr(i,j)*deg2rad)
        END DO
      END DO
# elif defined MERTZ_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f(i,j)=(4.0_r8*pi/86164.1_r8)*                                &
     &           SIN(latr(i,j)*deg2rad)
        END DO
      END DO
#else
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          angler(i,j)=0.0_r8
        END DO
      END DO
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          angler)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    angler)
#endif
!
!-----------------------------------------------------------------------
!  Compute Coriolis parameter (1/s) at RHO-points.
!-----------------------------------------------------------------------
!
#if defined BENCHMARK
      val1=2.0_r8*(2.0_r8*pi*366.25_r8/365.25_r8)/86400.0_r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          f(i,j)=val1*SIN(latr(i,j)*deg2rad)
        END DO
      END DO
#elif defined WEDDELL
      val1=10.4_r8/REAL(Lm(ng),r8)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          f(i,j)=2.0_r8*7.2E-05_r8*                                     &
     &           SIN((-79.0_r8+REAL(i-1,r8)*val1)*deg2rad)
        END DO
      END DO
# elif defined ICETEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f(i,j)=(4.0_r8*pi/86164.1_r8)*                                &
     &           SIN(latr(i,j)*deg2rad)
        END DO
      END DO
# elif defined MERTZ_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f(i,j)=(4.0_r8*pi/86164.1_r8)*                                &
     &           SIN(latr(i,j)*deg2rad)
        END DO
      END DO
#else
      val1=0.5_r8*Esize
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          f(i,j)=f0+beta*(yr(i,j)-val1)
        END DO
      END DO
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          f)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    f)
#endif
!
!-----------------------------------------------------------------------
!  Set bathymetry (meters; positive) at RHO-points.
!-----------------------------------------------------------------------
!
#if defined BENCHMARK
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=500.0_r8+1750.0_r8*(1.0+TANH((68.0_r8+latr(i,j))/dy))
        END DO
      END DO
#elif defined BL_TEST
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          val1=(xr(i,j)+500.0_r8)/15000.0_r8
          h(i,j)=14.0_r8+                                               &
     &           25.0_r8*(1.0_r8-EXP(-pi*xr(i,j)*1.0E-05_r8))-          &
     &           8.0_r8*EXP(-val1*val1)
        END DO
      END DO
#elif defined CANYON
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          val1=32000.0_r8-16000.0_r8*(SIN(pi*xr(i,j)/Xsize))**24
          h(i,j)=20.0_r8+0.5_r8*(depth-20.0_r8)*                        &
     &           (1.0_r8+TANH((yr(i,j)-val1)/10000.0_r8))
        END DO
      END DO
#elif defined ESTUARY_TEST
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=5.0_r8+(Xsize-xr(i,j))/Xsize*5.0_r8
        END DO
      END DO
#elif defined LAB_CANYON
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          r=0.35_r8+dx*(REAL(i-1,r8)+0.5_r8)
          theta=-pi+                                                    &
     &           0.5_r8*dth*((cff+1.0_r8)*(REAL(j-1,r8)+0.5_r8)+        &
     &                       (cff-1.0_r8)*(REAL(Mm(ng),r8)/twopi)*      &
     &                       SIN(dth*(REAL(j-1,r8)+0.5_r8)/             &
     &                           REAL(Mm(ng),r8)))
          val1=0.55_r8-0.15_r8*(COS(pi*theta*0.55_r8/0.2_r8)**2) !r_small
          val2=0.15_r8+0.15_r8*(COS(pi*theta*0.55_r8/0.2_r8)**2) !lambda
          IF (ABS(theta).ge.0.181818181818_r8) THEN
            IF (r.le.0.55_r8) THEN
              h(i,j)=0.025_r8                      ! shelf
            ELSE IF (r.ge.0.7_r8) THEN
              h(i,j)=0.125_r8                      ! deep
            ELSE
              h(i,j)=0.125_r8-0.1_r8*                                   &
     &               (COS(0.5_r8*pi*(r-0.55_r8)/0.15_r8)**2)
            END IF
          ELSE
            IF (r.le.val1) THEN
              h(i,j)=0.025_r8                      ! shelf
            ELSE IF (r.ge.0.7_r8) THEN
              h(i,j)=0.125_r8                      ! deep
            ELSE
              h(i,j)=0.125_r8-0.1_r8*                                   &
     &               (COS(0.5_r8*pi*(r-val1)/val2)**2)
            END IF
          END IF
        END DO
      END DO
#elif defined LAKE_SIGNELL
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=18.0_r8-16.0_r8*REAL(Mm(ng)-j,r8)/REAL(Mm(ng)-1,r8)
        END DO
      END DO
# elif defined MIXED_LAYER
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=50.0_r8
        END DO
      END DO
#elif defined OVERFLOW
      val1=200.0_r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=val1+0.5_r8*(depth-val1)*                              &
     &           (1.0_r8+TANH((yr(i,j)-100000.0_r8)/20000.0_r8))
        END DO
      END DO
#elif defined RIVERPLUME1
      DO j=JstrT,JendT
        DO i=IstrT,MIN(5,IendT)
          h(i,j)=15.0_r8
        END DO
        DO i=MAX(6,IstrT),IendT
          h(i,j)=depth+REAL(Lm(ng)-i,r8)*(15.0_r8-depth)/               &
     &                 REAL(Lm(ng)-6,r8)
        END DO
      END DO
#elif defined RIVERPLUME2
      DO j=JstrT,JendT
        DO i=IstrT,MIN(5,IendT)
          h(i,j)=15.0_r8
        END DO
        DO i=MAX(6,IstrT),IendT
          h(i,j)=depth+REAL(Lm(ng)-i,r8)*(15.0_r8-depth)/               &
     &                 REAL(Lm(ng)-6,r8)
        END DO
      END DO
#elif defined SEAMOUNT
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          val1=(xr(i,j)-0.5_r8*Xsize)/40000.0_r8
          val2=(yr(i,j)-0.5_r8*Esize)/40000.0_r8
          h(i,j)=depth-4500.0_r8*EXP(-(val1*val1+val2*val2))
        END DO
      END DO
#elif defined SED_TOY
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=20.0_r8
        END DO
      END DO
#elif defined SHOREFACE
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=11.75_r8-0.0125_r8*Xsize/REAL(Lm(ng)+1,r8)*REAL(i,r8)
        END DO
      END DO
#elif defined TEST_CHAN
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=10.0_r8+0.4040_r8*REAL(i,r8)/REAL(Lm(ng)+1,r8)
        END DO
      END DO
#elif defined UPWELLING
      IF (NSperiodic(ng)) THEN
        DO i=IstrT,IendT
          IF (i.le.Lm(ng)/2) THEN
            val1=REAL(i,r8)
          ELSE
            val1=REAL(Lm(ng)+1-i,r8)
          END IF
          val2=MIN(depth,84.5_r8+66.526_r8*TANH((val1-10.0_r8)/7.0_r8))
          DO j=JstrT,JendT
            h(i,j)=val2
          END DO
        END DO
      ELSE IF (EWperiodic(ng)) THEN
        DO j=JstrT,JendT
          IF (j.le.Mm(ng)/2) THEN
            val1=REAL(j,r8)
          ELSE
            val1=REAL(Mm(ng)+1-j,r8)
          END IF
          val2=MIN(depth,84.5_r8+66.526_r8*TANH((val1-10.0_r8)/7.0_r8))
          DO i=IstrT,IendT
            h(i,j)=val2
          END DO
        END DO
      END IF
#elif defined WEDDELL
      val1=98.80_r8
      val2=0.8270_r8
      DO k=-1,26
        xwrk(k)=REAL(k-1,r8)*15.0_r8*1000.0_r8
        hwrk(k)=375.0_r8
      END DO
      DO k=27,232
        zwrk=-2.0_r8+REAL(k-1,r8)*0.020_r8
        xwrk(k)=(520.0_r8+val1+zwrk*val1+                               &
     &           val1*val2*LOG(COSH(zwrk)))*1000.0_r8
        hwrk(k)=-75.0_r8+2198.0_r8*(1.0_r8+val2*TANH(zwrk))
      END DO
      DO k=233,235
        xwrk(k)=(850.0_r8+REAL(k-228,r8)*50.0_r8)*1000.0_r8
        hwrk(k)=4000.0_r8
      END DO
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=375.0_r8
          DO k=1,234
            IF ((xwrk(k).le.xr(i,1)).and.(xr(i,1).lt.xwrk(k+1))) THEN
               cff=1.0_r8/(xwrk(k+1)-xwrk(k))
               h(i,j)=cff*(xwrk(k+1)-xr(i,j))*hwrk(k  )+                &
     &                cff*(xr(i,j)-xwrk(k  ))*hwrk(k+1)
            END IF
          END DO
        END DO
      END DO
#elif defined WINDBASIN
      DO i=IstrT,IendT
        ival=INT(0.03_r8*REAL(Lm(ng)+1,r8))
        IF (i.lt.ival) THEN
          val1=1.0_r8-(REAL((i+1)-ival,r8)/REAL(ival,r8))**2
        ELSE IF ((Lm(ng)+1-i).lt.ival) THEN
          val1=1.0_r8-(REAL((Lm(ng)+1-i)-ival,r8)/REAL(ival,r8))**2
        ELSE
          val1=1.0_r8
        END IF
        DO j=JstrT,JendT
         val2=2.0_r8*REAL(j-(Mm(ng)+1)/2,r8)/REAL(Mm(ng)+1,r8)
         h(i,j)=depth*(0.08_r8+0.92_r8*val1*(1.0_r8-val2*val2))
        END DO
      END DO
# elif defined ICETEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          h(i,j)=900.0_r8
        END DO
      END DO
# elif defined ICESHELF2D
      DO j=JstrR,JendR
        DO i=IstrR,IendR
         h(i,j)=20.0_r8+REAL(j,r8)*(depth/Esize)*(Esize/REAL(Mm(ng),r8))
       write(6,*) h(i,j)
        END DO
      END DO
# elif defined MERTZ_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          h(i,j)=500.0_r8
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (j.le.22) THEN
            h(i,j)=500_r8
          ELSE IF (j.le.27) THEN
            h(i,j)=500.0_r8+(2000.0_r8/5.0_r8)*REAL(j-22,r8)
          ELSE
            h(i,j)=2500.0_r8
          END IF
        END DO
      END DO
#else
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=depth
        END DO
      END DO
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          h)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    h)
#endif
!
! Determine minimum depth: first, determine minimum values of depth
! within each subdomain, then determine global minimum by comparing
! these subdomain minima.
!
      my_min=h(IstrT,JstrT)
      my_max=h(IstrT,JstrT)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          my_min=MIN(my_min,h(i,j))
          my_max=MAX(my_max,h(i,j))
        END DO
      END DO
#ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
#endif
!$OMP CRITICAL (H_RANGE)
      IF (tile_count.eq.0) THEN
        hmin(ng)=my_min
        hmax(ng)=my_max
      ELSE
        hmin(ng)=MIN(hmin(ng),my_min)
        hmax(ng)=MAX(hmax(ng),my_max)
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#ifdef DISTRIBUTE
        buffer(1)=hmin(ng)
        buffer(2)=hmax(ng)
        op_handle(1)='MIN'
        op_handle(2)='MAX'
        CALL mp_reduce (ng, model, 2, buffer, op_handle)
        hmin(ng)=buffer(1)
        hmax(ng)=buffer(2)
#endif
      END IF
!$OMP END CRITICAL (H_RANGE)
#ifdef ICESHELF
!
!-----------------------------------------------------------------------
!  Set depth of ice shelf (meters; negative) at RHO-points.
!-----------------------------------------------------------------------
!
# ifdef WEDDELL
      val1=340.0_r8
      val2=val1/16.0_r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (i.gt.20) THEN
            zice(i,j)=0.0_r8
          ELSE IF (i.gt.4) THEN
            zice(i,j)=-val1+REAL(i-1,r8)*val2
          ELSE
            zice(i,j)=-val1
          END IF
        END DO
      END DO
#  elif defined ICETEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (j.eq.0) THEN
            zice(i,j)=-700.0_r8
          ELSE IF (j.le.41) THEN
            zice(i,j)=-700.0_r8+(500.0_r8/40.0_r8)*REAL(j-1,r8)
          ELSE
            zice(i,j)=0.0_r8
          END IF
!          IF (j.le.10) THEN
!            zice(i,j)=-700.0_r8
!          ELSE IF (j.lt.41) THEN
!            zice(i,j)=-700.0_r8+(500.0_r8/30.0_r8)*REAL(j-11,r8)
!          ELSE
!            zice(i,j)=-200.0_r8
!          END IF
        END DO
      END DO
#   elif defined ICESHELF2D
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (j.le.60) THEN
            zice(i,j)=-h(i,j)+20_r8
          ELSE 
         zice(i,j)=-(h(i,60)                                           &
     &             -atan(REAL(j-59,r8)/10)*(h(i,60)-300_r8))           &
     &             + 20_r8    
!         write(6,*) zice(i,j) 
          END IF
        END DO
      END DO
#  elif defined MERTZ_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (j.le.16) THEN
            zice(i,j)=-300.0_r8
!          ELSE IF (j.le.41) THEN
!            zice(i,j)=-700.0_r8+(500.0_r8/40.0_r8)*REAL(j-1,r8)
          ELSE
            zice(i,j)=0.0_r8
          END IF
        END DO
      END DO
# else
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          zice(i,j)=0.0_r8
        END DO
      END DO
# endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zice)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zice)
# endif
#endif

      RETURN
      END SUBROUTINE ana_grid_tile
