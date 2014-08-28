#include "cppdefs.h"
      MODULE step3d_t_mod
#if !defined TS_FIXED && (defined NONLINEAR && defined SOLVE3D)
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine time-steps tracer equations.  Notice that advective    !
!  and diffusive terms are time-stepped differently. It applies the    !
!  corrector time-step for horizontal/vertical advection,  vertical    !
!  diffusion, nudging if necessary, and lateral boundary conditions.   !
!                                                                      !
!  Notice that at input the tracer arrays have:                        !
!                                                                      !
!    t(:,:,:,nnew,:)   m Tunits  n+1     horizontal/vertical diffusion !
!                                        terms plus source/sink terms  !
!                                        (biology, sediment), if any   !
!                                                                      !
!    t(:,:,:,3   ,:)   Tunits    n+1/2   advective terms and vertical  !
!                                        diffusion predictor step      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: step3d_t
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_t (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef CLIMATOLOGY
      USE mod_clima
# endif
# ifdef DIAGNOSTICS_TS
      USE mod_diags
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
# ifdef RM_SUPERCOOLING
      USE mod_forces
# endif
# ifdef TS_PSOURCE
      USE mod_sources
# endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 35)
# endif
      CALL step3d_t_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
# if defined TS_PSOURCE || defined Q_PSOURCE
     &                    Msrc(ng), Nsrc(ng),                           &
     &                    SOURCES(ng) % Isrc,                           &
     &                    SOURCES(ng) % Jsrc,                           &
     &                    SOURCES(ng) % Tsrc,                           &
# endif
# ifdef TS_PSOURCE
     &                    SOURCES(ng) % Dsrc,                           &
# endif
# ifdef Q_PSOURCE
     &                    SOURCES(ng) % Qsrc,                           &
# endif
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
# endif
# ifdef TS_MPDATA
#  ifdef WET_DRY
     &                    GRID(ng) % rmask_wet,                         &
     &                    GRID(ng) % umask_wet,                         &
     &                    GRID(ng) % vmask_wet,                         &
#  endif
     &                    GRID(ng) % omn,                               &
     &                    GRID(ng) % om_u,                              &
     &                    GRID(ng) % om_v,                              &
     &                    GRID(ng) % on_u,                              &
     &                    GRID(ng) % on_v,                              &
# endif
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % z_r,                               &
# if defined TCLM_NUDGING && defined TCLIMATOLOGY
     &                    CLIMA(ng) % Tnudgcof,                         &
     &                    CLIMA(ng) % tclm,                             &
# endif
     &                    MIXING(ng) % Akt,                             &
     &                    OCEAN(ng) % W,                                &
# if defined FLOATS && defined FLOAT_VWALK
     &                    MIXING(ng) % dAktdz,                          &
# endif
#  ifdef RM_SUPERCOOLING
     &                    FORCES(ng) % stflx,                           &
#  endif
# ifdef DIAGNOSTICS_TS
     &                    DIAGS(ng) % DiaTwrk,                          &
# endif
     &                    OCEAN(ng) % t)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 35)
# endif
      RETURN
      END SUBROUTINE step3d_t
!
!***********************************************************************
      SUBROUTINE step3d_t_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nstp, nnew,                       &
# if defined TS_PSOURCE || defined Q_PSOURCE
     &                          Msrc, Nsrc,                             &
     &                          Isrc, Jsrc, Tsrc,                       &
# endif
# ifdef TS_PSOURCE
     &                          Dsrc,                                   &
# endif
# ifdef Q_PSOURCE
     &                          Qsrc,                                   &
# endif
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef TS_MPDATA
#  ifdef WET_DRY
     &                          rmask_wet, umask_wet, vmask_wet,        &
#  endif
     &                          omn, om_u, om_v, on_u, on_v,            &
# endif
     &                          pm, pn,                                 &
     &                          Hz, Huon, Hvom,                         &
     &                          z_r,                                    &
# if defined TCLM_NUDGING && defined TCLIMATOLOGY
     &                          Tnudgcof, tclm,                         &
# endif
     &                          Akt,                                    &
     &                          W,                                      &
# if defined FLOATS && defined FLOAT_VWALK
     &                          dAktdz,                                 &
# endif
# ifdef RM_SUPERCOOLING
     &                          stflx,                                  &
# endif
# ifdef DIAGNOSTICS_TS
     &                          DiaTwrk,                                &
# endif
     &                          t)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
# ifdef NESTING
      USE mod_nesting
# endif
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# ifdef DISTRIBUTE
#  if defined FLOATS && defined FLOAT_VWALK
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
      USE mp_exchange_mod, ONLY : mp_exchange4d
# endif
# ifdef TS_MPDATA
      USE mpdata_adiff_mod
# endif
# ifdef NESTING
      USE nesting_mod, ONLY : bry_fluxes
# endif
      USE t3dbc_mod, ONLY : t3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
# if defined TS_PSOURCE || defined Q_PSOURCE
      integer, intent(in) :: Msrc, Nsrc
# endif
!
# ifdef ASSUMED_SHAPE
#  if defined TS_PSOURCE || defined Q_PSOURCE
      integer, intent(in) :: Isrc(:)
      integer, intent(in) :: Jsrc(:)

      real(r8), intent(in) :: Tsrc(:,:,:)
#  endif
#  ifdef TS_PSOURCE
      real(r8), intent(in) :: Dsrc(:)
#  endif
#  ifdef Q_PSOURCE
      real(r8), intent(in) :: Qsrc(:,:)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef TS_MPDATA
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  if defined TCLM_NUDGING && defined TCLIMATOLOGY
      real(r8), intent(in) :: Tnudgcof(LBi:,LBj:,:)
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#  ifdef SUN
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  else
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
#  endif
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
#  ifdef RM_SUPERCOOLING
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
#  endif
#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SUN
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  else
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#  endif
#  if defined FLOATS && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:,LBj:,:)
#  endif

# else

#  if defined TS_PSOURCE || defined Q_PSOURCE
      integer, intent(in) :: Isrc(Msrc)
      integer, intent(in) :: Jsrc(Msrc)

      real(r8), intent(in) :: Tsrc(Msrc,N(ng),NT(ng))
#  endif
#  ifdef TS_PSOURCE
      real(r8), intent(in) :: Dsrc(Msrc)
#  endif
#  ifdef Q_PSOURCE
      real(r8), intent(in) :: Qsrc(Msrc,N(ng))
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef TS_MPDATA
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  if defined TCLM_NUDGING && defined TCLIMATOLOGY
      real(r8), intent(in) :: Tnudgcof(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: tclm(LBi:UBi,LBj:UBj,N(ng),NTCLM(ng))
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef RM_SUPERCOOLING
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
     &                                   NDT)
#  endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

#  if defined FLOATS && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:UBi,LBj:UBj,N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
# ifdef NESTING
      integer :: ILB, IUB, JLB, JUB
      integer :: dg, cr, rg
# endif
      integer :: i, ic, is, itrc, j, k, ltrc
# ifdef DIAGNOSTICS_TS
      integer :: idiag
# endif
      real(r8), parameter :: eps = 1.0E-16_r8
      real(r8), parameter :: a = -0.057_r8
      real(r8), parameter :: b = 0.0939_r8
      real(r8), parameter :: c = 7.61e-4_r8
      real(r8) :: cff, cff1, cff2, cff3, theta_f, S0_Cp_dT_Lf, delta_T
      real(r8) :: Tm
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
# ifdef TS_MPDATA
#  ifdef DIAGNOSTICS_TS
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,3) :: Dhadv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng),NT(ng)) :: Dvadv
#  endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng),NT(ng)) :: Ta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Ua
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Va
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: Wa
# endif

# include "set_bounds.h"

# ifdef NESTING
!
!  Notice that the trace flux boundary arrays are dimensioned with the
!  global dimensions of grid to facilitate processing.
!
      ILB=BOUNDS(ng)%LBi(-1)
      IUB=BOUNDS(ng)%UBi(-1)
      JLB=BOUNDS(ng)%LBj(-1)
      JUB=BOUNDS(ng)%UBj(-1)
# endif
!
!-----------------------------------------------------------------------
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!
!  Compute inverse thickness.
!
# ifdef TS_MPDATA
#  define I_RANGE Istrm2,Iendp2
#  define J_RANGE Jstrm2,Jendp2
# else
#  define I_RANGE Istr,Iend
#  define J_RANGE Jstr,Jend
# endif

      DO k=1,N(ng)
        DO j=J_RANGE
          DO i=I_RANGE
            oHz(i,j,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
      END DO
# undef I_RANGE
# undef J_RANGE
# ifdef TS_MPDATA
!
!  The MPDATA algorithm requires a three-point footprint, so exchange
!  boundary data on t(:,:,:,nnew,:) so other processes computed earlier
!  (horizontal diffusion, biology, or sediment) are accounted.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END DO
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nnew,:))
#  endif
# endif
!
!  Compute horizontal tracer advection fluxes.
!
      T_LOOP : DO itrc=1,NT(ng)
        K_LOOP : DO k=1,N(ng)

# if defined TS_C2HADVECTION
!
!  Second-order, centered differences horizontal advective fluxes.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=Huon(i,j,k)*                                      &
     &                0.5_r8*(t(i-1,j,k,3,itrc)+                        &
     &                        t(i  ,j,k,3,itrc))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=Hvom(i,j,k)*                                      &
     &                0.5_r8*(t(i,j-1,k,3,itrc)+                        &
     &                        t(i,j  ,k,3,itrc))
            END DO
          END DO

# elif defined TS_MPDATA
!
!  First-order, upstream differences horizontal advective fluxes.
!
          DO j=JstrVm2,Jendp2i
            DO i=IstrUm2,Iendp3
              cff1=MAX(Huon(i,j,k),0.0_r8)
              cff2=MIN(Huon(i,j,k),0.0_r8)
              FX(i,j)=cff1*t(i-1,j,k,3,itrc)+                           &
     &                cff2*t(i  ,j,k,3,itrc)
            END DO
          END DO
          DO j=JstrVm2,Jendp3
            DO i=IstrUm2,Iendp2i
              cff1=MAX(Hvom(i,j,k),0.0_r8)
              cff2=MIN(Hvom(i,j,k),0.0_r8)
              FE(i,j)=cff1*t(i,j-1,k,3,itrc)+                           &
     &                cff2*t(i,j  ,k,3,itrc)
            END DO
          END DO

# else
!
#  if defined TS_U3HADVECTION
!  Third-order, uptream-biased horizontal advective fluxes.
#  elif defined TS_A4HADVECTION
!  Fourth-order, Akima horizontal advective fluxes.
#  else
!  Fourth-order, centered differences horizontal advective fluxes.
#  endif
!
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
              FX(i,j)=t(i  ,j,k,3,itrc)-                                &
     &                t(i-1,j,k,3,itrc)
#  ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#  endif
            END DO
          END DO
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Istr-1,j)=FX(Istr,j)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Iend+2,j)=FX(Iend+1,j)
              END DO
            END IF
          END IF
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
#  if defined TS_U3HADVECTION
              curv(i,j)=FX(i+1,j)-FX(i,j)
#  elif defined TS_A4HADVECTION
              cff=2.0_r8*FX(i+1,j)*FX(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#  else
              grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
#  ifdef TS_U3HADVECTION
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,3,itrc)+                               &
     &                 t(i  ,j,k,3,itrc))-                              &
     &                cff1*(curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+        &
     &                      curv(i  ,j)*MIN(Huon(i,j,k),0.0_r8))
#  else
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,3,itrc)+                               &
     &                 t(i  ,j,k,3,itrc)-                               &
     &                 cff2*(grad(i  ,j)-                               &
     &                       grad(i-1,j)))
#  endif
            END DO
          END DO
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
              FE(i,j)=t(i,j  ,k,3,itrc)-                                &
     &                t(i,j-1,k,3,itrc)
#  ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#  endif
            END DO
          END DO
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jstr-1)=FE(i,Jstr)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jend+2)=FE(i,Jend+1)
              END DO
            END IF
          END IF
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
#  if defined TS_U3HADVECTION
              curv(i,j)=FE(i,j+1)-FE(i,j)
#  elif defined TS_A4HADVECTION
              cff=2.0_r8*FE(i,j+1)*FE(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#  else
              grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
#  ifdef TS_U3HADVECTION
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,3,itrc)+                               &
     &                 t(i,j  ,k,3,itrc))-                              &
     &                cff1*(curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+        &
     &                      curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8))
#  else
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,3,itrc)+                               &
     &                 t(i,j  ,k,3,itrc)-                               &
     &                 cff2*(grad(i,j  )-                               &
     &                       grad(i,j-1)))
#  endif
            END DO
          END DO
# endif
# if defined TS_PSOURCE && !defined Q_PSOURCE
!
!  Apply tracers point sources to the horizontal advection terms.
!
          DO is=1,Nsrc
            i=Isrc(is)
            j=Jsrc(is)
            IF (INT(Dsrc(is)).eq.0) THEN
#  ifdef TS_MPDATA
              IF (((IstrUm2.le.i).and.(i.le.Iendp3)).and.               &
     &            ((JstrVm2.le.j).and.(j.le.Jendp2i))) THEN
#  else
              IF (((Istr.le.i).and.(i.le.Iend+1)).and.                  &
     &            ((Jstr.le.j).and.(j.le.Jend))) THEN
#  endif
                IF (LtracerSrc(itrc,ng)) THEN
                  FX(i,j)=Huon(i,j,k)*Tsrc(is,k,itrc)
#  ifdef MASKING
                ELSE
                  IF ((rmask(i  ,j).eq.0.0_r8).and.                     &
     &                (rmask(i-1,j).eq.1.0_r8)) THEN
                    FX(i,j)=Huon(i,j,k)*t(i-1,j,k,3,itrc)
                  ELSE IF ((rmask(i  ,j).eq.1.0_r8).and.                &
     &                     (rmask(i-1,j).eq.0.0_r8)) THEN
                    FX(i,j)=Huon(i,j,k)*t(i  ,j,k,3,itrc)
                  END IF
#  endif
                END IF
              END IF
            ELSE IF (INT(Dsrc(is)).eq.1) THEN
#  ifdef TS_MPDATA
              IF (((IstrUm2.le.i).and.(i.le.Iendp2i)).and.              &
     &            ((JstrVm2.le.j).and.(j.le.Jendp3))) THEN
#  else
              IF (((Istr.le.i).and.(i.le.Iend)).and.                    &
     &            ((Jstr.le.j).and.(j.le.Jend+1))) THEN
#  endif
                IF (LtracerSrc(itrc,ng)) THEN
                  FE(i,j)=Hvom(i,j,k)*Tsrc(is,k,itrc)
#  ifdef MASKING
                ELSE
                  IF ((rmask(i,j  ).eq.0.0_r8).and.                     &
     &                (rmask(i,j-1).eq.1.0_r8)) THEN
                    FE(i,j)=Hvom(i,j,k)*t(i,j-1,k,3,itrc)
                  ELSE IF ((rmask(i,j  ).eq.1.0_r8).and.                &
     &                     (rmask(i,j-1).eq.0.0_r8)) THEN
                    FE(i,j)=Hvom(i,j,k)*t(i,j  ,k,3,itrc)
                  END IF
#  endif
                END IF
              END IF
            END IF
          END DO
# endif

# ifdef NESTING
!
!  If refinement grids, extract tracer horizontal advection fluxes
!  (Hz*u*T/n, Hz*v*T/m) at the grid contact boundary (physical
!  domain period) to be used in two-way nesting.
!
         IF (RefinedGrid(ng)) THEN
           DO cr=1,Ncontact
             dg=Rcontact(cr)%donor_grid
             rg=Rcontact(cr)%receiver_grid
             IF (ng.eq.rg) THEN
               CALL bry_fluxes (dg, rg, cr, iNLM, tile,                 &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          ILB, IUB, JLB, JUB,                     &
     &                          dt(ng), FX, FE,                         &
     &                          BRY_CONTACT(iwest, cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(ieast, cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(isouth,cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(inorth,cr)%Tflux(:,k,itrc))
             END IF
           END DO
         END IF
# endif
!
# ifdef TS_MPDATA
!  Time-step horizontal advection for intermediate diffusive tracer, Ta.
!  Advective fluxes have units of Tunits m3/s.  The new tracer has
!  units of m Tunits.
# else
!  Time-step horizontal advection term.  Advective fluxes have units
!  of Tunits m3/s.  The new tracer has units of m Tunits.
# endif
!
# ifdef TS_MPDATA
#  define I_RANGE IstrUm2,Iendp2i
#  define J_RANGE JstrVm2,Jendp2i
# else
#  define I_RANGE Istr,Iend
#  define J_RANGE Jstr,Jend
# endif
          DO j=J_RANGE
            DO i=I_RANGE
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j)-FX(i,j))
              cff2=cff*(FE(i,j+1)-FE(i,j))
              cff3=cff1+cff2
# ifdef TS_MPDATA
              Ta(i,j,k,itrc)=t(i,j,k,nnew,itrc)-cff3
# else
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff3
# endif
# ifdef DIAGNOSTICS_TS
#  ifdef TS_MPDATA
              Dhadv(i,j,iTxadv)=-cff1
              Dhadv(i,j,iTyadv)=-cff2
              Dhadv(i,j,iThadv)=-cff3
#  else
              DiaTwrk(i,j,k,itrc,iTxadv)=-cff1
              DiaTwrk(i,j,k,itrc,iTyadv)=-cff2
              DiaTwrk(i,j,k,itrc,iThadv)=-cff3
#  endif
# endif
            END DO
          END DO
# if defined DIAGNOSTICS_TS && defined TS_MPDATA
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DiaTwrk(i,j,k,itrc,iTxadv)=Dhadv(i,j,iTxadv)
              DiaTwrk(i,j,k,itrc,iTyadv)=Dhadv(i,j,iTyadv)
              DiaTwrk(i,j,k,itrc,iThadv)=Dhadv(i,j,iThadv)
            END DO
          END DO
# endif
        END DO K_LOOP
      END DO T_LOOP
!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
!
      DO j=J_RANGE
        DO itrc=1,NT(ng)

# if defined TS_SVADVECTION
!
!  Build conservative parabolic splines for the vertical derivatives
!  "FC" of the tracer.  Then, the interfacial "FC" values are
!  converted to vertical advective flux.
!
          DO i=Istr,Iend
#  ifdef NEUMANN
            FC(i,0)=1.5_r8*t(i,j,1,3,itrc)
            CF(i,1)=0.5_r8
#  else
            FC(i,0)=2.0_r8*t(i,j,1,3,itrc)
            CF(i,1)=1.0_r8
#  endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                             &
     &                    Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
              CF(i,k+1)=cff*Hz(i,j,k)
              FC(i,k)=cff*(3.0_r8*(Hz(i,j,k  )*t(i,j,k+1,3,itrc)+       &
     &                             Hz(i,j,k+1)*t(i,j,k  ,3,itrc))-      &
     &                     Hz(i,j,k+1)*FC(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef NEUMANN
            FC(i,N(ng))=(3.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/     &
     &                  (2.0_r8-CF(i,N(ng)))
#  else
            FC(i,N(ng))=(2.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/     &
     &                  (1.0_r8-CF(i,N(ng)))
#  endif
          END DO
          DO k=N(ng)-1,0,-1
            DO i=Istr,Iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
              FC(i,k+1)=W(i,j,k+1)*FC(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
          END DO

# elif defined TS_A4VADVECTION
!
!  Fourth-order, Akima vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=t(i,j,k+1,3,itrc)-                                &
     &                t(i,j,k  ,3,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=FC(i,1)
            FC(i,N(ng))=FC(i,N(ng)-1)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=2.0_r8*FC(i,k)*FC(i,k-1)
              IF (cff.gt.eps) THEN
                CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
              ELSE
                CF(i,k)=0.0_r8
              END IF
            END DO
          END DO
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,3,itrc)+                        &
     &                        t(i,j,k+1,3,itrc)-                        &
     &                        cff1*(CF(i,k+1)-CF(i,k)))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
          END DO
# elif defined TS_C2VADVECTION
!
!  Second-order, central differences vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,3,itrc)+                        &
     &                        t(i,j,k+1,3,itrc))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
          END DO
# elif defined TS_MPDATA
!
!  First_order, upstream differences vertical advective flux.
!
          DO i=I_RANGE
            DO k=1,N(ng)-1
              cff1=MAX(W(i,j,k),0.0_r8)
              cff2=MIN(W(i,j,k),0.0_r8)
              FC(i,k)=cff1*t(i,j,k  ,3,itrc)+                           &
     &                cff2*t(i,j,k+1,3,itrc)
            END DO
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
          END DO
# else
!
!  Fourth-order, central differences vertical advective flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                (cff2*(t(i,j,k  ,3,itrc)+                         &
     &                       t(i,j,k+1,3,itrc))-                        &
     &                 cff3*(t(i,j,k-1,3,itrc)+                         &
     &                       t(i,j,k+2,3,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*2.0_r8*                                    &
     &              (cff2*t(i,j,1,3,itrc)-                              &
     &               cff3*t(i,j,2,3,itrc))
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,1)=W(i,j,1)*                                           &
     &              (cff1*t(i,j,1,3,itrc)+                              &
     &               cff2*t(i,j,2,3,itrc)-                              &
     &               cff3*t(i,j,3,3,itrc))
            FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                               &
     &                    (cff1*t(i,j,N(ng)  ,3,itrc)+                  &
     &                     cff2*t(i,j,N(ng)-1,3,itrc)-                  &
     &                     cff3*t(i,j,N(ng)-2,3,itrc))
            FC(i,N(ng))=0.0_r8
          END DO
# endif
!
# ifdef TS_MPDATA
!  Time-step vertical advection for intermediate diffusive tracer, Ta
!  (Tunits).
# else
#  ifdef SPLINES
!  Time-step vertical advection term (Tunits).
#  else
!  Time-step vertical advection term (m Tunits).
#  endif
# endif
# ifdef DIAGNOSTICS_TS
!  Convert units of tracer diagnostic terms to Tunits.
# endif
!
          DO i=I_RANGE
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
# ifdef Q_PSOURCE
!
!  Apply mass point sources - Volume influx.
!
          DO is=1,Nsrc
            IF (LtracerSrc(itrc,ng)) THEN
              i=Isrc(is)
#  ifdef TS_MPDATA
              IF (((IstrUm2.le.i).and.(i.le.Iendp2i)).and.              &
     &            ((JstrVm2.le.j).and.(j.le.Jendp2i))) THEN
#  else
              IF (((IstrR.le.i).and.(i.le.IendR)).and.                  &
     &            ((JstrR.le.j).and.(j.le.JendR))) THEN
#  endif
                IF (j.eq.Jsrc(is)) THEN
                  DO k=1,N(ng)-1
                    FC(i,k)=FC(i,k)+0.5_r8*                             &
     &                      (Qsrc(is,k  )*Tsrc(is,k  ,itrc)+            &
     &                       Qsrc(is,k+1)*Tsrc(is,k+1,itrc))
                  END DO
                END IF
              END IF
            END IF
          END DO
# endif
          DO k=1,N(ng)
            DO i=I_RANGE
              cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
# ifdef TS_MPDATA
              Ta(i,j,k,itrc)=(Ta(i,j,k,itrc)-cff1)*oHz(i,j,k)
#  ifdef DIAGNOSTICS_TS
              Dvadv(i,j,k,itrc)=-cff1
#  endif
# else
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
#  ifdef SPLINES
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*oHz(i,j,k)
#  endif
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvadv)=-cff1
              DO idiag=1,NDT
                DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*    &
     &                                    oHz(i,j,k)
              END DO
#  endif
# endif
            END DO
          END DO
        END DO
# undef I_RANGE
# undef J_RANGE
# ifdef TS_MPDATA
      END DO
!
!-----------------------------------------------------------------------
!  Compute anti-diffusive velocities to corrected advected tracers
!  using MPDATA recursive method.  Notice that pipelined J-loop ended.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        CALL mpdata_adiff_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#  ifdef MASKING
     &                          rmask, umask, vmask,                    &
#  endif
#  ifdef WET_DRY
     &                          rmask_wet, umask_wet, vmask_wet,        &
#  endif
     &                          pm, pn, omn, om_u, on_v,                &
     &                          z_r, oHz,                               &
     &                          Huon, Hvom, W,                          &
     &                          t(:,:,:,3,itrc),                        &
     &                          Ta(:,:,:,itrc),  Ua, Va, Wa)
!
!  Compute anti-diffusive corrected advection fluxes.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff1=MAX(Ua(i,j,k),0.0_r8)
              cff2=MIN(Ua(i,j,k),0.0_r8)
              FX(i,j)=(cff1*Ta(i-1,j,k,itrc)+                           &
     &                 cff2*Ta(i  ,j,k,itrc))*                          &
     &                0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff1=MAX(Va(i,j,k),0.0_r8)
              cff2=MIN(Va(i,j,k),0.0_r8)
              FE(i,j)=(cff1*Ta(i,j-1,k,itrc)+                           &
     &                 cff2*Ta(i,j  ,k,itrc))*                          &
     &                0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
            END DO
          END DO
!
!  Time-step corrected horizontal advection (Tunits m).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j)-FX(i,j))
              cff2=cff*(FE(i,j+1)-FE(i,j))
              cff3=cff1+cff2
              t(i,j,k,nnew,itrc)=Ta(i,j,k,itrc)*Hz(i,j,k)-cff3
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTxadv)=DiaTwrk(i,j,k,itrc,iTxadv)-    &
     &                                   cff1
              DiaTwrk(i,j,k,itrc,iTyadv)=DiaTwrk(i,j,k,itrc,iTyadv)-    &
     &                                   cff2
              DiaTwrk(i,j,k,itrc,iThadv)=DiaTwrk(i,j,k,itrc,iThadv)-    &
     &                                   cff3
#  endif
            END DO
          END DO
        END DO
!
!  Compute anti-diffusive corrected vertical advection flux.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=MAX(Wa(i,j,k),0.0_r8)
              cff2=MIN(Wa(i,j,k),0.0_r8)
              FC(i,k)=cff1*Ta(i,j,k  ,itrc)+                            &
     &                cff2*Ta(i,j,k+1,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Time-step corrected vertical advection (Tunits).
#  ifdef DIAGNOSTICS_TS
!  Convert units of tracer diagnostic terms to Tunits.
#  endif
!
          DO i=Istr,Iend
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvadv)=Dvadv(i,j,k,itrc)-             &
     &                                   cff1
              DO idiag=1,NDT
                DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*    &
     &                                    oHz(i,j,k)
              END DO
#  endif
            END DO
          END DO
        END DO
      END DO
!
!  Start pipelined J-loop.
!
      DO j=Jstr,Jend
# endif /* TS_MPDATA */
!
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
!
        DO itrc=1,NT(ng)
          ltrc=MIN(NAT,itrc)

# if defined SPLINES && !defined TS_MPDATA
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hz(i,j,k  )-                                 &
     &                dt(ng)*Akt(i,j,k-1,ltrc)*oHz(i,j,k  )
              CF(i,k)=cff1*Hz(i,j,k+1)-                                 &
     &                dt(ng)*Akt(i,j,k+1,ltrc)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                     &
     &                dt(ng)*Akt(i,j,k,ltrc)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(t(i,j,k+1,nnew,itrc)-t(i,j,k,nnew,itrc)-     &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*Akt(i,j,k,ltrc)
              cff1=dt(ng)*oHz(i,j,k)*(DC(i,k)-DC(i,k-1))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff1
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+    &
     &                                   cff1
#  endif
            END DO
          END DO
# else
!
!  Compute off-diagonal coefficients FC [lambda*dt*Akt/Hz] for the
!  implicit vertical diffusion terms at future time step, located
!  at horizontal RHO-points and vertical W-points.
!  Also set FC at the top and bottom levels.
!
          cff=-dt(ng)*lambda
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
              FC(i,k)=cff*cff1*Akt(i,j,k,ltrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
              DC(i,k)=t(i,j,k,nnew,itrc)
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_TS
             cff1=t(i,j,N(ng),nnew,itrc)*oHz(i,j,N(ng))
#  endif
             DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/     &
     &                   (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
             t(i,j,N(ng),nnew,itrc)=DC(i,N(ng))
#  ifdef DIAGNOSTICS_TS
             DiaTwrk(i,j,N(ng),itrc,iTvdif)=                            &
     &                              DiaTwrk(i,j,N(ng),itrc,iTvdif)+     &
     &                              t(i,j,N(ng),nnew,itrc)-cff1
#  endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_TS
              cff1=t(i,j,k,nnew,itrc)*oHz(i,j,k)
#  endif
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              t(i,j,k,nnew,itrc)=DC(i,k)
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+    &
     &                                   t(i,j,k,nnew,itrc)-cff1
#  endif
            END DO
          END DO
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions and, if appropriate, nudge
!  to tracer data and apply Land/Sea mask.
!-----------------------------------------------------------------------
!
# if defined TCLM_NUDGING && defined TCLIMATOLOGY
      ic=0
# endif
      DO itrc=1,NT(ng)
!
!  Set lateral boundary conditions.
!
        CALL t3dbc_tile (ng, tile, itrc,                                &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   t)

# if defined TCLM_NUDGING && defined TCLIMATOLOGY
!
!  Nudge towards tracer climatology.
!
        IF (LtracerCLM(itrc,ng)) THEN
          ic=ic+1
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+                  &
     &                             dt(ng)*Tnudgcof(i,j,itrc)*           &
     &                             (tclm(i,j,k,ic)-t(i,j,k,nnew,itrc))
              END DO
            END DO
          END DO
        END IF
# endif
!-----------------------------------------------------------------------
!     If RM_SUPERCOOLING is activated, limit temperatures to the freezing temperature.
!     Remove super cooling by basal refreezing or pseudo frazil ice production.  
!     Heat input is assumed to be entirely sourced from latent heat of ice production.
!-----------------------------------------------------------------------

#if defined RM_SUPERCOOLING
       IF (itrc.eq.itemp) THEN
         Do j=JstrR,JendR
           Do i=IstrR,IendR
             Do k=1,N(ng)
              ! calculate freezing point in potential temperature
!              theta_f = pot_freezt(t(i,j,k,nnew,isalt),-z_r(i,j,k),i,j,k)
              theta_f=a*t(i,j,k,nnew,isalt)+b+c*-z_r(i,j,k)
              ! check if temperature is below freezing point
        CALL potit(t(i,j,k,nnew,isalt),t(i,j,k,nnew,itemp),-z_r(i,j,k),0.0_r8,Tm,i,j)
              delta_T = theta_f - Tm !t(i,j,k,nnew,itemp) !itrc)
              If (delta_T.gt.0.0_r8) THEN
                 ! (1) compute heat to be produced by ice formation per cell
                 ! Q = Cp * (theta_f - T) * rho0 * Hz * 1/pn * 1/pm
                 ! [J] = [J/kg/K] * [K] * [kg/m^3] * [m] * [m] * [m]
                 
                 ! (2) latent heat equates to fresh water flux per grid cell
                 ! Q / dt = - wb * 1/pn * 1/pm * rho0 * Lf 
                 ! [J/s] = [m/s] * [m] * [m] * [kg/m^3] * [J/kg]
                 
                 ! combine (1) and (2)
                 ! wb = - Cp * (theta_f - T) * Hz / (Lf * dt)
                 ! [m/s] = [J/kg/K] * [K] * [m] / ( [J/kg] * [s] )
                 
                 ! (3) adjust cell salinity (linear combination)
                 ! S1 = S0 * (1 - wb / Hz * dt)  ==  S0 * (1 + Cp * (theta_f - T) / Lf) == S0 + S0*Cp*(theta_f-T)/Lf
                 ! [] = [] * [m/s] / [m] * [s] == [] * [J/kg/K] * [K] / [J/Kg]   
                 
                 ! (4) adjust surface salt flux variable (vertical integral of freezing)
                 ! Sflux = Sflux + (S1 - S0) / dt * Hz == Sflux  +  S0*Cp*(theta_f-T)/Lf / dt * Hz   
                 ! [m/s] = [m/s] + [] / [s] * [m] == [m/s] + [] * [J/kg/K] * [K] / [J/kg] / [s] * [m]
                 
                 ! (5) adjust potential temperature in cell to freezing point
                 ! T = theta_f
                 ! [C] = [C]
                 
                 ! (3)
                 S0_Cp_dT_Lf = (t(i,j,k,nnew,isalt) *3985.0_r8* delta_T)&
     &                       / 3.34e5_r8
                 t(i,j,k,nnew,isalt) = t(i,j,k,nnew,isalt) + S0_Cp_dT_Lf
                 
                 ! (4)
                 stflx(i,j,isalt)=stflx(i,j,isalt)+S0_Cp_dT_Lf/dt(ng)   &
     &                           / Hz(i,j,k)
                 
                 ! (5)
       t(i,j,k,nnew,itemp) = thetaa(t(i,j,k,nnew,isalt),theta_f,-z_r(i,j,k),0.0_r8)
!                 t(i,j,k,nnew,itemp) = t(i,j,k,nnew,itemp) + delta_T
     CALL potit(t(i,j,k,nnew,isalt),t(i,j,k,nnew,itemp),-z_r(i,j,k),0.0_r8,Tm,i,j)
                 delta_T = theta_f - Tm !t(i,j,k,nnew,itemp) !itrc)
                  write(6,*) delta_T
!                 write(6,*) 'restoring freez temp at ',i,j,k, ' is ',   &
!     &                        theta_f, ' ', delta_T ,'degC' 
               endif
             enddo
           enddo
         enddo
       endif
#endif


# ifdef MASKING
!
!  Apply Land/Sea mask.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*rmask(i,j)
            END DO
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_TS
!
!  Compute time-rate-of-change diagnostic term.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
     &                                   DiaTwrk(i,j,k,itrc,iTrate)
!!            DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
!!   &                                   t(i,j,k,nstp,itrc)
            END DO
          END DO
        END DO
# endif
!
!  Apply periodic boundary conditions.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END IF
      END DO
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nnew,:))
# endif
# if defined FLOATS && defined FLOAT_VWALK
!
!-----------------------------------------------------------------------
!  Compute vertical gradient in vertical T-diffusion coefficient for
!  floats random walk.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          DO k=1,N(ng)
            dAktdz(i,j,k)=(Akt(i,j,k,1)-Akt(i,j,k-1,1))/Hz(i,j,k)
          END DO
        END DO
      END DO
!
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          dAktdz)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    dAktdz)
#  endif
# endif
      RETURN
      END SUBROUTINE step3d_t_tile

! *********************************************************************
      SUBROUTINE potit(Sal,theta,Pres,RPres,Temp,i,j)
! *********************************************************************
! Calculates from the salinity (sal, psu), potential temperature 
! (theta, degC) and reference pressure (pres, dbar) the in-situ 
! temperaure (Temp_insitu, degC) related to the in-situ pressure 
! (rfpres, dbar) with the help of an iterative method.
      USE mod_kinds

      integer, intent(in)   :: i, j
      real(r8), intent(in)  :: Sal, Pres,theta
      real(r8), intent(out) :: Temp

      integer               :: ind
      real(r8)              :: tpmd, theta1, thetad, epsi, RPres

      data tpmd / 0.001 /

      epsi = 0.
      do ind=1,100
      Temp   = theta+epsi
      thetad  = thetaa(Sal,Temp,Pres,RPres)-theta
      IF(abs(thetad).lt.tpmd) return
       epsi = epsi-thetad
      ENDdo
      write(6,*) ' WARNING!',                                          &
      & ' in-situ temperature calculation has not converged!', i,j
      RETURN
      END SUBROUTINE potit
! *********************************************************************
      REAL FUNCTION thetaa(Sal,Temp,Pres,RPres)
! Calculates from the salinity (sal, psu), the in-situ temperature 
! (Temp, degC) and the in-situ pressure press, dbar) the potential 
! temperature (Theta, degC) converted to the reference pressure
! (RPres, dbar). A Runge-Kutta procedure of the fourth order is used.
!
! Check value: theta   =    36.89073  degC
!         given sal    =    40.0      psu
!               Temp   =    40.0      degC
!               pres   = 10000.000    dbar
!               rfpres =     0.000    dbar
      USE mod_kinds

      real(r8), intent(in) ::  Sal,Temp,Pres,RPres
      real(r8)             ::  p,t,dp,dt,q,ct2,ct3,cq2a,cq2b,cq3a,cq3b

      data ct2 ,ct3  /0.29289322 ,  1.707106781/
      data cq2a,cq2b /0.58578644 ,  0.121320344/
      data cq3a,cq3b /3.414213562, -4.121320344/

      p  = Pres
      t  = Temp
      dp = RPres-Pres
      dt = dp*dTemp(Sal,t,p)
      t  = t +0.5*dt
      q = dt
      p  = p +0.5*dp
      dt = dp*dTemp(Sal,t,p)
      t  = t + ct2*(dt-q)
      q  = cq2a*dt + cq2b*q
      dt = dp*dTemp(Sal,t,p)
      t  = t + ct3*(dt-q)
      q  = cq3a*dt + cq3b*q
      p  = RPres
      dt = dp*dTemp(Sal,t,p)

      thetaa = t + (dt-q-q)/6.0

      END FUNCTION thetaa

! *********************************************************************
! *********************************************************************
      REAL FUNCTION dTemp(Sal,Temp,Pres)
! Calculates from the salinity (Sal,psu), the in-situ Temperature
! (Temp, degC) and the in-situ pressure (Pres, dbar) the adiabatic 
! temperature gradient (dTemp, K Dbar^-1).
!
! Check values: dTemp  =     3.255976E-4 K dbar^-1
!          given Sal    =    40.0         psu
!                Temp   =    40.0         degC
!                Pres   = 10000.000       dbar
      USE mod_kinds

      real(r8), intent(in) :: Sal, Temp, Pres
      real(r8)             :: s0,a0,a1,a2,a3,b0,b1,c0,c1,c2,c3
      real(r8)             :: d0,d1,e0,e1,e2,ds

      data s0 /35.0D0/
      data a0,a1,a2,a3 /3.5803D-5, 8.5258D-6, -6.8360D-8, 6.6228D-10/
      data b0,b1       /1.8932D-6, -4.2393D-8/
      data c0,c1,c2,c3 /1.8741D-8, -6.7795D-10, 8.7330D-12, -5.4481D-14/
      data d0,d1       /-1.1351D-10, 2.7759D-12/
      data e0,e1,e2    /-4.6206D-13,  1.8676D-14, -2.1687D-16/

      ds = Sal-s0
      dTemp = ( ( (e2*Temp + e1)*Temp + e0 )*Pres                      &
      &      + ( (d1*Temp + d0)*ds                                     &
      &      + ( (c3*Temp + c2)*Temp + c1 )*Temp + c0 ) )*Pres         &
      &      + (b1*Temp + b0)*ds +  ( (a3*Temp + a2)*Temp + a1 )*Temp  &
      &      + a0
      RETURN
      END FUNCTION dTemp

#endif
      END MODULE step3d_t_mod
