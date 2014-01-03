romsiceshelf version from master branch adapted  for Mertz domain

** General code improvements:
09/12/2013 -- modifed from Dave's edits
Moved calculation of surface stresses before ustar, so ustar is calculated correctly, in ROMS/Nonlinear/IceShelf/iceshelf_vbc.F (see bug report -- 11/2013: Bug - Ustar calculated from sustr and svstr, before sustr and svstr calculated)

10/12/2013 -- modified from Dave's edits
Bad code in ROMS/Nonlinear/IceShelf/iceshelf_vbc.F that executes with ICETEST (l.59, 61, 129, 131 --> commented `#ifdef ICESHELF_MORPH' and `#endif' AND commented l.367 `write(6,*))

20/12/2013 -- GAMMAS calculation and ustar
in Nonlinear/Iceshelf/iceshelf_mod.F
Changed: small=8.08e-6_r8 instead of 1.0e-10_r8.
Because in Nonlinear/Iceshelf/iceshelf_vbc.F -- GAMMA_Mc1987, ustar.gt.small but GammaT needs to be less or equal to the molecular thermal diffusivity of sea water (Kt = 1.4e-7) which gives ustar~8.08e-6

03/01/2014 -- ANA_SEAICE
1- in Nonlinear/Iceshelf/iceshelf_mod.F
Changed the position of the definition of gamma, refsalt and temp_f for IF DEFINED ICESHELF_3EQN_VBC out from IF DEFINED ANA_SEAICE
2- Nonlinear/IceShelf/iceshelf_vbc.F
Moved #endif at l.530 to l.512 (befose ELSE which indicate the zice=0 case and then allows to consider ANA_SEAICE to define seasonal flux for the open ocean)

** Adapted for an idealised Mertz Glacier region
Application now call MERTZ_TEST
Include ANA_MASK with a band of land south or the domain allowing a fjord in the middle of the south border
Iceshelf: ice tongue with curved base 
Self break: using atan function to avoid strong angle

Resolution of 5 km and 31 vertical levels. 



