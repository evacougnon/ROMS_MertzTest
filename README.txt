romsiceshelf version from master branch adapted  for Mertz domain

** General code improvements:
09/12/2013 -- modifed from Dave's edits
Moved calculation of surface stresses before ustar, so ustar is calculated correctly, in ROMS/Nonlinear/IceShelf/iceshelf_vbc.F (see bug report -- 11/2013: Bug - Ustar calculated from sustr and svstr, before sustr and svstr calculated)

10/12/2013 -- modified from Dave's edits
Bad code in ROMS/Nonlinear/IceShelf/iceshelf_vbc.F that executes with ICETEST (l.59, 61, 129, 131 --> commented `#ifdef ICESHELF_MORPH' and `#endif' AND commented l.367 `write(6,*))

20/12/2013 -- GAMMAS calculation and ustar
in Nonlinear/Iceshelf/iceshelf_mod.F
Change: small=8.08e-6_r8 instead of 1.0e-10_r8.
Because in Nonlinear/Iceshelf/iceshelf_vbc.F -- GAMMA_Mc1987, ustar.gt.small but GammaT needs to be less or equal to the molecular thermal diffusivity of sea water (Kt = 1.4e-7) which gives ustar~8.08e-6

** Adapted for an idealised Mertz Glacier region
Application now call MERTZ_TEST
Include ANA_MASK with a band of land south or the domain allowing a fjord in the middle of the south border
Iceshelf: ice tongue with curved base 
Self break: using atan function to avoid strong angle

Resolution of 5 km and 31 vertical levels. 



