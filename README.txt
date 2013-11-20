ROMS iceshelfBen's version -- romsiceshelf

with typo corrections:

** 20/11/2013 -- modified from Dave's edits
in ROMS/External/ocean_icetest.in: Fixed ocean_icetest.in to have the correct CPP flag
in ROMS/Nonlinear/step2d_LF_AM3.h: Corrected CPP code in step2d routine
in ROMS/Nonlinear/IceShelf/iceshelf_vbc.F: l.514 There is an error in the calculation of the drag stresses, sustr and svstr, in the quadratic parameterisation. This is a coding error, where cff2 is used instead of cff1.

** 15/10/2013 -- modified by evacougnon:
in ROMS/Nonlinear/IceShelf/iceshelf_vbc.F (from l.197 to 204): 
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EW_PERIODIC=.TRUE.
#  else
      logical :: EW_PERIODIC=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NS_PERIODIC=.TRUE.
#  else
      logical :: NS_PERIODIC=.FALSE.
#  endif
# endif

** 18/09/2013 -- modified by evacougnon: 
in ROMS/Functionals/ 2 files (ana_initial.h -- l.733 and ana_smflux.h -- l.338): "define ICESHELF2D" replaced by "defined ICESHELFD"

