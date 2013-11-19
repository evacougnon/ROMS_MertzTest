ROMS iceshelfBen's version -- romsiceshelf

with typo corrections:

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

