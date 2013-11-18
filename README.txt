Ben's version -- romsiceshelf

with 2 typo corrections:

* in ROMS/Functionals/ 2 files (ana_initial.h and ana_smflux.h): defined ICESHELFD corrected

* in ROMS/Nonlinear/IceShelf/iceshelf_vbc.F: 

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

