romsiceshelf version from mertz branch adapted for a 2 dimentional simple idealised model
All started with the base of icetest scripts

** General code improvements:
09/01/2014 by updating the code to adapt it for mertz_test, a mistake has been included in ana_initial and ana_smflux for icetest experiment! defined ICETEST had been removed while it should be included! This is now corrected and ana_initial and ana_smflux include defined ICETEST, ICETEST2D and MERTZ_TEST.

09/01/2014: The previous mistake may have induced an instability and small(=8.08e-6) needs to be changed! small=1.0e-5 

09/01/2014: there is a jump in the surface salinity flux over the open ocean. iceshelf_vbc has been updated --> l.245 is now (+0.7_r8 added):   
sfcSalt = saltMax + sRateDec * (259.0_r8 - tyear) + 0.7_r8

** Adapted for icetest 2d 
10/01/2014: code updated to use the new experiment called 'ICETEST2D'
ocean_icetest2d.in
ana_grid
ana_pair (but undef in cpp!)
cppdefs
icetest2d
lmd_vmix.F
ana_initial
ana_smflux




