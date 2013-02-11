README
***
=========
precursor 69b2
Changes - the main reason we went to a new version number is that the KWD wind was changed to allow for a new parameter to set the velocity at the base of the wind to be a multiple of the sound speed rather than fixed to be equal to the sound speed. This means that a pf file with a KWD wind from 69b wont work in 69c. Here is a list of the changes.
knigge.c
line to read in the new parameter - geo.kn_v_zero
line to multiply vzero by the new parameter
matom.c
deleted declaration of rho_test, make had been complaining it was never used, and I like to not have any errors other than the ones I have made!
photon2d.c
added a couple of counters to count up electron scatters vs resonant scatters in a cell. This populates two new elements in the xplasma structure nscat_es and nscat_res. There is also a change to py_wind to allow them to be written out.
python.h
added the new parameter into the geo.kn structures and also the two new counters mentioned above into the plasma structure
resonate.c
commented out the declaratiomns of ftemp and beta - not used
shell_wind.c
removed unused variable dvtemp
wind_updates2d.c
removed now unused d2_agn parameter. This was used to compute the IP,but we do it better now.


