# lare2d

A version of Lare2d, set up to model the solar corona. Full details to be in a future paper, probably.

The original Lare2d readme is present, for making changes to the original code (/my modifications)

TO USE:

    Obtain a version of SDF from https://github.com/keithbennett/SDF/tree/632d84fb9b2e48f1b7bca95a531559c9f4bf1c99 or equivalent. The necessary folders are C, FORTRAN and utilities, I'm not sure what the others do. Place these inside the lare2d folder.

    Compile SDF etc. Navigate to 'C' and type 'build' then navigate to FORTRAN and type 'make compiler=gfortran' (or intel if that is the case, like on Hamilton). The makefile in the FORTRAN folder will need to be updated on the line 'mpif90 ?=' to reflect where mpi is installed. On the maths linux computers it is '/usr/lib64/openmpi/bin/mpif90'. On Hamilton it is just 'mpiifort'.

    If necessary, change the output directories in 'runlare.py' and './src/control.f90'. It is currently set to ./Data.

    To run simulations, run the python wrapper 'runlare.py'. This contains various parameters that can be altered (11 I think). These are reasonably self-explanatory I think. Before each run, the selected parameters will be printed (just before the big LARE2D sign), and it's worth checking that this has been done correctly.

    To plot the outputs (which are saved as .sdf files) use 'plotlare.py'. This plots various diagnostics, and can be run at the same time as lare itself (the plotter will wait for the code to catch up if necessary). Note if you do this too early that the vmin and vmax on the colormaps are fixed from the start, and so later snapshots may be well out of these bounds.

Have fun.
