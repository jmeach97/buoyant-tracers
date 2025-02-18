timecorrelated-passive:
	gfortran -I/usr/local/include -o timecorrelated_passive.exe source/interpolation.f90 source/runge_kutta/velocitymodule.f90 source/runge_kutta/particlemodule.f90 source/velocity_models/timecorrelated.f90 source/particle_models/passivetracer.f90 source/main_programs/main_timecorrelated.f90 -lfftw3 -O3
