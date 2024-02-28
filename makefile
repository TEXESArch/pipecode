fife:	fife.o fits.o pltsub.o fft.o funcs.o
	gfortran fife.o fits.o pltsub.o fft.o funcs.o -o fife -L/usr/local/pgplot -lpgplot -L/opt/X11/lib -lX11
fife.o:	fife.F90 mods.h
	gfortran -O -c -ffixed-line-length-80 -fmax-errors=16 fife.F90
pltsub.o:	pltsub.F90
	gfortran -O -c -ffixed-line-length-80 -fmax-errors=10 pltsub.F90
fft.o:	fft.F90
	gfortran -O -c -ffixed-line-length-80 -fmax-errors=10 fft.F90
fits.o:	fits.F90 mods.h
	gfortran -O -c -ffixed-line-length-80 -fmax-errors=10 fits.F90
funcs.o:	funcs.F90 mods.h
	gfortran -O -c -ffixed-line-length-80 -fmax-errors=10 funcs.F90
