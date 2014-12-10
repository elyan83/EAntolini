CFLAGS = -O2 -I cfitsio/
IMAGEPROC2FO = imageproc2f.o jfkdtree.o readimage.o printerror.o
imageproc2f : $(IMAGEPROC2FO)
	gcc -O2 -o imageproc2f $(IMAGEPROC2FO) -lm cfitsio/libcfitsio.a -I cfitsio/
clean : 
	rm *.o
	rm imageproc2f
