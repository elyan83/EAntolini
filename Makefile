.SUFFIXES: .pdf .tex
.tex.pdf:
	pdflatex  $< -interaction=nonstopmode
	pdflatex  $< -interaction=nonstopmode
PROPFILES = budget.pdf proposal.pdf research.pdf
all : imagetest2 imagecreate2 imageproc2 imlist readwritetest $(PROPFILES)
budget.pdf : rti.pdf
proposal.pdf : rti.pdf
research.pdf : rti.pdf
CFLAGS = -O2 -I cfitsio/
IMAGETEST2O = imagetest2.o jkdtree.o kdtree.o ran1.o gasdev.o
imagetest2 : $(IMAGETEST2O)
	gcc -O2 -o imagetest2 $(IMAGETEST2O) -lm
IMAGEPROC2O = imageproc2.o jkdtree.o kdtree.o readimage.o printerror.o
imageproc2 : $(IMAGEPROC2O)
	gcc -O2 -o imageproc2 $(IMAGEPROC2O) -lm cfitsio/libcfitsio.a -I cfitsio/
IMAGEPROC2FO = imageproc2f.o jfkdtree.o kdtree.o readimage.o printerror.o
imageproc2f : $(IMAGEPROC2FO)
	gcc -O2 -o imageproc2f $(IMAGEPROC2FO) -lm cfitsio/libcfitsio.a -I cfitsio/
IMAGEPROC2FBO = imageproc2fb.o jfkdtree.o kdtree.o readimage.o printerror.o
imageproc2fb : $(IMAGEPROC2FBO)
	gcc -O2 -o imageproc2fb $(IMAGEPROC2FBO) -lm cfitsio/libcfitsio.a -I cfitsio/
imlist : imlist.c
	gcc -O2 -o imlist imlist.c cfitsio/libcfitsio.a -I cfitsio/ -lm
READWRITETESTO = readwritetest.o writeimage.o readimage.o printerror.o
readwritetest : $(READWRITETESTO)
	gcc -O2 -o readwritetest $(READWRITETESTO) -lm cfitsio/libcfitsio.a -I cfitsio/
READWRITEFITSO = readwritefits.o
readwritefits : $(READWRITEFITSO)
	gcc -O2 -o readwritefits $(READWRITEFITSO) -lm cfitsio/libcfitsio.a -I cfitsio/
IMAGECREATE2O = imagecreate2.o jkdtree.o kdtree.o ran1.o gasdev.o writeimage.o printerror.o
imagecreate2 : $(IMAGECREATE2O)
	gcc -O2 -o imagecreate2 $(IMAGECREATE2O) -lm cfitsio/libcfitsio.a -I cfitsio/
# IMAGECREATEO = imagecreate.o kdtree.o ran1.o gasdev.o writeimage.o printerror.o
# imagecreate : $(IMAGECREATEO)
# 	gcc -O2 -o imagecreate $(IMAGECREATEO) -lm cfitsio/libcfitsio.a -I cfitsio/
# IMAGEPROCO = imageproc.o kdtree.o readimage.o printerror.o
# imageproc : $(IMAGEPROCO)
# 	gcc -O2 -o imageproc $(IMAGEPROCO) -lm cfitsio/libcfitsio.a -I cfitsio/
# IMAGETESTO = imagetest.o kdtree.o ran1.o gasdev.o
# imagetest : $(IMAGETESTO)
# 	gcc -O2 -o imagetest $(IMAGETESTO) -lm
PACKAGECONTENTS = imagetest2.c imagecreate2.c imageproc2.c Makefile \
	jkdtree.c jkdtree.h kdtree.c kdtree.h ran1.c gasdev.c readimage.c \
	imlist.c writeimage.c printerror.c readwritetest.c
ursi.zip : $(PACKAGECONTENTS)
	zip ursi.zip $(PACKAGECONTENTS)
clean : 
	rm *.o
	rm imlist imagetest imagetest2 imagecreate imageproc readwritetest
