#include "fitsio.h"

void printerror( int status);

int readimage( unsigned short **data, int *nx, int *ny, const char *filename )

    /************************************************************************/
    /* Read a FITS image and determine the minimum and maximum pixel values */
    /************************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status,  nfound, anynull;
    long naxes[2], fpixel, nbuffer, npixels, ii;

#define buffsize 1000
    float datamin, datamax, nullval, buffer[buffsize];

    status = 0;

    if ( fits_open_file(&fptr, filename, READONLY, &status) )
         printerror( status );

    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
         printerror( status );

    *nx      = naxes[0];
    *ny      = naxes[1];
    
    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    fpixel   = 1;
    nullval  = 0;                /* don't check for null values in the image */
    datamin  = 1.0E30;
    datamax  = -1.0E30;
    if ((*data=(unsigned short *) malloc(npixels*sizeof(unsigned short)))==NULL) {
      printf("Failed to allocate data in readimage() %s:%d\n",__FILE__,__LINE__);
      return -1;
    }

    /* Note that even though the FITS images contains unsigned integer */
    /* pixel values (or more accurately, signed integer pixels with    */
    /* a bias of 32768),  this routine is reading the values into a    */
    /* float array.   Cfitsio automatically performs the datatype      */
    /* conversion in cases like this.                                  */

    if ( fits_read_img(fptr, TUSHORT, fpixel, npixels, &nullval,
		       *data, &anynull, &status) )
      printerror( status );

    if ( fits_close_file(fptr, &status) )
         printerror( status );

    return 0;
}

int readfimage( float **data, int *nx, int *ny, const char *filename )

    /************************************************************************/
    /* Read a FITS image and determine the minimum and maximum pixel values */
    /************************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status,  nfound, anynull;
    long naxes[2], fpixel, nbuffer, npixels, ii;

#define buffsize 1000
    float datamin, datamax, nullval, buffer[buffsize];

    status = 0;

    if ( fits_open_file(&fptr, filename, READONLY, &status) )
         printerror( status );

    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
         printerror( status );

    *nx      = naxes[0];
    *ny      = naxes[1];
    
    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    fpixel   = 1;
    nullval  = 0;                /* don't check for null values in the image */
    datamin  = 1.0E30;
    datamax  = -1.0E30;
    if ((*data=(float *) malloc(npixels*sizeof(float)))==NULL) {
      printf("Failed to allocate data in readfimage() %s:%d\n",__FILE__,__LINE__);
      return -1;
    }

    /* Note that even though the FITS images contains unsigned integer */
    /* pixel values (or more accurately, signed integer pixels with    */
    /* a bias of 32768),  this routine is reading the values into a    */
    /* float array.   Cfitsio automatically performs the datatype      */
    /* conversion in cases like this.                                  */

    if ( fits_read_img(fptr, TFLOAT, fpixel, npixels, &nullval,
		       *data, &anynull, &status) )
      printerror( status );

    if ( fits_close_file(fptr, &status) )
         printerror( status );

    return 0;
}
