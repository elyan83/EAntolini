#include <sys/times.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "jfkdtree.h"
//#include "kdtree.h"


//how close the ration of the triangle have to be in order to find a peak
#define DISTCUT 1e-5
//How far one peak has to be away for another
#define INTERPEAKDISTANCE 5.0
#define MAXPIXELDISTANCE 5.0
//intehgrazione di 5 pixels
#define INTEGRATEREGION 5.0
#define SEARCHRADIUS 5.0
//if the pixel che chiamiamo picco e' piu' brillante dei 5 piu' vicini e' un raggio cosmico
#define HOTPIXELTHRESHOLD 2.0
#define NBRIGHT 80

//per non essere il rimore un pixel deve essere 3sigma oltre il rumore
#define PIXELNOISEFACTOR 3
#define PEAKTHRESHOLDFACTOR 8

typedef struct pair_data {
  double length;
  int i, j, k;
} pdata;


int idum=10000, npeak;
pix_t *catarray, *peakarray;
double glob_param[6];


int
dcomp(const void *a, const void *b) {
  if (*((double *) a) <*((double *) b)) return 1;
  if (*((double *) a) >*((double *) b)) return -1;
  return 0;
}

int
pcomp(const void *a, const void *b) {
  if ( ((pdata *) a)->length <((pdata *)b)->length) return 1;
if ( ((pdata *) a)->length >((pdata *)b)->length) return -1;
  return 0;
}

int
comppixt(const void *a, const void *b) {
  if (*((pix_t *) a)>*((pix_t *) b)) {
    return -1;
  } else {
    return 1;
  }
}

void
calctransform(pix_t *p1, pix_t *res,  int index1[], int index2[], int npoint, double param[])
{
  double s01, s02, s03, s11, s12, s13, s22, s23, s33;
  double d;
  int i;

  /* if there are fewer than three points, than the determinant will be zero. */
  if (npoint<3) {
    printf("npoint<3 calctransform() %s:%d\n",__FILE__,__LINE__);
  }
  s01=s02=s03=s11=s12=s13=s22=s23=s33=0;
  for (i=0;i<npoint;i++) {
    s01+= (double) p1[3*index1[i]+1]* (double) res[3*index2[i]];
    s02+= (double) p1[3*index1[i]+2]* (double) res[3*index2[i]];
    s03+= (double) res[3*index2[i]];
    s11+= (double) p1[3*index1[i]+1]* (double) p1[3*index1[i]+1];
    s12+= (double) p1[3*index1[i]+1]* (double) p1[3*index1[i]+2];
    s13+= (double) p1[3*index1[i]+1];
    s22+= (double) p1[3*index1[i]+2]* (double) p1[3*index1[i]+2];
    s23+= (double) p1[3*index1[i]+2];
    s33++;
  }
  d=(s13*s13*s22-2*s12*s13*s23+s11*s23*s23+s12*s12*s33-s11*s22*s33);
  param[0]=((s03*s13*s22 - s03*s12*s23 - s02*s13*s23 + s01*s23*s23  + s02*s12*s33 - s01*s22*s33)/d);
  param[1]=((-s03*s12*s13 + s02*s13*s13  + s03*s11*s23 - s01*s13*s23 - s02*s11*s33 + s01*s12*s33)/d);
  param[2]=((s03*s12*s12  - s02*s12*s13 - s03*s11*s22 + s01*s13*s22 + s02*s11*s23 - s01*s12*s23)/d);
}

void
triangleoutput(int ih, int jh, int kh, int iah, int jah, int kah) {
  pdata pa[3], pb[3];
  int i;
  int index1[3], index2[3];
  double param[3];

  /* calculate side lengths from short list and include points */
  pa[0].length=hypot((double) catarray[3*iah+1]-(double) catarray[3*jah+1],(double) catarray[3*iah+2]-(double) catarray[3*jah+2]);
  pa[0].i=iah;
  pa[0].j=jah;
  pa[0].k=kah;
  pa[1].length=hypot((double) catarray[3*jah+1]-(double) catarray[3*kah+1],(double) catarray[3*jah+2]-(double) catarray[3*kah+2]);
  pa[1].i=jah;
  pa[1].j=kah;
  pa[1].k=iah;
  pa[2].length=hypot((double) catarray[3*iah+1]-(double) catarray[3*kah+1],(double) catarray[3*iah+2]-(double) catarray[3*kah+2]);
  pa[2].i=iah;
  pa[2].j=kah;
  pa[2].k=jah;
  /* sort them and include the point data */
  qsort((void *) pa,3,sizeof(pdata),pcomp);

  /* calculate side lengths from short list and include points */
  pb[0].length=hypot((double) peakarray[3*ih+1]-(double) peakarray[3*jh+1],(double) peakarray[3*ih+2]-(double) peakarray[3*jh+2]);
  pb[0].i=ih;
  pb[0].j=jh;
  pb[0].k=kh;
  pb[1].length=hypot((double) peakarray[3*jh+1]-(double) peakarray[3*kh+1],(double) peakarray[3*jh+2]-(double) peakarray[3*kh+2]);
  pb[1].i=jh;
  pb[1].j=kh;
  pb[1].k=ih;
  pb[2].length=hypot((double) peakarray[3*ih+1]-(double) peakarray[3*kh+1],(double) peakarray[3*ih+2]-(double) peakarray[3*kh+2]);
  pb[2].i=ih;
  pb[2].j=kh;
  pb[2].k=jh;
  /* sort them and include the point data */
  qsort((void *) pb,3,sizeof(pdata),pcomp);

  /* determine point correspondences */
  for (i=0;i<3;i++) {
    index1[i]=pa[i].k;
    index2[i]=pb[i].k;
  }


  /* calculate forward transformation */
  calctransform(peakarray,catarray+1,index2,index1,3,glob_param); 
  calctransform(peakarray,catarray+2,index2,index1,3,glob_param+3);
#if 0
  /* output the pair information */
  for (i=0;i<3;i++) {
    printf("Pair #%d %3d - %3d %10.4f ; %3d - %3d %10.4f\n",i,pa[i].i,pa[i].j,pa[i].length,pb[i].i,pb[i].j,pb[i].length);
  }
  /* output the point correspondances */
  printf("{%d %d %d} -> {%d %d %d}\n",index1[0],index1[1],index1[2],index2[0],index2[1],index2[2]);
  printf("Transforms from 1->2\n");
  printf("x2 = %g x1 + %g y1 + %g\n",glob_param[0],glob_param[1],glob_param[2]);
  printf("y2 = %g x1 + %g y1 + %g\n",glob_param[3],glob_param[4],glob_param[5]);

  printf("Transforms from 2->1\n");
  calctransform(catarray,peakarray+1,index1,index2,3,param);
  printf("x1 = %g x2 + %g y2 + %g\n",param[0],param[1],param[2]);
  calctransform(catarray,peakarray+2,index1,index2,3,param);
  printf("y1 = %g x2 + %g y2 + %g\n",param[0],param[1],param[2]);
#endif
}

int
main(int argc, char *argv[]) {
  struct tms first_time, start_time, second_time;
  pix_t *pixel_array, *image, *stripe;
  double cpu_time_used, dumx, dumy;
  //struct kdtree *kd_triangle;
  //struct jkdtree *jkd, *jkd_peak, *jkd_cat;
  //struct jkdres *jres, *pix_res;
  //struct kdres *res;
  int i, j, k,npixel, nx, ny;
  double x, y, key[2], pix_key[2];
  pix_t *pix_ptr, *pix_done, *store_ptr, skylevel, peakthreshold, pixelnoisethreshold;
  pix_t *peak_ptr, *peak_done;
  int imagesize, imagenarray, pixelnarray, nstars;
  char *filename, *data;
  double la[3], ratioarray[2];
  double box_image_x[4];
  double box_image_y[4];
    
  int readfimage( float **data, int *nx, int *ny, const char *filename );

#if 0 
  times(&first_time);
    
  
  /* load the star catalogue */

  /* allocate a k-d tree to quickly find the stars in the catalogue */
  
  
  if ((jkd_cat = jkd_create(2))==NULL) {
    printf("Failed to allocate jkd_cat in main() %s:%d\n",__FILE__,__LINE__);
    return -1;
  }


  /* start to read in the catalogue from standard input */


  /* notice that you can just set nstars in the file to be smaller to leave 
     out some stars at the end -- these might be at the edge of the image 
     so they may be classified as hot pixels */
  
  scanf("%d stars\n",&nstars);

  if ((catarray=(pix_t *) malloc(nstars*3*sizeof(pix_t)))==NULL) {
    printf("Failed to allocate catarray in main() %s:%d\n",__FILE__,__LINE__);
    return -1;
  }


  
  pix_ptr=catarray;
  for (i=0;i<nstars;i++) {
    scanf("%g %g %g\n",pix_ptr,pix_ptr+1,pix_ptr+2);
    //printf("%g %g \n",pix_ptr[1],pix_ptr[2]);
    key[0]=pix_ptr[1];
    key[1]=pix_ptr[2];
      
      
      

    /* save the catalogue coordinates and the flux in the k-d tree - jkd_cat */
    
      /*if (jkd_insert(jkd_cat,key,*pix_ptr)) {
          printf("Failed to insert star in stellar catalogue in main() %s:%d\n",__FILE__,__LINE__);
          return -1;
        }*/
    pix_ptr+=3;
  }
  times(&second_time);

  cpu_time_used = ((double) (second_time.tms_utime - first_time.tms_utime)) / sysconf(_SC_CLK_TCK);

  printf("Load catalogue %ld %g\n",(long) (second_time.tms_utime - first_time.tms_utime),cpu_time_used);
    
  
  times(&first_time);
  qsort((void *) catarray,nstars,sizeof(pix_t)*3,comppixt);
  times(&second_time);

  cpu_time_used = ((double) (second_time.tms_utime - first_time.tms_utime)) / sysconf(_SC_CLK_TCK);

  printf("Sort random catalogue %ld %g\n",(long) (second_time.tms_utime - first_time.tms_utime),cpu_time_used);

  /* create triangle catalogue from the brightest NBRIGHT stars */
  //times(&first_time);

  /* create the kd-tree for the triangles*/
 // kd_triangle = kd_create(2);
  /* designate a function to deallocate the data */
 // kd_data_destructor(kd_triangle,free);

  /* build the tree */
 // for (i=0;i<NBRIGHT-2;i++) {
    /*    printf("%d %d %d\n",catarray[3*i],catarray[3*i+1],catarray[3*i+2]); */
   // for (j=i+1;j<NBRIGHT-1;j++) {
     // for (k=j+1;k<NBRIGHT;k++) {
	/* calculate side lengths */
    
    
	//la[0]=hypot((double) catarray[3*i+1]-(double) catarray[3*j+1],(double) catarray[3*i+2]-(double) catarray[3*j+2]);
	//la[1]=hypot((double) catarray[3*j+1]-(double) catarray[3*k+1],(double) catarray[3*j+2]-(double) catarray[3*k+2]);
	//la[2]=hypot((double) catarray[3*i+1]-(double) catarray[3*k+1],(double) catarray[3*i+2]-(double) catarray[3*k+2]);
	/* sort them */
	//qsort((void *) la,3,sizeof(double),dcomp);
	/* if the smallest side is less than a pixel, skip this triangle */
	//if (la[2]<1) break;
	/* calculate the side ratio */
	//ratioarray[0]=la[1]/la[0];
	//ratioarray[1]=la[2]/la[0];
	/* allocate an array to hold the points */
	//data=(char *) malloc(sizeof(char)*3);
	//data[0]=i; data[1]=j; data[2]=k;
	/* add it to the tree */
	/*if (kd_insert(kd_triangle, ratioarray, (void *) data)) {
	  printf("Failed to insert triangle in kd_triangle in main() %s:%d\n",__FILE__,__LINE__);
	  return -1;
	}
      }
    }
  }
  times(&second_time);

  cpu_time_used = ((double) (second_time.tms_utime - first_time.tms_utime)) / sysconf(_SC_CLK_TCK);

  printf("Build triangle catalogue %ld %g\n",(long) (second_time.tms_utime - first_time.tms_utime),cpu_time_used);*/

#endif
  
  /* here is where the image processing begins; */
  /* everything up to here is either for the simulation or */
  /* could be done beforehand */


  /* what is the file to load? */
  if (argc<2) {
    filename="joe.fits";
  } else {
    filename=argv[1];
  }

  times(&first_time);
      
  /* load the image */
  readfimage( &image, &nx, &ny, filename);

  times(&second_time);

  cpu_time_used = ((double) (second_time.tms_utime - first_time.tms_utime)) / sysconf(_SC_CLK_TCK);

  printf("Load image %ld %g\n",(long) (second_time.tms_utime - first_time.tms_utime),cpu_time_used);

  imagesize=nx; imagenarray=nx*ny; pixelnarray=3*imagenarray;

  printf("# nx= %d ny= %d %f %f\n",nx,ny,image[(imagesize/2)*imagesize+imagesize/2],image[42010]);

//#if 0
  /* output the image */
  /*{ FILE *out;
    if ((out=fopen("image_test.pgm","w"))==NULL) {
      printf("Failed to open image.pgm in main() %s:%d\n",__FILE__,__LINE__);
      return -1;
    }
    fprintf(out,"P2\n%d %d\n65535\n",imagesize,imagesize);
    for (i=0;i<imagenarray;i++) {*/
      /*      fprintf(out,"%d\n",(int) (sqrt(image[i]/65535.0)*65535));*/
    /*  fprintf(out,"%d\n",image[i]);
    }
    fclose(out);
  }*/
//#endif


  /* find sky level by looking at pixel value distribution in a single stripe */

  /* allocation the stripe */
  times(&first_time);
  start_time=first_time;

  if ((stripe=(pix_t *) malloc(imagesize*sizeof(pix_t)))==NULL) {
    printf("Failed to allocate stripe in main() %s:%d\n",__FILE__,__LINE__);
    return -1;
  }
 
  memcpy((void *) stripe, (void *) (image+imagenarray/2),imagesize*sizeof(pix_t));
  qsort((void *) stripe,imagesize,sizeof(pix_t),comppixt);
  printf("25: %f 50: %f 75: %f\n",stripe[imagesize/4],
	 stripe[imagesize/2],stripe[3*imagesize/4]);

  skylevel=stripe[imagesize/2];  
  pixelnoisethreshold=skylevel+PIXELNOISEFACTOR*sqrt((double) skylevel);
  peakthreshold=skylevel+PEAKTHRESHOLDFACTOR*sqrt((double) skylevel);

  /* free the stripe data */
  free((void *) stripe);

  times(&second_time);

  cpu_time_used = ((double) (second_time.tms_utime - first_time.tms_utime)) / sysconf(_SC_CLK_TCK);

  printf("-------------------------------\nFinding sky level %ld %g\n",(long) (second_time.tms_utime - first_time.tms_utime),cpu_time_used);


  /* create pixel array */

  if ((pixel_array=(pix_t *) malloc(pixelnarray*sizeof(pix_t)))==NULL) {
    printf("Failed to allocate pixel_array in main() %s:%d\n",__FILE__,__LINE__);
    return -1;
  }

  times(&first_time);
  store_ptr=pixel_array;
  pix_done=image+imagenarray;
  for (pix_ptr=image;pix_ptr<pix_done;pix_ptr++) {
    /* don't include the dim pixels -- dramatically speeds things up */
    if (*pix_ptr>pixelnoisethreshold) {
      *(store_ptr++)=*pix_ptr;
      *(store_ptr++)=(pix_ptr-image)%imagesize;
      *(store_ptr++)=(pix_ptr-image)/imagesize;
    }
      
  }

  /* free the image, we don't need it any more */
  free((void *) image);

  /* remember how many pixels stored */
  pix_done=store_ptr;
  npixel=(pix_done-pixel_array)/3;

  times(&second_time);

  cpu_time_used = ((double) (second_time.tms_utime - first_time.tms_utime)) / sysconf(_SC_CLK_TCK);


  printf("-------------------------------\nMaking the pixel array %ld %g\n",(long) (second_time.tms_utime - first_time.tms_utime),cpu_time_used);

  /* sort the pixel array */
  times(&first_time);
  qsort((void *) pixel_array,npixel,sizeof(pix_t)*3,comppixt);
  times(&second_time);
  
  
  cpu_time_used = ((double) (second_time.tms_utime - first_time.tms_utime)) / sysconf(_SC_CLK_TCK);

  printf("Sort pixel array %ld %g\n",(long) (second_time.tms_utime - first_time.tms_utime),cpu_time_used);

  times(&first_time);
    
    
  /* create the kd-tree for the pixel data*/
  if ((jkd = jkd_create(2))==NULL) {
    printf("Failed to allocate jkd in main() %s:%d\n",__FILE__,__LINE__);
    return -1;
  }
  /* create the kd-tree for the peak data*/
  if ((jkd_peak = jkd_create(2))==NULL) {
    printf("Failed to allocate jkd_peak in main() %s:%d\n",__FILE__,__LINE__);
    return -1;
  }

  /* allocate peak array for triangle matching */
  if ((peakarray=(pix_t *) malloc(NBRIGHT*3*sizeof(pix_t)))==NULL) {
    printf("Failed to allocate peakarray in main() %s:%d\n",__FILE__,__LINE__);
    return -1;
  }

  /* build the pixel and peak trees */
    
  peak_ptr=peakarray; peak_done=peakarray+NBRIGHT;
  for (pix_ptr=pixel_array;pix_ptr<pix_done;pix_ptr+=3) {
    key[0]=pix_ptr[1];
    key[1]=pix_ptr[2];
    /* is the pixel brighter than the PEAKTHRESHOLD? */
    if (*pix_ptr>peakthreshold) {
      /* find any peaks within the INTERPEAKDISTANCE */
      jres=jkd_nearest_range(jkd_peak,key,INTERPEAKDISTANCE);
      if (jkd_res_end(jres)) {
	/* find any brighter pixels within the MAXPIXELDISTANCE */
          jkd_res_free(jres);
          jres=jkd_nearest_range(jkd,key,MAXPIXELDISTANCE);
      if (jkd_res_end(jres)) {
          if (peak_ptr<peak_done) {
              memcpy((void *) peak_ptr, (void *) pix_ptr,sizeof(pix_t)*3);
              peak_ptr+=3;
          }
	  if (jkd_insert(jkd_peak,key, *pix_ptr)) {
	    printf("Failed to insert peak in main() %s:%d\n",__FILE__,__LINE__);
	    return -1;
	  }
	}
      }
      jkd_res_free(jres);
    }
    /* add the pixel to the pixel tree */
    if (jkd_insert(jkd,key,*pix_ptr)) {
      printf("Failed to insert pixel in main() %s:%d\n",__FILE__,__LINE__);
      return -1;
    }
  }
  npeak=(peak_ptr-peakarray)/3;

  /* free the pixel data */
  free((void *) pixel_array);

  times(&second_time);

  cpu_time_used = ((double) (second_time.tms_utime - first_time.tms_utime)) / sysconf(_SC_CLK_TCK);

  printf("Build the tree %ld %g\n",(long) (second_time.tms_utime - first_time.tms_utime),cpu_time_used);

#if 0
  /* now let's do the triangle matching to find what is the transformation between the peaks and the catalogue */
  //times(&first_time);
  
 // for (i=0;i<npeak-2;i++) {
    /*    printf("%d %d %d\n",peakarray[3*i],peakarray[3*i+1],peakarray[3*i+2]); */
    //for (j=i+1;j<npeak-1;j++) {
      //for (k=j+1;k<npeak;k++) {
	/* calculate side lengths */
	/*la[0]=hypot((double) peakarray[3*i+1]-(double) peakarray[3*j+1],(double) peakarray[3*i+2]-(double) peakarray[3*j+2]);
	la[1]=hypot((double) peakarray[3*j+1]-(double) peakarray[3*k+1],(double) peakarray[3*j+2]-(double) peakarray[3*k+2]);
	la[2]=hypot((double) peakarray[3*i+1]-(double) peakarray[3*k+1],(double) peakarray[3*i+2]-(double) peakarray[3*k+2]);*/
	/* sort them */
	//qsort((void *) la,3,sizeof(double),dcomp);
	
	/* if the smallest side is less than a pixel, skip this triangle */
	//if (la[2]<1) break;
	
	/* calculate the ratio */
	/*ratioarray[0]=la[1]/la[0];
	ratioarray[1]=la[2]/la[0];*/

	/* find all the triangles from the first (short) list that are within the DISTCUT */
	//res=kd_nearest_range(kd_triangle,ratioarray,DISTCUT);

	/* if there are some triangles, then tell us about them */
	/*if (kd_res_size(res)>0) {
	  while( !kd_res_end( res ) ) {
	    double diff, pos[2];*/
	    /* get the data and position of the current result item */
	    /*data = (char*) kd_res_item( res, pos );
	    diff=hypot(pos[0]-ratioarray[0],pos[1]-ratioarray[1]);
	    printf("# diff= %g\n",diff);
	    printf("# %d %d %d -> %d %d %d\n",i,j,k,data[0],data[1],data[2]);*/

	    /* output some information about the triangle and calculate the      */
	    /* affine transformation --- will be in the glob_param array         */

	   // triangleoutput(i,j,k,data[0],data[1],data[2]);

	    /* set the counters to npeak to exit the loop because we have found a */
	    /* matching triangle --- we may have to do better here (find several  */
	    /* triangles with similar affine transformations like in k-d match    */
	       
	    //i=j=k=npeak;
	    /* go to the next entry */
	    /*kd_res_next( res );
	  }
	}*/
	/* free the results structure */
	/*kd_res_free(res);
      }
    }
  }*/
  
  /* free the triangle tree, because we are done with it */
  //kd_free(kd_triangle);

  /* free the catalogue data array; only needed this to make the triangle catalogue  */
  /*free((void *) catarray);

  cpu_time_used = ((double) (second_time.tms_utime - first_time.tms_utime)) / sysconf(_SC_CLK_TCK);

  //printf("Triangle matching %ld %g\n",(long) (second_time.tms_utime - first_time.tms_utime),cpu_time_used);*/
#endif
  
 
    /* now let's do the aperture photometry on each peak and see if there are all in the catalogue */
  times(&first_time);

  /* get me all of the peaks */
 
    
  jres=jkd_nearest_range(jkd_peak,key,2.0*imagesize);
    
    
  /* let's go through the peaks */
  j=i=0;
    
  /* This file will contains the Magnitude and coordinates information of the matching stars found  */

 
  FILE *FpeakMagValues, *FMissingStarsValues, *FImageBox, *FPeakStatObj;
#if 0    
  FpeakMagValues      = fopen("FpeakMagValues.txt", "w");
    
  // This file will contains information of the Image stars which are not contained in the Catalog
  FMissingStarsValues = fopen("FMissingStarsValues.txt", "w");
    
  // This file contains the information about the Box Image : xBox yBox --> To be deleted
  FImageBox = fopen("FImageBox.txt", "w");
    
  // This file contains the information about the Pixels statistics of the objects in the image
  FPeakStatObj = fopen("PeakStatObj.txt", "w");
    
#else
  FpeakMagValues      = fopen("/Users/Elisa/c/Files/FpeakMagValues.txt", "w");
    
  // This file will contains information of the Image stars which are not contained in the Catalog
  FMissingStarsValues = fopen("/Users/Elisa/c/Files/FMissingStarsValues.txt", "w");
    
  // This file contains the information about the Box Image : xBox yBox
  FImageBox = fopen("/Users/Elisa/c/Files/FImageBox.txt", "w");
    
  // This file contains the information about the Pixels statistics of the objects in the image
  FPeakStatObj = fopen("/Users/Elisa/c/Files/PeakStatObj.txt", "w");
    
#endif
    
    
/* Output the corners of the image box */

/*
//first Corner (0,0)

    box_image_x[0] = 0*glob_param[0]+0*glob_param[1]+glob_param[2];
    box_image_y[0] = 0*glob_param[3]+0*glob_param[4]+glob_param[5];
//Second Corner (0,ny)
    box_image_x[1] = 0*glob_param[0]+ny*glob_param[1]+glob_param[2];
    box_image_y[1] = 0*glob_param[3]+ny*glob_param[4]+glob_param[5];
//Third Corner  (nx,ny)
    box_image_x[2] = nx*glob_param[0]+ny*glob_param[1]+glob_param[2];
    box_image_y[2] = nx*glob_param[3]+ny*glob_param[4]+glob_param[5];
//Fourth Corner (nx,0)
    box_image_x[3] = nx*glob_param[0]+0*glob_param[1]+glob_param[2];
    box_image_y[3] = nx*glob_param[3]+0*glob_param[4]+glob_param[5];
    
   // printf("Image_Box = (%g,%g)-(%g,%g)-(%g,%g)-(%g,%g) \n",box_image_x[0],box_image_y[0],box_image_x[1],box_image_y[1],box_image_x[2],box_image_y[2],box_image_x[3],box_image_y[3]);

    for (j=0;j<4;j++) {
        fprintf(FImageBox,"%g %g \n",box_image_x[j],box_image_y[j]);
    }*/
    
    double flux;
    double SumX, x_xp; /* Sum(f(x-xp)), x-xp */
    double SumY, y_yp; /* Sum(f(y-yp)), y-yp */
    double SumX2, x_xp2 ; /* Sum(f(x-xp)^2), (x_xp)^2 */
    double SumY2, y_yp2 ; /* Sum(f(y-yp)^2), (y_yp)^2 */
    double SumXY, x_y; /* Sum(f(x-xp)*(y-yp)), (x-xp)*(y_yp) */
    
    double xbar,ybar, x2bar, y2bar, xybar;
    double stat1, stat2, Size;
    
    
    while( !jkd_res_end( jres ) ) {
        pix_t *data, top_pixel;
        double catkey[2];
        double allpixel;
      
        i++;
    /* get the pixel value and position of the current peak item */
        top_pixel= jkd_res_item( jres, key );
    /* let's sum over the nearby pixels */
        allpixel=0;
        SumX = 0;
        SumY = 0;
        SumX2 = 0;
        SumY2 = 0;
        SumXY = 0;
        npixel=0;
      
    /* find the nearby pixels from the tree */
        pix_res=jkd_nearest_range(jkd,key,INTEGRATEREGION);
    /* loop over the nearby pixels from the tree */
        while (!jkd_res_end( pix_res ) ) {
            allpixel += flux =jkd_res_item( pix_res, pix_key);
            npixel++;
            jkd_res_next( pix_res );
        
            printf("%g \n",allpixel);
    
      /* Compute the Statistics */
        
            x_xp = pix_key[0]-key[0];
            y_yp = pix_key[1]-key[1];
            x_xp2 = x_xp * x_xp;
            y_yp2 = y_yp * y_yp;
            x_y  = (pix_key[0]-key[0])*(pix_key[1]-key[1]);*/
        
        
            printf("x-xp = %g \n",x_xp);
            printf("y-yp = %g \n",y_yp);
            
            SumX  += x_xp*flux;
            SumX2 += flux*x_xp2;
            SumY  += y_yp*flux;
            SumY2 += flux*y_yp2;
            SumXY += flux*x_y;
      
        
        
            }
      /* free the list of neighbouring pixels */
        jkd_res_free(pix_res);
    
      /*printf(" Sum(f(x-xp))= %g \n",SumX);
      printf(" Sum(f(y-yp))= %g \n",SumY);
      printf(" Sum(f(x-xp)^2)= %g \n",SumX2);
      printf(" Sum(f(y-yp)^2)= %g \n",SumY2);
      printf(" Sum(f((x-xp)*(y-yp)))= %g \n",SumXY);
      printf(" allpixels   = %g \n",allpixel);*/
      
      xbar = (SumX/allpixel);
      ybar = (SumY/allpixel);
      x2bar = (SumX2-(allpixel * xbar*xbar))/allpixel;
      y2bar = (SumY2-(allpixel * ybar*ybar))/allpixel;
      xybar = (SumXY - allpixel*(xbar * ybar))/allpixel;
      
    /* printf(" Sum(f(x-xp))/Sum(f) = %g \n",xbar);
     printf(" Sum(f(y-yp))/Sum(f) = %g \n",ybar);
     
      printf(" x2bar = %g \n",x2bar);
      printf(" y2bar = %g \n",y2bar);
      printf(" xybar = %g \n",xybar);*/
      
      if(npixel > 4){
          stat1 = (x2bar-y2bar)/(x2bar+y2bar);
          stat2 = (xybar/(x2bar+y2bar));
          Size =  (x2bar+y2bar);
          fprintf(FPeakStatObj, "%g %g %g %g %g %g %g %d \n",key[0],key[1],allpixel,top_pixel,Size,stat1,stat2,npixel);
      }
     // printf("%g %g %g \n",Size,stat1,stat2);
      
      
      


    /* find the position of the peak in the coordinate system of the catalogue */
    /* the affine transformation between the image coordinates (key[0], key[1]) */
    /* and the catalogue coordinates (catkey[0],catkey[1]) is contained in */
    /* glob_param[0..5] */
      
    /*catkey[0]=key[0]*glob_param[0]+key[1]*glob_param[1]+glob_param[2]; // x coord
    catkey[1]=key[0]*glob_param[3]+key[1]*glob_param[4]+glob_param[5]; // y coord*/
      
      

    
      
      
      
    
      
    //printf("%d %g %g %g %g %g %d\n",i,top_pixel,key[0],key[1],catkey[0],catkey[1],allpixel);

    /* look for the peak in the catalogue, using the catalogue coordinates*/
   // pix_res=jkd_nearest_range(jkd_cat,catkey,SEARCHRADIUS);
      
    /* is it missing from the catalogue? */
   /* if (jkd_res_end(pix_res)) {
      j++;
      if (allpixel<HOTPIXELTHRESHOLD*top_pixel) { 
          printf("Likely hot pixel %g %g %g %g\n",key[0],key[1],catkey[0],catkey[1]);
      } else {
	printf("Not found %g %g %g %g \n",key[0],key[1],catkey[0],catkey[1]);
          
	// allpixels = flux image, top_pixel = peak value, catkey = catalog coordinates, key = Image coordinates
          fprintf(FMissingStarsValues, "%g %g %g %g %g %g %g %g %g \n",allpixel,top_pixel,catkey[0],catkey[1],key[0],key[1],Size,stat1,stat2);
      }
    } else {
        pix_t val=jkd_res_item( pix_res, pix_key);*/
        
        /*print out :CatMag,ImageMag (of the pixels around the peak),PeakMag, xCatCoord,yCatCoord,xImageCoord, yImageCoord, pix_key[0] pix_key[1] coordinates from the catalogue*/
        
       /* fprintf(FpeakMagValues, "%g %g %g %g %g %g %g %g %g %g %g %g \n",val,allpixel,top_pixel,catkey[0],catkey[1],key[0],key[1],pix_key[0],pix_key[1],Size,stat1,stat2);
        
      
    }*/
    /* go to the next peak */
    //jkd_res_next( jres );
  //}

  /* close the output files */
 /* fclose(FpeakMagValues);
  fclose(FMissingStarsValues);*/

  /* free the peak result list */
  //jkd_res_free(jres);

  /* free the trees */
  /*jkd_free(jkd);
  jkd_free(jkd_cat);
  jkd_free(jkd_peak);

  times(&second_time);
  
  cpu_time_used = ((double) (second_time.tms_utime - first_time.tms_utime)) / sysconf(_SC_CLK_TCK);

  printf("Query the peak tree, photometry and search %ld %g\n",(long) (second_time.tms_utime - first_time.tms_utime),cpu_time_used);
  printf("Total time (found %d peaks; %d not in catalogue, %d inserted from catalogue) %ld %g\n",
	 i,j,nstars,(long) (second_time.tms_utime - start_time.tms_utime),((double) (second_time.tms_utime - start_time.tms_utime)) / sysconf(_SC_CLK_TCK));*/


  /* I'm done! */
  return 0;

}
