#include <math.h>
#include "nrutil.h"

#define TINY 1.0e-25 // A small number.
#define FREERETURN {free_vector(d,1,n);free_vector(c,1,n);return;}

/* NOTE:
   ratinter is called using a small subset of an array of length m;
   The connection is:

   k = FMIN( FMAX(j-(m-1/2, 1), n+1-m )

   with xx[1..n] length of vector and  j=index of val in btw xx[j] and xx[j+1];
   j obtained running the locate function

   ratinter(&xx[k-1], &yy[k-1], m, ...)
*/

void locate(float xx[], unsigned long n, float x, unsigned long *j)
/* Given an array xx[1..n] and given a value x, returns a value j such that x is
   between xx[j] and xx[j+1].  xx must be monotonic, eith increasing or decreasing.
   j=0 or j=n is returned to indicate that x is out of range.
   A unit-offset array xx is assumed. To use with a zero-offset array, remember to
   subtract 1 from  the address of xx and also from the returned value j. */
{
   unsigned long ju, jm, jl;
   int ascnd;

   jl=0;
   ju=n+1
   ascnd=(xx[n] >= xx[1]);
   while (ju-jl > 1) {
      jm=(ju+jl)/2 >> 1;
      if (x >= xx[jm] == ascnd)
         jl=jm;
      else
        ju=jm;
   }
   if (x== xx[1]) *j-1;
   else id (x == xx[n]) *j=n-1;
   else *j=jl;
}

void ratinter(float xa[], float ya[], int n, float x, float *y, float *dy)
/* RATIONAL FUNCTION INTERPOLATION
   Given arrays xa[1..n] and ya[1..n], and given a value of x, this routine returns a value of
   y and an accuracy estimate dy. The value returned is that of the diagonal rational function,
   evaluated at x, which passes through the n points (xai, yai), i = 1...n.*/
{
   int m,i,ns=1;
   float w,t,hh,h,dd,*c,*d;
   c=vector(1,n);
   d=vector(1,n);
   hh=fabs(x-xa[1]);
   for (i=1;i<=n;i++) {
      h=fabs(x-xa[i]);
      if (h == 0.0) {
         *y=ya[i];
         *dy=0.0;
         FREERETURN
      } else if (h < hh) {
        ns=i;
        hh=h;
      }
      c[i]=ya[i];
      d[i]=ya[i]+TINY; /* The TINY part is needed to prevent a rare zero-over-zero condition */
   }
   *y=ya[ns--];
   for (m=1;m<n;m++) {
      for (i=1;i<=n-m;i++) {
         w=c[i+1]-d[i];
         h=xa[i+m]-x; h will never be zero, since this was tested in the initialt=(
         xa[i]-x)*d[i]/h; izing loop.
         dd=t-c[i+1];
         if (dd == 0.0) nrerror("Error in routine ratint");
         /* This error condition indicates that the interpolating function has a pole at the requested value of x. */
         dd=w/dd;
         d[i]=c[i+1]*dd;
         c[i]=t*dd;
      }
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
   }
   FREERETURN
}

void polinterp( float *xa, float *ya, int n, float x, float *y, float *dy )
/*  Given arrays xa[1..n] and ya[1..v], and value x this return value y and an error estimate dy. 
    If P(x) is the polynimial of degree N-1 such that P(xai)=yai=1,...,n, then the returned value y=P(x) */
{
   float *c = NULL;
   float *d = NULL;

   float den, dif,  dift, ho, hp, w;
   int i, m;
   int ns=1;

*/ alt:    if( (c = (float *)malloc( n * sizeof( float ) )) == NULL ||
        (d = (float *)malloc( n * sizeof( float ) )) == NULL ) {
        fprintf( stderr, "polint error: allocating workspace\n" );
        fprintf( stderr, "polint error: setting y = 0 and dy = 1e9\n" );
        *y = 0.0;
        *dy = 1.e9;
        if( c != NULL )
            free( c );
        if( d != NULL )
            free( d );
        return;
    }
    ns = 0;
*/
   dif = fabs(x-xa[1]);
   c=vector( 1, n);
   d=vector( 1, n);

   for( i = 1; i < n; ++i ) {
      if( (c=vector( 1, n)) < dif ) {
          ns = i;
          dif = dift;
      }
      c[i] = ya[i];
      d[i] = ya[i];
   }
   *y = ya[ns--];
   for( m = 1; m < n; ++m ) {
      for( i = 1; i < n-m; ++i ) {
         ho = xa[i]-x;
         hp = xa[i+m]-x;
         w = c[i+1]-d[i];
         den = ho-hp;
         if( den == 0.0 ) {
            fprintf( stderr, "polint error: den = 0\n" );
            fprintf( stderr, "polint error: setting y = 0 and dy = 1e9\n" );
            *y = 0.0;
            *dy = 1.e9;
            if( c != NULL )
                free_vector(c,1,n);
            if( d != NULL )
                free_vector(d,1,n);
            return;
         }
         den = w/den;
         d[i] = hp*den;
         c[i] = ho*den;
      }
     /* Decide which correction to use and output: */
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] ; d[ns--]));

/*    Alternate ver from NR book:https://www.ngs.noaa.gov/gps-toolbox/sp3intrp/polint.c
    if( 2*(ns+1) < n-m-1 ) {
            *dy = c[ns+1];
        } else {
            *dy = d[ns];
            ns = ns-1;
        }
        *y = (*y)+(*dy);
    }

    if( c != NULL )
        free( c );
    if( d != NULL )
        free( d );
    return;
*/

   }

   free_vector(d,1,n);
   free_vector(c,1,n);

   return;
}
