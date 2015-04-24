#include <math.h>

double bessi0(double x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}
/* (C) Copr. 1986-92 Numerical Recipes Software +@%R$s21-. */

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

double bessi_c_(int *ordp, double *xp)
{
	double bessi0(double x);
	int i, j;
	double bi,bim,bip,tox;
	double ans;
	double x = *xp;
	int ord = *ordp;

/*	if (ord < 2) fprintf(stderr,"Order n less than 2 in bessi"); */
	if (x == 0.0)
		ans = 0.0;
	else {
		tox=2.0/fabs(x);
		bip=0.0;
		ans=0.0;
		bi=1.0;
		for (j=2*(ord+(int) sqrt(ACC*ord));j>0;j--) {
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if (fabs(bi) > BIGNO) {
				ans *= BIGNI;
				bi *= BIGNI;
				bip *= BIGNI;
			}
			if (j == ord) ans=bip;
		}
		ans *= bessi0(x)/bi;
		if (x < 0.0 && (ord & 1)) ans *= -1;
	}
	return ans;
}
#undef ACC
#undef BIGNO
#undef BIGNI
/* (C) Copr. 1986-92 Numerical Recipes Software +@%R$s21-. */
