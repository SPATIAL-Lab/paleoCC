#include <cmath>
#include "nr.h"
using namespace std;

extern DP dxsav;
extern int kmax,kount;
extern Vec_DP *xp_p;
extern Mat_DP *yp_p;
extern double exprod, lyso, calburial, ioratio, caldiss, c13diss, ioratio,
		orgburial, cpratio, swtemp, dwtemp, inject, sco3, dco3, initsco3, initdco3, fracsa,
		fpoc, fpocburial;

void NR::odeint(Vec_IO_DP &ystart, const DP x1, const DP x2, const DP eps,
	const DP h1, const DP hmin, int &nok, int &nbad,
	void derivs(const DP, Vec_I_DP &, Vec_IO_DP &),
	void rkqs(Vec_IO_DP &, Vec_IO_DP &, DP &, const DP, const DP,
	Vec_I_DP &, DP &, DP &, void (*)(const DP, Vec_I_DP &, Vec_IO_DP &)))
{
	const int MAXSTP=1e7;
	const DP TINY=1.0e-30;
	int i,nstp;
	DP xsav,x,hnext,hdid,h;

	int nvar=ystart.size();
	Vec_DP yscal(nvar),y(nvar),dydx(nvar);
	Vec_DP &xp=*xp_p;
	Mat_DP &yp=*yp_p;
	x=x1;
	h=SIGN(h1,x2-x1);
	nok = nbad = kount = 0;
	for (i=0;i<nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=0;nstp<MAXSTP;nstp++) {
		derivs(x,y,dydx);
		for (i=0;i<nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			for (i=0;i<nvar;i++) yp[i][kount]=y[i];
                        yp[nvar][kount] = inject;
                        yp[nvar+1][kount] = swtemp;
                        yp[nvar+2][kount] = dwtemp;
                        yp[nvar+3][kount] = exprod;
                        yp[nvar+4][kount] = lyso;
                        yp[nvar+5][kount] = calburial;
                        yp[nvar+6][kount] = caldiss;
                        yp[nvar+7][kount] = c13diss;
                        yp[nvar+8][kount] = orgburial * exprod + fpocburial * fpoc;
                        yp[nvar+9][kount] = orgburial * exprod / cpratio;
                        yp[nvar+10][kount] = sco3 - initsco3;
                        yp[nvar+11][kount] = dco3 - initdco3;
                        yp[nvar+12][kount] = fracsa;
       			xp[kount++]=x;
			xsav=x;
                        cout << x << "\n";
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs);
		if (hdid == h) ++nok; else ++nbad;
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=0;i<nvar;i++) ystart[i]=y[i];
			if (kmax != 0) {
				for (i=0;i<nvar;i++) yp[i][kount]=y[i];
				xp[kount++]=x;
			}
			return;
		}
		if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
		h=hnext;
	}
	nrerror("Too many steps in routine odeint");
}
