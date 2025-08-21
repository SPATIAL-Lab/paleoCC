#ifndef _NR_H_
#define _NR_H_
#include <fstream>
#include <complex>
#include "nrutil.h"
#include "nrtypes.h"
using namespace std;

namespace NR {

void odeint(Vec_IO_DP &ystart, const DP x1, const DP x2, const DP eps,
	const DP h1, const DP hmin, int &nok, int &nbad,
	void derivs(const DP, Vec_I_DP &, Vec_IO_DP &),
	void rkqs(Vec_IO_DP &, Vec_IO_DP &, DP &, const DP, const DP,
	Vec_I_DP &, DP &, DP &, void (*)(const DP, Vec_I_DP &, Vec_IO_DP &)));
void rkck(Vec_I_DP &y, Vec_I_DP &dydx, const DP x,
	const DP h, Vec_O_DP &yout, Vec_O_DP &yerr,
	void derivs(const DP, Vec_I_DP &, Vec_IO_DP &));
void rkqs(Vec_IO_DP &y, Vec_IO_DP &dydx, DP &x, const DP htry,
	const DP eps, Vec_I_DP &yscal, DP &hdid, DP &hnext,
	void derivs(const DP, Vec_I_DP &, Vec_IO_DP &));
}
#endif
