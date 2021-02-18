#include "real.h"
#include <math.h>
#include <stdio.h>
#include <omp.h>

#define ind(a, b, c) (c) + ((twonby2p1) * ((b) + (((ng) * (a)))))

int displace_2lpt(int np, int ng, real boxsize, real *psi1, real g1, real gdot1,real *psi2, real g2, real gdot2, int direction) {
	long np3 = ((long) np) * np * np;
	long ng3 = ((long) ng) * ng * ng;
	long twonby2p1 = 2 * (ng / 2 + 1);
	real dxp = boxsize / ((real) np); 
	real rscale = ((real) ng) / boxsize;
	real qoffset = 0.5 / rscale;
	if (direction == 0) {
		hh = xx;
		hv = vx;
	}
	if (direction == 1) {
		hh = yy;
		hv = vy;
	}
	if (direction == 2) {
		hh = zz;
		hv = vz;
	}
	#pragma omp parallel
	{
		long iq, jq, kq, ip, jp, kp, ip1, jp1, kp1;
		real xi, yi, zi, xo, yo, zo, psi1x, psi2x, qxx, qyy, qzz, qq;
		#pragma omp for 
		for (long i = 0; i < np3; i++) {
			iq = i / np / np;
			jq = (i / np) % np;
			kq = i % np;
			qxx = iq * dxp;
			qyy = jq * dxp;
			qzz = kq * dxp;
			switch (direction) {
				case 0:
					qq = qxx;
					break;
				case 1:
					qq = qyy;
					break;
				case 2:
					qq = qzz;
			}
			if (ng == np) {
				ip = iq;
				jp = jq;
				kp = kq;
			}
			else {
				ip = floor(qxx * rscale);
				jp = floor(qyy * rscale);
				kp = floor(qzz * rscale);
			}
			xi = 1. - (qxx * rscale - (real) ip);
			yi = 1. - (qyy * rscale - (real) jp);
			zi = 1. - (qzz * rscale - (real) kp);
			ip = (ip + ng) % ng;
			jp = (jp + ng) % ng;
			kp = (kp + ng) % ng;
			ip1 = (ip - 1 + ng) % ng;
			jp1 = (jp - 1 + ng) % ng;
			kp1 = (kp - 1 + ng) % ng;
			xo = 1 - xi;
			yo = 1 - yi;
			zo = 1 - zi;
			psi1x = xi * yi * zi * psi1[ind(ip, jp, kp)] +
					xo * yi * zi * psi1[ind(ip, jp, kp1)] +
					xi * yo * zi * psi1[ind(ip, jp1, kp)] +
					xo * yo * zi * psi1[ind(ip, jp1, kp1)] +
					xi * yi * zo * psi1[ind(ip1, jp, kp)] +
					xo * yi * zo * psi1[ind(ip1, jp, kp1)] +
					xi * yo * zo * psi1[ind(ip1, jp1, kp)] +
					xo * yo * zo * psi1[ind(ip1, jp1, kp1)];
			psi2x = xi * yi * zi * psi2[ind(ip, jp, kp)] +
					xo * yi * zi * psi2[ind(ip, jp, kp1)] +
					xi * yo * zi * psi2[ind(ip, jp1, kp)] +
					xo * yo * zi * psi2[ind(ip, jp1, kp1)] +
					xi * yi * zo * psi2[ind(ip1, jp, kp)] +
					xo * yi * zo * psi2[ind(ip1, jp, kp1)] +
					xi * yo * zo * psi2[ind(ip1, jp1, kp)] +
					xo * yo * zo * psi2[ind(ip1, jp1, kp1)];
			hh[i] = qq + qoffset + g1 * psi1x + g2 * psi2x;
			if (hh[i] < 0.)
				hh[i] += boxsize;
			else if (hh[i] > boxsize)
				hh[i] -= boxsize;
			hv[i] = gdot1 * psi1x + gdot2 * psi2x;
		}
	}
	return 0;
}
