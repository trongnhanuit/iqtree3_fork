/*
 *
 * lbfgsb_new.cpp
 * HAL_HAS
 *
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2014, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * All rights reserved. CSIRO is willing to grant you a license to HAL-HAS on the terms of the GNU General Public
 * License version 3 as published by the Free Software Foundation (http://www.gnu.org/licenses/gpl.html), except
 * where otherwise indicated for third party material.
 * The following additional terms apply under clause 7 of that license:
 * EXCEPT AS EXPRESSLY STATED IN THIS AGREEMENT AND TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, THE SOFTWARE
 * IS PROVIDED "AS-IS". CSIRO MAKES NO REPRESENTATIONS, WARRANTIES OR CONDITIONS OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO ANY REPRESENTATIONS, WARRANTIES OR CONDITIONS REGARDING THE CONTENTS OR ACCURACY
 * OF THE SOFTWARE, OR OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, THE ABSENCE
 * OF LATENT OR OTHER DEFECTS, OR THE PRESENCE OR ABSENCE OF ERRORS, WHETHER OR NOT DISCOVERABLE.
 * TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT SHALL CSIRO BE LIABLE ON ANY LEGAL THEORY (INCLUDING,
 * WITHOUT LIMITATION, IN AN ACTION FOR BREACH OF CONTRACT, NEGLIGENCE OR OTHERWISE) FOR ANY CLAIM, LOSS, DAMAGES
 * OR OTHER LIABILITY HOWSOEVER INCURRED.  WITHOUT LIMITING THE SCOPE OF THE PREVIOUS SENTENCE THE EXCLUSION OF
 * LIABILITY SHALL INCLUDE: LOSS OF PRODUCTION OR OPERATION TIME, LOSS, DAMAGE OR CORRUPTION OF DATA OR RECORDS;
 * OR LOSS OF ANTICIPATED SAVINGS, OPPORTUNITY, REVENUE, PROFIT OR GOODWILL, OR OTHER ECONOMIC LOSS; OR ANY SPECIAL,
 * INCIDENTAL, INDIRECT, CONSEQUENTIAL, PUNITIVE OR EXEMPLARY DAMAGES, ARISING OUT OF OR IN CONNECTION WITH THIS
 * AGREEMENT, ACCESS OF THE SOFTWARE OR ANY OTHER DEALINGS WITH THE SOFTWARE, EVEN IF CSIRO HAS BEEN ADVISED OF
 * THE POSSIBILITY OF SUCH CLAIM, LOSS, DAMAGES OR OTHER LIABILITY.
 * APPLICABLE LEGISLATION SUCH AS THE AUSTRALIAN CONSUMER LAW MAY APPLY REPRESENTATIONS, WARRANTIES, OR CONDITIONS,
 * OR IMPOSES OBLIGATIONS OR LIABILITY ON CSIRO THAT CANNOT BE EXCLUDED, RESTRICTED OR MODIFIED TO THE FULL EXTENT
 * SET OUT IN THE EXPRESS TERMS OF THIS CLAUSE ABOVE "CONSUMER GUARANTEES".  TO THE EXTENT THAT SUCH CONSUMER
 * GUARANTEES CONTINUE TO APPLY, THEN TO THE FULL EXTENT PERMITTED BY THE APPLICABLE LEGISLATION, THE LIABILITY
 * OF CSIRO UNDER THE RELEVANT CONSUMER GUARANTEE IS LIMITED (WHERE PERMITTED AT CSIRO’S OPTION) TO ONE OF FOLLOWING
 * REMEDIES OR SUBSTANTIALLY EQUIVALENT REMEDIES:
 * (a)               THE REPLACEMENT OF THE SOFTWARE, THE SUPPLY OF EQUIVALENT SOFTWARE, OR SUPPLYING RELEVANT
 *                   SERVICES AGAIN;
 * (b)               THE REPAIR OF THE SOFTWARE;
 * (c)               THE PAYMENT OF THE COST OF REPLACING THE SOFTWARE, OF ACQUIRING EQUIVALENT SOFTWARE, HAVING THE
 *                   RELEVANT SERVICES SUPPLIED AGAIN, OR HAVING THE SOFTWARE REPAIRED.
 * IN THIS CLAUSE, CSIRO INCLUDES ANY THIRD PARTY AUTHOR OR OWNER OF ANY PART OF THE SOFTWARE OR MATERIAL DISTRIBUTED
 * WITH IT.  CSIRO MAY ENFORCE ANY RIGHTS ON BEHALF OF THE RELEVANT THIRD PARTY.
 * Third Party Components
 * The following third party components are distributed with the Software.  You agree to comply with the license
 * terms for these components as part of accessing the Software.  Other third party software may also be identified
 * in separate files distributed with the Software.
 * ___________________________________________________________________
 * 
 * R : A Computer Language for Statistical Data Analysis version 3.0.1 (http://cran.r-project.org/src/base/R-3/R-3.0.1.tar.gz)
 * Copyright (C) 2000-2004 The R Core Team
 * This software is licensed under GNU GPL
 * 
 * JACOBI_EIGENVALUE.C (http://people.sc.fsu.edu/~jburkardt/c_src/jacobi_eigenvalue/jacobi_eigenvalue.c)
 * Copyright (C) 2003-2013 John Burkardt
 * This software is licensed under GNU LGPL (http://www.gnu.org/licenses/lgpl.html)
 * ___________________________________________________________________
 */


#include "lbfgsb_new.h"
#include <algorithm>
//using namespace std;

static int c__1 = 1;
static int c__11 = 11;

#if 0

// Function to access the L-BFGS-B function
// 1. int n : The number of the variables
// 2. double* x : initial values of the variables
// 3. double* l : lower bounds of the variables
// 4. int maxit : max # of iterations
// 5. void* ex  : the wrapped variables for objective function
// After the function is invoked, the values of x will be updated
void lbfgsb_R(int n, double* x, double* l, int maxit, void* ex) {
	int i;
	double Fmin;
	int fail;
	int fncount;
	int grcount;
	char msg[100];

	int m = 5;          // number of BFGS updates retained in the "L-BFGS-B" method. It defaults to 5.

	double *u = nullptr;   // upper bounds of the variables;

	int *nbd;           // 0: unbounded; 1: lower bounded; 2: both lower & upper; 3: upper bounded
	nbd = new int[n];
	for (i=0; i<n; i++)
		nbd[i] = 1;

	double factr = 1e7; // control the convergence of the "L-BFGS-B" method.
	// Convergence occurs when the reduction in the object is within this factor
	// of the machine tolerance.
	// Default is 1e7, that is a tolerance of about 1e-8

	double pgtol = 0;   // helps control the convergence of the "L-BFGS-B" method.
	// It is a tolerance on the projected gradient in the current search direction.
	// Default is zero, when the check is suppressed

	int trace = 0;      // non-negative integer.
	// If positive, tracing information on the progress of the optimization is produced.
	// Higher values may produce more tracing information.

	int nREPORT = 10;   // The frequency of reports for the "L-BFGS-B" methods if "trace" is positive.
	// Defaults to every 10 iterations.

/*#ifdef USE_OLD_PARAM
	lbfgsb(n, m, x, l, u, nbd, &Fmin, fn, gr1, &fail, ex,
			factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
#else*/
	lbfgsb(n, m, x, l, u, nbd, &Fmin, fn, gr2, &fail, ex,
			factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
//#endif

	delete[] nbd;
}

// Function to access the L-BFGS-B function
// 1. int n : The number of the variables
// 2. double* x : initial values of the variables
// 3. double* l : lower bounds of the variables
// 4. double* u : upper bounds of the variables
// 5. int maxit : max # of iterations
// 6. void* ex  : the wrapped variables for objective function
// After the function is invoked, the values of x will be updated
void lbfgsb_R2(int n, double* x, double* l, double* u, int maxit, void* ex) {
	int i;
	double Fmin;
	int fail;
	int fncount;
	int grcount;
	char msg[100];

	int m = 5;          // number of BFGS updates retained in the "L-BFGS-B" method. It defaults to 5.

	int *nbd;           // 0: unbounded; 1: lower bounded; 2: both lower & upper; 3: upper bounded
	nbd = new int[n];
	for (i=0; i<n; i++)
		nbd[i] = 2;

	double factr = 1e7; // control the convergence of the "L-BFGS-B" method.
	// Convergence occurs when the reduction in the object is within this factor
	// of the machine tolerance.
	// Default is 1e7, that is a tolerance of about 1e-8

	double pgtol = 0;   // helps control the convergence of the "L-BFGS-B" method.
	// It is a tolerance on the projected gradient in the current search direction.
	// Default is zero, when the check is suppressed

	int trace = 0;      // non-negative integer.
	// If positive, tracing information on the progress of the optimization is produced.
	// Higher values may produce more tracing information.

	int nREPORT = 10;   // The frequency of reports for the "L-BFGS-B" methods if "trace" is positive.
	// Defaults to every 10 iterations.

/*#ifdef USE_OLD_PARAM
	lbfgsb(n, m, x, l, u, nbd, &Fmin, fn, gr1, &fail, ex,
			factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
#else*/
	lbfgsb(n, m, x, l, u, nbd, &Fmin, fn, gr2, &fail, ex,
			factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
//#endif

	delete[] nbd;
}

#endif

// ========================================================= //
// FUNCTIONS converted from R v.3.0.1
// ========================================================= //

void lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
		double *Fmin, optimfn fminfn, optimgr fmingr, int *fail,
		void *ex, double factr, double pgtol,
		int *fncount, int *grcount, int maxit, char *msg,
		int trace, int nREPORT)
{
	char task[60];
	double f, *g, dsave[29], *wa;
	int tr = -1, iter = 0, *iwa, isave[44], lsave[4];

	/* shut up gcc -Wall in 4.6.x */

	for(int i = 0; i < 4; i++) lsave[i] = 0;

	if(n == 0) { /* not handled in setulb */
		*fncount = 1;
		*grcount = 0;
		*Fmin = fminfn(n, u, ex);
		strcpy(msg, "NOTHING TO DO");
		*fail = 0;
		return;
	}
	if (nREPORT <= 0) {
		cerr << "REPORT must be > 0 (method = \"L-BFGS-B\")" << endl;
		exit(1);
	}
	switch(trace) {
	case 2: tr = 0; break;
	case 3: tr = nREPORT; break;
	case 4: tr = 99; break;
	case 5: tr = 100; break;
	case 6: tr = 101; break;
	default: tr = -1; break;
	}

	*fail = 0;
	g = (double*) malloc (n * sizeof(double));
	/* this needs to be zeroed for snd in mainlb to be zeroed */
	wa = (double *) malloc((2*m*n+4*n+11*m*m+8*m) * sizeof(double));
	iwa = (int *) malloc(3*n * sizeof(int));
	strcpy(task, "START");
	while(1) {
		/* Main workhorse setulb() from ../appl/lbfgsb.c : */
		setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task,
				tr, lsave, isave, dsave);
		/*    Rprintf("in lbfgsb - %s\n", task);*/
		if (strncmp(task, "FG", 2) == 0) {
			f = fminfn(n, x, ex);
			if (!isfinite(f)) {
				cerr << "L-BFGS-B needs finite values of 'fn'" << endl;
				exit(1);
			}
			fmingr(n, x, g, ex);
		} else if (strncmp(task, "NEW_X", 5) == 0) {
			iter++;
			if(trace == 1 && (iter % nREPORT == 0)) {
				cout << "iter " << iter << " value " << f << endl;
			}
			if (iter > maxit) {
				*fail = 1;
				break;
			}
		} else if (strncmp(task, "WARN", 4) == 0) {
			*fail = 51;
			break;
		} else if (strncmp(task, "CONV", 4) == 0) {
			break;
		} else if (strncmp(task, "ERROR", 5) == 0) {
			*fail = 52;
			break;
		} else { /* some other condition that is not supposed to happen */
			*fail = 52;
			break;
		}
	}
	*Fmin = f;
	*fncount = *grcount = isave[33];
	if (trace) {
		cout << "final value " << *Fmin << endl;
		if (iter < maxit && *fail == 0)
			cout << "converged" << endl;
		else
			cout << "stopped after " << iter << " iterations\n";
	}
	strcpy(msg, task);
	free(g);
	free(wa);
	free(iwa);
}

void setulb(int n, int m, double *x, double *l, double *u, int *nbd,
		double *f, double *g, double factr, double *pgtol,
		double *wa, int * iwa, char *task, int iprint,
		int *lsave, int *isave, double *dsave)
{
	/*     ************

	 Subroutine setulb

	 This subroutine partitions the working arrays wa and iwa, and
	 then uses the limited memory BFGS method to solve the bound
	 constrained optimization problem by calling mainlb.
	 (The direct method will be used in the subspace minimization.)

	 n is an integer variable.
	 On entry n is the dimension of the problem.
	 On exit n is unchanged.

	 m is an integer variable.
	 On entry m is the maximum number of variable metric corrections
	 used to define the limited memory matrix.
	 On exit m is unchanged.

	 x is a double precision array of dimension n.
	 On entry x is an approximation to the solution.
	 On exit x is the current approximation.

	 l is a double precision array of dimension n.
	 On entry l is the lower bound on x.
	 On exit l is unchanged.

	 u is a double precision array of dimension n.
	 On entry u is the upper bound on x.
	 On exit u is unchanged.

	 nbd is an integer array of dimension n.
	 On entry nbd represents the type of bounds imposed on the
	 variables, and must be specified as follows:
	 nbd(i)=0 if x(i) is unbounded,
	 1 if x(i) has only a lower bound,
	 2 if x(i) has both lower and upper bounds, and
	 3 if x(i) has only an upper bound.
	 On exit nbd is unchanged.

	 f is a double precision variable.
	 On first entry f is unspecified.
	 On final exit f is the value of the function at x.

	 g is a double precision array of dimension n.
	 On first entry g is unspecified.
	 On final exit g is the value of the gradient at x.

	 factr is a double precision variable.
	 On entry factr >= 0 is specified by the user.    The iteration
	 will stop when

	 (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch

	 where epsmch is the machine precision, which is automatically
	 generated by the code. Typical values for factr: 1.d+12 for
	 low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
	 high accuracy.
	 On exit factr is unchanged.

	 pgtol is a double precision variable.
	 On entry pgtol >= 0 is specified by the user.    The iteration
	 will stop when

	 max{|proj g_i | i = 1, ..., n} <= pgtol

	 where pg_i is the ith component of the projected gradient.
	 On exit pgtol is unchanged.

	 wa is a double precision working array of length
	 (2mmax + 4)nmax + 11mmax^2 + 8mmax.

	 iwa is an integer working array of length 3nmax.

	 task is a working string of characters of length 60 indicating
	 the current job when entering and quitting this subroutine.

	 iprint is an integer variable that must be set by the user.
	 It controls the frequency and type of output generated:
	 iprint<0    no output is generated;
	 iprint=0    print only one line at the last iteration;
	 0<iprint<99 print also f and |proj g| every iprint iterations;
	 iprint=99   print details of every iteration except n-vectors;
	 iprint=100  print also the changes of active set and final x;
	 iprint>100  print details of every iteration including x and g;
	 When iprint > 0, the file iterate.dat will be created to
	 summarize the iteration.

	 csave is a working string of characters of length 60.

	 lsave is a logical working array of dimension 4.
	 On exit with 'task' = NEW_X, the following information is
	 available:
	 If lsave(1) = .true. then  the initial X has been replaced by
	 its projection in the feasible set;
	 If lsave(2) = .true. then  the problem is constrained;
	 If lsave(3) = .true. then  each variable has upper and lower
	 bounds;

	 isave is an integer working array of dimension 44.
	 On exit with 'task' = NEW_X, the following information is
	 available:
	 isave(22) = the total number of intervals explored in the
	 search of Cauchy points;
	 isave(26) = the total number of skipped BFGS updates before
	 the current iteration;
	 isave(30) = the number of current iteration;
	 isave(31) = the total number of BFGS updates prior the current
	 iteration;
	 isave(33) = the number of intervals explored in the search of
	 Cauchy point in the current iteration;
	 isave(34) = the total number of function and gradient
	 evaluations;
	 isave(36) = the number of function value or gradient
	 evaluations in the current iteration;
	 if isave(37) = 0  then the subspace argmin is within the box;
	 if isave(37) = 1  then the subspace argmin is beyond the box;
	 isave(38) = the number of free variables in the current
	 iteration;
	 isave(39) = the number of active constraints in the current
	 iteration;
	 n + 1 - isave(40) = the number of variables leaving the set of
	 active constraints in the current iteration;
	 isave(41) = the number of variables entering the set of active
	 constraints in the current iteration.

	 dsave is a double precision working array of dimension 29.
	 On exit with 'task' = NEW_X, the following information is
	 available:
	 dsave(1) = current 'theta' in the BFGS matrix;
	 dsave(2) = f(x) in the previous iteration;
	 dsave(3) = factr*epsmch;
	 dsave(4) = 2-norm of the line search direction vector;
	 dsave(5) = the machine precision epsmch generated by the code;
	 dsave(7) = the accumulated time spent on searching for
	 Cauchy points;
	 dsave(8) = the accumulated time spent on
	 subspace minimization;
	 dsave(9) = the accumulated time spent on line search;
	 dsave(11) = the slope of the line search function at
	 the current point of line search;
	 dsave(12) = the maximum relative step length imposed in
	 line search;
	 dsave(13) = the infinity norm of the projected gradient;
	 dsave(14) = the relative step length in the line search;
	 dsave(15) = the slope of the line search function at
	 the starting point of the line search;
	 dsave(16) = the square of the 2-norm of the line search
	 direction vector.

	 Subprograms called:

	 L-BFGS-B Library ... mainlb.


	 References:

	 [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
	 memory algorithm for bound constrained optimization'',
	 SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

	 [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
	 limited memory FORTRAN code for solving bound constrained
	 optimization problems'', Tech. Report, NAM-11, EECS Department,
	 Northwestern University, 1994.

	 (Postscript files of these papers are available via anonymous
	 ftp to ece.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)

	 [Aug 2000: via http://www.ece.nwu.edu/~nocedal/lbfgsb.html]

	 *    *  *

	 NEOS, November 1994. (Latest revision April 1997.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	char csave[60];

	/* Local variables */
	int lsnd, ld, lr, lt;
	int lz, lwa, lwn, lss, lws, lwt, lsy, lwy;

	/* make sure csave is initialized */
	csave[0] = '\0';

	/* Parameter adjustments */
	--wa;
	--isave;

	/* Function Body */
	if (strncmp(task, "START", 5) == 0) {
		isave[1] = m * n;
		isave[2] = m * m;
		isave[3] = m * m << 2;
		isave[4] = 1;
		isave[5] = isave[4] + isave[1];
		isave[6] = isave[5] + isave[1];
		isave[7] = isave[6] + isave[2];
		isave[8] = isave[7] + isave[2];
		isave[9] = isave[8];
		isave[10] = isave[9] + isave[2];
		isave[11] = isave[10] + isave[3];
		isave[12] = isave[11] + isave[3];
		isave[13] = isave[12] + n;
		isave[14] = isave[13] + n;
		isave[15] = isave[14] + n;
		isave[16] = isave[15] + n;
	}
	lws = isave[4];
	lwy = isave[5];
	lsy = isave[6];
	lss = isave[7];
	lwt = isave[9];
	lwn = isave[10];
	lsnd = isave[11];
	lz = isave[12];
	lr = isave[13];
	ld = isave[14];
	lt = isave[15];
	lwa = isave[16];
	mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol,
			&wa[lws], &wa[lwy], &wa[lsy],&wa[lss], &wa[lwt],&wa[lwn],
			&wa[lsnd], &wa[lz], &wa[lr], &wa[ld], &wa[lt], &wa[lwa],
			iwa, &iwa[n], &iwa[n << 1], task, iprint,
			csave, lsave, &isave[22], dsave);
	return;
} /* setulb */


void mainlb(int n, int m, double *x,
		double *l, double *u, int *nbd, double *f, double *g,
		double factr, double *pgtol, double *ws, double * wy,
		double *sy, double *ss, double *wt, double *wn,
		double *snd, double *z, double *r, double *d,
		double *t, double *wa, int *indx, int *iwhere,
		int *indx2, char *task, int iprint,
		char *csave, int *lsave, int *isave, double *dsave)
{
	/*     ************
	 Subroutine mainlb

	 This subroutine solves bound constrained optimization problems by
	 using the compact formula of the limited memory BFGS updates.

	 n is an integer variable.
	 On entry n is the number of variables.
	 On exit n is unchanged.

	 m is an integer variable.
	 On entry m is the maximum number of variable metric
	 corrections allowed in the limited memory matrix.
	 On exit m is unchanged.

	 x is a double precision array of dimension n.
	 On entry x is an approximation to the solution.
	 On exit x is the current approximation.

	 l is a double precision array of dimension n.
	 On entry l is the lower bound of x.
	 On exit l is unchanged.

	 u is a double precision array of dimension n.
	 On entry u is the upper bound of x.
	 On exit u is unchanged.

	 nbd is an integer array of dimension n.
	 On entry nbd represents the type of bounds imposed on the
	 variables, and must be specified as follows:
	 nbd(i)=0 if x(i) is unbounded,
	 1 if x(i) has only a lower bound,
	 2 if x(i) has both lower and upper bounds,
	 3 if x(i) has only an upper bound.
	 On exit nbd is unchanged.

	 f is a double precision variable.
	 On first entry f is unspecified.
	 On final exit f is the value of the function at x.

	 g is a double precision array of dimension n.
	 On first entry g is unspecified.
	 On final exit g is the value of the gradient at x.

	 factr is a double precision variable.
	 On entry factr >= 0 is specified by the user.    The iteration
	 will stop when

	 (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch

	 where epsmch is the machine precision, which is automatically
	 generated by the code.
	 On exit factr is unchanged.

	 pgtol is a double precision variable.
	 On entry pgtol >= 0 is specified by the user.    The iteration
	 will stop when

	 max{|proj g_i | i = 1, ..., n} <= pgtol

	 where pg_i is the ith component of the projected gradient.
	 On exit pgtol is unchanged.

	 ws, wy, sy, and wt are double precision working arrays used to
	 store the following information defining the limited memory
	 BFGS matrix:
	 ws, of dimension n x m, stores S, the matrix of s-vectors;
	 wy, of dimension n x m, stores Y, the matrix of y-vectors;
	 sy, of dimension m x m, stores S'Y;
	 ss, of dimension m x m, stores S'S;
	 wt, of dimension m x m, stores the Cholesky factorization
	 of (theta*S'S+LD^(-1)L'); see eq.
	 (2.26) in [3].

	 wn is a double precision working array of dimension 2m x 2m
	 used to store the LEL^T factorization of the indefinite matrix
	 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'    ]
	 [L_a -R_z       theta*S'AA'S ]

	 where       E = [-I  0]
	 [ 0  I]

	 snd is a double precision working array of dimension 2m x 2m
	 used to store the lower triangular part of
	 N = [Y' ZZ'Y      L_a'+R_z']
	 [L_a +R_z  S'AA'S   ]

	 z(n),r(n),d(n),t(n),wa(8*m) are double precision working arrays.
	 z is used at different times to store the Cauchy point and
	 the Newton point.


	 indx is an integer working array of dimension n.
	 In subroutine freev, indx is used to store the free and fixed
	 variables at the Generalized Cauchy Point (GCP).

	 iwhere is an integer working array of dimension n used to record
	 the status of the vector x for GCP computation.
	 iwhere(i)=0 or -3 if x(i) is free and has bounds,
	 1       if x(i) is fixed at l(i), and l(i) .ne. u(i)
	 2       if x(i) is fixed at u(i), and u(i) .ne. l(i)
	 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
	 -1       if x(i) is always free, i.e., no bounds on it.

	 indx2 is an integer working array of dimension n.
	 Within subroutine cauchy, indx2 corresponds to the array iorder.
	 In subroutine freev, a list of variables entering and leaving
	 the free set is stored in indx2, and it is passed on to
	 subroutine formk with this information.

	 task is a working string of characters of length 60 indicating
	 the current job when entering and leaving this subroutine.

	 iprint is an INTEGER variable that must be set by the user.
	 It controls the frequency and type of output generated:
	 iprint<0    no output is generated;
	 iprint=0    print only one line at the last iteration;
	 0<iprint<99 print also f and |proj g| every iprint iterations;
	 iprint=99   print details of every iteration except n-vectors;
	 iprint=100  print also the changes of active set and final x;
	 iprint>100  print details of every iteration including x and g;
	 When iprint > 0, the file iterate.dat will be created to
	 summarize the iteration.

	 csave is a working string of characters of length 60.

	 lsave is a logical working array of dimension 4.

	 isave is an integer working array of dimension 23.

	 dsave is a double precision working array of dimension 29.


	 Subprograms called

	 L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk,

	 errclb, prn1lb, prn2lb, prn3lb, active, projgr,

	 freev, cmprlb, matupd, formt.

	 Minpack2 Library ... timer, dpmeps.

	 Linpack Library ... dcopy, ddot.


	 References:

	 [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
	 memory algorithm for bound constrained optimization'',
	 SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

	 [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
	 Subroutines for Large Scale Bound Constrained Optimization''
	 Tech. Report, NAM-11, EECS Department, Northwestern University,
	 1994.

	 [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
	 Quasi-Newton Matrices and their use in Limited Memory Methods'',
	 Mathematical Programming 63 (1994), no. 4, pp. 129-156.

	 (Postscript files of these papers are available via anonymous
	 ftp to ece.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)

	 *  *     *
	 */

	/*
	 NEOS, November 1994. (Latest revision April 1997.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */


	/* System generated locals */
	int ws_offset=0, wy_offset=0, sy_offset=0, ss_offset=0, wt_offset=0,
			wn_offset=0, snd_offset=0, i__1;
	double d__1, d__2;

	/* Local variables */
	int head;
	double fold;
	int nact;
	double ddum;
	int info;
	int nfgv, ifun, iter, nint;
	char word[4]; /* allow for terminator */
	int i, iback, k = 0; /* -Wall */
	double gdold;
	int nfree;
	int boxed;
	int itail;
	double theta;
	double dnorm;
	int nskip, iword;
	double xstep = 0.0, stpmx; /* xstep is printed before being used */
	double gd, dr, rr;
	int ileave;
	int itfile;
	double cachyt, epsmch;
	int updatd;
	double sbtime;
	int prjctd;
	int iupdat;
	int cnstnd;
	double sbgnrm;
	int nenter;
	double lnscht;
	int nintol;
	double dtd;
	int col;
	double tol;
	int wrk;
	double stp, cpu1, cpu2;

	/* Parameter adjustments */
	--indx2;
	--iwhere;
	--indx;
	--t;
	--d;
	--r;
	--z;
	--g;
	--nbd;
	--u;
	--l;
	--x;
	--wa;
	--lsave;
	--isave;
	--dsave;

	/* Function Body */
	if (strncmp(task, "START", 5) == 0) {
		/*      Generate the current machine precision. */
		epsmch = DBL_EPSILON;
		fold = 0.;
		dnorm = 0.;
		cpu1 = 0.;
		gd = 0.;
		sbgnrm = 0.;
		stp = 0.;
		xstep = 0.;
		stpmx = 0.;
		gdold = 0.;
		dtd = 0.;
		/*      Initialize counters and scalars when task='START'. */
		/*         for the limited memory BFGS matrices: */
		col = 0;
		head = 1;
		theta = 1.;
		iupdat = 0;
		updatd = 0;
		iback = 0;
		itail = 0;
		ifun = 0;
		iword = 0;
		nact = 0;
		ileave = 0;
		nenter = 0;
		/*         for operation counts: */
		iter = 0;
		nfgv = 0;
		nint = 0;
		nintol = 0;
		nskip = 0;
		nfree = n;
		/*         for stopping tolerance: */
		tol = factr * epsmch;
		/*         for measuring running time: */
		cachyt = 0.;
		sbtime = 0.;
		lnscht = 0.;
		/*         'word' records the status of subspace solutions. */
		strcpy(word, "---");
		/*         'info' records the termination information. */
		info = 0;
		itfile = 0;
		/*      Check the input arguments for errors. */
		errclb(n, m, factr, &l[1], &u[1], &nbd[1], task, &info, &k);
		if (strncmp(task, "ERROR", 5) == 0) {
			prn3lb(n, x+1, f, task, iprint, info,
					iter, nfgv, nintol, nskip, nact, sbgnrm,
					nint, word, iback, stp, xstep, k);
			return;
		}

		prn1lb(n, m, l+1, u+1, x+1, iprint, epsmch);

		/*      Initialize iwhere & project x onto the feasible set. */
		active(n, &l[1], &u[1], &nbd[1], &x[1], &iwhere[1], iprint, &prjctd,
				&cnstnd, &boxed);
		/*      The end of the initialization. */
	} else {
		/*        restore local variables. */
		prjctd = lsave[1];
		cnstnd = lsave[2];
		boxed = lsave[3];
		updatd = lsave[4];

		nintol = isave[1];
		itfile = isave[3];
		iback = isave[4];
		nskip = isave[5];
		head = isave[6];
		col = isave[7];
		itail = isave[8];
		iter = isave[9];
		iupdat = isave[10];
		nint = isave[12];
		nfgv = isave[13];
		info = isave[14];
		ifun = isave[15];
		iword = isave[16];
		nfree = isave[17];
		nact = isave[18];
		ileave = isave[19];
		nenter = isave[20];

		theta = dsave[1];
		fold = dsave[2];
		tol = dsave[3];
		dnorm = dsave[4];
		epsmch = dsave[5];
		cpu1 = dsave[6];
		cachyt = dsave[7];
		sbtime = dsave[8];
		lnscht = dsave[9];
		gd = dsave[11];
		stpmx = dsave[12];
		sbgnrm = dsave[13];
		stp = dsave[14];
		gdold = dsave[15];
		dtd = dsave[16];
		/*    After returning from the driver go to the point where execution */
		/*    is to resume. */
		if (strncmp(task, "FG_LN", 5) == 0)    goto L666;
		if (strncmp(task, "NEW_X", 5) == 0)     goto L777;
		if (strncmp(task, "FG_ST", 5) == 0)     goto L111;

		if (strncmp(task, "STOP", 4) == 0) {
			if (strncmp(task + 6, "CPU", 3) == 0) {
				// restore the previous iterate.
				dcopy(&n, &t[1], &c__1, &x[1], &c__1);
				dcopy(&n, &r[1], &c__1, &g[1], &c__1);
				*f = fold;
			}
			goto L999;
		}
	}
	/*     Compute f0 and g0. */
	strcpy(task, "FG_START");
	/*        return to the driver to calculate f and g; reenter at 111. */
	goto L1000;
	L111:
	nfgv = 1;
	/*     Compute the infinity norm of the (-) projected gradient. */
	projgr(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);

	if (iprint >= 1)
		cout << "At iterate " << iter << " f= " << *f << " |proj g|= " << sbgnrm << endl;

	if (sbgnrm <= *pgtol) {
		/*                  terminate the algorithm. */
		strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
		goto L999;
	}
	/* ----------------- the beginning of the loop -------------------------- */
	L222:
	if (iprint >= 99)
		cout << "Iteration " << iter << endl;
	iword = -1;

	if (! cnstnd && col > 0) {
		/*                          skip the search for GCP. */
		dcopy(&n, &x[1], &c__1, &z[1], &c__1);
		wrk = updatd;
		nint = 0;
		goto L333;
	}
	/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

	/*     Compute the Generalized Cauchy Point (GCP). */

	/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
	timer(&cpu1);
	cauchy(n, &x[1], &l[1], &u[1], &nbd[1], &g[1], &indx2[1], &iwhere[1], &t[
																			 1], &d[1], &z[1], m, &wy[wy_offset], &ws[ws_offset], &sy[
																																	  sy_offset], &wt[wt_offset], &theta, &col, &head, &wa[1], &wa[(m
																																			  << 1) + 1], &wa[(m << 2) + 1], &wa[m * 6 + 1], &nint, iprint, &
																																			  sbgnrm, &info, &epsmch);
	if (info != 0) {
		/*       singular triangular system detected; refresh the lbfgs memory. */
		if (iprint >= 1) {
			cout << "Singular triangular system detected;" << endl;
			cout << "   refresh the lbfgs memory and restart the iteration." << endl;
		}
		info = 0;
		col = 0;
		head = 1;
		theta = 1.;
		iupdat = 0;
		updatd = 0;
		timer(&cpu2);
		cachyt = cachyt + cpu2 - cpu1;
		goto L222;
	}
	timer(&cpu2);
	cachyt = cachyt + cpu2 - cpu1;
	nintol += nint;
	/*     Count the entering and leaving variables for iter > 0; */
	/*     find the index set of free and active variables at the GCP. */
	freev(n, &nfree, &indx[1], &nenter, &ileave, &indx2[1], &iwhere[1], &
			wrk, &updatd, &cnstnd, iprint, &iter);
	nact = n - nfree;
	L333:
	/*     If there are no free variables or B=theta*I, then */
	/*                      skip the subspace minimization. */
	if (nfree == 0 || col == 0) {
		goto L555;
	}
	/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

	/*     Subspace minimization. */

	/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
	timer(&cpu1);
	/*     Form  the LEL^T factorization of the indefinite */
	/*     matrix       K = [-D -Y'ZZ'Y/theta     L_a'-R_z'    ] */
	/*               [L_a -R_z       theta*S'AA'S ] */
	/*     where       E = [-I  0] */
	/*               [ 0  I] */
	if (wrk) {
		formk(n, &nfree, &indx[1], &nenter, &ileave, &indx2[1], &iupdat, &
				updatd, &wn[wn_offset], &snd[snd_offset], m, &ws[ws_offset], &
				wy[wy_offset], &sy[sy_offset], &theta, &col, &head, &info);
	}
	if (info != 0) {
		/*        nonpositive definiteness in Cholesky factorization; */
		/*        refresh the lbfgs memory and restart the iteration. */
		if (iprint >= 0) {
			cout << "Nonpositive definiteness in Cholesky factorization in formk;" << endl;
			cout << "   refresh the lbfgs memory and restart the iteration." << endl;
		}
		info = 0;
		col = 0;
		head = 1;
		theta = 1.;
		iupdat = 0;
		updatd = 0;
		timer(&cpu2);
		sbtime = sbtime + cpu2 - cpu1;
		goto L222;
	}
	/*      compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x) */
	/*                             from 'cauchy'). */
	cmprlb(n, m, &x[1], &g[1], &ws[ws_offset], &wy[wy_offset], &sy[sy_offset]
																   , &wt[wt_offset], &z[1], &r[1], &wa[1], &indx[1], &theta, &
																   col, &head, &nfree, &cnstnd, &info);
	if (info != 0) {
		goto L444;
	}
	/*     call the direct method. */
	subsm(n, m, &nfree, &indx[1], &l[1], &u[1], &nbd[1], &z[1], &r[1], &
			ws[ws_offset], &wy[wy_offset], &theta, &col, &head, &iword, &wa[1]
																			, &wn[wn_offset], iprint, &info);
	L444:
	if (info != 0) {
		/*        singular triangular system detected; */
		/*        refresh the lbfgs memory and restart the iteration. */
		if (iprint >= 1) {
			cout << "Singular triangular system detected;" << endl;
			cout << "   refresh the lbfgs memory and restart the iteration." << endl;
		}
		info = 0;
		col = 0;
		head = 1;
		theta = 1.;
		iupdat = 0;
		updatd = 0;
		timer(&cpu2);
		sbtime = sbtime + cpu2 - cpu1;
		goto L222;
	}
	timer(&cpu2);
	sbtime = sbtime + cpu2 - cpu1;
	L555:
	/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

	/*     Line search and optimality tests. */

	/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
	/*     Generate the search direction d:=z-x. */
	i__1 = n;
	for (i = 1; i <= i__1; ++i) {
		d[i] = z[i] - x[i];
		/* L40: */
	}
	timer(&cpu1);
	L666:
	lnsrlb(n, &l[1], &u[1], &nbd[1], &x[1], f, &fold, &gd, &gdold, &g[1], &
			d[1], &r[1], &t[1], &z[1], &stp, &dnorm, &dtd, &xstep, &
			stpmx, &iter, &ifun, &iback, &nfgv, &info, task, &boxed, &cnstnd,
			csave, &isave[22], &dsave[17]);
	if (info != 0 || iback >= 20) {
		/*        restore the previous iterate. */
		dcopy(&n, &t[1], &c__1, &x[1], &c__1);
		dcopy(&n, &r[1], &c__1, &g[1], &c__1);
		*f = fold;
		if (col == 0) {
			/*           abnormal termination. */
			if (info == 0) {
				info = -9;
				/*          restore the actual number of f and g evaluations etc. */
				--nfgv;
				--ifun;
				--iback;
			}
			strcpy(task, "WARNING: ABNORMAL_TERMINATION_IN_LNSRCH");
			++iter;
			goto L999;
		} else {
			/*           refresh the lbfgs memory and restart the iteration. */
			if (iprint >= 1) {
				cout << "Bad direction in the line search;" << endl;
				cout << "   refresh the lbfgs memory and restart the iteration." << endl;
			}
			if (info == 0) {
				--nfgv;
			}
			info = 0;
			col = 0;
			head = 1;
			theta = 1.;
			iupdat = 0;
			updatd = 0;
			strcpy(task, "RESTART_FROM_LNSRCH");
			timer(&cpu2);
			lnscht = lnscht + cpu2 - cpu1;
			goto L222;
		}
	} else if (strncmp(task, "FG_LN", 5) == 0) {
		/*        return to the driver for calculating f and g; reenter at 666. */
		goto L1000;
	} else {
		/*        calculate and print out the quantities related to the new X. */
		timer(&cpu2);
		lnscht = lnscht + cpu2 - cpu1;
		++iter;
		/*      Compute the infinity norm of the projected (-)gradient. */
		projgr(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
		/*      Print iteration information. */
		prn2lb(n, x+1, f, g+1, iprint, iter, nfgv, nact,
				sbgnrm, nint, word, iword, iback, stp, xstep);
		goto L1000;
	}
	L777:
	/*     Test for termination. */
	if (sbgnrm <= *pgtol) {
		/*                  terminate the algorithm. */
		strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
		goto L999;
	}
	/* Computing MAX */
	d__1 = fabs(fold), d__2 = fabs(*f), d__1 = max(d__1,d__2);
	ddum = max((double)d__1,(double)1.);
	if (fold - *f <= tol * ddum) {
		/*                      terminate the algorithm. */
		strcpy(task, "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH");
		if (iback >= 10) info = -5;
		/*         i.e., to issue a warning if iback>10 in the line search. */
		goto L999;
	}
	/*     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's. */
	i__1 = n;
	for (i = 1; i <= i__1; ++i) {
		r[i] = g[i] - r[i];
		/* L42: */
	}
	rr = ddot(&n, &r[1], &c__1, &r[1], &c__1);
	if (stp == 1.) {
		dr = gd - gdold;
		ddum = -gdold;
	} else {
		dr = (gd - gdold) * stp;
		dscal(&n, &stp, &d[1], &c__1);
		ddum = -gdold * stp;
	}
	if (dr <= epsmch * ddum) {
		/*                  skip the L-BFGS update. */
		++nskip;
		updatd = 0;
		if (iprint >= 1) {
			cout << "ys=" << dr << "   -gs=" << ddum << ", BFGS update SKIPPED" << endl;
		}
		goto L888;
	}
	/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

	/*     Update the L-BFGS matrix. */

	/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
	updatd = 1;
	++iupdat;
	/*     Update matrices WS and WY and form the middle matrix in B. */
	matupd(n, m, &ws[ws_offset], &wy[wy_offset], &sy[sy_offset], &ss[
																	 ss_offset], &d[1], &r[1], &itail, &iupdat, &col, &head, &
																	 theta, &rr, &dr, &stp, &dtd);
	/*     Form the upper half of the pds T = theta*SS + L*D^(-1)*L'; */
	/*      Store T in the upper triangular of the array wt; */
	/*      Cholesky factorize T to J*J' with */
	/*         J' stored in the upper triangular of wt. */
	formt(m, &wt[wt_offset], &sy[sy_offset], &ss[ss_offset], &col, &theta, &
			info);
	if (info != 0) {
		/*        nonpositive definiteness in Cholesky factorization; */
		/*        refresh the lbfgs memory and restart the iteration. */
		if (iprint >= 0) {
			cout << "Nonpositive definiteness in Cholesky factorization in formk;" << endl;
			cout << "   refresh the lbfgs memory and restart the iteration." << endl;
		}
		info = 0;
		col = 0;
		head = 1;
		theta = 1.;
		iupdat = 0;
		updatd = 0;
		goto L222;
	}
	/*     Now the inverse of the middle matrix in B is */
	/*     [  D^(1/2)     O ] [ -D^(1/2)     D^(-1/2)*L' ] */
	/*     [ -L*D^(-1/2)     J ] [    0     J'         ] */
	L888:
	/* -------------------- the end of the loop ----------------------------- */
	goto L222;
	L999:
	L1000:
	/*     Save local variables. */
	lsave[1] = prjctd;
	lsave[2] = cnstnd;
	lsave[3] = boxed;
	lsave[4] = updatd;
	isave[1] = nintol;
	isave[3] = itfile;
	isave[4] = iback;
	isave[5] = nskip;
	isave[6] = head;
	isave[7] = col;
	isave[8] = itail;
	isave[9] = iter;
	isave[10] = iupdat;
	isave[12] = nint;
	isave[13] = nfgv;
	isave[14] = info;
	isave[15] = ifun;
	isave[16] = iword;
	isave[17] = nfree;
	isave[18] = nact;
	isave[19] = ileave;
	isave[20] = nenter;
	dsave[1] = theta;
	dsave[2] = fold;
	dsave[3] = tol;
	dsave[4] = dnorm;
	dsave[5] = epsmch;
	dsave[6] = cpu1;
	dsave[7] = cachyt;
	dsave[8] = sbtime;
	dsave[9] = lnscht;
	dsave[11] = gd;
	dsave[12] = stpmx;
	dsave[13] = sbgnrm;
	dsave[14] = stp;
	dsave[15] = gdold;
	dsave[16] = dtd;
	prn3lb(n, x+1, f, task, iprint, info,
			iter, nfgv, nintol, nskip, nact, sbgnrm,
			nint, word, iback, stp, xstep, k);
	return;
} /* mainlb */
/* ======================= The end of mainlb ============================= */

void errclb(int n, int m, double factr, double *l, double *u,
		int *nbd, char *task, int *info, int *k)
{
	/*    ************
	 Subroutine errclb

	 This subroutine checks the validity of the input data.

	 *     *  *

	 NEOS, November 1994. (Latest revision April 1997.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	/* Local variables */
	int i;

	/* Parameter adjustments */
	--nbd;
	--u;
	--l;

	/* Function Body */
	/*     Check the input arguments for errors. */
	if (n <= 0)
		strcpy(task, "ERROR: N .LE. 0");
	if (m <= 0)
		strcpy(task, "ERROR: M .LE. 0");
	if (factr < 0.)
		strcpy(task, "ERROR: FACTR .LT. 0");

	/*     Check the validity of the arrays nbd(i), u(i), and l(i). */
	for (i = 1; i <= n; ++i) {
		if (nbd[i] < 0 || nbd[i] > 3) {
			/*                             return */
			strcpy(task, "ERROR: INVALID NBD");
			*info = -6;
			*k = i;
		}
		if (nbd[i] == 2) {
			if (l[i] > u[i]) {
				/*                      return */
				strcpy(task, "ERROR: NO FEASIBLE SOLUTION");
				*info = -7;
				*k = i;
			}
		}
	}
	return;
} /* errclb */
/* ======================= The end of errclb ============================= */

void active(int n, double *l, double *u,
		int *nbd, double *x, int *iwhere, int iprint,
		int *prjctd, int *cnstnd, int *boxed)
{
	/*    ************

	 Subroutine active

	 This subroutine initializes iwhere and projects the initial x to
	 the feasible set if necessary.

	 iwhere is an integer array of dimension n.
	 On entry iwhere is unspecified.
	 On exit iwhere(i)=-1    if x(i) has no bounds
	 3    if l(i)=u(i)
	 0    otherwise.
	 In cauchy, iwhere is given finer gradations.


	 *     *  *

	 NEOS, November 1994. (Latest revision June 1996.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	/* Local variables */
	int nbdd, i;

	/* Parameter adjustments */
	--iwhere;
	--x;
	--nbd;
	--u;
	--l;

	/* Function Body */

	/*Initialize nbdd, prjctd, cnstnd and boxed. */
	nbdd = 0;
	*prjctd = 0;
	*cnstnd = 0;
	*boxed = 1;
	/*     Project the initial x to the easible set if necessary. */
	for (i = 1; i <= n; ++i) {
		if (nbd[i] > 0) {
			if (nbd[i] <= 2 && x[i] <= l[i]) {
				if (x[i] < l[i]) {
					*prjctd = 1;
					x[i] = l[i];
				}
				++nbdd;
			} else if (nbd[i] >= 2 && x[i] >= u[i]) {
				if (x[i] > u[i]) {
					*prjctd = 1;
					x[i] = u[i];
				}
				++nbdd;
			}
		}
	}

	/*     Initialize iwhere and assign values to cnstnd and boxed. */
	for (i = 1; i <= n; ++i) {
		if (nbd[i] != 2) {
			*boxed = 0;
		}
		if (nbd[i] == 0) {
			/*                  this variable is always free */
			iwhere[i] = -1;
			/*         otherwise set x(i)=mid(x(i), u(i), l(i)). */
		} else {
			*cnstnd = 1;
			if (nbd[i] == 2 && u[i] - l[i] <= 0.) {
				/*             this variable is always fixed */
				iwhere[i] = 3;
			} else {
				iwhere[i] = 0;
			}
		}
	}
	if (iprint >= 0) {
		if (*prjctd)
			cout << "The initial X is infeasible.  Restart with its projection." << endl;
		if (!*cnstnd)
			cout << "This problem is unconstrained." << endl;
	}
	if (iprint > 0)
		cout << "At X0, " << nbdd << " variables are exactly at the bounds" << endl;

	return;
} /* active */
/* ======================= The end of active ============================= */

void cauchy(int n, double *x, double *l, double *u, int *nbd,
		double *g, int *iorder, int * iwhere, double *t,
		double *d, double *xcp, int m,
		double *wy, double *ws, double *sy, double *wt,
		double *theta, int *col, int *head, double *p,
		double *c, double *wbp, double *v, int *nint,
		int iprint, double *sbgnrm, int *info, double * epsmch)
{
	/*     ************
	 Subroutine cauchy

	 For given x, l, u, g (with sbgnrm > 0), and a limited memory
	 BFGS matrix B defined in terms of matrices WY, WS, WT, and
	 scalars head, col, and theta, this subroutine computes the
	 generalized Cauchy point (GCP), defined as the first local
	 minimizer of the quadratic

	 Q(x + s) = g's + 1/2 s'Bs

	 along the projected gradient direction P(x-tg,l,u).
	 The routine returns the GCP in xcp.

	 n is an integer variable.
	 On entry n is the dimension of the problem.
	 On exit n is unchanged.

	 x is a double precision array of dimension n.
	 On entry x is the starting point for the GCP computation.
	 On exit x is unchanged.

	 l is a double precision array of dimension n.
	 On entry l is the lower bound of x.
	 On exit l is unchanged.

	 u is a double precision array of dimension n.
	 On entry u is the upper bound of x.
	 On exit u is unchanged.

	 nbd is an integer array of dimension n.
	 On entry nbd represents the type of bounds imposed on the
	 variables, and must be specified as follows:
	 nbd(i)=0 if x(i) is unbounded,
	 1 if x(i) has only a lower bound,
	 2 if x(i) has both lower and upper bounds, and
	 3 if x(i) has only an upper bound.
	 On exit nbd is unchanged.

	 g is a double precision array of dimension n.
	 On entry g is the gradient of f(x).  g must be a nonzero vector.
	 On exit g is unchanged.

	 iorder is an integer working array of dimension n.
	 iorder will be used to store the breakpoints in the piecewise
	 linear path and free variables encountered. On exit,
	 iorder(1),...,iorder(nleft) are indices of breakpoints
	 which have not been encountered;
	 iorder(nleft+1),...,iorder(nbreak) are indices of
	 encountered breakpoints; and
	 iorder(nfree),...,iorder(n) are indices of variables which
	 have no bound constraits along the search direction.

	 iwhere is an integer array of dimension n.
	 On entry iwhere indicates only the permanently fixed (iwhere=3)
	 or free (iwhere= -1) components of x.
	 On exit iwhere records the status of the current x variables.
	 iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
	 0   if x(i) is free and has bounds, and is moved
	 1   if x(i) is fixed at l(i), and l(i) .ne. u(i)
	 2   if x(i) is fixed at u(i), and u(i) .ne. l(i)
	 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
	 -1  if x(i) is always free, i.e., it has no bounds.

	 t is a double precision working array of dimension n.
	 t will be used to store the break points.

	 d is a double precision array of dimension n used to store
	 the Cauchy direction P(x-tg)-x.

	 xcp is a double precision array of dimension n used to return the
	 GCP on exit.

	 m is an integer variable.
	 On entry m is the maximum number of variable metric corrections
	 used to define the limited memory matrix.
	 On exit m is unchanged.

	 ws, wy, sy, and wt are double precision arrays.
	 On entry they store information that defines the
	 limited memory BFGS matrix:
	 ws(n,m) stores S, a set of s-vectors;
	 wy(n,m) stores Y, a set of y-vectors;
	 sy(m,m) stores S'Y;
	 wt(m,m) stores the
	 Cholesky factorization of (theta*S'S+LD^(-1)L').
	 On exit these arrays are unchanged.

	 theta is a double precision variable.
	 On entry theta is the scaling factor specifying B_0 = theta I.
	 On exit theta is unchanged.

	 col is an integer variable.
	 On entry col is the actual number of variable metric
	 corrections stored so far.
	 On exit col is unchanged.

	 head is an integer variable.
	 On entry head is the location of the first s-vector
	 (or y-vector) in S (or Y).
	 On exit col is unchanged.

	 p is a double precision working array of dimension 2m.
	 p will be used to store the vector p = W^(T)d.

	 c is a double precision working array of dimension 2m.
	 c will be used to store the vector c = W^(T)(xcp-x).

	 wbp is a double precision working array of dimension 2m.
	 wbp will be used to store the row of W corresponding
	 to a breakpoint.

	 v is a double precision working array of dimension 2m.

	 nint is an integer variable.
	 On exit nint records the number of quadratic segments explored
	 in searching for the GCP.

	 iprint is an INTEGER variable that must be set by the user.
	 It controls the frequency and type of output generated:
	 iprint<0    no output is generated;
	 iprint=0    print only one line at the last iteration;
	 0<iprint<99 print also f and |proj g| every iprint iterations;
	 iprint=99   print details of every iteration except n-vectors;
	 iprint=100  print also the changes of active set and final x;
	 iprint>100  print details of every iteration including x and g;
	 When iprint > 0, the file iterate.dat will be created to
	 summarize the iteration.

	 sbgnrm is a double precision variable.
	 On entry sbgnrm is the norm of the projected gradient at x.
	 On exit sbgnrm is unchanged.

	 info is an integer variable.
	 On entry info is 0.
	 On exit info = 0    for normal return,
	 = nonzero for abnormal return when the the system
	 used in routine bmv is singular.

	 Subprograms called:

	 L-BFGS-B Library ... hpsolb, bmv.

	 Linpack ... dscal dcopy, daxpy.


	 References:

	 [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
	 memory algorithm for bound constrained optimization'',
	 SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

	 [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
	 Subroutines for Large Scale Bound Constrained Optimization''
	 Tech. Report, NAM-11, EECS Department, Northwestern University, 1994.

	 (Postscript files of these papers are available via anonymous
	 ftp to ece.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)

	 *    *  *

	 NEOS, November 1994. (Latest revision April 1997.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	/* System generated locals */
	int wy_dim1, wy_offset, ws_dim1, ws_offset, sy_dim1, sy_offset,
	wt_dim1, wt_offset, i__2;
	double d__1;

	/* Local variables */
	double bkmin, dibp, dibp2, zibp, neggi, tsum;
	double f1, f2, f2_org__, dt, tj, tj0, tl= 0.0, tu=0.0, dtm, wmc, wmp, wmw;

	int i, j, ibp, iter, bnded, nfree, nleft, nbreak, ibkmin, pointr;
	int xlower, xupper, col2;

	/* Parameter adjustments */
	--xcp;
	--d;
	--t;
	--iwhere;
	--iorder;
	--g;
	--nbd;
	--u;
	--l;
	--x;
	--v;
	--wbp;
	--c;
	--p;
	wt_dim1 = m;    wt_offset = 1 + wt_dim1 * 1;    wt -= wt_offset;
	sy_dim1 = m;    sy_offset = 1 + sy_dim1 * 1;    sy -= sy_offset;
	ws_dim1 = n;    ws_offset = 1 + ws_dim1 * 1;    ws -= ws_offset;
	wy_dim1 = n;    wy_offset = 1 + wy_dim1 * 1;    wy -= wy_offset;

	/* Function Body */

	/*     Check the status of the variables, reset iwhere(i) if necessary;
	 *     compute the Cauchy direction d and the breakpoints t; initialize
	 *     the derivative f1 and the vector p = W'd (for theta = 1).
	 */

	if (*sbgnrm <= 0.) {
		if (iprint >= 0) cout << "Subgnorm = 0.  GCP = X.\n";
		dcopy(&n, &x[1], &c__1, &xcp[1], &c__1);
		return;
	}
	bnded = 1;
	nfree = n + 1;
	nbreak = 0;
	ibkmin = 0;
	bkmin = 0.;
	col2 = *col << 1;
	f1 = 0.;
	if (iprint >= 99)
		cout << "\n---------------- CAUCHY entered-------------------\n\n";

	/*     We set p to zero and build it up as we determine d. */
	for (i = 1; i <= col2; ++i)
		p[i] = 0.;

	/*     In the following loop we determine for each variable its bound */
	/*      status and its breakpoint, and update p accordingly. */
	/*      Smallest breakpoint is identified. */

	for (i = 1; i <= n; ++i) {
		neggi = -g[i];
		if (iwhere[i] != 3 && iwhere[i] != -1) {
			/*           if x(i) is not a constant and has bounds, */
			/*           compute the difference between x(i) and its bounds. */
			if (nbd[i] <= 2) {
				tl = x[i] - l[i];
			}
			if (nbd[i] >= 2) {
				tu = u[i] - x[i];
			}
			/*         If a variable is close enough to a bound */
			/*           we treat it as at bound. */
			xlower = nbd[i] <= 2 && tl <= 0.;
			xupper = nbd[i] >= 2 && tu <= 0.;
			/*        reset iwhere(i). */
			iwhere[i] = 0;
			if (xlower) {
				if (neggi <= 0.) {
					iwhere[i] = 1;
				}
			} else if (xupper) {
				if (neggi >= 0.) {
					iwhere[i] = 2;
				}
			} else {
				if (fabs(neggi) <= 0.) {
					iwhere[i] = -3;
				}
			}
		}
		pointr = *head;
		if (iwhere[i] != 0 && iwhere[i] != -1) {
			d[i] = 0.;
		} else {
			d[i] = neggi;
			f1 -= neggi * neggi;
			/*           calculate p := p - W'e_i* (g_i). */
			i__2 = *col;
			for (j = 1; j <= i__2; ++j) {
				p[j] += wy[i + pointr * wy_dim1] * neggi;
				p[*col + j] += ws[i + pointr * ws_dim1] * neggi;
				pointr = pointr % m + 1;
			}
			if (nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.) {
				/*                   x(i) + d(i) is bounded; compute t(i). */
				++nbreak;
				iorder[nbreak] = i;
				t[nbreak] = tl / (-neggi);
				if (nbreak == 1 || t[nbreak] < bkmin) {
					bkmin = t[nbreak];
					ibkmin = nbreak;
				}
			} else if (nbd[i] >= 2 && neggi > 0.) {
				/*                   x(i) + d(i) is bounded; compute t(i). */
				++nbreak;
				iorder[nbreak] = i;
				t[nbreak] = tu / neggi;
				if (nbreak == 1 || t[nbreak] < bkmin) {
					bkmin = t[nbreak];
					ibkmin = nbreak;
				}
			} else {/*          x(i) + d(i) is not bounded. */
				--nfree;
				iorder[nfree] = i;
				if (fabs(neggi) > 0.)
					bnded = 0;
			}
		}
		/* L50: */
	} /* for(i = 1:n) */

	/*     The indices of the nonzero components of d are now stored */
	/*     in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n). */
	/*     The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin. */
	if (*theta != 1.) {
		/*             complete the initialization of p for theta not= one. */
		dscal(col, theta, &p[*col + 1], &c__1);
	}
	/*     Initialize GCP xcp = x. */
	dcopy(&n, &x[1], &c__1, &xcp[1], &c__1);
	if (nbreak == 0 && nfree == n + 1) {
		/*            is a zero vector, return with the initial xcp as GCP. */
		if (iprint > 100) {
			cout << "Cauchy X =  ";
			for(i = 1; i <= n; i++) cout << xcp[i] << " ";
			cout << "\n";
		}
		return;
	}
	/*     Initialize c = W'(xcp - x) = 0. */
	for (j = 1; j <= col2; ++j)
		c[j] = 0.;

	/*     Initialize derivative f2. */
	f2 = -(*theta) * f1;
	f2_org__ = f2;
	if (*col > 0) {
		bmv(m, &sy[sy_offset], &wt[wt_offset], col, &p[1], &v[1], info);
		if (*info != 0) {
			return;
		}
		f2 -= ddot(&col2, &v[1], &c__1, &p[1], &c__1);
	}
	dtm = -f1 / f2;
	tsum = 0.;
	*nint = 1;
	if (iprint >= 99) cout << "There are " << nbreak << " breakpoints\n";

	/*     If there are no breakpoints, locate the GCP and return. */
	if (nbreak == 0) {
		goto L888;
	}
	nleft = nbreak;
	iter = 1;
	tj = 0.;
	/* ------------------- the beginning of the loop ------------------------- */
	L777:
	/*     Find the next smallest breakpoint; */
	/*     compute dt = t(nleft) - t(nleft + 1). */
	tj0 = tj;
	if (iter == 1) {
		/*       Since we already have the smallest breakpoint we need not do */
		/*       heapsort yet. Often only one breakpoint is used and the */
		/*       cost of heapsort is avoided. */
		tj = bkmin;
		ibp = iorder[ibkmin];
	} else {
		if (iter == 2) {
			/* Replace the already used smallest breakpoint with the */
			/* breakpoint numbered nbreak > nlast, before heapsort call. */
			if (ibkmin != nbreak) {
				t[ibkmin] = t[nbreak];
				iorder[ibkmin] = iorder[nbreak];
			}
		}
		/* Update heap structure of breakpoints */
		/* (if iter=2, initialize heap). */
		hpsolb(nleft, &t[1], &iorder[1], iter - 2);
		tj = t[nleft];
		ibp = iorder[nleft];
	}
	dt = tj - tj0;

	if (dt != 0 && iprint >=  100) {
		cout << "\nPiece    " << *nint << " f1, f2 at start point " << f1 << " " << f2 << "\n",
				cout << "Distance to the next break point =  " << dt << "\n";
		cout << "Distance to the stationary point =  " << dtm << "\n";
	}

	/*     If a minimizer is within this interval, */
	/*     locate the GCP and return. */
	if (dtm < dt) {
		goto L888;
	}
	/*     Otherwise fix one variable and */
	/*     reset the corresponding component of d to zero. */
	tsum += dt;
	--nleft;
	++iter;
	dibp = d[ibp];
	d[ibp] = 0.;
	if (dibp > 0.) {
		zibp = u[ibp] - x[ibp];
		xcp[ibp] = u[ibp];
		iwhere[ibp] = 2;
	} else {
		zibp = l[ibp] - x[ibp];
		xcp[ibp] = l[ibp];
		iwhere[ibp] = 1;
	}
	if (iprint >= 100) cout << "Variable  " << ibp << "  is fixed.\n";
	if (nleft == 0 && nbreak == n) {
		/*                           all n variables are fixed, */
		/*                          return with xcp as GCP. */
		dtm = dt;
		goto L999;
	}
	/*     Update the derivative information. */
	++(*nint);
	dibp2 = dibp * dibp;
	/*     Update f1 and f2. */
	/*      temporarily set f1 and f2 for col=0. */
	f1 += dt * f2 + dibp2 - *theta * dibp * zibp;
	f2 -= *theta * dibp2;
	if (*col > 0) {
		/*                update c = c + dt*p. */
		daxpy(&col2, &dt, &p[1], &c__1, &c[1], &c__1);
		/*         choose wbp, */
		/*         the row of W corresponding to the breakpoint encountered. */
		pointr = *head;
		for (j = 1; j <= *col; ++j) {
			wbp[j] = wy[ibp + pointr * wy_dim1];
			wbp[*col + j] = *theta * ws[ibp + pointr * ws_dim1];
			pointr = pointr % m + 1;
		}
		/*         compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'. */
		bmv(m, &sy[sy_offset], &wt[wt_offset], col, &wbp[1], &v[1], info);
		if (*info != 0) {
			return;
		}
		wmc = ddot(&col2,  &c[1], &c__1, &v[1], &c__1);
		wmp = ddot(&col2,  &p[1], &c__1, &v[1], &c__1);
		wmw = ddot(&col2,&wbp[1], &c__1, &v[1], &c__1);
		/*         update p = p - dibp*wbp. */
		d__1 = -dibp;
		daxpy(&col2, &d__1, &wbp[1], &c__1, &p[1], &c__1);
		/*         complete updating f1 and f2 while col > 0. */
		f1 += dibp * wmc;
		f2 += (2. * dibp * wmp - dibp2 * wmw);
	}
	if(f2 < (d__1 = *epsmch * f2_org__)) f2 = d__1;
	if (nleft > 0) {
		dtm = -f1 / f2;
		goto L777;
		/*           to repeat the loop for unsearched intervals. */
	} else if (bnded) {
		f1 = 0.;
		f2 = 0.;
		dtm = 0.;
	} else {
		dtm = -f1 / f2;
	}
	/* ------------------- the end of the loop ------------------------------- */
	L888:
	if (iprint >= 99) {
		cout << "\nGCP found in this segment\n";
		cout << "Piece    " << *nint << " f1, f2 at start point " << f1 << " " << f2 << "\n";
		cout << "Distance to the stationary point =  " << dtm << "\n";
	}

	if (dtm <= 0.) {
		dtm = 0.;
	}
	tsum += dtm;
	/*     Move free variables (i.e., the ones w/o breakpoints) and */
	/*     the variables whose breakpoints haven't been reached. */
	daxpy(&n, &tsum, &d[1], &c__1, &xcp[1], &c__1);
	L999:
	/*     Update c = c + dtm*p = W'(x^c - x) */
	/*     which will be used in computing r = Z'(B(x^c - x) + g). */
	if (*col > 0) {
		daxpy(&col2, &dtm, &p[1], &c__1, &c[1], &c__1);
	}
	if (iprint >= 100) {
		cout << "Cauchy X =  ";
		for(i = 1; i <= n; i++) cout << xcp[i] << " ";
		cout << "\n";
	}

	if (iprint >= 99)
		cout << "\n---------------- exit CAUCHY----------------------\n\n";
	return;
} /* cauchy */
/* ====================== The end of cauchy ============================== */

void freev(int n, int *nfree, int *indx,
		int *nenter, int *ileave, int *indx2, int *iwhere,
		int *wrk, int *updatd, int *cnstnd, int iprint,
		int *iter)
{
	/*    ************

	 Subroutine freev

	 This subroutine counts the entering and leaving variables when
	 iter > 0, and finds the index set of free and active variables
	 at the GCP.

	 cnstnd is a int variable indicating whether bounds are present

	 indx is an int array of dimension n
	 for i=1,...,nfree, indx(i) are the indices of free variables
	 for i=nfree+1,...,n, indx(i) are the indices of bound variables
	 On entry after the first iteration, indx gives
	 the free variables at the previous iteration.
	 On exit it gives the free variables based on the determination
	 in cauchy using the array iwhere.

	 indx2 is an int array of dimension n
	 On entry indx2 is unspecified.
	 On exit with iter>0, indx2 indicates which variables
	 have changed status since the previous iteration.
	 For i= 1,...,nenter, indx2(i) have changed from bound to free.
	 For i= ileave+1,...,n, indx2(i) have changed from free to bound.


	 *     *  *

	 NEOS, November 1994. (Latest revision June 1996.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	/* System generated locals */
	int i__1;

	/* Local variables */
	int iact, i, k;

	/* Parameter adjustments */
	--iwhere;
	--indx2;
	--indx;

	/* Function Body */
	*nenter = 0;
	*ileave = n + 1;
	if (*iter > 0 && *cnstnd) {/* count the entering and leaving variables. */
		i__1 = *nfree;
		for (i = 1; i <= i__1; ++i) {
			k = indx[i];
			if (iwhere[k] > 0) {
				--(*ileave);
				indx2[*ileave] = k;
				if (iprint >= 100)
					cout << "Variable " << k << " leaves the set of free variables\n";
			}
			/* L20: */
		}
		for (i = *nfree + 1; i <= n; ++i) {
			k = indx[i];
			if (iwhere[k] <= 0) {
				++(*nenter);
				indx2[*nenter] = k;
				if (iprint >= 100)
					cout << "Variable " << k << " enters the set of free variables\n";
			}
			/* L22: */
			if (iprint >= 100)
				cout << n + 1 - *ileave << " variables leave; " << *nenter << " variables enter\n";
		}
	}
	*wrk = *ileave < n + 1 || *nenter > 0 || *updatd;
	/*     Find the index set of free and active variables at the GCP. */
	*nfree = 0;
	iact = n + 1;
	for (i = 1; i <= n; ++i) {
		if (iwhere[i] <= 0) {
			++(*nfree);
			indx[*nfree] = i;
		} else {
			--iact;
			indx[iact] = i;
		}
	}
	if (iprint >= 99)
		cout << *nfree << "  variables are free at GCP on iteration " << *iter + 1 << endl;
	return;
} /* freev */
/* ======================= The end of freev ============================== */

void formk(int n, int *nsub, int *ind, int * nenter, int *ileave,
		int *indx2, int *iupdat, int * updatd, double *wn,
		double *wn1, int m, double *ws, double *wy, double *sy,
		double *theta, int *col, int *head, int *info)
{
	/*     ************

	 Subroutine formk

	 This subroutine forms  the LEL^T factorization of the indefinite

	 matrix     K = [-D -Y'ZZ'Y/theta       L_a'-R_z'  ]
	 [L_a -R_z         theta*S'AA'S ]
	 where E = [-I  0]
	 [ 0  I]
	 The matrix K can be shown to be equal to the matrix M^[-1]N
	 occurring in section 5.1 of [1], as well as to the matrix
	 Mbar^[-1] Nbar in section 5.3.

	 n is an integer variable.
	 On entry n is the dimension of the problem.
	 On exit n is unchanged.

	 nsub is an integer variable
	 On entry nsub is the number of subspace variables in free set.
	 On exit nsub is not changed.

	 ind is an integer array of dimension nsub.
	 On entry ind specifies the indices of subspace variables.
	 On exit ind is unchanged.

	 nenter is an integer variable.
	 On entry nenter is the number of variables entering the
	 free set.
	 On exit nenter is unchanged.

	 ileave is an integer variable.
	 On entry indx2(ileave),...,indx2(n) are the variables leaving
	 the free set.
	 On exit ileave is unchanged.

	 indx2 is an integer array of dimension n.
	 On entry indx2(1),...,indx2(nenter) are the variables entering
	 the free set, while indx2(ileave),...,indx2(n) are the
	 variables leaving the free set.
	 On exit indx2 is unchanged.
	 p
	 iupdat is an integer variable.
	 On entry iupdat is the total number of BFGS updates made so far.
	 On exit iupdat is unchanged.

	 updatd is a logical variable.
	 On entry 'updatd' is true if the L-BFGS matrix is updatd.
	 On exit 'updatd' is unchanged.

	 wn is a double precision array of dimension 2m x 2m.
	 On entry wn is unspecified.
	 On exit the upper triangle of wn stores the LEL^T factorization
	 of the 2*col x 2*col indefinite matrix
	 [-D -Y'ZZ'Y/theta       L_a'-R_z'  ]
	 [L_a -R_z         theta*S'AA'S ]

	 wn1 is a double precision array of dimension 2m x 2m.
	 On entry wn1 stores the lower triangular part of
	 [Y' ZZ'Y    L_a'+R_z']
	 [L_a+R_z    S'AA'S     ]
	 in the previous iteration.
	 On exit wn1 stores the corresponding updated matrices.
	 The purpose of wn1 is just to store these inner products
	 so they can be easily updated and inserted into wn.

	 m is an integer variable.
	 On entry m is the maximum number of variable metric corrections
	 used to define the limited memory matrix.
	 On exit m is unchanged.

	 ws, wy, sy, and wtyy are double precision arrays;
	 theta is a double precision variable;
	 col is an integer variable;
	 head is an integer variable.
	 On entry they store the information defining the
	 limited memory BFGS matrix:
	 ws(n,m) stores S, a set of s-vectors;
	 wy(n,m) stores Y, a set of y-vectors;
	 sy(m,m) stores S'Y;
	 wtyy(m,m) stores the Cholesky factorization
	 of (theta*S'S+LD^(-1)L')
	 theta is the scaling factor specifying B_0 = theta I;
	 col is the number of variable metric corrections stored;
	 head is the location of the 1st s- (or y-) vector in S (or Y).
	 On exit they are unchanged.

	 info is an integer variable.
	 On entry info is unspecified.
	 On exit info =  0 for normal return;
	 = -1 when the 1st Cholesky factorization failed;
	 = -2 when the 2st Cholesky factorization failed.

	 Subprograms called:

	 Linpack ... dcopy, dpofa, dtrsl.


	 References:
	 [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
	 memory algorithm for bound constrained optimization'',
	 SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

	 [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
	 limited memory FORTRAN code for solving bound constrained
	 optimization problems'', Tech. Report, NAM-11, EECS Department,
	 Northwestern University, 1994.

	 (Postscript files of these papers are available via anonymous
	 ftp to ece.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)

	 *  *     *

	 NEOS, November 1994. (Latest revision April 1997.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	/* System generated locals */
	int wn_dim1, wn_offset, wn1_dim1, wn1_offset, ws_dim1, ws_offset,
	wy_dim1, wy_offset, sy_dim1, sy_offset, i__1, i__2;

	/* Local variables */
	int dend, pend;
	int upcl;
	double temp1, temp2, temp3, temp4;
	int i, k;
	int ipntr, jpntr, k1, m2, dbegin, is, js, iy, jy, pbegin, is1, js1,
	col2;

	/* Parameter adjustments */
	--indx2;
	--ind;
	sy_dim1 = m;
	sy_offset = 1 + sy_dim1 * 1;
	sy -= sy_offset;
	wy_dim1 = n;
	wy_offset = 1 + wy_dim1 * 1;
	wy -= wy_offset;
	ws_dim1 = n;
	ws_offset = 1 + ws_dim1 * 1;
	ws -= ws_offset;
	wn1_dim1 = 2 * m;
	wn1_offset = 1 + wn1_dim1 * 1;
	wn1 -= wn1_offset;
	wn_dim1 = 2 * m;
	wn_offset = 1 + wn_dim1 * 1;
	wn -= wn_offset;

	/* Function Body */

	/*     Form the lower triangular part of */
	/*         WN1 = [Y' ZZ'Y      L_a'+R_z'] */
	/*               [L_a+R_z      S'AA'S   ] */
	/*      where L_a is the strictly lower triangular part of S'AA'Y */
	/*        R_z is the upper triangular part of S'ZZ'Y. */

	if (*updatd) {
		if (*iupdat > m) {/*        shift old part of WN1. */
			i__1 = m - 1;
			for (jy = 1; jy <= i__1; ++jy) {
				js = m + jy;
				i__2 = m - jy;
				dcopy(&i__2, &wn1[jy + 1 + (jy + 1)* wn1_dim1], &c__1,
						&wn1[jy + jy * wn1_dim1], &c__1);
				dcopy(&i__2, &wn1[js + 1 + (js + 1)* wn1_dim1], &c__1,
						&wn1[js + js * wn1_dim1], &c__1);
				i__2 = m - 1;
				dcopy(&i__2, &wn1[m + 2 + (jy + 1) * wn1_dim1], &c__1,
						&wn1[m + 1 + jy * wn1_dim1], &c__1);
				/* L10: */
			}
		}
		/*        put new rows in blocks (1,1), (2,1) and (2,2). */
		pbegin = 1;
		pend = *nsub;
		dbegin = *nsub + 1;
		dend = n;
		iy = *col;
		is = m + *col;
		ipntr = *head + *col - 1;
		if (ipntr > m) {
			ipntr -= m;
		}
		jpntr = *head;
		i__1 = *col;
		for (jy = 1; jy <= i__1; ++jy) {
			js = m + jy;
			temp1 = 0.;
			temp2 = 0.;
			temp3 = 0.;
			/*           compute element jy of row 'col' of Y'ZZ'Y */
			for (k = pbegin; k <= pend; ++k) {
				k1 = ind[k];
				temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
			}
			/*           compute elements jy of row 'col' of L_a and S'AA'S */
			for (k = dbegin; k <= dend; ++k) {
				k1 = ind[k];
				temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
				temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
			}
			wn1[iy + jy * wn1_dim1] = temp1;
			wn1[is + js * wn1_dim1] = temp2;
			wn1[is + jy * wn1_dim1] = temp3;
			jpntr = jpntr % m + 1;
			/* L20: */
		}
		/*        put new column in block (2,1). */
		jy = *col;
		jpntr = *head + *col - 1;
		if (jpntr > m) {
			jpntr -= m;
		}
		ipntr = *head;
		i__1 = *col;
		for (i = 1; i <= i__1; ++i) {
			is = m + i;
			temp3 = 0.;
			/*           compute element i of column 'col' of R_z */
			for (k = pbegin; k <= pend; ++k) {
				k1 = ind[k];
				temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
			}
			ipntr = ipntr % m + 1;
			wn1[is + jy * wn1_dim1] = temp3;
			/* L30: */
		}
		upcl = *col - 1;
	} else {
		upcl = *col;
	}
	/*     modify the old parts in blocks (1,1) and (2,2) due to changes */
	/*     in the set of free variables. */
	ipntr = *head;
	for (iy = 1; iy <= upcl; ++iy) {
		is = m + iy;
		jpntr = *head;
		for (jy = 1; jy <= iy; ++jy) {
			js = m + jy;
			temp1 = 0.;
			temp2 = 0.;
			temp3 = 0.;
			temp4 = 0.;
			for (k = 1; k <= *nenter; ++k) {
				k1 = indx2[k];
				temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
				temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
			}
			for (k = *ileave; k <= n; ++k) {
				k1 = indx2[k];
				temp3 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
				temp4 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
			}
			wn1[iy + jy * wn1_dim1] = wn1[iy + jy * wn1_dim1] + temp1 - temp3;
			wn1[is + js * wn1_dim1] = wn1[is + js * wn1_dim1] - temp2 + temp4;
			jpntr = jpntr % m + 1;
			/* L40: */
		}
		ipntr = ipntr % m + 1;
		/* L45: */
	}
	/*     modify the old parts in block (2,1). */
	ipntr = *head;
	for (is = m + 1; is <= m + upcl; ++is) {
		jpntr = *head;
		for (jy = 1; jy <= upcl; ++jy) {
			temp1 = 0.;
			temp3 = 0.;
			for (k = 1; k <= *nenter; ++k) {
				k1 = indx2[k];
				temp1 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
			}
			for (k = *ileave; k <= n; ++k) {
				k1 = indx2[k];
				temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
			}
			if (is <= jy + m) {
				wn1[is + jy * wn1_dim1] +=  temp1 - temp3;
			} else {
				wn1[is + jy * wn1_dim1] += -temp1 + temp3;
			}
			jpntr = jpntr % m + 1;
			/* L55: */
		}
		ipntr = ipntr % m + 1;
		/* L60: */
	}
	/*     Form the upper triangle of WN = [D+Y' ZZ'Y/theta      -L_a'+R_z' ] */
	/*                       [-L_a +R_z     S'AA'S*theta] */
	m2 = m << 1;
	i__1 = *col;
	for (iy = 1; iy <= i__1; ++iy) {
		is = *col + iy;
		is1 = m + iy;
		i__2 = iy;
		for (jy = 1; jy <= i__2; ++jy) {
			js = *col + jy;
			js1 = m + jy;
			wn[jy + iy * wn_dim1] = wn1[iy + jy * wn1_dim1] / *theta;
			wn[js + is * wn_dim1] = wn1[is1 + js1 * wn1_dim1] * *theta;
			/* L65: */
		}
		i__2 = iy - 1;
		for (jy = 1; jy <= i__2; ++jy) {
			wn[jy + is * wn_dim1] = -wn1[is1 + jy * wn1_dim1];
		}
		i__2 = *col;
		for (jy = iy; jy <= i__2; ++jy) {
			wn[jy + is * wn_dim1] = wn1[is1 + jy * wn1_dim1];
		}
		wn[iy + iy * wn_dim1] += sy[iy + iy * sy_dim1];
		/* L70: */
	}
	/*     Form the upper triangle of */
	/*        WN= [  LL'          L^-1(-L_a'+R_z')] */
	/*        [(-L_a +R_z)L'^-1   S'AA'S*theta  ] */
	/*      first Cholesky factor (1,1) block of wn to get LL' */
	/*                with L' stored in the upper triangle of wn. */
	dpofa(&wn[wn_offset], &m2, col, info);
	if (*info != 0) {
		*info = -1;
		return;
	}
	/*      then form L^-1(-L_a'+R_z') in the (1,2) block. */
	col2 = *col << 1;
	for (js = *col + 1; js <= col2; ++js) {
		dtrsl(&wn[wn_offset], &m2, col,
				&wn[js * wn_dim1 + 1], &c__11, info);
	}
	/*     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the */
	/*      upper triangle of (2,2) block of wn. */
	for (is = *col + 1; is <= col2; ++is) {
		for (js = is; js <= col2; ++js) {
			wn[is + js * wn_dim1] +=
					ddot(col, &wn[is * wn_dim1 + 1], &c__1,
							&wn[js * wn_dim1 + 1], &c__1);
		}
		/* L72: */
	}
	/*     Cholesky factorization of (2,2) block of wn. */
	dpofa(&wn[*col + 1 + (*col + 1) * wn_dim1], &m2, col, info);
	if (*info != 0) {
		*info = -2;
		return;
	}
	return;
} /* formk */
/* ======================= The end of formk ============================== */

void cmprlb(int n, int m, double *x,
		double *g, double *ws, double *wy, double *sy,
		double *wt, double *z, double *r, double *wa,
		int *indx, double *theta, int *col, int *head,
		int *nfree, int *cnstnd, int *info)
{
	/*    ************

	 Subroutine cmprlb

	 This subroutine computes r=-Z'B(xcp-xk)-Z'g by using
	 wa(2m+1)=W'(xcp-x) from subroutine cauchy.

	 Subprograms called:

	 L-BFGS-B Library ... bmv.


	 *     *  *

	 NEOS, November 1994. (Latest revision June 1996.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	/* System generated locals */
	int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset,
	wt_dim1, wt_offset, Col, n_f;

	/* Local variables */
	int i, j, k;
	double a1, a2;
	int pointr;

	/* Parameter adjustments */
	--indx;
	--r;
	--z;
	--g;
	--x;
	--wa;
	wt_dim1 = m;
	wt_offset = 1 + wt_dim1 * 1;
	wt -= wt_offset;
	sy_dim1 = m;
	sy_offset = 1 + sy_dim1 * 1;
	sy -= sy_offset;
	wy_dim1 = n;
	wy_offset = 1 + wy_dim1 * 1;
	wy -= wy_offset;
	ws_dim1 = n;
	ws_offset = 1 + ws_dim1 * 1;
	ws -= ws_offset;

	/* Function Body */
	Col = *col;
	if (! (*cnstnd) && Col > 0) {
		for (i = 1; i <= n; ++i)
			r[i] = -g[i];
	}
	else {
		n_f = *nfree;
		for (i = 1; i <= n_f; ++i) {
			k = indx[i];
			r[i] = -(*theta) * (z[k] - x[k]) - g[k];
		}
		bmv(m, &sy[sy_offset], &wt[wt_offset], col,
				&wa[(m << 1) + 1], &wa[1], info);
		if (*info != 0) {
			*info = -8;
			return;
		}
		pointr = *head;
		for (j = 1; j <= Col; ++j) {
			a1 = wa[j];
			a2 = *theta * wa[Col + j];
			for (i = 1; i <= n_f; ++i) {
				k = indx[i];
				r[i] += wy[k + pointr * wy_dim1] * a1 +
						ws[k + pointr * ws_dim1] * a2;
			}
			pointr = pointr % m + 1;
		}
	}
	return;
} /* cmprlb */
/* ======================= The end of cmprlb ============================= */


void subsm(int n, int m, int *nsub, int *ind,
		double *l, double *u, int *nbd, double *x,
		double *d, double *ws, double *wy, double *theta,
		int *col, int *head, int *iword, double *wv,
		double *wn, int iprint, int *info)
{
	/*    ************

	 Subroutine subsm

	 Given xcp, l, u, r, an index set that specifies
	 the active set at xcp, and an l-BFGS matrix B
	 (in terms of WY, WS, SY, WT, head, col, and theta),
	 this subroutine computes an approximate solution
	 of the subspace problem

	 (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)

	 subject to l <= x <= u
	 x_i = xcp_i   for all i in A(xcp)

	 along the subspace unconstrained Newton direction

	 d = -(Z'BZ)^(-1) r.

	 The formula for the Newton direction, given the L-BFGS matrix
	 and the Sherman-Morrison formula, is

	 d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.

	 where
	 K = [-D -Y'ZZ'Y/theta        L_a'-R_z'  ]
	 [L_a -R_z          theta*S'AA'S ]

	 Note that this procedure for computing d differs
	 from that described in [1]. One can show that the matrix K is
	 equal to the matrix M^[-1]N in that paper.

	 n is an integer variable.
	 On entry n is the dimension of the problem.
	 On exit n is unchanged.

	 m is an integer variable.
	 On entry m is the maximum number of variable metric corrections
	 used to define the limited memory matrix.
	 On exit m is unchanged.

	 nsub is an integer variable.
	 On entry nsub is the number of free variables.
	 On exit nsub is unchanged.

	 ind is an integer array of dimension nsub.
	 On entry ind specifies the coordinate indices of free variables.
	 On exit ind is unchanged.

	 l is a double precision array of dimension n.
	 On entry l is the lower bound of x.
	 On exit l is unchanged.

	 u is a double precision array of dimension n.
	 On entry u is the upper bound of x.
	 On exit u is unchanged.

	 nbd is a integer array of dimension n.
	 On entry nbd represents the type of bounds imposed on the
	 variables, and must be specified as follows:
	 nbd(i)=0 if x(i) is unbounded,
	 1 if x(i) has only a lower bound,
	 2 if x(i) has both lower and upper bounds, and
	 3 if x(i) has only an upper bound.
	 On exit nbd is unchanged.

	 x is a double precision array of dimension n.
	 On entry x specifies the Cauchy point xcp.
	 On exit x(i) is the minimizer of Q over the subspace of
	 free variables.

	 d is a double precision array of dimension n.
	 On entry d is the reduced gradient of Q at xcp.
	 On exit d is the Newton direction of Q.

	 ws and wy are double precision arrays;
	 theta is a double precision variable;
	 col is an integer variable;
	 head is an integer variable.
	 On entry they store the information defining the
	 limited memory BFGS matrix:
	 ws(n,m) stores S, a set of s-vectors;
	 wy(n,m) stores Y, a set of y-vectors;
	 theta is the scaling factor specifying B_0 = theta I;
	 col is the number of variable metric corrections stored;
	 head is the location of the 1st s- (or y-) vector in S (or Y).
	 On exit they are unchanged.

	 iword is an integer variable.
	 On entry iword is unspecified.
	 On exit iword specifies the status of the subspace solution.
	 iword = 0 if the solution is in the box,
	 1 if some bound is encountered.

	 wv is a double precision working array of dimension 2m.

	 wn is a double precision array of dimension 2m x 2m.
	 On entry the upper triangle of wn stores the LEL^T factorization
	 of the indefinite matrix

	 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
	 [L_a -R_z           theta*S'AA'S ]
	 where E = [-I  0]
	 [ 0  I]
	 On exit wn is unchanged.

	 iprint is an INTEGER variable that must be set by the user.
	 It controls the frequency and type of output generated:
	 iprint<0    no output is generated;
	 iprint=0    print only one line at the last iteration;
	 0<iprint<99 print also f and |proj g| every iprint iterations;
	 iprint=99   print details of every iteration except n-vectors;
	 iprint=100  print also the changes of active set and final x;
	 iprint>100  print details of every iteration including x and g;
	 When iprint > 0, the file iterate.dat will be created to
	 summarize the iteration.

	 info is an integer variable.
	 On entry info is unspecified.
	 On exit info = 0       for normal return,
	 = nonzero for abnormal return
	 when the matrix K is ill-conditioned.

	 Subprograms called:

	 Linpack dtrsl.


	 References:

	 [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
	 memory algorithm for bound constrained optimization'',
	 SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.



	 *  *  *

	 NEOS, November 1994. (Latest revision June 1996.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	/* System generated locals */
	int ws_offset, wn_dim1, wn_offset;

	/* Local variables */
	double alpha, dk, temp1, temp2;
	int i, j, k, m2, js, jy, pointr, ibd = 0, col2, ns;

	/* Parameter adjustments */
	--d;
	--u;
	--l;
	--x;
	--ind;
	--nbd;
	--wv;
	wn_dim1 = 2 * m;
	wn_offset = 1 + wn_dim1 * 1;
	wn -= wn_offset;
	/* ws[] and wy[] are both  [n x m ] :*/
	ws_offset = 1 + n * 1;
	ws -= ws_offset;
	wy -= ws_offset;

	ns = *nsub;
	if (ns <= 0)
		return;

	/*     Compute wv = W'Zd. */
	pointr = *head;
	for (i = 1; i <= *col; ++i) {
		temp1 = 0.;
		temp2 = 0.;
		for (j = 1; j <= ns; ++j) {
			k = ind[j];
			temp1 += wy[k + pointr * n] * d[j];
			temp2 += ws[k + pointr * n] * d[j];
		}
		wv[i] = temp1;
		wv[*col + i] = *theta * temp2;
		pointr = pointr % m + 1;
		/* L20: */
	}
	/*     Compute wv:=K^(-1)wv. */
	m2 = m << 1;
	col2 = *col << 1;
	dtrsl(&wn[wn_offset], &m2, &col2, &wv[1], &c__11, info);
	if (*info != 0) {
		return;
	}
	for (i = 1; i <= *col; ++i)
		wv[i] = -wv[i];

	dtrsl(&wn[wn_offset], &m2, &col2, &wv[1], &c__1, info);
	if (*info != 0) {
		return;
	}
	/*     Compute d = (1/theta)d + (1/theta**2)Z'W wv. */
	pointr = *head;
	for (jy = 1; jy <= *col; ++jy) {
		js = *col + jy;
		for (i = 1; i <= ns; ++i) {
			k = ind[i];
			d[i] += (wy[k + pointr * n] * wv[jy] / *theta +
					ws[k + pointr * n] * wv[js]);
		}
		pointr = pointr % m + 1;
		/* L40: */
	}

	for (i = 1; i <= ns; ++i)
		d[i] /= *theta;

	/*     Backtrack to the feasible region. */
	alpha = 1.;
	temp1 = alpha;
	for (i = 1; i <= ns; ++i) {
		k = ind[i];
		dk = d[i];
		if (nbd[k] != 0) {
			if (dk < 0. && nbd[k] <= 2) {
				temp2 = l[k] - x[k];
				if (temp2 >= 0.) {
					temp1 = 0.;
				} else if (dk * alpha < temp2) {
					temp1 = temp2 / dk;
				}
			} else if (dk > 0. && nbd[k] >= 2) {
				temp2 = u[k] - x[k];
				if (temp2 <= 0.) {
					temp1 = 0.;
				} else if (dk * alpha > temp2) {
					temp1 = temp2 / dk;
				}
			}
			if (temp1 < alpha) {
				alpha = temp1;
				ibd = i;
			}
		}
		/* L60: */
	}
	if (alpha < 1.) {
		dk = d[ibd];
		k = ind[ibd];
		if (dk > 0.) {
			x[k] = u[k];
			d[ibd] = 0.;
		} else if (dk < 0.) {
			x[k] = l[k];
			d[ibd] = 0.;
		}
	}
	for (i = 1; i <= ns; ++i)
		x[ind[i]] += alpha * d[i];

	*iword = (alpha < 1.) ? 1 : 0;

	return;
} /* subsm */
/* ====================== The end of subsm =============================== */

void lnsrlb(int n, double *l, double *u,
		int *nbd, double *x, double *f, double *fold,
		double *gd, double *gdold, double *g, double *d,
		double *r, double *t, double *z, double *stp,
		double *dnorm, double *dtd, double *xstep,
		double *stpmx, int *iter, int *ifun, int *iback, int *nfgv,
		int *info, char *task, int *boxed, int *cnstnd,
		char *csave, int *isave, double *dsave)
{
	/*     **********

	 Subroutine lnsrlb

	 This subroutine calls subroutine dcsrch from the Minpack2 library
	 to perform the line search.  Subroutine dscrch is safeguarded so
	 that all trial points lie within the feasible region.

	 Subprograms called:

	 Minpack2 Library ... dcsrch.

	 Linpack ... dtrsl, ddot.


	 *     *  *

	 NEOS, November 1994. (Latest revision June 1996.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 **********
	 */

	/* For dcsrch(): */
	const double stpmin = 0.;
	const double ftol = .001;
	const double gtol = .9;
	const double xtol = .1;

	/* System generated locals */
	double d1;

	/* Local variables */
	int i;
	double a1, a2;

	/* Parameter adjustments */
	--z;
	--t;
	--r;
	--d;
	--g;
	--x;
	--nbd;
	--u;
	--l;

	/* Function Body */
	if (strncmp(task, "FG_LN", 5) == 0) {
		goto L556;
	}
	*dtd = ddot(&n, &d[1], &c__1, &d[1], &c__1);
	*dnorm = sqrt(*dtd);
	/*     Determine the maximum step length. */
	*stpmx = 1e10;
	if (*cnstnd) {
		if (*iter == 0) {
			*stpmx = 1.;
		} else {
			for (i = 1; i <= n; ++i) {
				a1 = d[i];
				if (nbd[i] != 0) {
					if (a1 < 0. && nbd[i] <= 2) {
						a2 = l[i] - x[i];
						if (a2 >= 0.) {
							*stpmx = 0.;
						} else if (a1 * *stpmx < a2) {
							*stpmx = a2 / a1;
						}
					} else if (a1 > 0. && nbd[i] >= 2) {
						a2 = u[i] - x[i];
						if (a2 <= 0.) {
							*stpmx = 0.;
						} else if (a1 * *stpmx > a2) {
							*stpmx = a2 / a1;
						}
					}
				}
				/* L43: */
			}
		}
	}
	if (*iter == 0 && ! (*boxed)) {
		d1 = 1. / *dnorm;
		*stp = min(d1,*stpmx);
	} else {
		*stp = 1.;
	}
	dcopy(&n, &x[1], &c__1, &t[1], &c__1);
	dcopy(&n, &g[1], &c__1, &r[1], &c__1);
	*fold = *f;
	*ifun = 0;
	*iback = 0;
	strcpy(csave, "START");
	L556:
	*gd = ddot(&n, &g[1], &c__1, &d[1], &c__1);
	if (*ifun == 0) {
		*gdold = *gd;
		if (*gd >= 0.) {
			/*                 the directional derivative >=0. */
			/*                 Line search is impossible. */
			*info = -4;
			return;
		}
	}
	dcsrch(f, gd, stp,
			ftol, gtol, xtol,
			stpmin, *stpmx,
			csave, isave, dsave);
	*xstep = *stp * *dnorm;
	if (strncmp(csave, "CONV", 4) != 0 && strncmp(csave, "WARN", 4) != 0) {
		strcpy(task, "FG_LNSRCH");
		++(*ifun);
		++(*nfgv);
		*iback = *ifun - 1;
		if (*stp == 1.) {
			dcopy(&n, &z[1], &c__1, &x[1], &c__1);
		} else {
			for (i = 1; i <= n; ++i) {
				x[i] = *stp * d[i] + t[i];
			}
		}
	} else {
		strcpy(task, "NEW_X");
	}
	return;
} /* lnsrlb */
/* ======================= The end of lnsrlb ============================= */

void matupd(int n, int m, double *ws,
		double *wy, double *sy, double *ss, double *d,
		double *r, int *itail, int *iupdat, int *col,
		int *head, double *theta, double *rr, double *dr,
		double *stp, double *dtd)
{
	/*    ************

	 Subroutine matupd

	 This subroutine updates matrices WS and WY, and forms the
	 middle matrix in B.

	 Subprograms called:

	 Linpack ... dcopy, ddot.


	 *     *  *

	 NEOS, November 1994. (Latest revision June 1996.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	/* System generated locals */
	int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset,
	ss_dim1, ss_offset, i__1, i__2;

	/* Local variables */
	int j;
	int pointr;

	/* Parameter adjustments */
	--r;
	--d;
	ss_dim1 = m;
	ss_offset = 1 + ss_dim1 * 1;
	ss -= ss_offset;
	sy_dim1 = m;
	sy_offset = 1 + sy_dim1 * 1;
	sy -= sy_offset;
	wy_dim1 = n;
	wy_offset = 1 + wy_dim1 * 1;
	wy -= wy_offset;
	ws_dim1 = n;
	ws_offset = 1 + ws_dim1 * 1;
	ws -= ws_offset;

	/* Function Body */

	/*     Set pointers for matrices WS and WY. */
	if (*iupdat <= m) {
		*col = *iupdat;
		*itail = (*head + *iupdat - 2) % m + 1;
	} else {
		*itail = *itail % m + 1;
		*head = *head % m + 1;
	}
	/*     Update matrices WS and WY. */
	dcopy(&n, &d[1], &c__1, &ws[*itail * ws_dim1 + 1], &c__1);
	dcopy(&n, &r[1], &c__1, &wy[*itail * wy_dim1 + 1], &c__1);
	/*     Set theta=yy/ys. */
	*theta = *rr / *dr;
	/*     Form the middle matrix in B. */
	/*      update the upper triangle of SS, */
	/*                       and the lower triangle of SY: */
	if (*iupdat > m) {
		/*                move old information */
		i__1 = *col - 1;
		for (j = 1; j <= i__1; ++j) {
			dcopy(&j, &ss[(j + 1) * ss_dim1 + 2], &c__1,
					&ss[j * ss_dim1 + 1], &c__1);
			i__2 = *col - j;
			dcopy(&i__2, &sy[j + 1 + (j + 1) * sy_dim1], &c__1,
					&sy[j + j * sy_dim1], &c__1);
			/* L50: */
		}
	}
	/*      add new information: the last row of SY */
	/*                           and the last column of SS: */
	pointr = *head;
	i__1 = *col - 1;
	for (j = 1; j <= i__1; ++j) {
		sy[*col + j * sy_dim1] =
				ddot(&n, &d[1], &c__1, &wy[pointr * wy_dim1 + 1], &c__1);
		ss[j + *col * ss_dim1] =
				ddot(&n, &ws[pointr * ws_dim1 + 1], &c__1, &d[1], &c__1);
		pointr = pointr % m + 1;
		/* L51: */
	}
	if (*stp == 1.) {
		ss[*col + *col * ss_dim1] = *dtd;
	} else {
		ss[*col + *col * ss_dim1] = *stp * *stp * *dtd;
	}
	sy[*col + *col * sy_dim1] = *dr;
	return;
} /* matupd */
/* ======================= The end of matupd ============================= */

void formt(int m, double *wt, double *sy, double *ss,
		int *col, double *theta, int *info)
{
	/*     ************

	 Subroutine formt

	 This subroutine forms the upper half of the pos. def. and symm.
	 T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
	 of the array wt, and performs the Cholesky factorization of T
	 to produce J*J', with J' stored in the upper triangle of wt.

	 Subprograms called:

	 Linpack ... dpofa.


	 *  *     *

	 NEOS, November 1994. (Latest revision June 1996.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	/* System generated locals */
	int wt_dim1, wt_offset, sy_dim1, sy_offset, ss_dim1, ss_offset, i__1;

	/* Local variables */
	double ddum;
	int i, j, k;
	int k1;

	/* Parameter adjustments */
	ss_dim1 = m;
	ss_offset = 1 + ss_dim1 * 1;
	ss -= ss_offset;
	sy_dim1 = m;
	sy_offset = 1 + sy_dim1 * 1;
	sy -= sy_offset;
	wt_dim1 = m;
	wt_offset = 1 + wt_dim1 * 1;
	wt -= wt_offset;

	/* Function Body */

	/*     Form the upper half of  T = theta*SS + L*D^(-1)*L', */
	/*      store T in the upper triangle of the array wt. */
	i__1 = *col;
	for (j = 1; j <= i__1; ++j) {
		wt[j * wt_dim1 + 1] = *theta * ss[j * ss_dim1 + 1];
	}
	for (i = 2; i <= i__1; ++i) {
		for (j = i; j <= i__1; ++j) {
			k1 = min(i,j) - 1;
			ddum = 0.;
			for (k = 1; k <= k1; ++k) {
				ddum += sy[i + k * sy_dim1] * sy[j + k * sy_dim1] / sy[k +
																	   k * sy_dim1];
			}
			wt[i + j * wt_dim1] = ddum + *theta * ss[i + j * ss_dim1];
		}
		/* L55: */
	}
	/*     Cholesky factorize T to J*J' with */
	/*      J' stored in the upper triangle of wt. */
	dpofa(&wt[wt_offset], &m, col, info);
	if (*info != 0) {
		*info = -3;
	}
	return;
} /* formt */

/* ======================= The end of formt ============================== */

void bmv(int m, double *sy, double *wt,
		int *col, double *v, double *p, int *info)
{
	/*     ************

	 *     Subroutine bmv

	 *     This subroutine computes the product of the 2m x 2m middle matrix
	 *     in the compact L-BFGS formula of B and a 2m vector v;
	 *     it returns the product in p.

	 *     m is an integer variable.
	 *     On entry m is the maximum number of variable metric corrections
	 *       used to define the limited memory matrix.
	 *     On exit m is unchanged.

	 *     sy is a double precision array of dimension m x m.
	 *     On entry sy specifies the matrix S'Y.
	 *     On exit sy is unchanged.

	 *     wt is a double precision array of dimension m x m.
	 *     On entry wt specifies the upper triangular matrix J' which is
	 *       the Cholesky factor of (thetaS'S+LD^(-1)L').
	 *     On exit wt is unchanged.

	 *     col is an integer variable.
	 *     On entry col specifies the number of s-vectors (or y-vectors)
	 *       stored in the compact L-BFGS formula.
	 *     On exit col is unchanged.

	 *     v is a double precision array of dimension 2col.
	 *     On entry v specifies vector v.
	 *     On exit v is unchanged.

	 *     p is a double precision array of dimension 2col.
	 *     On entry p is unspecified.
	 *     On exit p is the product Mv.

	 *     info is an integer variable.
	 *     On entry info is unspecified.
	 *     On exit info = 0    for normal return,
	 *              = nonzero for abnormal return when the system
	 *                  to be solved by dtrsl is singular.

	 *     Subprograms called:

	 *     Linpack ... dtrsl.


	 *                 *    *  *

	 *     NEOS, November 1994. (Latest revision June 1996.)
	 *     Optimization Technology Center.
	 *     Argonne National Laboratory and Northwestern University.
	 *     Written by
	 *              Ciyou Zhu
	 *     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 *     ************
	 */

	/* System generated locals */
	int sy_dim1, sy_offset, wt_dim1, wt_offset, Col;


	/* Local variables */
	int i, k;
	int i2;
	double sum;

	/* Parameter adjustments */
	wt_dim1 = m;
	wt_offset = 1 + wt_dim1 * 1;
	wt -= wt_offset;
	sy_dim1 = m;
	sy_offset = 1 + sy_dim1 * 1;
	sy -= sy_offset;
	--p;
	--v;

	/* Function Body */
	if (*col == 0) {
		return;
	}
	/*    PART I: solve [     D^(1/2)      O ] [ p1 ] = [ v1 ]
	 *              [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].
	 *    solve Jp2=v2+LD^(-1)v1.
	 */
	Col = *col;
	p[*col + 1] = v[*col + 1];
	for (i = 2; i <= Col; ++i) {
		i2 = *col + i;
		sum = 0.;
		for (k = 1; k <= i - 1; ++k) {
			sum += sy[i + k * sy_dim1] * v[k] / sy[k + k * sy_dim1];
		}
		p[i2] = v[i2] + sum;
		/* L20: */
	}
	/*     Solve the triangular system */
	dtrsl(&wt[wt_offset], &m, col, &p[*col + 1], &c__11, info);
	if (*info != 0) {
		return;
	}
	/*     solve D^(1/2)p1=v1. */
	for (i = 1; i <= Col; ++i) {
		p[i] = v[i] / sqrt(sy[i + i * sy_dim1]);
	}

	/*    PART II: solve [ -D^(1/2)   D^(-1/2)*L'     ] [ p1 ] = [ p1 ]
	 *               [  0        J'         ] [ p2 ]   [ p2 ].
	 *    solve J^Tp2=p2.
	 */
	dtrsl(&wt[wt_offset], &m, col, &p[*col + 1], &c__1, info);
	if (*info != 0) {
		return;
	}
	/*     compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2) */
	/*           =-D^(-1/2)p1 + D^(-1)L'p2. */
	for (i = 1; i <= Col; ++i) {
		p[i] = -p[i] / sqrt(sy[i + i * sy_dim1]);
	}
	for (i = 1; i <= Col; ++i) {
		sum = 0.;
		for (k = i + 1; k <= Col; ++k) {
			sum += sy[k + i * sy_dim1] * p[*col + k] / sy[i + i * sy_dim1];
		}
		p[i] += sum;
		/* L60: */
	}
	return;
} /* bmv */
/* ======================== The end of bmv =============================== */

void hpsolb(int n, double *t, int *iorder, int iheap)
{
	/*    ************

	 Subroutine hpsolb

	 This subroutine sorts out the least element of t, and puts the
	 remaining elements of t in a heap.

	 n is an int variable.
	 On entry n is the dimension of the arrays t and iorder.
	 On exit n is unchanged.

	 t is a double precision array of dimension n.
	 On entry t stores the elements to be sorted,
	 On exit t(n) stores the least elements of t, and t(1) to t(n-1)
	 stores the remaining elements in the form of a heap.

	 iorder is an int array of dimension n.
	 On entry iorder(i) is the index of t(i).
	 On exit iorder(i) is still the index of t(i), but iorder may be
	 permuted in accordance with t.

	 iheap is an int variable specifying the task.
	 On entry iheap should be set as follows:
	 iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,
	 iheap .ne. 0 if otherwise.
	 On exit iheap is unchanged.


	 References:
	 Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.

	 *     *  *

	 NEOS, November 1994. (Latest revision June 1996.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */

	/* Local variables */
	double ddum;
	int i, j, k, indxin, indxou;
	double out;

	/* Parameter adjustments */
	--iorder;
	--t;

	/* Function Body */
	if (iheap == 0) {
		/*      Rearrange the elements t(1) to t(n) to form a heap. */
		for (k = 2; k <= n; ++k) {
			ddum = t[k];
			indxin = iorder[k];
			/*         Add ddum to the heap. */
			i = k;
			h_loop:
			if (i > 1) {
				j = i / 2;
				if (ddum < t[j]) {
					t[i] = t[j];
					iorder[i] = iorder[j];
					i = j;
					goto h_loop;
				}
			}
			t[i] = ddum;
			iorder[i] = indxin;
			/* L20: */
		}
	}
	/*     Assign to 'out' the value of t(1), the least member of the heap, */
	/*      and rearrange the remaining members to form a heap as */
	/*      elements 1 to n-1 of t. */
	if (n > 1) {
		i = 1;
		out = t[1];
		indxou = iorder[1];
		ddum = t[n];
		indxin = iorder[n];
		/*      Restore the heap */
		Loop:
		j = i + i;
		if (j <= n - 1) {
			if (t[j + 1] < t[j]) {
				++j;
			}
			if (t[j] < ddum) {
				t[i] = t[j];
				iorder[i] = iorder[j];
				i = j;
				goto Loop;
			}
		}
		t[i] = ddum;
		iorder[i] = indxin;
		/*     Put the least member in t(n). */
		t[n] = out;
		iorder[n] = indxou;
	}
	return;
} /* hpsolb */
/* ====================== The end of hpsolb ============================== */

void dcsrch(double *f, double *g, double *stp,
		/*Chgd: the next five are no longer pointers:*/
		double ftol, double gtol, double xtol,
		double stpmin, double stpmax,
		char *task, int *isave, double *dsave)
{
	/*    **********

	 Subroutine dcsrch

	 This subroutine finds a step that satisfies a sufficient
	 decrease condition and a curvature condition.

	 Each call of the subroutine updates an interval with
	 endpoints stx and sty. The interval is initially chosen
	 so that it contains a minimizer of the modified function

	 psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).

	 If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
	 interval is chosen so that it contains a minimizer of f.

	 The algorithm is designed to find a step that satisfies
	 the sufficient decrease condition

	 f(stp) <= f(0) + ftol*stp*f'(0),

	 and the curvature condition

	 abs(f'(stp)) <= gtol*abs(f'(0)).

	 If ftol is less than gtol and if, for example, the function
	 is bounded below, then there is always a step which satisfies
	 both conditions.

	 If no step can be found that satisfies both conditions, then
	 the algorithm stops with a warning. In this case stp only
	 satisfies the sufficient decrease condition.

	 A typical invocation of dcsrch has the following outline:

	 task = 'START'
	 10 continue
	 call dcsrch( ... )
	 if (task .eq. 'FG') then
	 Evaluate the function and the gradient at stp
	 goto 10
	 end if

	 NOTE: The user must no alter work arrays between calls.

	 The subroutine statement is

	 subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
	 task,isave,dsave)
	 where

	 f is a double precision variable.
	 On initial entry f is the value of the function at 0.
	 On subsequent entries f is the value of the
	 function at stp.
	 On exit f is the value of the function at stp.

	 g is a double precision variable.
	 On initial entry g is the derivative of the function at 0.
	 On subsequent entries g is the derivative of the
	 function at stp.
	 On exit g is the derivative of the function at stp.

	 stp is a double precision variable.
	 On entry stp is the current estimate of a satisfactory
	 step. On initial entry, a positive initial estimate
	 must be provided.
	 On exit stp is the current estimate of a satisfactory step
	 if task = 'FG'. If task = 'CONV' then stp satisfies
	 the sufficient decrease and curvature condition.

	 ftol is a double precision variable.
	 On entry ftol specifies a nonnegative tolerance for the
	 sufficient decrease condition.
	 On exit ftol is unchanged.

	 gtol is a double precision variable.
	 On entry gtol specifies a nonnegative tolerance for the
	 curvature condition.
	 On exit gtol is unchanged.

	 xtol is a double precision variable.
	 On entry xtol specifies a nonnegative relative tolerance
	 for an acceptable step. The subroutine exits with a
	 warning if the relative difference between sty and stx
	 is less than xtol.
	 On exit xtol is unchanged.

	 stpmin is a double precision variable.
	 On entry stpmin is a nonnegative lower bound for the step.
	 On exit stpmin is unchanged.

	 stpmax is a double precision variable.
	 On entry stpmax is a nonnegative upper bound for the step.
	 On exit stpmax is unchanged.

	 task is a character variable of length at least 60.
	 On initial entry task must be set to 'START'.
	 On exit task indicates the required action:

	 If task(1:2) = 'FG' then evaluate the function and
	 derivative at stp and call dcsrch again.

	 If task(1:4) = 'CONV' then the search is successful.

	 If task(1:4) = 'WARN' then the subroutine is not able
	 to satisfy the convergence conditions. The exit value of
	 stp contains the best point found during the search.

	 If task(1:5) = 'ERROR' then there is an error in the
	 input arguments.

	 On exit with convergence, a warning or an error, the
	 variable task contains additional information.

	 isave is an integer work array of dimension 2.

	 dsave is a double precision work array of dimension 13.

	 Subprograms called

	 MINPACK-2 ... dcstep


	 MINPACK-1 Project. June 1983.
	 Argonne National Laboratory.
	 Jorge J. More' and David J. Thuente.

	 MINPACK-2 Project. October 1993.
	 Argonne National Laboratory and University of Minnesota.
	 Brett M. Averick, Richard G. Carter, and Jorge J. More'.

	 **********
	 */

	/* Local variables */
	int stage;
	double finit, ginit, width, ftest, gtest, stmin, stmax, width1, fm,
	gm, fx, fy, gx, gy;
	int brackt;
	double fxm, fym, gxm, gym, stx, sty;

	/* Parameter adjustments */
	--dsave;
	--isave;

	/* Function Body */

	/*     Initialization block. */
	if (strncmp(task, "START", 5) == 0) {
		/*      Check the input arguments for errors. */
		if (*stp < stpmin)    strcpy(task, "ERROR: STP .LT. STPMIN");
		if (*stp > stpmax)    strcpy(task, "ERROR: STP .GT. STPMAX");
		if (*g >= 0.)        strcpy(task, "ERROR: INITIAL G .GE. ZERO");
		if (ftol < 0.)        strcpy(task, "ERROR: FTOL .LT. ZERO");
		if (gtol < 0.)        strcpy(task, "ERROR: GTOL .LT. ZERO");
		if (xtol < 0.)        strcpy(task, "ERROR: XTOL .LT. ZERO");
		if (stpmin < 0.)    strcpy(task, "ERROR: STPMIN .LT. ZERO");
		if (stpmax < stpmin)    strcpy(task, "ERROR: STPMAX .LT. STPMIN");

		/*      Exit if there are errors on input. */
		if (strncmp(task, "ERROR", 5) == 0) {
			return;
		}
		/*      Initialize local variables. */
		brackt = 0;
		stage = 1;
		finit = *f;
		ginit = *g;
		gtest = ftol * ginit;
		width = stpmax - stpmin;
		width1 = width / .5;
		/*      The variables stx, fx, gx contain the values of the step, */
		/*      function, and derivative at the best step. */
		/*      The variables sty, fy, gy contain the value of the step, */
		/*      function, and derivative at sty. */
		/*      The variables stp, f, g contain the values of the step, */
		/*      function, and derivative at stp. */
		stx = 0.;    fx = finit;    gx = ginit;
		sty = 0.;    fy = finit;    gy = ginit;
		stmin = 0.;
		stmax = *stp + *stp * 4.;
		strcpy(task, "FG");
		goto L1000;
	} else {
		/*      Restore local variables. */
		if (isave[1] == 1) {
			brackt = 1;
		} else {
			brackt = 0;
		}
		stage = isave[2];
		ginit = dsave[1];
		gtest = dsave[2];
		gx = dsave[3];
		gy = dsave[4];
		finit = dsave[5];
		fx = dsave[6];
		fy = dsave[7];
		stx = dsave[8];
		sty = dsave[9];
		stmin = dsave[10];
		stmax = dsave[11];
		width = dsave[12];
		width1 = dsave[13];
	}
	/*      If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
	/*      algorithm enters the second stage. */
	ftest = finit + *stp * gtest;
	if (stage == 1 && *f <= ftest && *g >= 0.) {
		stage = 2;
	}
	/*    Test for warnings. */
	if (brackt && (*stp <= stmin || *stp >= stmax))
		strcpy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS");
	if (brackt && stmax - stmin <= xtol * stmax)
		strcpy(task, "WARNING: XTOL TEST SATISFIED");
	if (*stp == stpmax && *f <= ftest && *g <= gtest)
		strcpy(task, "WARNING: STP = STPMAX");
	if (*stp == stpmin && (*f > ftest || *g >= gtest))
		strcpy(task, "WARNING: STP = STPMIN");
	/*    Test for convergence. */
	if (*f <= ftest && fabs(*g) <= gtol * (-ginit))
		strcpy(task, "CONVERGENCE");
	/*    Test for termination. */
	if (strncmp(task, "WARN", 4) == 0 || strncmp(task, "CONV", 4) == 0)
		goto L1000;

	/*     A modified function is used to predict the step during the */
	/*     first stage if a lower function value has been obtained but */
	/*     the decrease is not sufficient. */
	if (stage == 1 && *f <= fx && *f > ftest) {
		/*      Define the modified function and derivative values. */
		fm = *f - *stp * gtest;
		fxm = fx - stx * gtest;
		fym = fy - sty * gtest;
		gm = *g - gtest;
		gxm = gx - gtest;
		gym = gy - gtest;
		/*      Call dcstep to update stx, sty, and to compute the new step. */
		dcstep(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &
				stmin, &stmax);
		/*      Reset the function and derivative values for f. */
		fx = fxm + stx * gtest;
		fy = fym + sty * gtest;
		gx = gxm + gtest;
		gy = gym + gtest;
	} else {
		/*     Call dcstep to update stx, sty, and to compute the new step. */
		dcstep(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &
				stmax);
	}
	/*     Decide if a bisection step is needed. */
	if (brackt) {
		if (fabs(sty - stx) >= width1 * .66) {
			*stp = stx + (sty - stx) * .5;
		}
		width1 = width;
		width = fabs(sty - stx);
	}
	/*     Set the minimum and maximum steps allowed for stp. */
	if (brackt) {
		stmin = min(stx,sty);
		stmax = max(stx,sty);
	} else {
		stmin = *stp + (*stp - stx) * 1.1;
		stmax = *stp + (*stp - stx) * 4.;
	}
	/*     Force the step to be within the bounds stpmax and stpmin. */
	if(*stp < stpmin) *stp = stpmin;
	if(*stp > stpmax) *stp = stpmax;

	/*     If further progress is not possible, let stp be the best */
	/*     point obtained during the search. */
	if ((brackt && (*stp <= stmin || *stp >= stmax)) ||
			(brackt && (stmax - stmin <= xtol * stmax))) {
		*stp = stx;
	}
	/*     Obtain another function and derivative. */
	strcpy(task, "FG");
	L1000:
	/*     Save local variables. */
	if (brackt) {
		isave[1] = 1;
	} else {
		isave[1] = 0;
	}
	isave[2] = stage;
	dsave[1] = ginit;
	dsave[2] = gtest;
	dsave[3] = gx;
	dsave[4] = gy;
	dsave[5] = finit;
	dsave[6] = fx;
	dsave[7] = fy;
	dsave[8] = stx;
	dsave[9] = sty;
	dsave[10] = stmin;
	dsave[11] = stmax;
	dsave[12] = width;
	dsave[13] = width1;
	return;
} /* dcsrch */
/* ====================== The end of dcsrch ============================== */


void dcstep(double *stx, double *fx, double *dx,
		double *sty, double *fy, double *dy, double *stp,
		double *fp, double *dp, int *brackt, double *stpmin,
		double *stpmax)
{
	/*    **********

	 Subroutine dcstep

	 This subroutine computes a safeguarded step for a search
	 procedure and updates an interval that contains a step that
	 satisfies a sufficient decrease and a curvature condition.

	 The parameter stx contains the step with the least function
	 value. If brackt is set to .true. then a minimizer has
	 been bracketed in an interval with endpoints stx and sty.
	 The parameter stp contains the current step.
	 The subroutine assumes that if brackt is set to .true. then

	 min(stx,sty) < stp < max(stx,sty),

	 and that the derivative at stx is negative in the direction
	 of the step.

	 The subroutine statement is

	 subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
	 stpmin,stpmax)

	 where

	 stx is a double precision variable.
	 On entry stx is the best step obtained so far and is an
	 endpoint of the interval that contains the minimizer.
	 On exit stx is the updated best step.

	 fx is a double precision variable.
	 On entry fx is the function at stx.
	 On exit fx is the function at stx.

	 dx is a double precision variable.
	 On entry dx is the derivative of the function at
	 stx. The derivative must be negative in the direction of
	 the step, that is, dx and stp - stx must have opposite
	 signs.
	 On exit dx is the derivative of the function at stx.

	 sty is a double precision variable.
	 On entry sty is the second endpoint of the interval that
	 contains the minimizer.
	 On exit sty is the updated endpoint of the interval that
	 contains the minimizer.

	 fy is a double precision variable.
	 On entry fy is the function at sty.
	 On exit fy is the function at sty.

	 dy is a double precision variable.
	 On entry dy is the derivative of the function at sty.
	 On exit dy is the derivative of the function at the exit sty.

	 stp is a double precision variable.
	 On entry stp is the current step. If brackt is set to .true.
	 then on input stp must be between stx and sty.
	 On exit stp is a new trial step.

	 fp is a double precision variable.
	 On entry fp is the function at stp
	 On exit fp is unchanged.

	 dp is a double precision variable.
	 On entry dp is the the derivative of the function at stp.
	 On exit dp is unchanged.

	 brackt is an logical variable.
	 On entry brackt specifies if a minimizer has been bracketed.
	 Initially brackt must be set to .false.
	 On exit brackt specifies if a minimizer has been bracketed.
	 When a minimizer is bracketed brackt is set to .true.

	 stpmin is a double precision variable.
	 On entry stpmin is a lower bound for the step.
	 On exit stpmin is unchanged.

	 stpmax is a double precision variable.
	 On entry stpmax is an upper bound for the step.
	 On exit stpmax is unchanged.

	 MINPACK-1 Project. June 1983
	 Argonne National Laboratory.
	 Jorge J. More' and David J. Thuente.

	 MINPACK-2 Project. October 1993.
	 Argonne National Laboratory and University of Minnesota.
	 Brett M. Averick and Jorge J. More'.

	 **********
	 */

	/* System generated locals */
	double d__1, d__2;

	/* Local variables */
	double sgnd, stpc, stpf, stpq, p, q, gamm, r__, s, theta;

	sgnd = *dp * (*dx / fabs(*dx));
	/*     First case: A higher function value. The minimum is bracketed. */
	/*     If the cubic step is closer to stx than the quadratic step, the */
	/*     cubic step is taken, otherwise the average of the cubic and */
	/*     quadratic steps is taken. */
	if (*fp > *fx) {
		theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
		/* Computing MAX */
		d__1 = fabs(theta), d__2 = fabs(*dx),
				d__1 = max(d__1,d__2), d__2 = fabs(*dp);
		s = max(d__1,d__2);
		/* Computing 2nd power */
		d__1 = theta / s;
		gamm = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
		if (*stp < *stx) {
			gamm = -gamm;
		}
		p = gamm - *dx + theta;
		q = gamm - *dx + gamm + *dp;
		r__ = p / q;
		stpc = *stx + r__ * (*stp - *stx);
		stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2. * (*stp
				- *stx);
		if (fabs(stpc - *stx) < fabs(stpq - *stx)) {
			stpf = stpc;
		} else {
			stpf = stpc + (stpq - stpc) / 2.;
		}
		*brackt = 1;
		/*     Second case: A lower function value and derivatives of opposite */
		/*     sign. The minimum is bracketed. If the cubic step is farther from */
		/*     stp than the secant step, the cubic step is taken, otherwise the */
		/*     secant step is taken. */
	} else if (sgnd < 0.) {
		theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
		/* Computing MAX */
		d__1 = fabs(theta), d__2 = fabs(*dx),
				d__1 = max(d__1,d__2), d__2 = fabs(*dp);
		s = max(d__1,d__2);
		/* Computing 2nd power */
		d__1 = theta / s;
		gamm = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
		if (*stp > *stx) {
			gamm = -gamm;
		}
		p = gamm - *dp + theta;
		q = gamm - *dp + gamm + *dx;
		r__ = p / q;
		stpc = *stp + r__ * (*stx - *stp);
		stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
		if (fabs(stpc - *stp) > fabs(stpq - *stp)) {
			stpf = stpc;
		} else {
			stpf = stpq;
		}
		*brackt = 1;
		/*     Third case: A lower function value, derivatives of the same sign, */
		/*     and the magnitude of the derivative decreases. */
	} else if (fabs(*dp) < fabs(*dx)) {
		/*      The cubic step is computed only if the cubic tends to infinity */
		/*      in the direction of the step or if the minimum of the cubic */
		/*      is beyond stp. Otherwise the cubic step is defined to be the */
		/*      secant step. */
		theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
		/* Computing MAX */
		d__1 = fabs(theta), d__2 = fabs(*dx),
				d__1 = max(d__1,d__2), d__2 = fabs(*dp);
		s = max(d__1,d__2);
		/*      The case gamm = 0 only arises if the cubic does not tend */
		/*      to infinity in the direction of the step. */
		/* Computing MAX */
		/* Computing 2nd power */
		d__1 = theta / s;
		d__1 = d__1 * d__1 - *dx / s * (*dp / s);
		gamm = d__1 < 0 ? 0. : s * sqrt(d__1);
		if (*stp > *stx) {
			gamm = -gamm;
		}
		p = gamm - *dp + theta;
		q = gamm + (*dx - *dp) + gamm;
		r__ = p / q;
		if (r__ < 0. && gamm != 0.) {
			stpc = *stp + r__ * (*stx - *stp);
		} else if (*stp > *stx) {
			stpc = *stpmax;
		} else {
			stpc = *stpmin;
		}
		stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
		if (*brackt) {
			/*         A minimizer has been bracketed. If the cubic step is */
			/*         closer to stp than the secant step, the cubic step is */
			/*         taken, otherwise the secant step is taken. */
			if (fabs(stpc - *stp) < fabs(stpq - *stp)) {
				stpf = stpc;
			} else {
				stpf = stpq;
			}
			d__1 = *stp + (*sty - *stp) * .66;
			if (*stp > *stx) {
				stpf = min(d__1,stpf);
			} else {
				stpf = max(d__1,stpf);
			}
		} else {
			/*         A minimizer has not been bracketed. If the cubic step is */
			/*         farther from stp than the secant step, the cubic step is */
			/*         taken, otherwise the secant step is taken. */
			if (fabs(stpc - *stp) > fabs(stpq - *stp)) {
				stpf = stpc;
			} else {
				stpf = stpq;
			}
			stpf = min(*stpmax,stpf);
			stpf = max(*stpmin,stpf);
		}
		/*     Fourth case: A lower function value, derivatives of the */
		/*     same sign, and the magnitude of the derivative does not */
		/*     decrease. If the minimum is not bracketed, the step is either */
		/*     stpmin or stpmax, otherwise the cubic step is taken. */
	} else {
		if (*brackt) {
			theta = (*fp - *fy) * 3. / (*sty - *stp) + *dy + *dp;
			/* Computing MAX */
			d__1 = fabs(theta), d__2 = fabs(*dy), d__1 = max(d__1,d__2), d__2 =
					fabs(*dp);
			s = max(d__1,d__2);
			/* Computing 2nd power */
			d__1 = theta / s;
			gamm = s * sqrt(d__1 * d__1 - *dy / s * (*dp / s));
			if (*stp > *sty) {
				gamm = -gamm;
			}
			p = gamm - *dp + theta;
			q = gamm - *dp + gamm + *dy;
			r__ = p / q;
			stpc = *stp + r__ * (*sty - *stp);
			stpf = stpc;
		} else if (*stp > *stx) {
			stpf = *stpmax;
		} else {
			stpf = *stpmin;
		}
	}
	/*     Update the interval which contains a minimizer. */
	if (*fp > *fx) {
		*sty = *stp;
		*fy = *fp;
		*dy = *dp;
	} else {
		if (sgnd < 0.) {
			*sty = *stx;
			*fy = *fx;
			*dy = *dx;
		}
		*stx = *stp;
		*fx = *fp;
		*dx = *dp;
	}
	/*     Compute the new step. */
	*stp = stpf;
	return;
} /* dcstep */
/* ====================== The end of dcstep ============================== */


void prn3lb(int n, double *x, double *f, char *task, int iprint,
		int info, int iter, int nfgv, int nintol, int nskip,
		int nact, double sbgnrm, int nint,
		char *word, int iback, double stp, double xstep,
		int k)
{
	if(strncmp(task, "CONV", 4) == 0) {
		if (iprint >= 0) {
			cout << endl;
			cout << "iterations " << iter << endl;
			cout << "function evaluations " << nfgv << endl;
			cout << "segments explored during Cauchy searches " << nintol << endl;
			cout << "BFGS updates skipped " << nskip << endl;
			cout << "active bounds at final generalized Cauchy point " << nact << endl;
			cout << "norm of the final projected gradient " << sbgnrm << endl;
			cout << "inal function value " << *f << endl;
			cout << endl;
		}
		if (iprint >= 100) pvector((char*)"X =", x, n);
		if (iprint >= 1)
			cout << "F = " << *f << endl;
	}
	if (iprint >= 0) {
		switch(info) {
		case -1:
			cout << "Matrix in 1st Cholesky factorization in formk is not Pos. Def.";
			break;
		case -2:
			cout << "Matrix in 2st Cholesky factorization in formk is not Pos. Def.";
			break;
		case -3:
			cout << "Matrix in the Cholesky factorization in formt is not Pos. Def.";
			break;
		case -4:
			cout << "Derivative >= 0, backtracking line search impossible.";
			break;
		case -5:
			cout << "l(" << k << ") > u(" << k << ").  No feasible solution";
			break;
		case -6:
			cout << "Input nbd(" << k << ") is invalid";
			break;
		case -7:
			cout << "Warning:  more than 10 function and gradient evaluations" << endl;
			cout << "   in the last line search" << endl;
			break;
		case -8:
			cout << "The triangular system is singular." << endl;
			break;
		case -9:
			cout << "Line search cannot locate an adequate point after 20 function" << endl;
			cout << "and gradient evaluations" << endl;
			break;
		default:
			break;
		}
	}
}

void prn1lb(int n, int m, double *l, double *u, double *x,
		int iprint, double epsmch)
{
	if (iprint >=  0) {
		cout << "N = " << n << ", M = " << m << " machine precision = " << epsmch << endl;
		if (iprint >= 100){
			pvector((char*)"L =", l, n);
			pvector((char*)"X0 =",x, n);
			pvector((char*)"U =", u, n);
		}
	}
}

void projgr(int n, double *l, double *u,
		int *nbd, double *x, double *g, double *sbgnrm)
{
	/*    ************

	 Subroutine projgr

	 This subroutine computes the infinity norm of the projected gradient.


	 *     *  *

	 NEOS, November 1994. (Latest revision April 1997.)
	 Optimization Technology Center.
	 Argonne National Laboratory and Northwestern University.
	 Written by
	 Ciyou Zhu
	 in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

	 ************
	 */
	int i;
	double gi, d__1;

	*sbgnrm = 0.;
	for (i = 0; i < n; ++i) {
		gi = g[i];
		if (nbd[i] != 0) {
			if (gi < 0.) {
				if (nbd[i] >= 2) {
					if(gi < (d__1 = x[i] - u[i]))
						gi = d__1;
				}
			} else {
				if (nbd[i] <= 2) {
					if(gi > (d__1 = x[i] - l[i]))
						gi = d__1;
				}
			}
		}
		if(*sbgnrm < (d__1 = fabs(gi))) *sbgnrm = d__1;
	}
	return;
} /* projgr */

void prn2lb(int n, double *x, double *f, double *g, int iprint,
		int iter, int nfgv, int nact, double sbgnrm,
		int nint, char *word, int iword, int iback,
		double stp, double xstep)
{
	if (iprint >=  99) {
		cout << "LINE SEARCH " << iback << " times; norm of step = " << xstep << "\n";
		if (iprint > 100) {
			pvector((char*)"X =", x, n);
			pvector((char*)"G =", g, n);
		}
	} else if (iprint > 0 && iter%iprint == 0) {
		cout << "At iterate " << iter << "  f = " << *f << "  |proj g|=  " << sbgnrm << "\n";
	}
}

void pvector(char *title, double *x, int n)
{
	int i;
	cout << title;
	for (i = 0; i < n; i++) cout << x[i] << " ";
	cout << endl;
}


int dcopy(int *n, double *dx, int *incx,
		double *dy, int *incy)
{


	/* System generated locals */
	// int i__1;

	/* Local variables */
	int i, m, ix, iy, mp1;


	/*     copies a vector, x, to a vector, y.
	 uses unrolled loops for increments equal to one.
	 jack dongarra, linpack, 3/11/78.
	 modified 12/3/93, array(1) declarations changed to array(*)



	 Parameter adjustments
	 Function Body */


	if (*n <= 0) {
		return 0;
	}
	if (*incx == 1 && *incy == 1) {
		goto L20;
	}

	/*        code for unequal increments or equal increments
	 not equal to 1 */

	ix = 1;
	iy = 1;
	if (*incx < 0) {
		ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
		iy = (-(*n) + 1) * *incy + 1;
	}
	// i__1 = *n;
	for (i = 1; i <= *n; ++i) {
		DY(iy) = DX(ix);
		ix += *incx;
		iy += *incy;
		/* L10: */
	}
	return 0;

	/*        code for both increments equal to 1


	 clean-up loop */

	L20:
	m = *n % 7;
	if (m == 0) {
		goto L40;
	}
	// i__1 = m;
	for (i = 1; i <= m; ++i) {
		DY(i) = DX(i);
		/* L30: */
	}
	if (*n < 7) {
		return 0;
	}
	L40:
	mp1 = m + 1;
	// i__1 = *n;
	for (i = mp1; i <= *n; i += 7) {
		DY(i) = DX(i);
		DY(i + 1) = DX(i + 1);
		DY(i + 2) = DX(i + 2);
		DY(i + 3) = DX(i + 3);
		DY(i + 4) = DX(i + 4);
		DY(i + 5) = DX(i + 5);
		DY(i + 6) = DX(i + 6);
		/* L50: */
	}
	return 0;
} /* dcopy */

void timer(double * ttime)
{
	*ttime = 0.0;
}

int dscal(int *n, double *da, double *dx,
		int *incx)
{


	/* System generated locals */
	// int i__1, i__2;

	/* Local variables */
	int i, m, nincx, mp1;


	/*     scales a vector by a constant.
	 uses unrolled loops for increment equal to one.
	 jack dongarra, linpack, 3/11/78.
	 modified 3/93 to return if incx .le. 0.
	 modified 12/3/93, array(1) declarations changed to array(*)



	 Parameter adjustments
	 Function Body */


	if (*n <= 0 || *incx <= 0) {
		return 0;
	}
	if (*incx == 1) {
		goto L20;
	}

	/*        code for increment not equal to 1 */

	nincx = *n * *incx;
	// i__1 = nincx;
	// i__2 = *incx;
	for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
		DX(i) = *da * DX(i);
		/* L10: */
	}
	return 0;

	/*        code for increment equal to 1


	 clean-up loop */

	L20:
	m = *n % 5;
	if (m == 0) {
		goto L40;
	}
	// i__2 = m;
	for (i = 1; i <= m; ++i) {
		DX(i) = *da * DX(i);
		/* L30: */
	}
	if (*n < 5) {
		return 0;
	}
	L40:
	mp1 = m + 1;
	// i__2 = *n;
	for (i = mp1; i <= *n; i += 5) {
		DX(i) = *da * DX(i);
		DX(i + 1) = *da * DX(i + 1);
		DX(i + 2) = *da * DX(i + 2);
		DX(i + 3) = *da * DX(i + 3);
		DX(i + 4) = *da * DX(i + 4);
		/* L50: */
	}
	return 0;
} /* dscal */

double ddot(int *n, double *dx, int *incx, double *dy,
		int *incy)
{


	/* System generated locals */
	// int i__1;
	double ret_val;

	/* Local variables */
	int i, m;
	double dtemp;
	int ix, iy, mp1;


	/*     forms the dot product of two vectors.
	 uses unrolled loops for increments equal to one.
	 jack dongarra, linpack, 3/11/78.
	 modified 12/3/93, array(1) declarations changed to array(*)



	 Parameter adjustments
	 Function Body */


	ret_val = 0.;
	dtemp = 0.;
	if (*n <= 0) {
		return ret_val;
	}
	if (*incx == 1 && *incy == 1) {
		goto L20;
	}

	/*        code for unequal increments or equal increments
	 not equal to 1 */

	ix = 1;
	iy = 1;
	if (*incx < 0) {
		ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
		iy = (-(*n) + 1) * *incy + 1;
	}
	// i__1 = *n;
	for (i = 1; i <= *n; ++i) {
		dtemp += DX(ix) * DY(iy);
		ix += *incx;
		iy += *incy;
		/* L10: */
	}
	ret_val = dtemp;
	return ret_val;

	/*        code for both increments equal to 1


	 clean-up loop */

	L20:
	m = *n % 5;
	if (m == 0) {
		goto L40;
	}
	// i__1 = m;
	for (i = 1; i <= m; ++i) {
		dtemp += DX(i) * DY(i);
		/* L30: */
	}
	if (*n < 5) {
		goto L60;
	}
	L40:
	mp1 = m + 1;
	// i__1 = *n;
	for (i = mp1; i <= *n; i += 5) {
		dtemp = dtemp + DX(i) * DY(i) + DX(i + 1) * DY(i + 1) + DX(i + 2) *
				DY(i + 2) + DX(i + 3) * DY(i + 3) + DX(i + 4) * DY(i + 4);
		/* L50: */
	}
	L60:
	ret_val = dtemp;
	return ret_val;
} /* ddot */

int daxpy(int *n, double *da, double *dx,
		int *incx, double *dy, int *incy)
{


	/* System generated locals */
	// int i__1;

	/* Local variables */
	int i, m, ix, iy, mp1;


	/*     constant times a vector plus a vector.
	 uses unrolled loops for increments equal to one.
	 jack dongarra, linpack, 3/11/78.
	 modified 12/3/93, array(1) declarations changed to array(*)



	 Parameter adjustments
	 Function Body */    

	if (*n <= 0) {
		return 0;
	}
	if (*da == 0.) {
		return 0;
	}
	if (*incx == 1 && *incy == 1) {
		goto L20;
	}

	/*        code for unequal increments or equal increments
	 not equal to 1 */

	ix = 1;
	iy = 1;
	if (*incx < 0) {
		ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
		iy = (-(*n) + 1) * *incy + 1;
	}
	// i__1 = *n;
	for (i = 1; i <= *n; ++i) {
		DY(iy) += *da * DX(ix);
		ix += *incx;
		iy += *incy;
		/* L10: */
	}
	return 0;

	/*        code for both increments equal to 1


	 clean-up loop */

	L20:
	m = *n % 4;
	if (m == 0) {
		goto L40;
	}
	// i__1 = m;
	for (i = 1; i <= m; ++i) {
		DY(i) += *da * DX(i);
		/* L30: */
	}
	if (*n < 4) {
		return 0;
	}
	L40:
	mp1 = m + 1;
	// i__1 = *n;
	for (i = mp1; i <= *n; i += 4) {
		DY(i) += *da * DX(i);
		DY(i + 1) += *da * DX(i + 1);
		DY(i + 2) += *da * DX(i + 2);
		DY(i + 3) += *da * DX(i + 3);
		/* L50: */
	}
	return 0;
} /* daxpy */

int dpofa(double *a, int *lda, int *n, int *info)
{
	/* System generated locals */
	int a_dim1, a_offset, i__1, i__2, i__3;

	/* Local variables */
	int j, k;
	double s, t;
	int jm1;

	/*<       integer lda,n,info >*/
	/*<       double precision a(lda,1) >*/

	/*     dpofa factors a double precision symmetric positive definite */
	/*     matrix. */

	/*     dpofa is usually called by dpoco, but it can be called */
	/*     directly with a saving in time if  rcond  is not needed. */
	/*     (time for dpoco) = (1 + 18/n)*(time for dpofa) . */

	/*     on entry */

	/*        a       double precision(lda, n) */
	/*                the symmetric matrix to be factored.  only the */
	/*                diagonal and upper triangle are used. */

	/*        lda     integer */
	/*                the leading dimension of the array  a . */

	/*        n       integer */
	/*                the order of the matrix  a . */

	/*     on return */

	/*        a       an upper triangular matrix  r  so that  a = trans(r)*r */
	/*                where  trans(r)  is the transpose. */
	/*                the strict lower triangle is unaltered. */
	/*                if  info .ne. 0 , the factorization is not complete. */

	/*        info    integer */
	/*                = 0  for normal return. */
	/*                = k  signals an error condition.  the leading minor */
	/*                     of order  k  is not positive definite. */

	/*     linpack.  this version dated 08/14/78 . */
	/*     cleve moler, university of new mexico, argonne national lab. */

	/*     subroutines and functions */

	/*     blas ddot */
	/*     fortran dsqrt */

	/*     internal variables */

	/*<       double precision ddot,t >*/
	/*<       double precision s >*/
	/*<       integer j,jm1,k >*/
	/*     begin block with ...exits to 40 */


	/*<          do 30 j = 1, n >*/
	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		/*<             info = j >*/
		*info = j;
		/*<             s = 0.0d0 >*/
		s = 0.;
		/*<             jm1 = j - 1 >*/
		jm1 = j - 1;
		/*<             if (jm1 .lt. 1) go to 20 >*/
		if (jm1 < 1) {
			goto L20;
		}
		/*<             do 10 k = 1, jm1 >*/
		i__2 = jm1;
		for (k = 1; k <= i__2; ++k) {
			/*<                t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1) >*/
			i__3 = k - 1;
			t = a[k + j * a_dim1] - ddot(&i__3, &a[k * a_dim1 + 1], &c__1, &
					a[j * a_dim1 + 1], &c__1);
			/*<                t = t/a(k,k) >*/
			t /= a[k + k * a_dim1];
			/*<                a(k,j) = t >*/
			a[k + j * a_dim1] = t;
			/*<                s = s + t*t >*/
			s += t * t;
			/*<    10       continue >*/
			/* L10: */
		}
		/*<    20       continue >*/
		L20:
		/*<             s = a(j,j) - s >*/
		s = a[j + j * a_dim1] - s;
		/*     ......exit */
		/*<             if (s .le. 0.0d0) go to 40 >*/
		if (s <= 0.) {
			goto L40;
		}
		/*<             a(j,j) = dsqrt(s) >*/
		a[j + j * a_dim1] = sqrt(s);
		/*<    30    continue >*/
		/* L30: */
	}
	/*<          info = 0 >*/
	*info = 0;
	/*<    40 continue >*/
	L40:
	/*<       return >*/
	return 0;
	/*<       end >*/
} /* dpofa */

int dtrsl(double *t, int *ldt, int *n,
		double *b, int *job, int *info)
{
	/* System generated locals */
	int t_dim1, t_offset, i__1, i__2;

	/* Local variables */
	int j, jj, case__;
	double temp;


	/*     dtrsl solves systems of the form */

	/*                   t * x = b */
	/*     or */
	/*                   trans(t) * x = b */

	/*     where t is a triangular matrix of order n. here trans(t) */
	/*     denotes the transpose of the matrix t. */

	/*     on entry */

	/*         t         double precision(ldt,n) */
	/*                   t contains the matrix of the system. the zero */
	/*                   elements of the matrix are not referenced, and */
	/*                   the corresponding elements of the array can be */
	/*                   used to store other information. */

	/*         ldt       integer */
	/*                   ldt is the leading dimension of the array t. */

	/*         n         integer */
	/*                   n is the order of the system. */

	/*         b         double precision(n). */
	/*                   b contains the right hand side of the system. */

	/*         job       integer */
	/*                   job specifies what kind of system is to be solved. */
	/*                   if job is */

	/*                        00   solve t*x=b, t lower triangular, */
	/*                        01   solve t*x=b, t upper triangular, */
	/*                        10   solve trans(t)*x=b, t lower triangular, */
	/*                        11   solve trans(t)*x=b, t upper triangular. */

	/*     on return */

	/*         b         b contains the solution, if info .eq. 0. */
	/*                   otherwise b is unaltered. */

	/*         info      integer */
	/*                   info contains zero if the system is nonsingular. */
	/*                   otherwise info contains the index of */
	/*                   the first zero diagonal element of t. */

	/*     linpack. this version dated 08/14/78 . */
	/*     g. w. stewart, university of maryland, argonne national lab. */

	/*     subroutines and functions */

	/*     blas daxpy,ddot */
	/*     fortran mod */

	/*     internal variables */


	/*     begin block permitting ...exits to 150 */

	/*        check for zero diagonal elements. */

	/* Parameter adjustments */
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;
	--b;

	/* Function Body */
	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
		/*     ......exit */
		if (t[*info + *info * t_dim1] == 0.) {
			goto L150;
		}
		/* L10: */
	}
	*info = 0;

	/*        determine the task and go to it. */

	case__ = 1;
	if (*job % 10 != 0) {
		case__ = 2;
	}
	if (*job % 100 / 10 != 0) {
		case__ += 2;
	}
	switch (case__) {
	case 1:  goto L20;
	case 2:  goto L50;
	case 3:  goto L80;
	case 4:  goto L110;
	}

	/*        solve t*x=b for t lower triangular */

	L20:
	b[1] /= t[t_dim1 + 1];
	if (*n < 2) {
		goto L40;
	}
	i__1 = *n;
	for (j = 2; j <= i__1; ++j) {
		temp = -b[j - 1];
		i__2 = *n - j + 1;
		daxpy(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &c__1);
		b[j] /= t[j + j * t_dim1];
		/* L30: */
	}
	L40:
	goto L140;

	/*        solve t*x=b for t upper triangular. */

	L50:
	b[*n] /= t[*n + *n * t_dim1];
	if (*n < 2) {
		goto L70;
	}
	i__1 = *n;
	for (jj = 2; jj <= i__1; ++jj) {
		j = *n - jj + 1;
		temp = -b[j + 1];
		daxpy(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
		b[j] /= t[j + j * t_dim1];
		/* L60: */
	}
	L70:
	goto L140;

	/*        solve trans(t)*x=b for t lower triangular. */

	L80:
	b[*n] /= t[*n + *n * t_dim1];
	if (*n < 2) {
		goto L100;
	}
	i__1 = *n;
	for (jj = 2; jj <= i__1; ++jj) {
		j = *n - jj + 1;
		i__2 = jj - 1;
		b[j] -= ddot(&i__2, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
		b[j] /= t[j + j * t_dim1];
		/* L90: */
	}
	L100:
	goto L140;

	/*        solve trans(t)*x=b for t upper triangular. */

	L110:
	b[1] /= t[t_dim1 + 1];
	if (*n < 2) {
		goto L130;
	}
	i__1 = *n;
	for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		b[j] -= ddot(&i__2, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
		b[j] /= t[j + j * t_dim1];
		/* L120: */
	}
	L130:
	L140:
	L150:
	return 0;
} /* dtrsl */
