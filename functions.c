#include "header.h"
#include <gsl/gsl_sf_lambert.h>
#include "gauss_legendre.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>

// FUNCTIONS RELATED TO THE TOTAL SOJOURN TIME DISTRIBUTION W
// ==========================================================

// Exact evaluation of the alpha level set of the total sojourn time distribution
double Wait(double x, void* params) {
	// This function is used to evaluate the closed form expression for f^[a] given in Equation (49):
	// lambda_[a1] = f^[a] (lambda_[a2])

	struct Wait_Parmeters* p;
	p = (struct Wait_Parmeters*)params;
	double mu1 = p->mu1;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double T = p->T;
	double l1;

	double maxf;
	if (fabs(mu1 - mu2) > pow(10, -11)) {
		maxf = (1 - alpha) - mu1 * exp(-mu2 * T) / (mu1 - mu2) - mu2 * exp(-mu1 * T) / (mu2 - mu1);
	}
	else {
		maxf = 1 - alpha - exp(-mu1 * T) - mu1 * T * exp(-mu1 * T);
	}
	if (maxf > 0) {
		// Abcisa of the singularity point (Equation 48)
		double lsing;
		lsing = (1 / T) * (1 + mu2 * T + gsl_sf_lambert_Wm1((alpha - 1) / exp(1)));

		if ((x - lsing < -pow(10, -10)) && (((mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) *
			exp((1 - alpha) * (mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) >= -1 / exp(1) + pow(10, -10))) {
			l1 = mu1 + (1 - alpha) * (mu2 - x) / (exp(-(mu2 - x) * T) - (1 - alpha)) - (1 / T) * gsl_sf_lambert_Wm1(((mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) *
				exp((1 - alpha) * (mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))));
		}
		else if ((x - lsing > pow(10, -10)) && (((mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) *
			exp((1 - alpha) * (mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) >= -1 / exp(1) + pow(10, -10))) {
			l1 = mu1 + (1 - alpha) * (mu2 - x) / (exp(-(mu2 - x) * T) - (1 - alpha)) - (1 / T) * gsl_sf_lambert_W0(((mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) *
				exp((1 - alpha) * (mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))));
		}
		else {
			l1 = (1 / T) * (1 + mu1 * T + gsl_sf_lambert_Wm1((alpha - 1) / exp(1)));
		}
	}
	else {
		l1 = 0;
	}
	return l1;
}

// Exact evaluation of the derivative of the alpha level set of the total sojourn time distribution
double WaitDeriv(double x, void* params) {
	// This function is used to evaluate the closed form expression for the derivative of f^[a] given in Equation (50):

	struct Wait_Parmeters* p;
	p = (struct Wait_Parmeters*)params;
	double mu1 = p->mu1;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double T = p->T;
	double l1d;

	double maxf;
	if (fabs(mu1 - mu2) > pow(10, -11)) {
		maxf = (1 - alpha) - mu1 * exp(-mu2 * T) / (mu1 - mu2) - mu2 * exp(-mu1 * T) / (mu2 - mu1);
	}
	else {
		maxf = 1 - alpha - exp(-mu1 * T) - mu1 * T * exp(-mu1 * T);
	}
	if (maxf > 0) {
		// Abcisa of the singularity point (Equation 48)
		double lsing;
		lsing = (1 / T) * (1 + mu2 * T + gsl_sf_lambert_Wm1((alpha - 1) / exp(1)));

		if ((x - lsing < -pow(10, -10)) && (((mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) *
			exp((1 - alpha) * (mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) >= -1 / exp(1) + pow(10, -10))) {
			l1d = (pow((1 - alpha), 2) - (1 - alpha) * (1 + (mu2 - x) * T) * exp(-(mu2 - x) * T)) / pow(exp(-(mu2 - x) * T) - (1 - alpha), 2) +
				gsl_sf_lambert_Wm1(((mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) * exp((1 - alpha) * (mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha)))) /
				(1 + gsl_sf_lambert_Wm1(((mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) * exp((1 - alpha) * (mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))))) *
				(((1 - alpha) - (1 + (mu2 - x) * T) * exp(-(mu2 - x) * T)) * (exp(-(mu2 - x) * T) - (1 - alpha) * (1 - (mu2 - x) * T))) /
				(-T * (mu2 - x) * pow(exp(-(mu2 - x) * T) - (1 - alpha), 2));
		}
		else if ((x - lsing > pow(10, -10)) && (((mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) *
			exp((1 - alpha) * (mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) >= -1 / exp(1) + pow(10, -10))) {
			l1d = (pow((1 - alpha), 2) - (1 - alpha) * (1 + (mu2 - x) * T) * exp(-(mu2 - x) * T)) / pow(exp(-(mu2 - x) * T) - (1 - alpha), 2) +
				gsl_sf_lambert_W0(((mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) * exp((1 - alpha) * (mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha)))) /
				(1 + gsl_sf_lambert_W0(((mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))) * exp((1 - alpha) * (mu2 - x) * T / (exp(-(mu2 - x) * T) - (1 - alpha))))) *
				(((1 - alpha) - (1 + (mu2 - x) * T) * exp(-(mu2 - x) * T)) * (exp(-(mu2 - x) * T) - (1 - alpha) * (1 - (mu2 - x) * T))) /
				(-T * (mu2 - x) * pow(exp(-(mu2 - x) * T) - (1 - alpha), 2));
		}
		else {
			l1d = -1;
		}
	}
	else {
		l1d = 0;
	}
	return l1d;
}


// Evaluation of the total sojourn time distribution
double Service(double l1, double l2, void* params) {
	// This function evaluates the total sojourn time distribution given in Equation (6)

	struct Wait_Parmeters* p;
	p = (struct Wait_Parmeters*)params;
	double mu1 = p->mu1;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double T = p->T;
	double slevel;
	if (fabs((mu1 - l1) - (mu2 - l2)) > pow(10, -9)) {
		slevel = 1 - (mu1 - l1) * exp(-(mu2 - l2) * T) / ((mu1 - l1) - (mu2 - l2)) -
			(mu2 - l2) * exp(-(mu1 - l1) * T) / ((mu2 - l2) - (mu1 - l1));
	}
	else {
		slevel = 1 - exp(-(mu1 - l1) * T) - (mu1 - l1) * T * exp(-(mu1 - l1) * T);
	}
	return slevel;
}

// FUNCTIONS RELATED TO THE TOTAL SERVICE TIME DISTRIBUTION S
// ==========================================================

// Integrand of the convolution integral
double Integrand(double x, void* params) {
	// This function refers to the integrand of the convolution integral given in Equation (12)

	struct Integrand_Parameters* p;
	p = (struct Integrand_Parameters*)params;
	double lambda1 = p->lambda1;	// Total flow at hub a1
	double lambda2 = p->lambda2;	// Total flow at hub a2
	double mu1 = p->mu1;			// Processing rate of hub a1
	double mu2 = p->mu2;			// Processing rate of hub a2
	double T = p->T;				// Maximum service time requirement: tau
	double a = p->a;				// Scale parameter of the gamma density function
	double b = p->b;				// Shape parameter for the gamma density function
	double val = 0;
	if (fabs((mu1 - lambda1) - (mu2 - lambda2)) <= pow(10, -9)) {
		val = (1 - exp(-(mu1 - lambda1) * (T - x)) - (mu1 - lambda1) * (T - x) * exp(-(mu1 - lambda1) * (T - x))) * gsl_ran_gamma_pdf(x, a, b);
	}
	else {
		val = (1 - (mu1 - lambda1) * exp(-(mu2 - lambda2) * (T - x)) / ((mu1 - lambda1) - (mu2 - lambda2)) - (mu2 - lambda2) * exp(-(mu1 - lambda1) * (T - x)) / ((mu2 - lambda2) - (mu1 - lambda1))) *	// hypoexponential part
			gsl_ran_gamma_pdf(x, a, b);
	}
	return val;
}

// Numerical method for computing the confolution integral
double Convolve(int n, double (*Integrand)(double, void*), void* params) {
	// This function uses the Gauss-Legendre Quadrature method to compute the convolution integral in Equation (12)

	struct Convolve_Parameters* p;
	p = (struct Convolve_Parameters*)params;
	double lambda1 = p->lambda1;
	double lambda2 = p->lambda2;
	double mu1 = p->mu1;
	double mu2 = p->mu2;
	double T = p->T;
	double a = p->a;
	double b = p->b;
	double I = 0;
	struct Integrand_Parameters paramsI = { lambda1, lambda2, mu1, mu2, T, a, b };
	I = gauss_legendre(n, Integrand, &paramsI, 0, T);
	return I;
}

// For the evaluation of the alpha level set of the total service time distribution
double Service1Dist(double x, void* params) {
	struct Service1Dist_Parmeters* p;
	p = (struct Service1Dist_Parmeters*)params;
	double lambda2 = p->lambda2;
	double mu1 = p->mu1;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double T = p->T;
	double a = p->a;
	double b = p->b;
	int n = p->n;
	double conv;
	struct Convolve_Parameters paramsC = { x, lambda2, mu1, mu2, T, a, b };
	conv = Convolve(n, Integrand, &paramsC) - alpha;
	return conv;
}

// Root finding method for evaluating the alpha level set of the total service time distribution
double feval(double x, void* params) {
	struct feval_Parameters* p;
	p = (struct feval_Parameters*)params;
	double mu1 = p->mu1;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double T = p->T;
	double a = p->a;
	double b = p->b;
	int n = p->n;

	double l1 = 0, seval = 0;
	struct Convolve_Parameters paramsC = { 0, 0, mu1, mu2, T, a, b };
	seval = Convolve(n, Integrand, &paramsC);
	if (seval >= alpha) {
		int max_iter = 100;

		// Solver type
		const gsl_root_fsolver_type* Type3;
		gsl_root_fsolver* s3;
		Type3 = gsl_root_fsolver_brent;     // Use of the Brendt's algorithm
		s3 = gsl_root_fsolver_alloc(Type3);
		gsl_function F3;

		// Solver parameters
		double x_lo = -mu1;  // Lower limit
		double x_hi = mu1;  // Upper limit

		// Set up the function and the solver
		struct Service1Dist_Parmeters paramsS1 = { x, mu1, mu2, alpha, T, a, b, n };
		F3.function = &Service1Dist;
		F3.params = &paramsS1;
		gsl_root_fsolver_set(s3, &F3, x_lo, x_hi);

		int status = GSL_CONTINUE;
		for (int i = 1; i < max_iter && status == GSL_CONTINUE; i++) {
			// Iterate one step of the solver
			status = gsl_root_fsolver_iterate(s3);
			if (status != GSL_SUCCESS)
				break;

			// Get the solver current best solution and bounds
			l1 = gsl_root_fsolver_root(s3);
			x_lo = gsl_root_fsolver_x_lower(s3);
			x_hi = gsl_root_fsolver_x_upper(s3);
			status = gsl_root_test_interval(x_lo, x_hi, 0, 1e-6);
		}
		gsl_root_fsolver_free(s3);
	}
	return l1;
}

// Numerical method for evaluating the derivative of the alpha level set for the total service time distribution
double fderiv(double x, void* params) {
	struct feval_Parameters* p;
	p = (struct feval_Parameters*)params;
	double mu1 = p->mu1;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double T = p->T;
	double a = p->a;
	double b = p->b;
	int n = p->n;

	gsl_function F;
	double fder, abserr;

	struct feval_Parameters paramsS1 = { mu1, mu2, alpha, T, a, b, n };
	F.function = &feval;
	F.params = &paramsS1;

	gsl_deriv_central(&F, x, 1e-4, &fder, &abserr);
	double f1, f2, e = 1e-2;
	if (fder > 0) {
		f1 = feval(x - e, &paramsS1);
		f2 = feval(x + e, &paramsS1);
		fder = (f2 - f1) / ((x + e) - (x - e));
	}
	if (fder > 0) {
		e = 1e-1;
		f1 = feval(x - e, &paramsS1);
		f2 = feval(x + e, &paramsS1);
		fder = (f2 - f1) / ((x + e) - (x - e));
	}
	if (fder > 0) {
		e = 1;
		f1 = feval(x - e, &paramsS1);
		f2 = feval(x + e, &paramsS1);
		fder = (f2 - f1) / ((x + e) - (x - e));
	}
	if (fder > 0) {
		e = 2;
		f1 = feval(x - e, &paramsS1);
		f2 = feval(x + e, &paramsS1);
		fder = (f2 - f1) / ((x + e) - (x - e));
	}
	return fder;
}

// NUMERICAL EVALUATION OF THE SINGLUARITY POINT OF THE ALPHA LEVEL SET OF THE TOTAL SERVICE TIME DISTRIBUTION
// ===========================================================================================================

// Integrand of the convolution integral for computing the singularity point
double IntegrandSing(double x, void* params) {
	struct IntegrandSing_Parameters* p;
	p = (struct IntegrandSing_Parameters*)params;
	double lambda2 = p->lambda2;
	double mu2 = p->mu2;
	double T = p->T;
	double a = p->a;
	double b = p->b;
	double val = 0;
	val = (1 - exp(-(mu2 - lambda2) * (T - x)) - (mu2 - lambda2) * (T - x) * exp(-(mu2 - lambda2) * (T - x))) * gsl_ran_gamma_pdf(x, a, b);
	return val;
}

// Numerical evaluation of the convolution integrand for the singularity
double ConvolveSing(int n, double (*Integrand)(double, void*), void* params) {
	struct ConvolveSing_Parameters* p;
	p = (struct ConvolveSing_Parameters*)params;
	double lambda2 = p->lambda2;
	double mu2 = p->mu2;
	double T = p->T;
	double a = p->a;
	double b = p->b;
	double I = 0, error = 0;
	struct IntegrandSing_Parameters paramsI = { lambda2, mu2, T, a, b };
	I = gauss_legendre(n, IntegrandSing, &paramsI, 0, T);
	return I;
}

// Function used for finding the root
double SingDist(double x, void* params) {
	struct SingDist_Parmeters* p;
	p = (struct SingDist_Parmeters*)params;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double T = p->T;
	double a = p->a;
	double b = p->b;
	int n = p->n;
	double conv;
	struct ConvolveSing_Parameters paramsC = { x, mu2, T, a, b };
	conv = ConvolveSing(n, IntegrandSing, &paramsC) - alpha;
	return conv;
}

double Singeval(void* params) {
	struct Singeval_Parameters* p;
	p = (struct Singeval_Parameters*)params;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double T = p->T;
	double a = p->a;
	double b = p->b;
	int n = p->n;
	double l1 = 0;
	int max_iter = 100;

	// Solver type
	const gsl_root_fsolver_type* Type3;
	gsl_root_fsolver* s3;
	Type3 = gsl_root_fsolver_brent;     // Use of the Brendt's algorithm
	s3 = gsl_root_fsolver_alloc(Type3);
	gsl_function F3;

	// Solver parameters
	double x_lo = -mu2;  // Lower limit
	double x_hi = mu2;  // Upper limit

	// Set up the function and the solver
	struct SingDist_Parmeters paramsC = { mu2, alpha, T, a, b, n };
	F3.function = &SingDist;
	F3.params = &paramsC;
	gsl_root_fsolver_set(s3, &F3, x_lo, x_hi);

	int status = GSL_CONTINUE;
	for (int i = 1; i < max_iter && status == GSL_CONTINUE; i++) {
		// Iterate one step of the solver
		status = gsl_root_fsolver_iterate(s3);
		if (status != GSL_SUCCESS)
			break;

		// Get the solver current best solution and bounds
		l1 = gsl_root_fsolver_root(s3);
		x_lo = gsl_root_fsolver_x_lower(s3);
		x_hi = gsl_root_fsolver_x_upper(s3);
		status = gsl_root_test_interval(x_lo, x_hi, 0, 1e-6);
	}
	gsl_root_fsolver_free(s3);
	return l1;
}

// FUNCTIONS FOR ALGORITHM 2
// =========================

double TauEval(double x, void* params) {
	struct Tau_Params* p;
	p = (struct Tau_Params*)params;
	double lambda1 = p->lambda1;
	double lambda2 = p->lambda2;
	double mu1 = p->mu1;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double result;
	result = -(1 - alpha) * (mu1 - lambda1 - (mu2 - lambda2)) + (mu1 - lambda1) * exp(-(mu2 - lambda2) * x)
		- (mu2 - lambda2) * exp(-(mu1 - lambda1) * x);
	return result;
}

double TauRoot(void* params) {
	struct Tau_Params2* p;
	p = (struct Tau_Params2*)params;
	double lambda1 = p->lambda1;
	double lambda2 = p->lambda2;
	double mu1 = p->mu1;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double tauref = p->tauref;
	double l1 = 0;
	int max_iter = 100;

	// Solver type
	const gsl_root_fsolver_type* Type3;
	gsl_root_fsolver* s3;
	Type3 = gsl_root_fsolver_brent;     // Use of the Brendt's algorithm
	s3 = gsl_root_fsolver_alloc(Type3);
	gsl_function F3;

	// Solver parameters
	double x_lo = 0;  // Lower limit
	double x_hi = tauref;  // Upper limit

	// Set up the function and the solver
	struct Tau_Params paramsC = { lambda1, lambda2, mu1, mu2, alpha };
	F3.function = &TauEval;
	F3.params = &paramsC;
	gsl_root_fsolver_set(s3, &F3, x_lo, x_hi);

	int status = GSL_CONTINUE;
	for (int i = 1; i < max_iter && status == GSL_CONTINUE; i++) {
		// Iterate one step of the solver
		status = gsl_root_fsolver_iterate(s3);
		if (status != GSL_SUCCESS)
			break;

		// Get the solver current best solution and bounds
		l1 = gsl_root_fsolver_root(s3);
		x_lo = gsl_root_fsolver_x_lower(s3);
		x_hi = gsl_root_fsolver_x_upper(s3);
		status = gsl_root_test_interval(x_lo, x_hi, 0, 1e-6);
	}
	gsl_root_fsolver_free(s3);
	return l1;
}

double TauSel(void* params) { // Algorithm 2
	struct TauSelParams* p;
	p = (struct TauSelParams*)params;
	double mu1 = p->mu1;
	double mu2 = p->mu2;
	double alpha = p->alpha;
	double tauref = p->tauref;
	double apar = p->apar;
	double bpar = p->bpar;
	double lmax1 = p->lmax1;
	double lmax2 = p->lmax2;

	double tau = 0, lsing1, lsing2;

	struct Singeval_Parameters paramssing = { mu2, alpha, tauref, apar, bpar, 40 };
	lsing2 = Singeval(&paramssing);
	lsing1 = mu1 - mu2 + lsing2;
	tau = (1 + gsl_sf_lambert_Wm1((alpha - 1) / exp(1))) / (lsing1 - mu1);

	struct Wait_Parmeters paramsapprox = { mu1, mu2, alpha, tau};
	struct Wait_Parmeters paramsapprox2 = { mu2, mu1, alpha, tau};
	if ((Wait(0, &paramsapprox) < lmax1) || (Wait(0, &paramsapprox2) < lmax2)) {
		if (lmax1 < lmax2) {
			struct Tau_Params2 paramst1 = { lmax1, 0, mu1, mu2, alpha, tauref };
			tau = TauRoot(&paramst1);
		}
		else {
			struct Tau_Params2 paramst1 = { 0, lmax2, mu1, mu2, alpha, tauref };
			tau = TauRoot(&paramst1);
		}
	}

	return tau;
}


// OTHER FUNCTIONS
// ===============

// Sorting algorithm
void quicksort(double* number, int* idx, int first, int last) {
	int i, j, pivot, temp2;
	double temp;
	if (first < last) {
		pivot = first;
		i = first;
		j = last;
		while (i < j) {
			while (number[i] <= number[pivot] && i < last)
				i++;
			while (number[j] > number[pivot])
				j--;
			if (i < j) {
				temp = number[i];
				number[i] = number[j];
				number[j] = temp;

				temp2 = idx[i];
				idx[i] = idx[j];
				idx[j] = temp2;
			}
		}
		temp = number[pivot];
		number[pivot] = number[j];
		number[j] = temp;

		temp2 = idx[pivot];
		idx[pivot] = idx[j];
		idx[j] = temp2;

		quicksort(number, idx, first, j - 1);
		quicksort(number, idx, j + 1, last);
	}
}