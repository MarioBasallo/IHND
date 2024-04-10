#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include <errno.h>
#include<ilcplex/cplex.h>

void modelM1(void);
void modelM2(void);
void modelM3(void);
void free_memory(void);
void free_and_null(char**);
void initialize_memory(void);
void read_ap(const char*);
void quicksort(double*, int*, int, int);
FILE* open_file(const char*, const char*);

// Create vectors and matrices
int** create_int_matrix(int, int);
double** create_double_matrix(int, int);
double*** create_double_3Dmatrix(int, int, int);
double*** create_int_3Dmatrix(int, int, int);
int* create_int_vector(int);
double* create_double_vector(int);
void i_vector(int** vector, int n, char* s);
void d_vector(double** vector, int n, char* s);
void c_vector(char** vector, int n, char* s);

// Structure for node coordinates
typedef struct {
	double x, y;
} coordinate;

typedef struct xinfo {
	int* pos;
	int* node2;
	int* first;
	int dim;
}xinfo;

typedef struct bitf {
	int* arcpos;	// Arc positions of candidate arcs
	int dim;		// Domension of the ser Ak
	struct xinfo* node;
} bitf;

typedef struct crit2 {
	int* arc2;
	int* com;
	int dim;

	int* ord;

	int* pos;
	double* tau;
	int dimys;
} crit2;

typedef struct crit1 {
	int* arc1;	// Arc positions of candidate arcs
	struct crit2* arc2;
} crit1;

typedef struct COMMOD {
	int* i;
	int* j;
	int  dim;
} COMMOD;

// For storing the feasible hubs
typedef struct FHubs {
	int* hub;	// Feasible hub
	int* pos;	// Position in the ordered list of feasible hubs
	int  dim;	// Number of feasible hubs
} FHubs;


// Cut structure and cut callback function
static int CPXPUBLIC mycutcallback(CPXCENVptr, void*, int, void*, int*);
struct cutinfo {
	CPXLPptr lp;
	int     numcols;
	int     num;
	double* x;
};
typedef struct cutinfo CUTINFO, * CUTINFOptr;


// For storing the feasible inter-hub arcs
typedef struct FHubArcs {
	int* arcpos;
	int** pos;
	int  dim;
} FHubArcs;


// Functions related to the total sojourn time distribution
struct Wait_Parmeters {
	double mu1, mu2, alpha, T;
};
double Wait(double, void*);
double WaitDeriv(double, void*);
double Service(double, double, void*);


// Functions related to the total service time distributions
struct Integrand_Parameters {
	double lambda1, lambda2, mu1, mu2, T, a, b;
};
double Integrand(double x, void* params);
struct Convolve_Parameters {
	double lambda1, lambda2, mu1, mu2, T, a, b;
};
double Convolve(int n, double (*f)(double, void*), void* params);
struct Service1Dist_Parmeters {
	double lambda2, mu1, mu2, alpha, T, a, b, n;
};
double Service1Dist(double x, void* params);
struct feval_Parameters {
	double mu1, mu2, alpha, T, a, b, n;
};
double feval(double x, void* params);
double fderiv(double x, void* params);



// Computation of the singularity point of the alpha level set of the total service time distribution
// ==================================================================================================
struct IntegrandSing_Parameters {
	double lambda2, mu2, T, a, b;
};
double IntegrandSing(double x, void* params);
struct ConvolveSing_Parameters {
	double lambda2, mu2, T, a, b;
};
double ConvolveSing(int n, double (*f)(double, void*), void* params);
double Singeval(void* params);
struct Singeval_Parameters {
	double mu2, alpha, T, a, b, n;
};
struct SingDist_Parmeters {
	double mu2, alpha, T, a, b, n;
};
double SingDist(double x, void* params);


// Algorithm 2
// ===========
double TauEval(double x, void* params);
struct Tau_Params {
	double lambda1, lambda2, mu1, mu2, alpha;
};
double TauRoot(void* params);
struct Tau_Params2 {
	double lambda1, lambda2, mu1, mu2, alpha, tauref;
};

double TauSel(void* params);
struct TauSelParams {
	double mu1, mu2, alpha, tauref, apar, bpar, lmax1, lmax2;
};
