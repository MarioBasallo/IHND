#include "header.h"
#include <gsl/gsl_deriv.h>
#include "gauss_legendre.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_lambert.h>


extern double** c2, ** c, * f, * mu, ** W;
extern double** c_c, ** c_t, ** c_d;
extern int N;
extern COMMOD comm;
extern coordinate* pts;
extern bitf* Ak;
extern xinfo* node;
extern crit1* arc1;
extern crit2* arc2;
extern FHubs FeasHubs;

extern double rq, eta, apar, bpar, vel;
extern double alpha, muhat, * tau;
extern double*** lmaxk, * lmax;
extern int* FeasArc, * FeasArcPos;
extern int svars, * nbFeas, ** nCom;
extern int xvars, * xa1, * xa2, * xk;
extern int yvars, * ypos;

void read_ap(const char* name) {
	FILE* in;
	in = open_file(name, "r");
	fscanf(in, "%d", &N);

	comm.dim = N * N;
	initialize_memory();

	// Read node coordinates
	for (int p = 0; p < N; p++) {
		if (fscanf(in, "%lf %lf", &pts[p].x, &pts[p].y) != 2) {
			fprintf(stderr, "ERROR: Can't read coordinates.");
			exit(1);
		}
	}

	// Read the flow matrix
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (fscanf(in, "%lf", &W[i][j]) != 1) {
				fprintf(stderr, "ERROR: Can't flow matrix\n");
				exit(1);
			}
		}
	}

	// Read fixed costs
	for (int i = 0; i < N; i++) {
		if (fscanf(in, "%lf", &f[i]) != 1) {
			fprintf(stderr, "ERROR: Can't read fixed costs.\n");
			exit(1);
		}
		f[i] = f[i] / 2;
	}

	// Read hub capacities
	for (int i = 0; i < N; i++) {
		if (fscanf(in, "%lf", &mu[i]) != 1) {
			fprintf(stderr, "ERROR: Can't read capacities\n");
			exit(1);
		}
	}

	// Compute distance and cost matrices
	double dbar = 0, v = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			c[i][j] = floor(((sqrt(pow((pts[i].x - pts[j].x), 2) + pow((pts[i].y - pts[j].y), 2))) / 1000) * pow(10, 6)) * pow(10, -6);
			c_c[i][j] = 3 * c[i][j];
			c_t[i][j] = 0.75 * c[i][j];
			c_d[i][j] = 2 * c[i][j];
			c2[i][j] = c[i][j];		// Keep this distance matrix for computing transportation time related parameters related
			c[i][j] = 7 * c[i][j];
			dbar += c2[i][j];
		}
	}
	dbar = dbar / (N * (N - 1));
	free(pts);

	// Define the set of commodities K
	comm.dim = 0;
	for (int i = 0; i < N - 1; i++) {
		for (int j = i + 1; j < N; j++) { // Two successive commodities with strictly positive demand are of the form k1 = (i,j) and k2 = (j,i)
			if ((W[i][j] > 0) && (i != j)) {
				comm.i[comm.dim] = i;
				comm.j[comm.dim] = j;
				comm.dim++;
			}
			if ((W[j][i] > 0) && (i != j)) {
				comm.j[comm.dim] = i;
				comm.i[comm.dim] = j;
				comm.dim++;
			}
		}
	}

	// Compute speed of hub and spoke vehicles (Section 7.1)
	muhat = 0;
	for (int i = 0; i < N; i++) {
		muhat += mu[i];
	}
	muhat = muhat / N;
	vel = 0.2 * muhat * dbar;

	// Compute the service time requirement for each commodity (Equation 41)
	tau = (double*)calloc(comm.dim, sizeof(double));
	for (int k = 0; k < comm.dim; k++) {
		apar = 1 / pow(0.5, 2);		// Constant shape parameter of the gamma density function
		bpar = (c2[comm.i[k]][comm.j[k]] / vel) / apar;		// Scale parameter of the gamma density function
		tau[k] = rq * gsl_cdf_gamma_Pinv(0.7, apar, bpar);	
	}

	// Compute upper bounds for flow variables with respect to each commodity and hub arc
	int count, nbkp = 40, * candidate_arc;
	double seval = 0;
	xvars = 0;
	Ak = (bitf*)calloc(comm.dim, sizeof(bitf));
	candidate_arc = create_int_vector(N * N);
	lmaxk = create_double_3Dmatrix(comm.dim, N, N);
	int akpos = 0;
	for (int k = 0; k < comm.dim; k++) {
		count = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				// Preprocessing (Equation 1)
				if ((i != j) && (c_c[comm.i[k]][i] + c_t[i][j] + c_d[j][comm.j[k]] < c_c[comm.i[k]][j] + c_t[j][i] + c_d[i][comm.j[k]]) && (c_c[comm.i[k]][i] + c_t[i][j] + c_d[j][comm.j[k]] < c[comm.i[k]][comm.j[k]])) {
					apar = 1 / pow(0.5, 2);
					bpar = (c2[comm.i[k]][i] + eta * c2[i][j] + c2[j][comm.j[k]]) / (apar * vel);
					struct Convolve_Parameters paramsC = { W[comm.i[k]][comm.j[k]], W[comm.i[k]][comm.j[k]], mu[i], mu[j], tau[k], apar, bpar };
					seval = Convolve(nbkp, Integrand, &paramsC);

					// Feasibility check for hub arc (Equation 16)
					if (seval >= alpha + pow(10, -6)) {

						struct feval_Parameters params = { mu[i], mu[j], alpha, tau[k], apar, bpar, nbkp };
						lmaxk[k][i][j] = feval(0, &params);

						struct feval_Parameters params2 = { mu[j], mu[i], alpha, tau[k], apar, bpar, nbkp };
						lmaxk[k][j][i] = feval(0, &params2);

						candidate_arc[count] = i * N + j;
						count++;

						akpos++;
					}
				}
			}
		}
		Ak[k].arcpos = (int*)calloc(count, sizeof(int));
		Ak[k].dim = 0;
		for (int ii = 0; ii < count; ii++) {
			Ak[k].arcpos[Ak[k].dim++] = candidate_arc[ii];
		}
		xvars += count;
	}
	free(candidate_arc);

	// Feasible (undirected) hub arcs
	FeasArc = create_int_vector(N * (N - 1) / 2);
	FeasArcPos = create_int_vector((N - 2) * N + (N - 1));
	svars = 0;
	FeasHubs.dim = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			count = 0;
			for (int k = 0; k < comm.dim; k++) {
				if (lmax[i] < lmaxk[k][i][j]) {
					lmax[i] = lmaxk[k][i][j];
				}
				if ((lmaxk[k][i][j] > 0)) {
					count++;
				}
			}
			if ((count > 0) && (j > i)) {
				FeasArc[svars] = i * N + j;
				FeasArcPos[i * N + j] = svars;
				svars++;
			}
		}
		if (lmax[i] > 0) {
			FeasHubs.dim++;
		}
	}



	// Feasible hubs
	FeasHubs.hub = (FHubs*)calloc(FeasHubs.dim, sizeof(FHubs));
	FeasHubs.pos = (FHubs*)calloc(N, sizeof(FHubs));
	arc1 = (crit1*)calloc(N, sizeof(crit1));
	count = 0;
	for (int i = 0; i < N; i++) {
		arc1[i].arc2 = (crit2*)calloc(N, sizeof(crit2));
		if ((lmax[i] > 0)) {
			FeasHubs.hub[count] = i;
			FeasHubs.pos[i] = count;
			count++;
		}
	}

	int a1, a2;
	nbFeas = create_int_vector(N);
	for (int k = 0; k < comm.dim; k++) {
		Ak[k].node = (xinfo*)calloc(N, sizeof(xinfo));
		for (int a = 0; a < Ak[k].dim; a++) {
			a1 = (int)floor(Ak[k].arcpos[a] / N);
			a2 = Ak[k].arcpos[a] - N * a1;
			nbFeas[a1]++;
			nbFeas[a2]++;

			Ak[k].node[a1].dim++;
			Ak[k].node[a2].dim++;
		}
	}

	akpos = 0;
	int* countN;
	xa1 = create_int_vector(xvars);
	xa2 = create_int_vector(xvars);
	xk = create_int_vector(xvars);
	for (int k = 0; k < comm.dim; k++) {
		for (int i = 0; i < N; i++) {
			Ak[k].node[i].pos = (int*)calloc(Ak[k].node[i].dim, sizeof(int));
			Ak[k].node[i].node2 = (int*)calloc(N, sizeof(int));
			Ak[k].node[i].first = (int*)calloc(N, sizeof(int));
		}
		countN = create_int_vector(N);
		for (int a = 0; a < Ak[k].dim; a++) {
			a1 = (int)floor(Ak[k].arcpos[a] / N);
			a2 = Ak[k].arcpos[a] - N * a1;

			Ak[k].node[a1].pos[countN[a1]] = akpos;
			Ak[k].node[a2].pos[countN[a2]] = akpos;

			Ak[k].node[min(a1, a2)].node2[max(a1, a2)] = akpos;
			Ak[k].node[min(a1, a2)].first[max(a1, a2)] = a1;

			countN[a1]++;
			countN[a2]++;

			xa1[akpos] = a1;
			xa2[akpos] = a2;
			xk[akpos] = k;
			akpos++;
		}
		free(countN);
	}

	// Ordering of commodities by criticality
	int* candidate_com;
	double* comm_crit, tauka;
	candidate_com = create_int_vector(comm.dim);
	comm_crit = create_double_vector(comm.dim);
	yvars = 0;
	int i, j;
	for (int a = 0; a < svars; a++) {
		i = (int)floor(FeasArc[a] / N);
		j = FeasArc[a] - N * i;
		count = 0;
		arc1[i].arc2[j].dimys = 0;
		for (int k = 0; k < comm.dim; k++) {
			if (lmaxk[k][i][j] > 0) {
				candidate_com[count] = k;
				comm_crit[count] = 1 - lmaxk[k][i][j] * lmaxk[k][j][i] / (mu[i] * mu[j]);
				count++;
			}
		}
		if (count > 0) {
			quicksort(comm_crit, candidate_com, 0, count - 1);
			arc1[i].arc2[j].com = (int*)calloc(count, sizeof(int));
			arc1[i].arc2[j].dim = count;
			for (int kk = 0; kk < count; kk++) {
				arc1[i].arc2[j].com[kk] = candidate_com[kk];
			}
		}
	}
	free(candidate_com);
	free(comm_crit);

	// Ordering of commodities by the homothetic projection
	yvars = 0;
	ypos = create_int_vector(xvars);
	double* delta, * taus;
	int* post;
	delta = create_double_vector(xvars);
	taus = create_double_vector(xvars);
	post = create_double_vector(xvars);
	for (int a = 0; a < svars; a++) {
		a1 = (int)floor(FeasArc[a] / N);
		a2 = FeasArc[a] - N * a1;
		arc1[a1].arc2[a2].dimys = 0;
		int k = 0, sel;
		double lmax1 = 0, lmax2 = 0, lsing1, lsing2, tauka = 0, tauaux = 0;
		while (k < comm.dim) {
			sel = 0;
			if ((comm.i[k] == comm.j[k + 1]) && (comm.j[k] == comm.i[k + 1]) && 
				(xa1[Ak[k].node[a1].node2[a2]] == xa2[Ak[k + 1].node[a1].node2[a2]]) && 
				(xa1[Ak[k + 1].node[a1].node2[a2]] == xa2[Ak[k].node[a1].node2[a2]]) && (k < comm.dim - 1)) {
				if (lmaxk[k][a1][a2] > 0) {
					sel++;
					tauaux = tau[k];
					lmax1 = lmaxk[k][a1][a2];
					lmax2 = lmaxk[k][a2][a1];
					apar = 1 / pow(0.5, 2);
					if (Ak[k].node[a1].first[a2] == a1) {
						bpar = (c2[comm.i[k]][a1] + eta * c2[a1][a2] + c2[a2][comm.j[k]]) / (apar * vel);
					}
					else {
						bpar = (c2[comm.i[k]][a2] + eta * c2[a1][a2] + c2[a1][comm.j[k]]) / (apar * vel);
					}
					akpos = Ak[k].node[a1].node2[a2];
					ypos[akpos] = yvars;
				}
				if (lmaxk[k + 1][a1][a2] > 0) {
					if (sel == 0) {
						sel++;
						tauaux = tau[k + 1];
						lmax1 = lmaxk[k + 1][a1][a2];
						lmax2 = lmaxk[k + 1][a2][a1];
						apar = 1 / pow(0.5, 2);
						if (Ak[k + 1].node[a1].first[a2] == a1) {
							bpar = (c2[comm.j[k]][a1] + eta * c2[a1][a2] + c2[a2][comm.i[k]]) / (apar * vel);
						}
						else {
							bpar = (c2[comm.j[k]][a2] + eta * c2[a1][a2] + c2[a1][comm.i[k]]) / (apar * vel);
						}
					}
					akpos = Ak[k + 1].node[a1].node2[a2];
					ypos[akpos] = yvars;
				}
				if ((lmaxk[k][a1][a2] > 0) || (lmaxk[k + 1][a1][a2] > 0)) {
					struct TauSelParams tauselpar = { mu[a1], mu[a2], alpha, tauaux, 1 / pow(0.5, 2), bpar, lmax1, lmax2 };
					tauka = TauSel(&tauselpar);
					lsing1 = (1 / tauka) * (1 + mu[a1] * tauka + gsl_sf_lambert_Wm1((alpha - 1) / exp(1)));
					lsing2 = (1 / tauka) * (1 + mu[a2] * tauka + gsl_sf_lambert_Wm1((alpha - 1) / exp(1)));

					delta[arc1[a1].arc2[a2].dimys] = mu[a1] - lsing1 + mu[a2] - lsing2;
					taus[yvars] = tauka;
					post[arc1[a1].arc2[a2].dimys] = yvars;

					arc1[a1].arc2[a2].dimys++;
					yvars++;
				}
				k = k + 2;
			}
			else {
				if (lmaxk[k][a1][a2] > 0) {
					sel++;
					tauaux = tau[k];
					lmax1 = lmaxk[k][a1][a2];
					lmax2 = lmaxk[k][a2][a1];
					apar = 1 / pow(0.5, 2);
					if (Ak[k].node[a1].first[a2] == a1) {
						bpar = (c2[comm.i[k]][a1] + eta * c2[a1][a2] + c2[a2][comm.j[k]]) / (apar * vel);
					}
					else {
						bpar = (c2[comm.i[k]][a2] + eta * c2[a1][a2] + c2[a1][comm.j[k]]) / (apar * vel);
					}
					struct TauSelParams tauselpar = { mu[a1], mu[a2], alpha, tauaux, 1 / pow(0.5, 2), bpar, lmax1, lmax2 };
					tauka = TauSel(&tauselpar);
					lsing1 = (1 / tauka) * (1 + mu[a1] * tauka + gsl_sf_lambert_Wm1((alpha - 1) / exp(1)));
					lsing2 = (1 / tauka) * (1 + mu[a2] * tauka + gsl_sf_lambert_Wm1((alpha - 1) / exp(1)));

					delta[arc1[a1].arc2[a2].dimys] = mu[a1] - lsing1 + mu[a2] - lsing2;
					taus[yvars] = tauka;
					post[arc1[a1].arc2[a2].dimys] = yvars;
					akpos = Ak[k].node[a1].node2[a2];
					ypos[akpos] = yvars;
					arc1[a1].arc2[a2].dimys++;
					yvars++;
				}
				k = k + 1;
			}
		}
		// Homothetic ordering of commodities
		quicksort(delta, post, 0, arc1[a1].arc2[a2].dimys - 1);
		arc1[a1].arc2[a2].pos = (int*)calloc(arc1[a1].arc2[a2].dimys, sizeof(int));
		arc1[a1].arc2[a2].tau = (double*)calloc(arc1[a1].arc2[a2].dimys, sizeof(double));
		arc1[a1].arc2[a2].ord = (int*)calloc(yvars, sizeof(int));
		for (int i = 0; i < arc1[a1].arc2[a2].dimys; i++) {
			arc1[a1].arc2[a2].pos[i] = post[i];
			arc1[a1].arc2[a2].tau[i] = taus[post[i]];
			arc1[a1].arc2[a2].ord[post[i]] = i + 1;
		}
	}
	free(delta);
	free(taus);
	free(post);

	int l, m;
	nCom = create_int_matrix(N, N);
	for (int i = 0; i < FeasHubs.dim; i++) {
		for (int j = 0; j < FeasHubs.dim; j++) {
			a1 = FeasHubs.hub[i];
			a2 = FeasHubs.hub[j];
			for (int kk = 0; kk < comm.dim; kk++) { 
				for (int aa = 0; aa < Ak[kk].dim; aa++) {
					l = (int)floor(Ak[kk].arcpos[aa] / N);
					m = Ak[kk].arcpos[aa] - N * l;
					if ((l == a1) || (m == a1) || (l == a2) || (m == a2)) {
						nCom[a1][a2]++;
					}
				}
			}
		}
	}
	fclose(in);
}


void initialize_memory(void) {
	c = create_double_matrix(N, N);   // Matrix for unit transport costs
	c2 = create_double_matrix(N, N);   // Matrix for unit transport costs
	c_c = create_double_matrix(N, N); // Matrix for the unit collection costs
	c_t = create_double_matrix(N, N); // Matrix for the unit transfer costs
	c_d = create_double_matrix(N, N); // Matrix for the unit distribution costs
	f = create_double_vector(N);    // Matrix for fixed costs
	W = create_double_matrix(N, N);    // Matrix for demands
	mu = create_double_vector(N);	   // Processing rate of outbound units
	lmax = create_double_vector(N);
	pts = (coordinate*)malloc((N) * sizeof(coordinate));  // Pointer for node coordinates
	comm.i = (int*)calloc(comm.dim, sizeof(int));   // Pointer for the origin node in a commodity
	comm.j = (int*)calloc(comm.dim, sizeof(int));   // Pointer for the destination node in a commodity
}

void free_memory(void) {

	for (int k = 0; k < comm.dim; k++) {
		for (int i = 0; i < N; i++) {
			free(lmaxk[k][i]);
			free(Ak[k].node[i].pos);
			free(Ak[k].node[i].node2);
			free(Ak[k].node[i].first);
		}
		free(lmaxk[k]);
		free(Ak[k].node);
	}
	free(lmaxk);

	for (int i = 0; i < N; i++) {
		free(W[i]);
		free(c[i]);
		free(c2[i]);
		free(c_c[i]);
		free(c_t[i]);
		free(c_d[i]);
	}
	free(W);
	free(f);
	free(mu);
	free(c);
	free(c2);
	free(c_c);
	free(c_t);
	free(c_d);

	free(lmax);
	free(tau);
	free(FeasHubs.hub);
	free(FeasHubs.pos);
	free(comm.i);
	free(comm.j);

	for (int i = 0; i < comm.dim; i++) {
		free(Ak[i].arcpos);
	}
	free(Ak);
	free(nbFeas);

	for (int i = 0; i < N - 1; i++) {
		for (int j = i + 1; j < N; j++) {
			free(arc1[i].arc2[j].com);
			free(arc1[i].arc2[j].pos);
			free(arc1[i].arc2[j].tau);
			free(arc1[i].arc2[j].ord);
		}
		free(arc1[i].arc2);
	}
	free(arc1);
	free(xa1);
	free(xa2);
	free(xk);
	free(ypos);
	free(FeasArc);
	free(FeasArcPos);
}

// Openning a file
FILE* open_file(const char* name, const char* mode) {
	FILE* file;
	if ((file = fopen(name, mode)) == NULL) {
		printf("\nError: Unable to read the file.\n");
		exit(8);
	}
}

// Create vectors and matrices
double** create_double_matrix(int rows, int cols) {
	double** ptr;
	if ((ptr = (double**)calloc(rows, sizeof(double*))) == NULL) {
		printf("\nError: Not enough memory.\n");
		exit(8);
	}
	for (int i = 0; i < rows; i++) {
		ptr[i] = create_double_vector(cols);
	}
	return ptr;
}

// Create vectors and matrices
double*** create_double_3Dmatrix(int idx, int rows, int cols) {
	double*** ptr;
	if ((ptr = (double***)calloc(idx, sizeof(double**))) == NULL) {
		printf("\nError: Not enough memory.\n");
		exit(8);
	}
	for (int i = 0; i < idx; i++) {
		ptr[i] = create_double_matrix(rows, cols);
	}
	return ptr;
}

double*** create_int_3Dmatrix(int idx, int rows, int cols) {
	double*** ptr;
	if ((ptr = (int***)calloc(idx, sizeof(int**))) == NULL) {
		printf("\nError: Not enough memory.\n");
		exit(8);
	}
	for (int i = 0; i < idx; i++) {
		ptr[i] = create_double_matrix(rows, cols);
	}
	return ptr;
}


double* create_double_vector(int dim) {
	double* ptr;
	if ((ptr = (double*)calloc(dim, sizeof(double))) == NULL) {
		printf("\nError: Not enough memory.\n");
		return ptr;
	}
}

int** create_int_matrix(int rows, int cols) {
	int** ptr;
	if ((ptr = (int**)calloc(rows, sizeof(int*))) == NULL) {
		printf("\nError: Not enough memory.\n");
		exit(8);
	}
	for (int i = 0; i < rows; i++) {
		ptr[i] = create_int_vector(cols);
	}
	return ptr;
}

int* create_int_vector(int dim) {
	int* ptr;
	if ((ptr = (int*)calloc(dim, sizeof(int))) == NULL) {
		printf("\nError: Not enough memory.\n");
		return ptr;
	}
}

void i_vector(int** vector, int n, char* s) {
	if ((*vector = (int*)calloc(n, sizeof(int))) == NULL)
		//error(s);
		printf("Error \n");
	return;
}

void d_vector(double** vector, int n, char* s) {
	if ((*vector = (double*)calloc(n, sizeof(double))) == NULL)
		// error(s);
		printf("Error \n");
	return;
}

void c_vector(char** vector, int n, char* s) {
	if ((*vector = (char*)calloc(n, sizeof(char))) == NULL)
		//error(s);
		printf("Error \n");
	return;
}