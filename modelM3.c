#include "header.h"
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_deriv.h>
#include "gauss_legendre.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

extern double** c2, ** c, * f, * mu, ** W;						// Transport costs, fixed costs, Capacity levels, Demand
extern double** c_c, ** c_t, ** c_d;
extern int N;
extern COMMOD comm;
extern coordinate* pts;					// Node coordinares
extern bitf* Ak;						// Hub arcs structure for each commodity
extern int* FeasArc, * FeasArcPos;
extern xinfo* node;
extern crit1* arc1;
extern crit2* arc2;
extern FHubs FeasHubs;
FILE* CutInfo,* CutInfo2;

extern double rq, eta, apar, bpar, vel;
extern double alpha, muhat, * tau, *** lmaxk, * lmax, cuttime;
extern char instancia[10];
extern int svars, * nbFeas, ** nCom;
extern int xvars, * xa1, * xa2, * xk;
extern int yvars, * ypos;
int*** cutCount;
extern char instancia[10];

int* xcuts;
extern int counter;
int nccuts;
int ncuts;

#define EPS 1e-6
void modelM3(void) {

	clock_t start, end;
	double cputime;

	//Variables required for cplex
	CPXLPptr  lp = NULL;// Estructura de datos para almacenar un problema en cplex ......
	CPXENVptr env = NULL;// Entorno cplex ...............................................
	CPXENVptr copy = NULL;// Entorno cplex ..............................................
	int numcols;		// Número de variables ..........................................
	int numrows;		// Número de restricciones sin contar las cotas  ................
	int numnz;			// Número de elementos no nulos de la matriz ....................
	double* obj;		// Coeficientes de la función objetivo ..........................
	double* rhs;		// Términos independientes de las restricciones .................
	char* sense;		// Sentido de las restricciones (<=: 'L', =:'E', >=:'G') ........
	int* matbeg;		// Índice del primer coeficiente no nulo de cada columna ........
	int* matind;		// Fila a la que corresponde cada elemento no nulo  .............
	double* matval;		// Valores de los coeficientes no nulos en las restricciones.....
	double* lb;			// Cotas inferiores de las variables ............................
	double* ub;			// Cotas superiores de las variables ............................
	char* ctype;		// Tipos de variable ('C', 'I', 'B') solo si hay enteras ........
	double* x;			// Vector para solucion (double, aunque el problema sea entero) .
	int status = 0;		// Para recuperar el estatus de la optimizacion .................
	char probname[16];	// Nombre del problema para cplex ...............................

	start = clock();

	// Initialize CPLEX environment and problem
	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}
	strcpy(probname, "HLPPSC");
	lp = CPXcreateprob(env, &status, probname);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}
	status = CPXchgobjsen(env, lp, CPX_MIN);


	// Create columns
	// ==============

	// Define z[h] variables 
	int colind = 0;
	numcols = FeasHubs.dim;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	for (int i = 0; i < FeasHubs.dim; i++) {
		obj[colind] = f[FeasHubs.hub[i]];
		lb[colind] = 0;
		ub[colind] = 1;
		ctype[colind] = 'B';
		colind++;
	}
	status = CPXnewcols(env, lp, colind, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);

	// Define x[k][a] variables
	colind = 0;
	numcols = xvars;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	int l, m;
	for (int k = 0; k < comm.dim; k++) {
		for (int a = 0; a < Ak[k].dim; a++) {
			l = (int)floor(Ak[k].arcpos[a] / N);
			m = Ak[k].arcpos[a] - N * l;
			obj[colind] = W[comm.i[k]][comm.j[k]] * (c_c[comm.i[k]][l] +
				c_t[l][m] +
				c_d[m][comm.j[k]]) - W[comm.i[k]][comm.j[k]] * c[comm.i[k]][comm.j[k]];
			lb[colind] = 0;
			ub[colind] = 1;
			ctype[colind] = 'C';
			colind++;
		}
	}
	status = CPXnewcols(env, lp, colind, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);

	// Define v[k][a] variables
	colind = 0;
	numcols = yvars;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	for (int i = 0; i < yvars; i++) {
		obj[colind] = 0;
		lb[colind] = 0;
		ub[colind] = 1;
		ctype[colind] = 'B';
		colind++;
	}
	status = CPXnewcols(env, lp, colind, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);

	// Define y[k][a] variables
	colind = 0;
	numcols = yvars;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	for (int i = 0; i < yvars; i++) {
		obj[colind] = 0;
		lb[colind] = 0;
		ub[colind] = 1;
		ctype[colind] = 'B';
		colind++;
	}
	status = CPXnewcols(env, lp, colind, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);



	// CONSTRAINTS
	// ===========

	// Add constraint sum{a in A[k]} x[a,k] <= 1, forall k in K
	int rowind;
	int nzind;
	int akpos = 0;
	for (int k = 0; k < comm.dim; k++) {
		rowind = 0;
		nzind = 0;
		numnz = Ak[k].dim;
		d_vector(&rhs, 1, "open_cplex:2");
		c_vector(&sense, 1, "open_cplex:3");
		i_vector(&matbeg, 1, "open_cplex:4");
		i_vector(&matind, numnz, "open_cplex:6");
		d_vector(&matval, numnz, "open_cplex:7");
		sense[0] = 'L';
		rhs[0] = 1;
		matbeg[0] = nzind;
		for (int a = 0; a < Ak[k].dim; a++) {
			matind[nzind] = akpos + FeasHubs.dim;
			matval[nzind] = 1;
			akpos++;
			nzind++;
		}
		status = CPXaddrows(env, lp, 0, 1, nzind, rhs, sense, matbeg, matind, matval, NULL, NULL);
		free(matbeg);
		free(matind);
		free(matval);
		free(sense);
		free(rhs);
	}

	// Add constraint -z[h] + sum{a in A[k] | i in a} x[a,k] <= 0, forall i in H, k in K
	rowind = 0;
	nzind = 0;
	numrows = comm.dim * FeasHubs.dim;
	numnz = FeasHubs.dim * comm.dim + 2 * xvars;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");
	for (int i = 0; i < FeasHubs.dim; i++) {
		akpos = 0;
		for (int k = 0; k < comm.dim; k++) {
			sense[rowind] = 'L';
			rhs[rowind] = 0;
			matbeg[rowind] = nzind;
			matind[nzind] = i;
			matval[nzind] = -1;
			nzind++;
			for (int a = 0; a < Ak[k].dim; a++) {
				l = (int)floor(Ak[k].arcpos[a] / N);
				m = Ak[k].arcpos[a] - N * l;
				if ((l == FeasHubs.hub[i]) || (m == FeasHubs.hub[i])) {
					matind[nzind] = akpos + FeasHubs.dim;
					matval[nzind] = 1;
					nzind++;
				}
				akpos++;
			}
			rowind++;
		}
	}
	status = CPXaddrows(env, lp, 0, rowind, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	// Add constraints sum{k in K}sum{a in Ak, h in a} w[k] * x[k][a] <= Lambda[h]*z[h]
	int a1, a2;
	rowind = 0;
	nzind = 0;
	numrows = FeasHubs.dim;
	numnz = xvars * 2 + FeasHubs.dim;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");
	for (int i = 0; i < FeasHubs.dim; i++) {
		a1 = FeasHubs.hub[i];
		sense[rowind] = 'L';
		rhs[rowind] = 0;
		matbeg[rowind] = nzind;
		matind[nzind] = i;
		matval[nzind] = -lmax[a1];
		nzind++;
		akpos = 0;
		for (int k = 0; k < comm.dim; k++) {
			for (int a = 0; a < Ak[k].dim; a++) {
				l = (int)floor(Ak[k].arcpos[a] / N);
				m = Ak[k].arcpos[a] - N * l;
				if (l == a1) {
					matind[nzind] = akpos + FeasHubs.dim;
					matval[nzind] = W[comm.i[k]][comm.j[k]];
					nzind++;
				}
				if (m == a1) {
					matind[nzind] = akpos + FeasHubs.dim;
					matval[nzind] = W[comm.i[k]][comm.j[k]];
					nzind++;
				}
				akpos++;
			}
		}
		rowind++;
	}
	status = CPXaddrows(env, lp, 0, rowind, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	// Add constraints x[k][a] <= y[k][a]
	rowind = 0;
	nzind = 0;
	numrows = xvars;
	numnz = xvars * 2;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");
	akpos = 0;
	for (int k = 0; k < comm.dim; k++) {
		for (int a = 0; a < Ak[k].dim; a++) {
			sense[rowind] = 'L';
			rhs[rowind] = 0;
			matbeg[rowind] = nzind;
			matind[nzind] = akpos + FeasHubs.dim;
			matval[nzind] = 1;
			nzind++;
			matind[nzind] = ypos[akpos] + FeasHubs.dim + xvars;
			matval[nzind] = -1;
			nzind++;
			rowind++;
			akpos++;
		}
	}
	status = CPXaddrows(env, lp, 0, rowind, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	// Incremental constraints v[k][a] <= v[k-1][a]
	rowind = 0;
	nzind = 0;
	numrows = yvars;
	numnz = yvars * 2 - svars;
	int akpos2;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");
	akpos = 0;
	for (int a = 0; a < svars; a++) {
		a1 = (int)floor(FeasArc[a] / N);
		a2 = FeasArc[a] - N * a1;
		for (int k = 0; k < arc1[a1].arc2[a2].dimys; k++) {
			sense[rowind] = 'L';
			rhs[rowind] = 0;
			matbeg[rowind] = nzind;

			if (k == 0) {
				matind[nzind] = arc1[a1].arc2[a2].pos[k] + FeasHubs.dim + xvars;
				matval[nzind] = -1;
				nzind++;
			}
			else {
				matind[nzind] = arc1[a1].arc2[a2].pos[k] + FeasHubs.dim + xvars;
				matval[nzind] = 1;
				nzind++;

				matind[nzind] = arc1[a1].arc2[a2].pos[k - 1] + FeasHubs.dim + xvars;
				matval[nzind] = -1;
				nzind++;
			}
			rowind++;
		}
	}
	status = CPXaddrows(env, lp, 0, rowind, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	// Add constraints v[0][a] <= z[m]
	rowind = 0;
	nzind = 0;
	numrows = 2 * svars;
	numnz = 4 * svars;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");
	akpos = 0;
	for (int a = 0; a < svars; a++) {
		a1 = (int)floor(FeasArc[a] / N);
		a2 = FeasArc[a] - N * a1;
		int dimy;
		dimy = arc1[a1].arc2[a2].dimys;

		sense[rowind] = 'L';
		rhs[rowind] = 0;

		matbeg[rowind] = nzind;
		matind[nzind] = FeasHubs.pos[a1];
		matval[nzind] = -1;
		nzind++;

		matind[nzind] = arc1[a1].arc2[a2].pos[dimy - 1] + FeasHubs.dim + xvars;
		matval[nzind] = 1;
		nzind++;
		rowind++;

		sense[rowind] = 'L';
		rhs[rowind] = 0;

		matbeg[rowind] = nzind;
		matind[nzind] = FeasHubs.pos[a2];
		matval[nzind] = -1;
		nzind++;

		matind[nzind] = arc1[a1].arc2[a2].pos[dimy - 1] + FeasHubs.dim + xvars;
		matval[nzind] = 1;
		nzind++;
		rowind++;
	}
	status = CPXaddrows(env, lp, 0, rowind, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	// Tuning CPLEX parameters for branch and bound algorithm
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);		// Output display 
	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);		// different levels of output display
	CPXsetdblparam(env, CPX_PARAM_TILIM, 86400);		// Execution time limit
	CPXsetintparam(env, CPX_PARAM_THREADS, 1);			// Number of threads to use
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.000001);		// e-optimal solution (% gap)
	CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0000001);	// Integer precission

	// Branching priority
	int cnt = FeasHubs.dim + yvars;
	int* indices;
	int* priority;
	indices = create_int_vector(cnt);
	priority = create_int_vector(cnt);
	for (int i = 0; i < FeasHubs.dim; i++) {
		indices[i] = i;
		priority[i] = 10;
	}
	akpos = 0;
	for (int i = 0; i < yvars; i++) {
		indices[i + FeasHubs.dim] = i + FeasHubs.dim + xvars;
		priority[i + FeasHubs.dim] = 1;
	}
	status = CPXcopyorder(env, lp, cnt, indices, priority, NULL);
	free(indices);
	free(priority);

	//CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF);
	CPXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_DUALONLY);


	CUTINFO cutinfo;
	cutinfo.x = NULL;
	int cur_numcols;
	cur_numcols = CPXgetnumcols(env, lp);
	cutinfo.lp = lp;
	cutinfo.numcols = cur_numcols;
	cutinfo.x = (double*)malloc(cur_numcols * sizeof(double));

	// Set up to use MIP callback 
	status = CPXsetlazyconstraintcallbackfunc(env, mycutcallback, &cutinfo);


	// Creation of a file to store detailed information of the fine perspective cuts added during the optimization process
	char charvar[100];
	sprintf(charvar, "Convexity_check_CPC-r%.0f-alpha%2.0f", rq, 100 * alpha);
	CutInfo = open_file(strcat(strcat(charvar, "-"), instancia), "a+");
	fprintf(CutInfo, "Parameter r: %.1lf = \n\n", rq);
	fprintf(CutInfo, "a1\ta2\tRay intercept with vertical axis\tDerivative of tilde{f}^[ka]\n");

	char charvar2[100];
	sprintf(charvar2, "Convexity_check_FPC-r%.0f-alpha%2.0f", rq, 100 * alpha);
	CutInfo2 = open_file(strcat(strcat(charvar2, "-"), instancia), "a+");
	fprintf(CutInfo2, "Parameter r: %.1lf = \n\n", rq);
	fprintf(CutInfo2, "k\ta1\ta2\tlambda_[a2]\tDerivative of f^[ka]\tmu_[a1]\tmu_[a2]\tAlpha\tTau_[k]\n");

	double sum = 0;
	for (int k = 0; k < comm.dim; k++) {
		sum += W[comm.i[k]][comm.j[k]] * c[comm.i[k]][comm.j[k]];
	}

	ncuts = 0;
	cutCount = create_int_3Dmatrix(comm.dim, N, N);
	xcuts = create_int_vector(xvars);
	numcols = CPXgetnumcols(env, lp);
	d_vector(&x, numcols, "open_cplex:0");
	if (FeasHubs.dim > 0) {
		CPXmipopt(env, lp);

		printf("Status: ");
		int i = CPXgetstat(env, lp);
		if (i == 101)
			printf("\tOptimal solution\n\n");
		else if (i == 102)
			printf("\te-optimal\n\n");
		else if (i == 107)
			printf("\tTime\n\n");
		else if (i == 110)
			printf("\tNo integer solution (%d)\n\n", i);
		else if (i == 103)
			printf("\tInteger infeasible (%d)\n\n", i);
		else
			printf("\tNot known (%d)\n\n", i);

		// Print solution
		// ==============
		double BestUB, BestLB;
		CPXgetmipobjval(env, lp, &BestUB);
		CPXgetbestobjval(env, lp, &BestLB);

		// Bounds and branch and bound nodes
		printf("Upper bound: %f\n", BestUB + sum);
		printf("Lower bound: %f\n\n", BestLB + sum);
		int nodecount;
		nodecount = CPXgetnodecnt(env, lp);
		printf("Number of BB nodes : %ld  \n\n", nodecount);

		// Hub locations
		double InvCost = 0;
		int CountHubs = 0;
		CPXgetmipx(env, lp, x, 0, numcols - 1);
		printf("Optimal Hubs: ");
		for (int i = 0; i < FeasHubs.dim; i++) {
			if (x[i] > 0.5) {
				InvCost += f[FeasHubs.hub[i]];
				printf("%d ", FeasHubs.hub[i] + 1);
				CountHubs++;
			}
		}
		printf("\n\n");

		end = clock();
		cputime = (double)(end - start) / CLOCKS_PER_SEC;
		printf("Time: %.2f\n\n", cputime);

		int countarcs2 = 0, count = 0;
		double totalcom = 0, commviol = 0;
		double la1, la2, seval, minslev = 1, maxslev = 0, averslev = 0;
		for (int a = 0; a < svars; a++) {
			a1 = (int)floor(FeasArc[a] / N);
			a2 = FeasArc[a] - N * a1;
			la1 = 0;
			la2 = 0;
			for (int k = 0; k < comm.dim; k++) {
				for (int aa = 0; aa < Ak[k].node[a1].dim; aa++) {
					la1 += W[comm.i[k]][comm.j[k]] * x[Ak[k].node[a1].pos[aa] + FeasHubs.dim];
				}
				for (int aa = 0; aa < Ak[k].node[a2].dim; aa++) {
					la2 += W[comm.i[k]][comm.j[k]] * x[Ak[k].node[a2].pos[aa] + FeasHubs.dim];
				}
			}

			int com, a1x, a2x, sel = 0;
			for (int k = 0; k < arc1[a1].arc2[a2].dim; k++) {
				com = arc1[a1].arc2[a2].com[k];
				akpos = Ak[arc1[a1].arc2[a2].com[k]].node[a1].node2[a2];
				a1x = xa1[akpos];
				a2x = xa2[akpos];
				if ((x[akpos + FeasHubs.dim] > 0.0)) {
					apar = 1 / pow(0.5, 2);
					bpar = (c2[comm.i[com]][a1x] + eta * c2[a1x][a2x] + c2[a2x][comm.j[com]]) / (apar * vel);
					if ((a1x == a1)) {
						struct Convolve_Parameters paramsC = { la1, la2, mu[a1], mu[a2], tau[com], apar, bpar };
						seval = Convolve(40, Integrand, &paramsC);
					}
					else {
						struct Convolve_Parameters paramsC = { la2, la1, mu[a2], mu[a1], tau[com], apar, bpar };
						seval = Convolve(40, Integrand, &paramsC);
					}
					sel++;
					count++;
					averslev = ((count - 1) * averslev + seval) / count;
					if (minslev > seval)
						minslev = seval;
					if (maxslev < seval)
						maxslev = seval;

					if (seval < alpha - pow(10, -6)) {
						commviol = commviol + 1;
					}
					totalcom = totalcom + 1;
				}
			}
		}

		printf("Minimum network service level: %.6lf \n", minslev);
		printf("Average network service level: %.6lf \n", averslev);
		printf("Maximum network service level: %.6lf \n\n", maxslev);

		printf("--------------------------------------------------------\n\n");
	}
	else {
		printf("Instance: %s\nSolution status: Noe feasible hubs\nTotal cost: %lf\n", instancia, sum);
	}

	// Close the problem
	if (lp != NULL) {
		status = CPXfreeprob(env, &lp);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env != NULL) {
		status = CPXcloseCPLEX(&env);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}
	free(xcuts);
	free(x);
	for (int i = 0; i < comm.dim; i++) {
		for (int j = 0; j < N; j++) {
			free(cutCount[i][j]);
		}
		free(cutCount[i]);
	}
	free(cutCount);
	free_and_null((char**)&cutinfo.x);
	fclose(CutInfo);
	fclose(CutInfo2);
}


// ========================
// Set lazy constraints
// ========================
static int CPXPUBLIC
mycutcallback(CPXCENVptr env,
	void* cbdata,
	int wherefrom,
	void* cbhandle,
	int* useraction_p) {

	clock_t start, end;
	start = clock();

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;

	int      numcols = cutinfo->numcols;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* matindc;
	double* matvalc;

	int status = 0;
	int nzind;
	double objval;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	status = CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, numcols - 1);
	if (status != 0)
		return status;

	int ccounter = 0;
	double epsilon = pow(10, -1);

	int a1, a2, count = 0, akpos, akpos2 = 0, l, m, k;
	double la1, la2, l1, l2, fa1, fa2, fka1, lmax1, lmax2, L1, L2, delta1, delta2, seval;
	for (int a = 0; a < svars; a++) {
		a1 = (int)floor(FeasArc[a] / N);
		a2 = FeasArc[a] - N * a1;
		if ((x[FeasHubs.pos[a1]] > 0.5) && (x[FeasHubs.pos[a2]] > 0.5)) {
			count = 0;
			la1 = 0;
			la2 = 0;
			for (int k = 0; k < comm.dim; k++) {
				for (int aa = 0; aa < Ak[k].node[a1].dim; aa++) {
					la1 += W[comm.i[k]][comm.j[k]] * x[Ak[k].node[a1].pos[aa] + FeasHubs.dim];
				}
				for (int aa = 0; aa < Ak[k].node[a2].dim; aa++) {
					la2 += W[comm.i[k]][comm.j[k]] * x[Ak[k].node[a2].pos[aa] + FeasHubs.dim];
				}
			}

			struct Wait_Parmeters plim1 = { mu[a1], mu[a2], alpha, arc1[a1].arc2[a2].tau[0] };
			struct Wait_Parmeters plim2 = { mu[a2], mu[a1], alpha, arc1[a1].arc2[a2].tau[0] };
			L1 = fmax(lmax[a1], Wait(0, &plim1));
			L2 = fmax(lmax[a2], Wait(0, &plim2));

			for (int k = arc1[a1].arc2[a2].dimys - 1; k >= 0; k--) {

				if (x[arc1[a1].arc2[a2].pos[k] + FeasHubs.dim + xvars] > 0.5) {
					struct Wait_Parmeters paramsapprox1 = { mu[a1], mu[a2], alpha, arc1[a1].arc2[a2].tau[k] };
					struct Wait_Parmeters paramsapprox2 = { mu[a2], mu[a1], alpha, arc1[a1].arc2[a2].tau[k] };

					lmax1 = Wait(0, &paramsapprox1);
					lmax2 = Wait(0, &paramsapprox2);
					if ((la1 > lmax1) && (Service(la1, la2, &paramsapprox1) < alpha - epsilon)) {
						l1 = lmax1;
						l2 = 0;
						fka1 = WaitDeriv(l2, &paramsapprox1);

						fprintf(CutInfo, "%d \t %d \t %.12lf \t %.12lf \n", a1 + 1, a2 + 1, mu[a1] - mu[a2] * (mu[a1] - l1) / (mu[a2] - l2), fka1);

						nzind = 0;
						matindc = (int*)malloc((nCom[a1][a2] + arc1[a1].arc2[a2].dimys + 2) * sizeof(int));
						matvalc = (double*)malloc((nCom[a1][a2] + arc1[a1].arc2[a2].dimys + 2) * sizeof(double));

						matindc[nzind] = FeasHubs.pos[a1];
						matvalc[nzind] = -L1;
						nzind++;

						matindc[nzind] = FeasHubs.pos[a2];
						matvalc[nzind] = -L2 * fabs(fka1);
						nzind++;

						akpos = 0;
						for (int kk = 0; kk < comm.dim; kk++) {
							for (int aa = 0; aa < Ak[kk].dim; aa++) {
								l = (int)floor(Ak[kk].arcpos[aa] / N);
								m = Ak[kk].arcpos[aa] - N * l;
								if (((l == a1) || (m == a1)) && ((l == a2) || (m == a2))) {
									matindc[nzind] = akpos + FeasHubs.dim;
									matvalc[nzind] = (1 + fabs(fka1)) * W[comm.i[kk]][comm.j[kk]];
									nzind++;
								}
								if (((l == a1) || (m == a1)) && ((l != a2) && (m != a2))) {
									matindc[nzind] = akpos + FeasHubs.dim;
									matvalc[nzind] = W[comm.i[kk]][comm.j[kk]];
									nzind++;
								}
								if (((l != a1) && (m != a1)) && ((l == a2) || (m == a2))) {
									matindc[nzind] = akpos + FeasHubs.dim;
									matvalc[nzind] = W[comm.i[kk]][comm.j[kk]] * fabs(fka1);
									nzind++;
								}
								akpos++;
							}
						}
						for (int kk = 0; kk < arc1[a1].arc2[a2].dimys; kk++) {
							if (kk == 0) {
								delta1 = L1 - (mu[a1] + (l1 - mu[a1]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]) +
									fabs(fka1) * (L2 - (mu[a2] + (l2 - mu[a2]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]));
								matindc[nzind] = arc1[a1].arc2[a2].pos[kk] + FeasHubs.dim + xvars;
								matvalc[nzind] = delta1;
								nzind++;
							}
							else {
								delta1 = L1 - (mu[a1] + (l1 - mu[a1]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk - 1]) +
									fabs(fka1) * (L2 - (mu[a2] + (l2 - mu[a2]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk - 1]));
								delta2 = L1 - (mu[a1] + (l1 - mu[a1]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]) +
									fabs(fka1) * (L2 - (mu[a2] + (l2 - mu[a2]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]));
								if (delta2 < delta1) {
									printf("Error on delta");
								}
								matindc[nzind] = arc1[a1].arc2[a2].pos[kk] + FeasHubs.dim + xvars;
								matvalc[nzind] = delta2 - delta1;
								nzind++;
							}
						}
						status = CPXcutcallbackadd(env, cbdata, wherefrom, nCom[a1][a2] + arc1[a1].arc2[a2].dimys + 2, 0, 'L', matindc, matvalc, CPX_USECUT_FORCE);
						free(matindc);
						free(matvalc);

						nccuts++;
						ccounter++;
						count++;
					}
					if ((la2 > lmax2) && (Service(la1, la2, &paramsapprox1) < alpha - epsilon)) {
						l1 = 0;
						l2 = lmax2;
						fka1 = WaitDeriv(l2, &paramsapprox1);

						fprintf(CutInfo, "%d \t %d \t %.12lf \t %.12lf \n", a1 + 1, a2 + 1, mu[a1] - mu[a2] * (mu[a1] - l1) / (mu[a2] - l2), fka1);

						nzind = 0;
						matindc = (int*)malloc((nCom[a1][a2] + arc1[a1].arc2[a2].dimys + 2) * sizeof(int));
						matvalc = (double*)malloc((nCom[a1][a2] + arc1[a1].arc2[a2].dimys + 2) * sizeof(double));

						matindc[nzind] = FeasHubs.pos[a1];
						matvalc[nzind] = -L1;
						nzind++;

						matindc[nzind] = FeasHubs.pos[a2];
						matvalc[nzind] = -L2 * fabs(fka1);
						nzind++;

						akpos = 0;
						for (int kk = 0; kk < comm.dim; kk++) {
							for (int aa = 0; aa < Ak[kk].dim; aa++) {
								l = (int)floor(Ak[kk].arcpos[aa] / N);
								m = Ak[kk].arcpos[aa] - N * l;
								if (((l == a1) || (m == a1)) && ((l == a2) || (m == a2))) {
									matindc[nzind] = akpos + FeasHubs.dim;
									matvalc[nzind] = (1 + fabs(fka1)) * W[comm.i[kk]][comm.j[kk]];
									nzind++;
								}
								if (((l == a1) || (m == a1)) && ((l != a2) && (m != a2))) {
									matindc[nzind] = akpos + FeasHubs.dim;
									matvalc[nzind] = W[comm.i[kk]][comm.j[kk]];
									nzind++;
								}
								if (((l != a1) && (m != a1)) && ((l == a2) || (m == a2))) {
									matindc[nzind] = akpos + FeasHubs.dim;
									matvalc[nzind] = W[comm.i[kk]][comm.j[kk]] * fabs(fka1);
									nzind++;
								}
								akpos++;
							}
						}
						for (int kk = 0; kk < arc1[a1].arc2[a2].dimys; kk++) {
							if (kk == 0) {
								delta1 = L1 - (mu[a1] + (l1 - mu[a1]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]) +
									fabs(fka1) * (lmax[a2] - (mu[a2] + (l2 - mu[a2]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]));
								matindc[nzind] = arc1[a1].arc2[a2].pos[kk] + FeasHubs.dim + xvars;
								matvalc[nzind] = delta1;
								nzind++;
							}
							else {
								delta1 = L1 - (mu[a1] + (l1 - mu[a1]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk - 1]) +
									fabs(fka1) * (L2 - (mu[a2] + (l2 - mu[a2]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk - 1]));
								delta2 = L1 - (mu[a1] + (l1 - mu[a1]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]) +
									fabs(fka1) * (L2 - (mu[a2] + (l2 - mu[a2]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]));
								if (delta2 < delta1) {
									printf("Error on delta");
								}
								matindc[nzind] = arc1[a1].arc2[a2].pos[kk] + FeasHubs.dim + xvars;
								matvalc[nzind] = delta2 - delta1;
								nzind++;
							}
						}
						status = CPXcutcallbackadd(env, cbdata, wherefrom, nCom[a1][a2] + arc1[a1].arc2[a2].dimys + 2, 0, 'L', matindc, matvalc, CPX_USECUT_FORCE);
						free(matindc);
						free(matvalc);

						nccuts++;
						ccounter++;
						count++;
					}

					if ((Service(la1, la2, &paramsapprox1) < alpha - epsilon) && (la1 <= lmax1 + pow(10, -10)) && (la2 <= lmax2 + pow(10, -10))) {
						fa1 = Wait(la2, &paramsapprox1);
						fa2 = Wait(la1, &paramsapprox2);
						l1 = fa1;
						l2 = la2;
						fka1 = WaitDeriv(l2, &paramsapprox1);

						fprintf(CutInfo, "%d \t %d \t %.12lf \t %.12lf \n", a1 + 1, a2 + 1, mu[a1] - mu[a2] * (mu[a1] - l1) / (mu[a2] - l2), fka1);

						nzind = 0;
						matindc = (int*)malloc((nCom[a1][a2] + arc1[a1].arc2[a2].dimys + 2) * sizeof(int));
						matvalc = (double*)malloc((nCom[a1][a2] + arc1[a1].arc2[a2].dimys + 2) * sizeof(double));

						matindc[nzind] = FeasHubs.pos[a1];
						matvalc[nzind] = -L1;
						nzind++;

						matindc[nzind] = FeasHubs.pos[a2];
						matvalc[nzind] = -L2 * fabs(fka1);
						nzind++;

						akpos = 0;
						for (int kk = 0; kk < comm.dim; kk++) {
							for (int aa = 0; aa < Ak[kk].dim; aa++) {
								l = (int)floor(Ak[kk].arcpos[aa] / N);
								m = Ak[kk].arcpos[aa] - N * l;
								if (((l == a1) || (m == a1)) && ((l == a2) || (m == a2))) {
									matindc[nzind] = akpos + FeasHubs.dim;
									matvalc[nzind] = (1 + fabs(fka1)) * W[comm.i[kk]][comm.j[kk]];
									nzind++;
								}
								if (((l == a1) || (m == a1)) && ((l != a2) && (m != a2))) {
									matindc[nzind] = akpos + FeasHubs.dim;
									matvalc[nzind] = W[comm.i[kk]][comm.j[kk]];
									nzind++;
								}
								if (((l != a1) && (m != a1)) && ((l == a2) || (m == a2))) {
									matindc[nzind] = akpos + FeasHubs.dim;
									matvalc[nzind] = W[comm.i[kk]][comm.j[kk]] * fabs(fka1);
									nzind++;
								}
								akpos++;
							}
						}
						for (int kk = 0; kk < arc1[a1].arc2[a2].dimys; kk++) {
							if (kk == 0) {
								delta1 = L1 - (mu[a1] + (l1 - mu[a1]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]) +
									fabs(fka1) * (lmax[a2] - (mu[a2] + (l2 - mu[a2]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]));
								matindc[nzind] = arc1[a1].arc2[a2].pos[kk] + FeasHubs.dim + xvars;
								matvalc[nzind] = delta1;
								nzind++;
							}
							else {
								delta1 = L1 - (mu[a1] + (l1 - mu[a1]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk - 1]) +
									fabs(fka1) * (L2 - (mu[a2] + (l2 - mu[a2]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk - 1]));
								delta2 = L1 - (mu[a1] + (l1 - mu[a1]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]) +
									fabs(fka1) * (L2 - (mu[a2] + (l2 - mu[a2]) * arc1[a1].arc2[a2].tau[k] / arc1[a1].arc2[a2].tau[kk]));
								if ((delta2 < delta1) || (delta1 < 0) || (delta2 < 0)) {
									printf("Error on delta");
								}
								matindc[nzind] = arc1[a1].arc2[a2].pos[kk] + FeasHubs.dim + xvars;
								matvalc[nzind] = delta2 - delta1;
								nzind++;
							}
						}
						status = CPXcutcallbackadd(env, cbdata, wherefrom, nCom[a1][a2] + arc1[a1].arc2[a2].dimys + 2, 0, 'L', matindc, matvalc, CPX_USECUT_FORCE);
						free(matindc);
						free(matvalc);

						nccuts++;
						ccounter++;
						count++;
					}
				}
				if (count > 0) {
					break;
				}
			}
		}
	}

	// Introduce fine perspective cuts
	double fka;
	if (ccounter == 0) {
		for (int i = 0; i < FeasHubs.dim - 1; i++) {
			for (int j = i + 1; j < FeasHubs.dim; j++) {
				if ((x[i] > 0.5) && (x[j] > 0.5)) {
					count = 0;
					for (int com = arc1[FeasHubs.hub[i]].arc2[FeasHubs.hub[j]].dim - 1; com >= 0; com--) {
						k = arc1[FeasHubs.hub[i]].arc2[FeasHubs.hub[j]].com[com];
						akpos = Ak[k].node[FeasHubs.hub[i]].node2[FeasHubs.hub[j]];
						a1 = Ak[k].node[FeasHubs.hub[i]].first[FeasHubs.hub[j]];
						if (a1 == FeasHubs.hub[i]) {
							a2 = FeasHubs.hub[j];
						}
						else {
							a2 = FeasHubs.hub[i];
						}
						la1 = 0;
						la2 = 0;
						for (int kk = 0; kk < comm.dim; kk++) {
							for (int aa = 0; aa < Ak[kk].node[a1].dim; aa++) {
								la1 += W[comm.i[kk]][comm.j[kk]] * x[Ak[kk].node[a1].pos[aa] + FeasHubs.dim];
							}
							for (int aa = 0; aa < Ak[kk].node[a2].dim; aa++) {
								la2 += W[comm.i[kk]][comm.j[kk]] * x[Ak[kk].node[a2].pos[aa] + FeasHubs.dim];
							}
						}

						if ((x[akpos + FeasHubs.dim] > pow(10, -6))) {
							apar = 1 / pow(0.5, 2);
							bpar = (c2[comm.i[k]][a1] + eta * c2[a1][a2] + c2[a2][comm.j[k]]) / (apar * vel);
							struct Convolve_Parameters paramsC = { la1, la2, mu[a1], mu[a2], tau[k], apar, bpar };
							seval = Convolve(40, Integrand, &paramsC);

							if ((cutCount[k][a1][a2] == 0) && (seval < alpha - pow(10, -6))) {
								matindc = (int*)malloc((nbFeas[a1] + 2) * sizeof(int));
								matvalc = (double*)malloc((nbFeas[a1] + 2) * sizeof(double));
								nzind = 0;
								matindc[nzind] = FeasHubs.pos[a1];
								matvalc[nzind] = -lmax[a1];
								nzind++;
								for (int kk = 0; kk < comm.dim; kk++) {
									for (int aa = 0; aa < Ak[kk].node[a1].dim; aa++) {
										matindc[nzind] = Ak[kk].node[a1].pos[aa] + FeasHubs.dim;
										matvalc[nzind] = W[comm.i[kk]][comm.j[kk]];
										nzind++;
									}
								}
								matindc[nzind] = ypos[akpos] + FeasHubs.dim + xvars + yvars;
								matvalc[nzind] = lmax[a1] - lmaxk[k][a1][a2];
								nzind++;
								status = CPXcutcallbackadd(env, cbdata, wherefrom, (nbFeas[a1] + 2), 0, 'L', matindc, matvalc, CPX_USECUT_FORCE);
								free(matindc);
								free(matvalc);

								ncuts++;
								cutCount[k][a1][a2]++;
								count++;
							}
							if ((cutCount[k][a2][a1] == 0) && (seval < alpha - pow(10, -6))) {
								// Upper bound
								matindc = (int*)malloc((nbFeas[a2] + 2) * sizeof(int));
								matvalc = (double*)malloc((nbFeas[a2] + 2) * sizeof(double));
								nzind = 0;
								matindc[nzind] = FeasHubs.pos[a2];
								matvalc[nzind] = -lmax[a2];
								nzind++;
								for (int kk = 0; kk < comm.dim; kk++) {
									for (int aa = 0; aa < Ak[kk].node[a2].dim; aa++) {
										matindc[nzind] = Ak[kk].node[a2].pos[aa] + FeasHubs.dim;
										matvalc[nzind] = W[comm.i[kk]][comm.j[kk]];
										nzind++;
									}
								}
								matindc[nzind] = ypos[akpos] + FeasHubs.dim + xvars + yvars;
								matvalc[nzind] = lmax[a2] - lmaxk[k][a2][a1];
								nzind++;
								status = CPXcutcallbackadd(env, cbdata, wherefrom, (nbFeas[a2] + 2), 0, 'L', matindc, matvalc, CPX_USECUT_FORCE);
								free(matindc);
								free(matvalc);

								ncuts++;
								cutCount[k][a2][a1]++;
								count++;
							}
							if ((seval < alpha - pow(10, -6)) && (la1 <= lmaxk[k][a1][a2] + pow(10, -10)) && (la2 <= lmaxk[k][a2][a1] + pow(10, -10))) {

								struct feval_Parameters params = { mu[a1], mu[a2], alpha, tau[k], apar, bpar, 40 };
								fka = feval(la2, &params);
								fka1 = fderiv(la2, &params);
								// Print cut information
								//printf("%d \t %d \t %d \t %d \t %.12lf \t %lf - %lf \t %lf - %lf \t %lf \n", k, a1, a2, akpos, seval, la1, lmaxk[k][a1][a2], la2, lmaxk[k][a2][a1], x[akpos + FeasHubs.dim + xvars]);
								fprintf(CutInfo2, "%d \t %d \t %d \t %.12lf \t %.12lf \t %.12lf \t %.12lf \t %.2lf \t %.12e \n", k + 1, a1 + 1, a2 + 1, la2, fka1, mu[a1], mu[a2], alpha, tau[k]);

								nzind = 0;
								matindc = (int*)malloc((nCom[a1][a2] + 3) * sizeof(int));
								matvalc = (double*)malloc((nCom[a1][a2] + 3) * sizeof(double));

								matindc[nzind] = FeasHubs.pos[a1];
								matvalc[nzind] = -lmax[a1];
								nzind++;

								matindc[nzind] = FeasHubs.pos[a2];
								matvalc[nzind] = -lmax[a2] * ((-1) * fka1);
								nzind++;

								int akpos2 = 0, l, m;
								for (int kk = 0; kk < comm.dim; kk++) {
									for (int aa = 0; aa < Ak[kk].dim; aa++) {
										l = (int)floor(Ak[kk].arcpos[aa] / N);
										m = Ak[kk].arcpos[aa] - N * l;
										if (akpos == akpos2) {
											matindc[nzind] = akpos + FeasHubs.dim;
											matvalc[nzind] = W[comm.i[k]][comm.j[k]] * (1 + ((-1) * fka1));
											nzind++;
										}
										if (((l == a1) || (m == a1)) && ((l == a2) || (m == a2)) && (kk != k)) {
											matindc[nzind] = akpos2 + FeasHubs.dim;
											matvalc[nzind] = (1 + fka1 * (-1)) * W[comm.i[kk]][comm.j[kk]];
											nzind++;
										}
										if (((l == a1) || (m == a1)) && ((l != a2) && (m != a2))) {
											matindc[nzind] = akpos2 + FeasHubs.dim;
											matvalc[nzind] = W[comm.i[kk]][comm.j[kk]];
											nzind++;
										}
										if (((l != a1) && (m != a1)) && ((l == a2) || (m == a2))) {
											matindc[nzind] = akpos2 + FeasHubs.dim;
											matvalc[nzind] = W[comm.i[kk]][comm.j[kk]] * fka1 * (-1);
											nzind++;
										}
										akpos2++;
									}
								}
								matindc[nzind] = ypos[akpos] + FeasHubs.dim + xvars + yvars;
								matvalc[nzind] = lmax[a1] - fka + (lmax[a2] - la2) * ((-1) * fka1);
								nzind++;
								status = CPXcutcallbackadd(env, cbdata, wherefrom, nCom[a1][a2] + 3, 0, 'L', matindc, matvalc, CPX_USECUT_FORCE);
								free(matindc);
								free(matvalc);

								ncuts++;
								count++;
							}
						}
						if (count > 0) {
							if (xcuts[akpos] == 0) {
								matindc = (int*)malloc(2 * sizeof(int));
								matvalc = (double*)malloc(2 * sizeof(double));
								nzind = 0;
								matindc[nzind] = akpos + FeasHubs.dim;
								matvalc[nzind] = 1;
								nzind++;
								matindc[nzind] = ypos[akpos] + FeasHubs.dim + xvars + yvars;
								matvalc[nzind] = -1;
								nzind++;
								status = CPXcutcallbackadd(env, cbdata, wherefrom, 2, 0, 'L', matindc, matvalc, CPX_USECUT_FORCE);
								free(matindc);
								free(matvalc);
								xcuts[akpos]++;
							}
							break;
						}
					}
				}
			}
		}
	}
	end = clock();
	cuttime += (double)(end - start) / CLOCKS_PER_SEC;

	if (status != 0)
		return status;
	*useraction_p = CPX_CALLBACK_SET;
	return 0;
}