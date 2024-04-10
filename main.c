#include "header.h"

double** c2, ** c, * f, * mu, ** W;
double** c_c, ** c_t, ** c_d;
int N;
COMMOD comm;
coordinate* pts;
bitf* Ak;
xinfo* node;
crit1* arc1;
crit2* arc2;
FHubs FeasHubs;

double apar, bpar, vel;
char instancia[10];
double alpha, cuttime;
double muhat, * tau;
double*** lmaxk, * lmax;
int svars, * nbFeas, ** nCom;
int xvars,* xa1, * xa2, * xk;
int yvars, * ypos;
int* FeasArc, * FeasArcPos;
int counter = 0;

double alphas[5] = { 0.8, 0.85, 0.9, 0.95, 0.99 };	// Service level requirements
double rq = 3.0;	// Definition of parameter r (Equation 41)
double eta = 1.5;	// Time correction factor for the inter-hub transport (Equation 42)

// Main function
int main() {
	FILE* ini;
	int num_inst;

	ini = open_file("Problem_Instances_ap.txt", "r");
	fscanf(ini, "%d", &num_inst);

	for (int i = 0; i < num_inst; i++) {

		fscanf(ini, "%s", &instancia);

		for (int j = 0; j < 5; j++) {
			alpha = alphas[j];
			cuttime = 0;

			printf("Instance: %.*s\nService level: %.2lf\nParameter r: %.1lf\n\n", 4, instancia, 100 * alpha, rq);

			read_ap(instancia);

			//modelM1();
			//modelM2();
			modelM3();

			free_memory();

			counter++;
		}
	}
	fclose(ini);
}