#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "riemann.c"


/*
  C adaptation from C++ code written by Richard J. Gonsalves.
*/ 


double L = 4.0;                    // length of shock tube
double gama = 1.4;               // ratio of specific heats
int N = 5000;                     // number of grid points

double CFL = 0.4;                // Courant-Friedrichs-Lewy number
//double nu = 0.0;                 // artificial viscosity coefficient

double **U = NULL;                      // solution with 3 components
//double **newU = NULL;                   // new solution
double **F = NULL;                      // flux with 3 components
//double *vol = NULL;                     // for Roe solver

double h;                        // lattice spacing
double tau;                      // time step
double c;

int step;

void allocate();
double cMax();
void initialize();
void boundaryConditions(double **U);
void upwindGodunovStep();
void solve(double tMax);


int main()
{
  //solve(LaxFriedrichsStep, 1.0, "LaxFriedrichs", 5);
  solve( 1.0);
   return 0; 
}


void allocate() {
  int j;

  U = malloc(N * sizeof(double *));
  // newU = malloc(N * sizeof(double *));
  F = malloc(N * sizeof(double *));
  // vol = malloc(N * sizeof(double *));
  for (j = 0; j < N; j++) {
    U[j] = malloc(3 * sizeof(double));
    //newU[j] = malloc(3 * sizeof(double));
    F[j] = malloc(3 * sizeof(double));
  }
}

void initialize() {
  int j;
  double rho,p,u,e;
  allocate();
  h = 1.0 * L / (N - 1);

  for (j = 0; j < N; j++) {
    rho = 1; p = 1; u = 0;
    if (j > N / 2){
      rho = 0.125;
      p = 0.1;
    }
    e = p / (gama - 1) + rho * u * u / 2.0;
    U[j][0] = rho;
    U[j][1] = rho * u;
    U[j][2] = e;
    //vol[j] = 1;

    
    F[j][0] = rho*u;
    F[j][1] = rho*u*u+p;
    F[j][2] = u*(e+p);
  }
  tau = CFL * h / cMax();
  step = 0;
}


double cMax() {
    double uMax = 0;
    double rho, u, p, c;
    int i;
    for (i = 0; i < N; i++) {
      if (U[i][0] == 0)
	continue;

      rho = U[i][0];
      u = U[i][1] / rho;
      p = (U[i][2] - rho * u * u / 2) * (gama - 1);
      c = sqrt(gama * fabs(p) / rho);

      if (uMax < (c + fabs(u)))
	uMax = c + fabs(u);
    }    
    return uMax;
}


void boundaryConditions(double **U) {

    // reflection boundary conditions at the tube ends
  U[0][0] = U[1][0];
  U[0][1] = -U[1][1];
  U[0][2] = U[1][2];
  U[N - 1][0] = U[N - 2][0];
  U[N - 1][1] = -U[N - 2][1];
  U[N - 1][2] = U[N - 2][2];
}

void upwindGodunovStep() {
  int i, j;
    // find fluxes using Riemann solver
  for (j = 0; j < N - 1; j++){
        Riemann(U[j], U[j + 1], F[j]);
	//printf(" %f %f %f \n" , F[j][0] , F[j][1] , F[j][2]);
  }
    // update U
  for (j = 1; j < N - 1; j++){
    for (i = 0; i < 3; i++){
      U[j][i] -= tau / h * (F[j][i] - F[j - 1][i]);
    }
  }
}





void solve( double tMax)
{
    initialize();

    double t = 0.0;
 
 
    FILE *out;
    
    int j;
    
    double rho, u, e, P;

    tau = CFL * h / cMax();
 

      
      

      while (t < tMax) {
	boundaryConditions(U);
	tau = CFL * h / cMax();
	upwindGodunovStep();
	t += tau;

      }
      
      out = fopen("UpWindGoudonovfinal.dat" , "w");
      
      for (j =0 ; j< N; j++) {

	rho = U[j][0];
	u = U[j][1]/rho;
	e = U[j][2];
	P = (gama-1.0)*(e-rho*u*u/2.0);
	fprintf(out, "%d\t%f\t%f\t%f\t%f\n", j, rho, u, e, P);
	
      }

    

      fclose(out);   
    
}


