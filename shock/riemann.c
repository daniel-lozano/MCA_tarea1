
/*
  Translated from a C++ Riemann solver written by Richard
  J. Gonsalves,  which in turn is rewrite of a Riemann solver
  found in Laney's textbook Computational Gasdynamics (Cambridge
  University Press.)     
 */ 
double min(double x,double y) { 
  if(x>y){
    return y;
  }
  if (x<y){
    return x;
  }

}

 double fg(double x) { //inline double fg(double x) {
    const double gamma = 1.4;
    const double g2 = (gamma + 1) / (2 * gamma);
    return (x-1) / sqrt(g2 * (x - 1) + 1);
}



 
/*
Translated from a C++ Riemann solver written by Richard
J. Gonsalves, which in turn is rewrite of a Riemann solver
found in Laney's textbook Computational Gasdynamics (Cambridge
University Press.)
*/
#ifndef RIEMANN_H_INCLUDED
#define RIEMANN_H_INCLUDED


void Riemann(double *U4, double *U1, double *F) {


  const double gamma = 1.4;
  const double g1 = (gamma-1) / (2*gamma);
  const double g2 = (gamma+1) / (2*gamma);
  const double g3 = (gamma+1) / (gamma-1);
  const double tol = 1e-10;


  int i;


  double x, y, z, fz;
  double a1, a2, a3, a4;
  double p1, p2, p3, p4;
  double u1, u2, u3, u4;
  double s1, s2, s3, s4;
  double rho1, rho2, rho3, rho4;
// compute primitive variables


 rho1 = U1[0];
 u1 = U1[1]/rho1;
 p1 = (U1[2]-rho1*u1*u1/2)*(gamma-1);
 a1 = sqrt(gamma*p1/rho1);
 rho4 = U4[0];
 u4 = U4[1]/rho4;
 p4 = (U4[2]-rho4*u4*u4/2)*(gamma-1);
 a4 = sqrt(gamma*p4/rho4);

 double pprima = pow(p4/p1, g1);

 double du = u4 - u1;
// apply the bisection method to find p2
 x = 0.05*p4/p1; //x=p2/p1
 y = 0.5*p4/p1; //y=p2/p1
 while(y-x > tol){
   z = x+(y-x)/2.0;
   fz = pprima
     - pow(z, g1) / (1 + g1 * (gamma * du - a1 * fg(z)) / a4);
   if (fz >0) {
     x=z;
   }
   else{
     y=z;
   }
 }
// Compute shock
 p2 = p1 * x;
 u2 = u1 + (a1/gamma)*((x-1)/sqrt(g2*(x-1)+1));
 a2 = a1 * sqrt(x*(g3+x)/(1+g3*x));
 rho2 = gamma*p2/(a2*a2);
//s = u1 + a1*sqrt(g2*(x-1)+1);
 s1 = u1 + a1*sqrt(g2*(x-1)+1);
// Compute contact
 p3 = p2;
 u3 = u2;
 a3 = (u4+2.0*a4/(gamma-1)-u3)*(gamma-1)/2.0;
 rho3 = gamma*p3/(a3*a3);
 s2 = u2;
//printf("%f,%f,%f,%f\n",p1,p2,p3,p4);
// Compute expansion
 s3 = u3 - a3;
 s4 = u4 - a4;
// Compute fluxes
 double f1, f2, f3;
 double a, u, p, rho;
// double e1, e2, e3, e4;
// e1 = p1/(gamma-1)+rho1*u1*u1/2.0;
// e2 = p2/(gamma-1)+rho2*u2*u2/2.0;
// e3 = p3/(gamma-1)+rho3*u3*u3/2.0;
// e4 = p4/(gamma-1)+rho4*u4*u4/2.0;
 if(s4 > 0) {
   f1 = rho4*u4;
   f2 = rho4*u4*u4 + p4;
   f3 = .5*rho4*u4*u4*u4 + rho4*a4*a4*u4/(gamma-1.);
 } else if (s3 > 0) {
   u = ((gamma-1.)*u4+2.*a4)/(gamma+1.);
   a = u;
   p = p4*pow(a/a4, 2.*gamma/(gamma-1.));
   if (a < 0 || p < 0) {
     printf("Negative a or p in Riemann");
   }
   rho = gamma*p/(a*a);
   f1 = rho*u;
   f2 = rho*u*u + p;
   f3 = .5*rho*u*u*u + rho*a*a*u/(gamma-1.);
 } else if (s2 > 0) {
   f1 = rho3*u3;
   f2 = rho3*u3*u3 + p3;
   f3 = .5*rho3*u3*u3*u3 + rho3*a3*a3*u3/(gamma-1.);
 } else if (s1 > 0) {
   f1 = rho2*u2;
   f2 = rho2*u2*u2 + p2;
   f3 = .5*rho2*u2*u2*u2 + rho2*a2*a2*u2/(gamma-1.);
 } else {
   f1 = rho1*u1;
   f2 = rho1*u1*u1 + p1;
   f3 = .5*rho1*u1*u1*u1 + rho1*a1*a1*u1/(gamma-1.);
 }
 F[0] = f1;
 F[1] = f2;
 F[2] = f3;
}
#endif /* RIEMANN_H_INCLUDED */

 



  /* 
    const double gamma = 1.4;
    const double g1 = (gamma - 1) / (2 * gamma);
    const double g2 = (gamma + 1) / (2 * gamma);
    const double g3 = (gamma + 1) / (gamma - 1);
    const double tol = 1e-10;
    double f1, f2, f3, a, u, rho;

    // compute primitive variables  4 -> L  y 1 -> R 
    double rho1 = U1[0];
    double u1 = U1[1] / rho1;
    double p1 = (U1[2] - rho1 * u1 * u1 / 2) * (gamma - 1);
    double h1 = (U1[2] + p1)/rho1;

    double rho4 = U4[0];
    double u4 = U4[1] / rho4;
    double p4 = (U4[2] - rho4 * u4 * u4 / 2) * (gamma - 1);
    double h4 = (U4[2] + p4)/rho4;

    // switch states if necessary so high pressure is on left
    bool revflag = false;
    if (p4 < p1) {
        double swap = p1; p1 = p4; p4 = swap;
        swap = u1; u1 = -u4; u4 = -swap;
        swap = rho1; rho1 = rho4; rho4 = swap;
        revflag = true;
    }

    double a1 = sqrt(gamma * p1 / rho1);
    double a4 = sqrt(gamma * p4 / rho4);
    double p = pow(p4/p1, g1);

    double du = u4 - u1;
    double dp = p4 - p1;
    double drho = rho4 - rho1 ;

    //compute r l 
    
    double rhol = rho4;
    double ul = u4;
    double pl = p4; 
    
    //compute Roe-averaged 
    
    double rhorl = sqrt(rho1*rho4);
    double url = (sqrt(rho4)*u4 + sqrt(rho1)*u1)/(sqrt(rho4) + sqrt(rho1));
    double hrl = (sqrt(rho4)*h4 + sqrt(rho1)*h1)/(sqrt(rho4) + sqrt(rho1));
    double arl = sqrt((gamma - 1 )*(hrl - 0.5*url*url));
    
    // compute Roe-averga wave speed
    
    double lambda1= url;
    double lambda2 = url + arl;
    double lambda3 = url - arl;
    
    //compute wave strengths
    
    double v1 = drho - (dp/(arl*arl));
    double v2 = du + (dp/(rhorl*arl));
    double v3 = du - (dp/(rhorl*arl));
    
    //compute solution
    
    f1 = rhol*ul + min(0,lambda1)*v1 + (0.5*rhorl/arl)*(min(0,lambda2)*v2 - min(0,lambda3)*v3);
    f2 = rhol*ul*ul+ pl+url*min(0,lambda1)*v1 + (0.5*rhorl/arl)*(lambda2*min(0,lambda2)*v2 -lambda3* min(0,lambda3)*v3);
    f3 =  rhol*h4*ul +  (0.5*url*url)*min(0,lambda1)*v1  + (0.5*rhorl/arl)*((hrl+arl*url)*min(0,lambda2)*v2 - (hrl-arl*url)*min(0,lambda3)*v3);


    F[0] = f1;
    F[1] = f2;
    F[2] = f3;
}

  */
    /*
    // apply the secant method
    // initial guesses
    double x = 0.05 * p4 / p1;
    double y = 0.5 * p4 / p1;
    double fx = p - pow(x, g1) / (1 + g1 * (gamma * du - a1 * fg(x)) / a4);
    double fy = p - pow(y, g1) / (1 + g1 * (gamma * du - a1 * fg(y)) / a4);
    bool converge = false;

    for (int i = 1; i <= 20; i++) {

        double z = y - fy * (y - x) / (fy - fx);
        double fz = p - pow(z, g1) / (1 + g1 * (gamma * du - a1 * fg(z)) / a4);

        if (abs(fz) < tol && abs(z - y) < tol) {
            
            break;
        }

        x = y;
        fx = fy;
        y = z;
        fy = fz;
    }




    // Compute shock
    double p2 = p1 * x;
    double u2 = u1 + a1 * fg(x) / gamma;
    //     u2 = u4 + 2.*a4*(1.-(x**g1)/p)/(gamma-1.)
    double a2 = a1 * sqrt(x * (g3 + x) / (1 + g3 * x));
    double rho2 = gamma * p2 / (a2 * a2);
    double s1 = u1 + a1 * sqrt(g2 *(x - 1) + 1);
    //     s1 = (rho1*u1 - rho2*u2)/(rho1-rho2)

    // Compute contact
    double p3 = p2;
    double u3 = u2;
    double a3 = a4 + 0.5 * (gamma - 1) * (u4 - u3);
    double s2 = u2;
    double rho3 = gamma * p3/(a3 * a3);

    // Compute expansion
    double s3 = u3 - a3;
    double s4 = u4 - a4;

    // Compute fluxes
    double f1, f2, f3, a, u, rho;
    if (revflag) {
        if (s4 > 0) {
            f1 = -rho4 * u4;
            f2 = rho4 * u4 * u4 + p4;
            f3 = -0.5 * rho4 * u4 * u4 * u4
                 - rho4 * a4 * a4 * u4 / (gamma - 1);
        } else if (s3 > 0) {
            u = (-(gamma-1.)*u4+2.*a4)/(gamma+1.);
            a = u;
            p = p4*pow(a/a4, 2.*gamma/(gamma-1.));
            
            rho = gamma*p/(a*a);
            f1 = -rho*u;
            f2 = rho*u*u + p ;
            f3 = -.5*rho*u*u*u - rho*a*a*u/(gamma-1.);
        } else if (s2 > 0) {
            f1 = -rho3*u3;
            f2 = rho3*u3*u3 + p3;
            f3 =  -.5*rho3*u3*u3*u3 - rho3*a3*a3*u3/(gamma-1.);
        } else if (s1 > 0) {
            f1 = -rho2*u2;
            f2 = rho2*u2*u2 + p2;
            f3 = -.5*rho2*u2*u2*u2 - rho2*a2*a2*u2/(gamma-1.);
        } else {
            f1 = -rho1*u1;
            f2 = rho1*u1*u1 + p1;
            f3 = -.5*rho1*u1*u1*u1 - rho1*a1*a1*u1/(gamma-1.);
        }
    } else {
        if(s4 > 0) {
            f1 = rho4*u4;
            f2 = rho4*u4*u4 + p4;
            f3 = .5*rho4*u4*u4*u4 + rho4*a4*a4*u4/(gamma-1.);
        } else if (s3 > 0) {
            u = ((gamma-1.)*u4+2.*a4)/(gamma+1.);
            a = u;
            p = p4*pow(a/a4, 2.*gamma/(gamma-1.));
            
            rho = gamma*p/(a*a);
            f1 = rho*u;
            f2 = rho*u*u + p;
            f3 = .5*rho*u*u*u + rho*a*a*u/(gamma-1.);
        } else if (s2 > 0) {
            f1 = rho3*u3;
            f2 = rho3*u3*u3 + p3;
            f3 =  .5*rho3*u3*u3*u3 + rho3*a3*a3*u3/(gamma-1.);
        } else if (s1 > 0) {
            f1 = rho2*u2;
            f2 = rho2*u2*u2 + p2;
            f3 = .5*rho2*u2*u2*u2 + rho2*a2*a2*u2/(gamma-1.);
        } else {
            f1 = rho1*u1;
            f2 = rho1*u1*u1 + p1;
            f3 = .5*rho1*u1*u1*u1 + rho1*a1*a1*u1/(gamma-1.);
        }
    }
*/
 
