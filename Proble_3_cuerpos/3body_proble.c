#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>


void RK4(double dt,double *q1, double *q3, double *p1, double *p3,double e);
void leapfrog_step(double dt, double t, double *q1, double *q3,double *p1, double *p3,double e);
double q1prime(double dt, double q1, double q3, double p1, double p3, double e);
double q3prime(double dt, double q1, double q3, double p1, double p3, double e);
double p1prime(double dt, double q1, double q3, double p1, double p3, double e);
double p3prime(double dt, double q1, double q3, double p1, double p3, double e);

int  main(int argc, char **argv){

  double Tm=2800;
  double dt=0.006;
  double e=1.0;

  if(!argv[1]){
    fprintf(stderr,"introduzca un parametro de entrada\n");
    exit(1);
  }
  double condicion=atof(argv[1]);

  //Condiciones iniciales

  double T=0;
  double P=drand48();
  double Q=drand48();
  double q10= condicion;
  double q30=Q;
  double p10=0;
  double p30=P;

  double Q10= condicion;
  double Q30=Q;
  double P10=0;
  double P30=P;

  FILE *f;
  f=fopen("evolucion.dat","w");

  //escribiendo primera coordenada, condicion inicial

  fprintf(f,"%e %e %e %e %e %e %e %e %e \n",T,q10,q30,p10,p30,Q10,Q30,P10,P30);

  while(T<Tm){

    //avanza tiempo

    T=T+dt;
    RK4(dt,&q10,&q30,&p10,&p30,e);
    leapfrog_step(dt,T,&Q10,&Q30,&P10,&P30,e);
    fprintf(f,"%e %e %e %e %e %e %e %e %e \n",T,q10,q30,p10,p30,Q10,Q30,P10,P30);
   
  }
  fclose(f);
  return 0;
}

void RK4(double dt,double *q1, double *q3, double *p1, double *p3,double e){
  double k1,k2,k3,k4;//para q1
  double l1,l2,l3,l4;//para p1
  double m1,m2,m3,m4;//para q3
  double n1,n2,n3,n4;//para p3
  double q1in=*q1;
  double q3in=*q3;
  double p1in=*p1;
  double p3in=*p3;

  //paso 1
  k1=q1prime(dt,q1in,q3in,p1in,p3in,e);
  l1=p1prime(dt,q1in,q3in,p1in,p3in,e);
  m1=q3prime(dt,q1in,q3in,p1in,p3in,e);
  n1=p3prime(dt,q1in,q3in,p1in,p3in,e);

  //paso 2
  k2=q1prime(dt,q1in+k1*dt/2,q3in+m1*dt/2,p1in+l1*dt/2 ,p3in+n1*dt/2,e);
  l2=p1prime(dt,q1in+k1*dt/2,q3in+m1*dt/2,p1in+l1*dt/2 ,p3in+n1*dt/2,e);
  m2=q1prime(dt,q1in+k1*dt/2,q3in+m1*dt/2,p1in+l1*dt/2 ,p3in+n1*dt/2,e);
  n2=p1prime(dt,q1in+k1*dt/2,q3in+m1*dt/2,p1in+l1*dt/2 ,p3in+n1*dt/2,e);

  //paso 3
  k3=q1prime(dt,q1in+k2*dt/2,q3in+m2*dt/2,p1in+l2*dt/2 ,p3in+n2*dt/2,e);
  l3=p1prime(dt,q1in+k2*dt/2,q3in+m2*dt/2,p1in+l2*dt/2 ,p3in+n2*dt/2,e);
  m3=q1prime(dt,q1in+k2*dt/2,q3in+m2*dt/2,p1in+l2*dt/2 ,p3in+n2*dt/2,e);
  n3=p1prime(dt,q1in+k2*dt/2,q3in+m2*dt/2,p1in+l2*dt/2 ,p3in+n2*dt/2,e);  

  //paso 4
  k4=q1prime(dt,q1in+k3*dt,q3in+m3*dt,p1in+l3*dt ,p3in+n3*dt,e);
  l4=p1prime(dt,q1in+k3*dt,q3in+m3*dt,p1in+l3*dt ,p3in+n3*dt,e);
  m4=q1prime(dt,q1in+k3*dt,q3in+m3*dt,p1in+l3*dt ,p3in+n3*dt,e);
  n4=p1prime(dt,q1in+k3*dt,q3in+m3*dt,p1in+l3*dt ,p3in+n3*dt,e); 

  q1in+=(k1+2*k2+2*k3+k4)*dt/6.0;
  p1in+=(l1+2*l2+2*l3+l4)*dt/6.0;
  q3in+=(m1+2*m2+2*m3+m4)*dt/6.0;
  p3in+=(n1+2*n2+2*n3+n4)*dt/6.0;

  *q1=q1in;
  *q3=q3in;
  *p1=p1in;
  *p3=p3in;
  
}

void leapfrog_step(double dt, double t, double *q1, double *q3,double *p1, double *p3,double e){
  double q1in;
  double q3in;
  double p1in;
  double p3in;
  double a0_2=pow(2,-1/3)/(2*(2-pow(2,1/3)));
  double a1_2=1/(2*(2-pow(2,1/3)));
    
  q1in = *q1;
  p1in = *p1;
  q3in = *q3;
  p3in = *p3;
  
  /*kick*/
  p1in+=a1_2*p1prime(dt,q1in,q3in,p1in,p3in,e)* dt;
  /*drift*/
  q1in+=2*a1_2* p1in * dt;
  /*kick*/
  p1in+= a1_2* p1prime(dt,q1in,q3in,p1in,p3in,e)* dt;

  /*kick*/
  p1in+=a0_2*p1prime(dt,q1in,q3in,p1in,p3in,e)* dt;
  /*drift*/
  q1in+=2*a0_2* p1in * dt;
  /*kick*/
  p1in+= a0_2* p1prime(dt,q1in,q3in,p1in,p3in,e)* dt;

   /*kick*/
  p1in+=a1_2*p1prime(dt,q1in,q3in,p1in,p3in,e)* dt;
  /*drift*/
  q1in+=2*a1_2* p1in * dt;
  /*kick*/
  p1in+= a1_2* p1prime(dt,q1in,q3in,p1in,p3in,e)* dt;

  //integra la particula pequen\~a

  /*kick*/
  p3in+=a1_2*p3prime(dt,q1in,q3in,p1in,p3in,e)* dt;
  /*drift*/
  q3in+=2*a1_2* p3in * dt;
  /*kick*/
  p3in+= a1_2* p3prime(dt,q1in,q3in,p1in,p3in,e)* dt;
  
  /*kick*/
  p3in+=a0_2*p3prime(dt,q1in,q3in,p1in,p3in,e)* dt;
  /*drift*/
  q3in+=2*a0_2* p3in * dt;
  /*kick*/
  p3in+= a0_2* p3prime(dt,q1in,q3in,p1in,p3in,e)* dt;
  
  /*kick*/
  p3in+=a1_2*p3prime(dt,q1in,q3in,p1in,p3in,e)* dt;
  /*drift*/
  q3in+=2*a1_2* p3in * dt;
  /*kick*/
  p3in+= a1_2* p3prime(dt,q1in,q3in,p1in,p3in,e)* dt;
  


 

    

  *q1=q1in;
  *p1=p1in;
  *q3=q3in;
  *p3=p3in;
}




double q1prime(double dt,double q1, double q3, double p1, double p3, double e){

  return p1;

}
double q3prime(double dt,double q1, double q3, double p1, double p3, double e){
  return p3;
}

double p1prime(double dt, double q1, double q3, double p1, double p3, double e){
  return -2.0*q1/pow(4*q1*q1+e*e,3/2);

}
double p3prime(double dt, double q1, double q3, double p1, double p3, double e){

  return (q1-q3)/pow((q1-q3)*(q1-q3)+e*e/4,3/2) - (q1+q3)/pow((q1+q3)*(q1+q3)+e*e/4,3/2);

}
