// main.cpp
// implements a finite horizon LQR to control a lander system

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <cstdio> // for gnuplot


#include "kalman.h" // for state estimates
#include "DARE.h" // controller


#define P_SIGMA_H  0.01
#define P_SIGMA_G 0.1

#define M_SIGMA 1.


double g,hight,v,t;

double dt;

// from fall.cpp
void fall_init(double x0, double gravity)
{
	g=gravity;
        hight= x0;
        v=0;

}


// falling with noise injected and control input u
double falling(double noise_p,double noise_g,double u)
{

  hight=hight+dt*v + dt*noise_p;
  v=v-dt*(g+noise_g)+dt*u;

  return(hight);

}



// PLOT with gnuplot
void plot(void)
{
  FILE* gp = popen("gnuplot", "w");
  fprintf(gp, "set terminal png\n");

  fprintf(gp, "set output 'landing.png'\n");
  fprintf(gp, "plot 'DATA' u 1:6 w l t 'height', '' u 1:7 w l t 'velocity'\n"); // using, with lines, title

  fprintf(gp, "set output 'gravity.png'\n");
  fprintf(gp, "plot 'DATA' u 1:(-$5) w l t 'estimated gravity'\n");

  pclose(gp);
}


int main(int argc, char* argv[]) {

  int n = 3; // Number of states
  int m = 1; // Number of measurements

  // SOME NOISE
  std::default_random_engine generator;
  std::normal_distribution<double> p_noise_h(0.0,P_SIGMA_H*P_SIGMA_H);
  std::normal_distribution<double> p_noise_g(0.0,P_SIGMA_G*P_SIGMA_G);
  std::normal_distribution<double> m_noise(0.0,M_SIGMA*M_SIGMA);


  dt = 1.0/30; // Time step

  Eigen::MatrixXd A(n, n); // System dynamics matrix
  Eigen::MatrixXd B(n, m); // Output matrix /// NOT USED IN JUST FALLING

  Eigen::MatrixXd C(m, n); // Output matrix
  Eigen::MatrixXd Q(n, n); // Process noise covariance
  Eigen::MatrixXd R(m, m); // Measurement noise covariance
  Eigen::MatrixXd P(n, n); // Estimate error covariance



  A << 1.,dt,0,  0,1,dt,  0,0,1;
  B << 0.,dt,0;   // thrust changes the velocity: v += dt*u 

  C << 1,0,0;

  // Reasonable covariance matrices
  Q <<P_SIGMA_H*P_SIGMA_H  ,0.,0, 0,0,0, 0,0,P_SIGMA_G*P_SIGMA_G ;
  R << M_SIGMA*M_SIGMA;

  P <<  1,0,0,  0,1,0, 0,0,1;

  std::cout << "A: \n" << A << std::endl;
  std::cout << "C: \n" << C << std::endl;
  std::cout << "Q: \n" << Q << std::endl;
  std::cout << "R: \n" << R << std::endl;
  std::cout << "P: \n" << P << std::endl;


  KalmanFilter kf(dt,A,B, C, Q, R, P);


  Eigen::VectorXd x0(n);
  // third is initial guess of gravity
  x0 << 110,0,0;
  //x0 << 100,0,2; ///  WOrks even with wrong sign
  kf.init(dt,x0);

  // CONTROLLER (finite-horizon LQR) same A and B as KF.
  // thrust is cheap (small Rc), so the controller can brake hard enough
  int horizon = 450; // long enough to watch system converge 
  Eigen::MatrixXd Qc(n, n); // state cost
  Eigen::MatrixXd Rc(1, 1);   // control cost
  Qc.setZero();
  Rc << 1e-4;

  // Construct the  SYSTEM
  DARE contr(A, B, Qc, Rc, horizon);

  contr.init();

  fall_init(100,5.5); // real plant conditions

  FILE* log;
  log=fopen("DATA","w");


  double t = 0;
  int num=100;
  Eigen::VectorXd y(m);


  y << falling(0.,0.,0.);

  double u = 0; // thrust
  int i = 0; // step number, picks the gain K_i

  // from fall.cpp
  //printf("%f\t %f\t%f\t%f\t%f \n",t, y(0),  kf.state()[0],kf.state()[1],kf.state()[2]);
  fprintf(log,"%f\t %f\t%f\t%f\t%f\t%f\t%f\t%f \n",t, y(0),  kf.state()[0],kf.state()[1],kf.state()[2],
          hight, v, u); // using true h, v and u

  while(hight>0. && i < horizon-1){ // true height (y is noisy) to stop at ground or end of horizon
    t += dt;
    kf.update(y,u);  

    // CONTROL on the ESTIMATE: u = -K_i * x_hat   
    u = -(contr.getK(i) * kf.state())(0);
    if (u < 0) u = 0;  // limit the lander to only brake (not truly optimal)
    i++;

    //printf("%f\t %f\t%f\t%f\t%f \n",t, y(0),  kf.state()[0],kf.state()[1],kf.state()[2]);
    fprintf(log,"%f\t %f\t%f\t%f\t%f\t%f\t%f\t%f \n",t, y(0),  kf.state()[0],kf.state()[1],kf.state()[2],
            hight, v, u);                             
    y << falling(p_noise_h(generator),p_noise_g(generator),u); 
    y(0)+=m_noise(generator);
  }

  fclose(log);
  std::cout << "Touchdown: t= " << t << "  h= " << hight << "  v= " << v
            << "  g_est= " << -kf.state()[2] << std::endl;
  plot();

  return 0;
}
