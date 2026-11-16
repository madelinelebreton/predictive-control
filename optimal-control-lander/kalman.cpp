// kalman.cpp
// standard linear Kalman filter, with variance P and gain K
// observer for lander: estimate height, speed, and gravity from noisy height measurements


#include <iostream>
#include <stdexcept>

#include "kalman.h"

// constructor, model and noise matrices
KalmanFilter::KalmanFilter(
    double dt,
    const Eigen::MatrixXd& A,
    const Eigen::MatrixXd& B,
    const Eigen::MatrixXd& C,
    const Eigen::MatrixXd& Q,
    const Eigen::MatrixXd& R,
    const Eigen::MatrixXd& P)
    // model uncertainty Q, sensor noise R
  : A(A), B(B), C(C), Q(Q), R(R), P0(P),
    m(C.rows()), n(A.rows()), dt(dt), initialized(false),
    I(n, n), x_hat(n), x_hat_new(n)
{
  I.setIdentity();
}

KalmanFilter::KalmanFilter() {}

// starting guess x0 and uncertainty P0 (smaller means trust x0 more)
void KalmanFilter::init(double t0, const Eigen::VectorXd& x0) {
  x_hat = x0;
  P = P0;
  this->t0 = t0;
  t = t0;
  initialized = true;
}

void KalmanFilter::init() {
  x_hat.setZero();
  P = P0;
  t0 = 0;
  t = t0;
  initialized = true;
}

// the filtering
void KalmanFilter::update(const Eigen::VectorXd& y,double u) {

  if(!initialized)
    throw std::runtime_error("Filter is not initialized!");

  x_hat_new = A * x_hat +B*u; // predict the state
  P = A*P*A.transpose() + Q; // predict uncertainty
  K = P*C.transpose()*(C*P*C.transpose() + R).inverse(); // compute the gain
  x_hat_new += K * (y - C*x_hat_new); // correct the state
  P = (I - K*C)*P; // correct the uncertainty. P should settle to a constant before start controlling
  x_hat = x_hat_new;

  t += dt;
}

// add time step to A and dt
void KalmanFilter::update(const Eigen::VectorXd& y, double dt, const Eigen::MatrixXd A) {

  this->A = A;
  this->dt = dt;
  update(y,0);
}
