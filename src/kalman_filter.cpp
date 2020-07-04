#include "kalman_filter.h"
#define PI 3.14159265

#include <iostream>
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &Hj_in, MatrixXd &R_in,
                        MatrixXd &R_ekf_in, MatrixXd &Q_in)
{
  cout << "In KalmanFilter::Init" << endl;
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  
  Hj_ = Hj_in;
  R_ = R_in;
  
  R_ekf_ = R_ekf_in;
  Q_ = Q_in;
  
  I_ = Eigen::MatrixXd::Identity(4,4);
}

void KalmanFilter::Predict() 
{
  // predict the state
  x_ = F_*x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_*P_*Ft + Q_;
}

// Kalman filter used for Lidar measurements
void KalmanFilter::Update(const VectorXd &z) 
{
  // Update the state using Kalman Filter equations
  VectorXd y = z - H_* x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_*P_*Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_*Ht*Si;

  // New state
  x_ = x_ + (K * y);
  P_ = (I_ - (K * H_))*P_;
}

// Extended Kalman Filter used for Radar measurements
void KalmanFilter::UpdateEKF(const VectorXd &z) 
{
  float px, py, vx, vy;
  px = x_[0];
  py = x_[1];
  vx = x_[2];
  vy = x_[3];

  // If rho == 0, skip the update step to avoid dividing by zero.
  // This is crude but should be fairly robust on our data set.
  if(px == 0. && py == 0.)
    return;

  Hj_ = tools.CalculateJacobian(x_);

  // for radar measurement, we have h(x)
  VectorXd hofx(3);

  float rho = sqrt(px*px + py*py);
  float theta = atan2(py, px);
  float rho_dot = (px*vx+py*vy)/rho;

  hofx << rho, theta, rho_dot;

  // Update the state using xtended Kalman Filter equations
  VectorXd y = z - hofx;

  if(y[1] > PI)
    y[1] -= 2.f * PI;
  if(y[1] < -PI)
    y[1] += 2.f * PI;

  MatrixXd Hjt = Hj_.transpose();
  MatrixXd S = (Hj_* P_ * Hjt) + R_ekf_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Hjt * Si;

  // Compute new state
  x_ = x_ + (K * y);
  P_ = (I_ - (K * Hj_)) * P_;
}
