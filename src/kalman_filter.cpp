#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;

  // defined
  
  // define identity matrix here ...
  I = Eigen::MatrixXd::Identity(4,4);
  

}

void KalmanFilter::Predict()
{
  // predict state
  x_ = F_*x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_*P_*Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z)
{
  // Update the state using Kalman Filter equations

  VectorXd y = z - H_*x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_*P_*Ht + R_;

  MatrixXd Si = S.inverse();

  // Calculate Kalman gain
  MatrixXd K = P_*Ht*Si;
  
  // Calculate new state and covariance matrix
  x_ = x_ + K*y;
  P_ = (I - K*H_) *P_;

  std::cout<<"Update for Lidar executed!!!"<<std::endl;

}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  // update the state by using Extended Kalman Filter equations

  float pos_x = x_[0];
  float pos_y = x_[1];
  float vel_x = x_[2];
  float vel_y = x_[3];

  // If rho == 0, avoid division by zero.

  if( pos_x == 0. && pos_y == 0. )
    return;

  Tools tools;

  MatrixXd Hj = tools.CalculateJacobian(x_);

  // vector h(x') definition
  VectorXd hof_x(3);

  float rho = sqrt( pos_x*pos_x + pos_y*pos_y );
  float angle = atan2(pos_y, pos_x);
  float rho_dot = (pos_x*vel_x + pos_y*vel_y)/rho;

  hof_x << rho, angle, rho_dot; 

  // Update the state using Extended Kalman Filter equations
  VectorXd y = z - hof_x;

  // if the angle exceeds +180 or -180 degrees, add +/- 2PI
  if(y[1] > M_PI)
    y[1] -= 2.f*M_PI;
  if(y[1] < -M_PI)
    y[1] += 2.f*M_PI;


  // Calculations
  MatrixXd Hjt = Hj.transpose();
  MatrixXd S = Hj*P_*Hjt + R_;
  MatrixXd Si = S.inverse();
  
  // Calculate Kalman Gain
  MatrixXd K = P_*Hjt*Si;

  // Compute new state
  x_ = x_ + (K*y);
  P_ = (I - K*Hj)*P_;

  std::cout<<"Update EKF for Radar executed!!!"<<std::endl;

}
