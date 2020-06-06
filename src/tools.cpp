#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth)
{
   // TODO: Calculate the RMSE here.
   VectorXd rmse(4);
   rmse << 0,0,0,0;
   
   // the estimation vector size should not be zero
   // the estimation vector size should equal ground truth vector size
   if(estimations.size() == 0 || estimations.size() != ground_truth.size())
   {
      cout<<"Invalid size of estimations and ground truth vectors!!!"<<endl;
      return rmse;
   }
   
   for (int i=0; i < estimations.size(); ++i) {
      VectorXd residual = estimations[i] - ground_truth[i];
      residual = residual.array()*residual.array();
      rmse += residual;
  }
  
  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
    MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // TODO: YOUR CODE HERE 

  // check division by zero
  if (px == 0.0 & py == 0.0)
  {
      cout<<"Error! - Division by zero!"<<endl;
      return;
  }
  else
  {
  // compute the Jacobian matrix
  float rho = (px*px + py*py);
  
  float elem11 = px/(sqrt(rho));
  float elem12 = py/(sqrt(rho));
  
  float elem21 = (-1)*(py/rho);
  float elem22 = (px/rho);
  
  float elem31 = (py*(vx*py - vy*px))/(sqrt(rho*rho*rho));
  float elem32 = (px*(vy*px - vx*py))/(sqrt(rho*rho*rho));
  float elem33 = px/sqrt(rho);
  float elem34 = py/sqrt(rho);
  
  Hj<< elem11, elem12, 0, 0,
        elem21, elem22, 0, 0,
        elem31, elem32, elem33, elem34;
  }
  return Hj;

}
