#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;


//Constructor
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_ = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // Finish initializing the FusionEKF
  
  // Set the process and measurement noises
  noise_ax = 9.0;
  noise_ay = 9.0;

}

// Destructor
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  
  // Initialization
  
  if (!is_initialized_) {

    // Initialize the state ekf_.x_ with the first measurement.
    // TODO: Create the covariance matrix.

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {  
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      float rho = measurement_pack.raw_measurements_[0];
      float angle = measurement_pack.raw_measurements_[1];
      ekf_.x_ << rho*cos(angle), rho*sin(angle), 0.0, 0.0;
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0, 0.0;

    }

    // get the previous timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // initialize the state covariance matrix
    // since we do not any knowledge of the velocities before, we will
    // set their covariance to a large value, something like 100 or 1000

    // object covariance matrix
    MatrixXd P(4,4);
    P << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 100, 0,
         0, 0, 0, 100;

    // state transition matrix
    MatrixXd F(4,4);
    F << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;

    // initialize laser measurement matrix

    H_ << 1, 0, 0, 0,
          0, 1, 0, 0;

    // measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
              0, 0.0225;

    // measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
    
    // initialize with first state vector
    MatrixXd Q(4,4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      ekf_.Init(ekf_.x_, P, F, H_, R_radar_, Q);
    }
    else
    {
      ekf_.Init(ekf_.x_, P, F, H_, R_laser_, Q);
    }
    

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction Step
  */
 
  // TODO: Update the state transition matrix F according to the new elapsed time.

  float dt;

  dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;

  previous_timestamp_ = measurement_pack.timestamp_;

  // check if time difference is greater than 0
  if (dt > 0)
  {
    // motion model updated
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    // update covariance matrices

    // START HERE
    float dt2 = dt * dt;
    float dt3 = dt * dt * dt;
    float dt4 = dt * dt * dt * dt;

    ekf_.Q_ << 0.25*dt4*noise_ax, 0.5*dt3*noise_ax, 0,
         0, 0.25*dt4*noise_ay, 0, 0.5*dt3*noise_ay,
         0.5*dt3*noise_ax, 0, dt2*noise_ax, 0,
         0, 0.5*dt3*noise_ay, 0, dt2*noise_ay;

    ekf_.Predict();

  }

  // Use the sensor type to perform the update step.
  // Update the state and covariance matrices.

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
    // For radar update, we use Extended Kalman Filter
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else
  {
    // Laser updates
    // For Lidar update, we use Kalman Filter
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
