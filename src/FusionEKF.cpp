#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::cout;
using std::endl;


// Constructor

FusionEKF::FusionEKF() {

  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices

  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // 1measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  // Finish initializing the FusionEKF.

  // Set the process and measurement noises
  noise_ax = 9.0;
  noise_ay = 9.0;

}

// Destructor

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  // Initialization //

  if (!is_initialized_) 
  {
    // Initialize the state ekf_.x_ with the first measurement.
    // create empty state vector

    VectorXd x(4);

    // object covariance matrix
    // since we do not any knowledge of the velocities before, we will
    // set their covariance to a large value, something like 100 or 1000

    MatrixXd P(4, 4);
    P << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1000, 0,
         0, 0, 0, 1000;

    if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
    {
      // Convert radar from polar to cartesian coordinates and initialize state.

      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rhot_dot = measurement_pack.raw_measurements_[2];

      // assign measurements
      x<<rho*cos(phi), rho*sin(phi), 0, 0;
      cout<<"First measurements are from radar!"<<endl;
    }
    else if(measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
    {
      // Initialize state.
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      cout<<"First measurements are from LiDAR!"<<endl;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // Initialize transition matrix
    MatrixXd F(4, 4);
    F << 1, 0, 0, 0,
	       0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;

    // Initialize measurement matrix for laser measurements
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    // initialize empty Q matrix
    MatrixXd Q(4,4);

    // Initialize ekf_ with the first state vector, 
    ekf_.Init(x, P, F, H_laser_, Hj_, R_laser_, R_radar_, Q);

    // done initializing, no need to predict or update
    // set initialized to true
    is_initialized_ = true;
    return;

  }


  // Prediction Step //

  // Update the state transition matrix F according to the new elapsed time.
  // Compute the time from the previous measurement in seconds (we divide as time is in microseconds)

  float dt;
  
  dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;
  
  // Update the motion model matrix for a timestep dt.

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // update covariance matrices
  double dt2 = dt * dt;
  double dt3 = dt * dt * dt;
  double dt4 = dt * dt * dt * dt;
  
  ekf_.Q_ << (dt4/4)*noise_ax, 0, (dt3/2)*noise_ax, 0,
            0, (dt4/4)*(noise_ay), 0, (dt3/2)*(noise_ay),
            (dt3/2)*(noise_ax), 0, (dt2*noise_ax), 0,
            0, (dt3/2)*(noise_ay), 0, dt2*noise_ay;

  ekf_.Predict();

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 
  else 
  {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout<< "x_ = " << ekf_.x_ <<endl;
  cout<< "P_ = " << ekf_.P_ <<endl;

}

