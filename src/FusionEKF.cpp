#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  is_first_predict_= true;

 // initializing matrices------------------|
  
  //measurement matrix - laser
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1,0,0,0,
              0,1,0,0;

  //measurement matrix - radar
  Hj_ = MatrixXd(3, 4);
  Hj_.setZero();


  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // preparing the for the initial transition matrix F_
  
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;



  // to initial  state covariance matrix 
  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;
         //  high incertitude (1000), on carthesian speed,  in case the first measuremnt comes from a radar 
         //( since Rho dot does not convert exactly in vx vy , we initialized it at vx, vy at zero)
         // this incertitude will diminish through iterations as the filter converges

  // process civariance matrix place holder   
  Q_ = MatrixXd(4, 4);
  Q_.setZero();

  // preparing the for the initial state vector 

  x_ = VectorXd(4);
  x_ << 1, 1, 1, 1;  


   //* TODO: Set the process and measurement noises 


}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) 

  {
  


    // to initial transition matrix F_, 
  
    //ekf_.F_ = MatrixXd(4, 4);
    /*ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;
    */

     

    // TODO: Initialize the state ekf_.x_ with the first measurement.
    //cout << "EKF: " << endl;
    /*
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    */



    // vector x_ to store the first measurement.
    // this vector x_ will be pass as an argument of the Kalman Filter init function to initialize the state vector 
    x_ = VectorXd(4);
    x_ << 1, 1, 1, 1;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 

    {
      // Convert radar measurement from polar to cartesian coordinates .
      // initialize a vector with carthesian coordinate ( 2D inital position ) and zero carthesian velocity

      x_ << measurement_pack.raw_measurements_(0)* cos(measurement_pack.raw_measurements_(1)), 
            measurement_pack.raw_measurements_(0)* sin(measurement_pack.raw_measurements_(1)), 
            0,
            0;
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 

    {
      // use laser measurement carthesian coordinates (2D inital position) and zero velocity ( not measured by Laser )

      x_ << measurement_pack.raw_measurements_(0), 
            measurement_pack.raw_measurements_(1), 
            0, 
            0;
    }

    ekf_.Init(x_, P_, F_, Hj_, R_radar_, Q_);
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction: just pass F, Q in init !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   */

    
    //Compute the time elapsed between the current and previous measurements
    
    float dt = (measurement_pack.timestamp_ - previous_timestamp_)/ 1000000.0 ; // feed is in microsecondes 
    previous_timestamp_ = measurement_pack.timestamp_;

    //Update the state transition matrix F according to the new elapsed time.

    F_(0, 2) = dt;
    F_(1, 3) = dt;
    

    //  process noise is set at 9 for this project 
    int noise_ax = 9;
    int noise_ay = 9;


    // set the process covariance matrix Q
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
           0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
           dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
           0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


  if (is_first_predict_)
    {
    ekf_.Init(ekf_.x_, P_, F_, Hj_, R_radar_, Q_); // F_, Hj_, R_radar_, Q_ are only initilaized as placeholders, 
                                              // but not used in preduction 
    is_first_predict_=false;
    }

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {
    //  Radar updates
    // Hj_ needs to be calculated at everypass
    Hj_= tools.CalculateJacobian(ekf_.x_); 
    ekf_.Init(ekf_.x_, ekf_.P_, F_, Hj_, R_radar_, Q_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } 
  else 
  {
    // TODO: Laser updates
    ekf_.Init(ekf_.x_, ekf_.P_, F_, H_laser_, R_laser_, Q_);
    ekf_.Update(measurement_pack.raw_measurements_);
    //cmake .. && make
  // ./ExtendedKF


  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}





