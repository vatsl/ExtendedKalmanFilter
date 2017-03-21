#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0.0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
                0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

    // measurement function matrix - laser
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    // state covariance matrix P
    // diagonal elements only
    // low values for high certainity and high values for low certainity
    MatrixXd P_ = MatrixXd(4, 4);
    P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

    // state transition matrix
    MatrixXd F_ = MatrixXd(4, 4);
    //the initial transition matrix F_
    F_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

    // new covariance matrix based on noise vector
    MatrixXd Q_ = MatrixXd(4, 4);
    Q_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          1, 0, 1, 0,
          0, 1, 0, 1;

  /**
    * Initializing the FusionEKF.
    * Set the process and measurement noises
  */

    // noise measurements
    noise_ax = 9.0;
    noise_ay = 9.0;

    // Will be overwritten by tools.CalculateJacobian()
//    Hj_ << 1, 0, 0, 0,
//           0, 1, 0, 0,
//           0, 0, 1, 0;


    VectorXd x_ = VectorXd(4);
    x_ << 1, 1, 1, 1;
    ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
      // first measurement
      cout << "EKF: " << endl;
      ekf_.x_ = VectorXd(4);
      ekf_.x_ << 1, 1, 1, 1;

      // all variables for position and velocity unknown
      float px = 0.0;
      float py = 0.0;
      float vx = 0.0;
      float vy = 0.0;

      if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
         Convert radar from polar to cartesian coordinates and initialize state.
      */
          float rho = measurement_pack.raw_measurements_[0];
          float phi = measurement_pack.raw_measurements_[1];
          float rho_dot = measurement_pack.raw_measurements_[2];

          px = rho * cos(phi);
          py = rho * sin(phi);
          //vx = rho_dot * cos(phi);
          //vy = rho_dot * sin(phi);

      }

      else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
         Initialize state.
      */
          px = measurement_pack.raw_measurements_[0];
          py = measurement_pack.raw_measurements_[1];
          //vx = 0;
          //vy = 0;

      }

      ekf_.x_ << px, py, vx, vy;
      previous_timestamp_ = measurement_pack.timestamp_;
      // done initializing, no need to predict or update
      is_initialized_ = true;
      cout << "EKF Initialized " << endl;
      return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

    // dt - time difference in seconds
    float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    float dt2 = dt*dt;
    float dt3 = (dt2*dt);
    float dt4 = (dt3*dt);

    dt3 = dt3/2;
    dt4 = dt4/4;

    // noise for accelration
    noise_ax = 9;
    noise_ay = 9;

    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    // new covariance matrix
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt4*noise_ax, 0, dt3*noise_ax, 0,
               0, dt4*noise_ay, 0, dt3*noise_ay,
               dt3*noise_ax, 0, dt2*noise_ax, 0,
               0, dt3*noise_ay, 0, dt2*noise_ay;

    ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates
      try {
          ekf_.R_ = R_radar_;
          ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
          ekf_.UpdateEKF(measurement_pack.raw_measurements_);
      } catch(...){
          // Ignore if errors - mostly during the Jacobian calculations
          return;
      }

  } else {
      // Laser updates
      ekf_.R_ = R_laser_;
      ekf_.H_ = H_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
