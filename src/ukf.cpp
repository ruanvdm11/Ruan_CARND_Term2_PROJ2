#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;
  time_us_ = 0.0;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3-n_x_;
  weights_ = VectorXd(2*n_aug_+1);
  P_ << 1, 0, 0, 0, 0,
  		0, 1, 0, 0, 0,
  		0, 0, 1, 0, 0,
  		0, 0, 0, 1, 0,
  		0, 0, 0, 0, 1;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for(int i = 1; i < 2*n_aug_+1; i++){
  	double weight_ = 0.5/(n_aug_+lambda_);
  	weights_(i) = weight_;
  }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	if(is_initialized_ == false){
		std::cout<<"Initializing x state"<<std::endl;
		if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
			//
			float roh = meas_package.raw_measurements_[0];
      		float phi = meas_package.raw_measurements_[1];
      		float rohdot = meas_package.raw_measurements_[2];
      		phi = atan2(sin(phi),cos(phi));
      		x_ << roh*cos(phi), roh*sin(phi), 0, 0, 0;
      		std::cout<<"x_initial = "<<x_<<std::endl;
      		use_radar_ = true;

			//
		}
		else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
			//
			//x_ <<
			std::cout<<"Initializing with laser"<<std::endl;
			use_laser_ = true;
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0., 0., 0.;
			std::cout<<"x_initial = "<<x_<<std::endl;

			//
		}
		time_us_ = meas_package.timestamp_;

		is_initialized_ = true;
		return;
		
	}

	float delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
	time_us_ = meas_package.timestamp_;
	Prediction(delta_t);
	if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
		UpdateRadar(meas_package);
	}
	else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
		UpdateLidar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
	///////////////// GENERATE SIGMA POINTS /////////////////

	///////////////// Augmentation Sigma PTS /////////////////
	VectorXd x_aug_(n_aug_);
	MatrixXd P_aug_(n_aug_,n_aug_);
	MatrixXd Xsig_aug_(n_aug_, 2*n_aug_+1);
	x_aug_.head(n_x_) = x_;
	x_aug_(5) = 0;
	x_aug_(6) = 0;

	P_aug_.fill(0.0);
	P_aug_.topLeftCorner(n_x_,n_x_) = P_;
	P_aug_(5,5) = std_a_*std_a_;
	P_aug_(6,6) = std_yawdd_*std_yawdd_;

	MatrixXd L = P_aug_.llt().matrixL();

	Xsig_aug_.col(0) = x_aug_;
	for(int j = 0; j < n_aug_; j++){
		Xsig_aug_.col(j+1)        = x_aug_ + sqrt(lambda_+n_aug_) * L.col(j);
		Xsig_aug_.col(j+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L.col(j);
	}
	///////////////// Predicting Sigma PTS /////////////////
	
	for(int k = 0; k < 2*n_aug_+1; k++){
		// Values for easier readability
		double p_x_ = Xsig_aug_(0,k);
		double p_y_ = Xsig_aug_(1,k);
		double v_ = Xsig_aug_(2,k);
		double yaw_ = Xsig_aug_(3,k);
		double yawd_ = Xsig_aug_(4,k);
		double nu_a_ = Xsig_aug_(5,k);
		double nu_yawdd_ = Xsig_aug_(6,k);

		double px_p, py_p;
		
		// Avoid division by zero
		if(fabs(yawd_)>0.001){
			px_p = p_x_ + v_/yawd_*(sin(yaw_+yawd_*delta_t) - sin(yaw_));
			py_p = p_y_ + v_/yawd_*(cos(yaw_)-cos(yaw_+yawd_*delta_t));
		} else{
			px_p = p_x_ + v_*delta_t*cos(yaw_);
			py_p = p_y_ + v_*delta_t*sin(yaw_);
		}

		double v_p = v_;
		double yaw_p = yaw_ + yawd_*delta_t;
		double yawd_p = yawd_;

		// Add noise
		px_p = px_p + 0.5*nu_a_*delta_t*delta_t*cos(yaw_);
		py_p = py_p + 0.5*nu_a_*delta_t*delta_t*sin(yaw_);
		v_p = v_p + nu_a_*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd_*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd_*delta_t;

		Xsig_pred_(0,k) = px_p;
		Xsig_pred_(1,k) = py_p;
		Xsig_pred_(2,k) = v_p;
		Xsig_pred_(3,k) = yaw_p;
		Xsig_pred_(4,k) = yawd_p;
	}
	///////////////// Predicting Mean Cova /////////////////
	// Predicted Mean state
	x_.fill(0.0);
	for(int k = 0; k < 2*n_aug_+1; k++){
		x_ = x_ + weights_(k)*Xsig_pred_.col(k);
	}
	//Predicted state Covariance matrix
	P_.fill(0.0);
	for(int i = 0; i < 2*n_aug_+1; i++){
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		x_diff(3) = atan2(sin(x_diff(3)),cos(x_diff(3)));
		P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
	}

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	use_radar_ = false;
	int n_z_ = 2;
	MatrixXd Z_sig = MatrixXd(n_z_, 2*n_aug_+1);


	for(int i = 0; i < 2*n_aug_+1; i++){
		double p_xl = Xsig_pred_(0,i);
		double p_yl = Xsig_pred_(1,i);
		
		Z_sig(0,i) = p_xl;
		Z_sig(1,i) = p_yl;
	}
	// Mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);
	z_pred.fill(0.0);
	for(int i = 0; i < 2*n_aug_+1; i++){
		z_pred = z_pred + weights_(i)*Z_sig.col(i);
	}
	// Measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_,n_z_);
	S.fill(0.0);
	for(int i = 0; i < 2*n_aug_+1;i++){
		VectorXd z_diff = Z_sig.col(i) - z_pred;

		//while (z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
		//while (z_diff(1)<M_PI) z_diff(1)+=2.*M_PI;

		S = S + weights_(i)*z_diff*z_diff.transpose();
	}

	MatrixXd R = MatrixXd(n_z_,n_z_);
	R << std_laspx_*std_laspx_, 0,
		 0, std_laspy_*std_laspy_;

	S += R;

	///////////////// UPDATING /////////////////

	MatrixXd Tc = MatrixXd(n_x_,n_z_);

	Tc.fill(0.0);
	for(int i = 0; i < 2*n_aug_+1; i++){
		
		VectorXd z_diff = Z_sig.col(i) - z_pred;
		z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		Tc += weights_(i)*x_diff*z_diff.transpose();
	}
	VectorXd z = VectorXd(n_z_);

	z = meas_package.raw_measurements_;

	VectorXd z_diff = z - z_pred;
	MatrixXd K = Tc*S.inverse();

	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	std::cout<<"x_ = "<<x_<<std::endl;
	std::cout<<"P_ = "<<P_<<std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	use_radar_ = false;
	int n_z_ = 3;
	MatrixXd Z_sig = MatrixXd(n_z_, 2*n_aug_+1);


	for(int i = 0; i < 2*n_aug_+1; i++){
		double p_x = Xsig_pred_(0,i);
		double p_y = Xsig_pred_(1,i);
		double v_ = Xsig_pred_(2,i);
		double yaw_ = Xsig_pred_(3,i);

		double v1 = cos(yaw_)*v_;
		double v2 = sin(yaw_)*v_;

		Z_sig(0,i) = sqrt(p_x*p_x + p_y*p_y);
		Z_sig(1,i) = atan2(p_y,p_x);
		Z_sig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
	}
	// Mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);
	z_pred.fill(0.0);
	for(int i = 0; i < 2*n_aug_+1; i++){
		z_pred = z_pred + weights_(i)*Z_sig.col(i);
	}
	z_pred(1) = atan2(sin(z_pred(1)),cos(z_pred(1)));
	// Measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_,n_z_);
	S.fill(0.0);
	for(int i = 0; i < 2*n_aug_+1;i++){
		VectorXd z_diff = Z_sig.col(i) - z_pred;

		z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));

		S = S + weights_(i)*z_diff*z_diff.transpose();
	}
	MatrixXd R = MatrixXd(n_z_,n_z_);
	R << std_radr_*std_radr_, 0, 0,
		 0, std_radphi_*std_radphi_, 0,
		 0, 0, std_radrd_*std_radrd_;
	S += R;

	///////////////// UPDATING /////////////////

	MatrixXd Tc = MatrixXd(n_x_,n_z_);

	Tc.fill(0.0);
	for(int i = 0; i < 2*n_aug_+1; i++){
		
		VectorXd z_diff = Z_sig.col(i) - z_pred;
		z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));


		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		x_diff(3) = atan2(sin(x_diff(3)),cos(x_diff(3)));

		Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
	}
	VectorXd z = VectorXd(n_z_);

	z = meas_package.raw_measurements_;

	VectorXd z_diff = z - z_pred;

	z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));

	MatrixXd K = Tc*S.inverse();

	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	std::cout<<"x_ = "<<x_<<std::endl;
	std::cout<<"P_ = "<<P_<<std::endl;
}
