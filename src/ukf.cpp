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
	is_initialized = false;
	use_laser = true;  // if this is false, laser measurements will be ignored (except during init)
	use_radar = true;  // if this is false, radar measurements will be ignored (except during init)

	x = VectorXd(5);  // initial state vector
	P = MatrixXd(5, 5);  // initial covariance matrix

	// TODO: Adjust these
	std_a = 30;  // Process noise standard deviation longitudinal acceleration in m/s^2
	std_yawdd = 30;  // Process noise standard deviation yaw acceleration in rad/s^2

	// Measurement noise values below should not be changed

	std_las_px = 0.15;  // Laser measurement noise standard deviation position1 in m
	std_las_py = 0.15;  // Laser measurement noise standard deviation position2 in m

	std_rad_r = 0.3;  // Radar measurement noise standard deviation radius in m
	std_rad_phi = 0.03;  // Radar measurement noise standard deviation angle in rad
	std_rad_rd = 0.3;  // Radar measurement noise standard deviation radius change in m/s

	/**
	TODO:
	Complete the initialization. See ukf.h for other member properties.
	Hint: one or more values initialized above might be wildly off...
	*/
}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} meas_package The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
	TODO:
	Complete this function! Make sure you switch between lidar and radar measurements.
	*/

	is_initialized = true;
}

// void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
// 	// create sigma point matrix
// 	MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);
//
// 	// calculate square root of P
// 	MatrixXd A = P.llt().matrixL();
// 	float c = sqrt(lambda + n_x);
// 	MatrixXd M = c*A;
//
// 	// first row is just x
// 	Xsig.col(0) = x;
// 	// set sigma points as columns of matrix Xsig
// 	for(int i = 1; i < n_x + 1; i++) {
// 		Xsig.col(i) = x + M.col(i - 1);
// 	}
// 	for(int i = n_x + 1; i < 2 * n_x + 1; i++) {
// 		Xsig.col(i) = x - M.col(i - (n_x + 1));
// 	}
//
// 	*Xsig_out = Xsig;
// }

void UKF::GenerateAugmentedSigmaPoints(MatrixXd* Xsig_out) {
	VectorXd x_aug = VectorXd(7);  // augmented mean vector
	MatrixXd P_aug = MatrixXd(7, 7);  // augmented state covariance
	MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);  // sigma point matrix

	// create augmented mean state
	x_aug.head(n_x) = x;
	x_aug(n_x) = 0;
	x_aug(n_x+1) = 0;

	// create augmented covariance matrix
	P_aug.fill(0);
	P_aug.topLeftCorner(n_x, n_x) = P;

	MatrixXd Q = MatrixXd(2, 2);
	Q << std_a * std_a, 0,
		 0, std_yawdd * std_yawdd;
	P_aug.bottomRightCorner(2, 2) = Q;

	// create square root matrix
	MatrixXd M = P_aug.llt().matrixL();
	M *= sqrt(lambda + n_aug);

	// create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for(int i = 0; i < n_aug; i++) {
		Xsig_aug.col(i + 1) = x_aug + M.col(i);
		Xsig_aug.col(n_aug + i + 1) = x_aug - M.col(i);
	}

	*Xsig_out = Xsig_aug;
}

VectorXd UKF::ProcessModelPrediction(VectorXd Xsig_aug, float Dt) {
	float x = Xsig_aug(0);
	float y = Xsig_aug(1);
	float vel = Xsig_aug(2);
	float yaw = Xsig_aug(3);
	float yaw_rate = Xsig_aug(4);
	float nu_a = Xsig_aug(5);
	float nu_p = Xsig_aug(6);

	// Prediction = current (x_k) + change (x_p) + noise (x_n)

	VectorXd x_k = Xsig_aug.head(5);

	VectorXd x_p = VectorXd(5);
	if(abs(yaw_rate) < 0.001) {  // Avoid division by zero
		x_p <<  vel * cos(yaw) * Dt,
				vel * sin(yaw) * Dt,
				0,
				0,
				0;
	}
	else {
		x_p <<  (vel / yaw_rate) * (sin(yaw + yaw_rate * Dt) - sin(yaw)),
				(vel / yaw_rate) * (-cos(yaw + yaw_rate * Dt) + cos(yaw)),
				0,
				yaw_rate * Dt,
				0;
	}

	float c1 = Dt * Dt / 2;
	VectorXd x_n = VectorXd(5);
	x_n <<  c1 * cos(yaw) * nu_a,
			c1 * sin(yaw) * nu_a,
			nu_a * Dt,
			c1 * nu_p,
			nu_p * Dt;

	VectorXd Xsig_pred = x_k + x_p + x_n;
	return Xsig_pred;
}

MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig_aug, float Dt) {
	// pass each augmented sigma point through the process model
	for(int i = 0; i < Xsig_aug.cols(); i++) {
		Xsig_pred.col(i) = ProcessModelPrediction(Xsig_aug.col(i), Dt);
	}

	*Xsig_out = Xsig_pred;
}

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} - Dt the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Prediction(double Dt) {
	/**
	TODO:
	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/

	// Generate (Augmented) Sigma points


	// Predict Sigma points


	// Predict mean and covariance
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

	// Predict measurement


	// Update state

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

	// Predict measurement


	// Update state
}
