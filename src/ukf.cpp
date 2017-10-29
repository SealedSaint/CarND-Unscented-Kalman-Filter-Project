#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

UKF::UKF() {
	is_initialized = false;
	use_laser = false;  // if this is false, laser measurements will be ignored (except during init)
	use_radar = true;  // if this is false, radar measurements will be ignored (except during init)

	n_x = 5;
	n_aug = 7;
	n_z_rad = 3;
	n_z_las = 2;

	x = VectorXd(n_x);  // initial state vector
	P = MatrixXd(n_x, n_x);  // initial covariance matrix

	lambda = 3 - n_aug;  // spreading parameter

	// Set the weights
	weights = VectorXd(2 * n_aug + 1);
	for(int i = 0; i < weights.size(); i++) {
		if(i == 0) weights(i) = lambda / double(lambda + n_aug);
		else weights(i) = 1 / double(2 * (lambda + n_aug));
	}

	// TODO: Adjust these
	std_a = 2;  // Process noise standard deviation longitudinal acceleration in m/s^2
	std_yawdd = M_PI / 8;  // Process noise standard deviation yaw acceleration in rad/s^2

	/**** Measurement noise values below should not be changed ****/
	std_las_px = 0.15;  // Laser measurement noise standard deviation position1 in m
	std_las_py = 0.15;  // Laser measurement noise standard deviation position2 in m
	std_rad_r = 0.3;  // Radar measurement noise standard deviation radius in m
	std_rad_phi = 0.03;  // Radar measurement noise standard deviation angle in rad
	std_rad_rd = 0.3;  // Radar measurement noise standard deviation radius change in m/s
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage measurement_package) {
	if(!is_initialized) {  // Set x and P with our first measurement
		cout << "Initializing..." << endl;
		if(measurement_package.sensor_type == MeasurementPackage::LASER) {
			cout << "First measurement is: LASER" << endl;
			x << measurement_package.raw_measurements(0),
				 measurement_package.raw_measurements(1),
				 0,
				 0,
				 0;
		}
		else if(measurement_package.sensor_type == MeasurementPackage::RADAR) {
			cout << "First measurement is: RADAR" << endl;
			double r = measurement_package.raw_measurements(0);
			double radians = measurement_package.raw_measurements(1);
			x << r * cos(radians),
				 r * sin(radians),
				 0,
				 0,
				 0;
		}

		P.setIdentity(n_x, n_x);

		// Set the timestamp
		last_time = measurement_package.timestamp;

		is_initialized = true;
		cout << "Initialization complete!" << endl;
		return;
	}

	// Is initialized, so do our normal predict -> update flow

	// Check for skips
	if(measurement_package.sensor_type == MeasurementPackage::LASER && !use_laser) return;
	else if(measurement_package.sensor_type == MeasurementPackage::RADAR && !use_radar) return;

	double Dt = (measurement_package.timestamp - last_time) / 1000000.0;
	cout << Dt << endl;
	last_time = measurement_package.timestamp;

	// PREDICT

	// cout << "Beginning Prediction..." << endl;
	MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
	VectorXd x_pred = VectorXd(n_x);
	MatrixXd P_pred = MatrixXd(n_x, n_x);
	Predict(Dt, &Xsig_pred, &x_pred, &P_pred);  // This sets the three values above

	// UPDATE

	// cout << "Beginning Update..." << endl;
	VectorXd z = measurement_package.raw_measurements;
	if(measurement_package.sensor_type == MeasurementPackage::LASER) {
		// cout << "Updating for: LASER" << endl;
		n_z = n_z_las;
	}
	else if(measurement_package.sensor_type == MeasurementPackage::RADAR) {
		// cout << "Updating for: RADAR" << endl;
		n_z = n_z_rad;
		UpdateRadar(Xsig_pred, x_pred, P_pred, z);
	}
}

// PREDICTION STEP

MatrixXd UKF::GenerateAugmentedSigmaPoints() {
	VectorXd x_aug = VectorXd(n_aug);  // augmented mean vector
	MatrixXd P_aug = MatrixXd(n_aug, n_aug);  // augmented state covariance
	MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);  // sigma point matrix

	// create augmented mean state
	x_aug.head(n_x) = x;
	x_aug(n_x) = 0;
	x_aug(n_x + 1) = 0;

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

	return Xsig_aug;
}

VectorXd UKF::ProcessModelPrediction(VectorXd x_sig_aug, double Dt) {
	double x = x_sig_aug(0);
	double y = x_sig_aug(1);
	double vel = x_sig_aug(2);
	double yaw = x_sig_aug(3);
	double yaw_rate = x_sig_aug(4);
	double nu_a = x_sig_aug(5);
	double nu_p = x_sig_aug(6);

	// Prediction = current (x_k) + change (x_p) + noise (x_n)

	VectorXd x_k = x_sig_aug.head(n_x);

	VectorXd x_p = VectorXd(n_x);
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
	VectorXd x_n = VectorXd(n_x);
	x_n <<  c1 * cos(yaw) * nu_a,
			c1 * sin(yaw) * nu_a,
			nu_a * Dt,
			c1 * nu_p,
			nu_p * Dt;

	VectorXd x_sig_pred = x_k + x_p + x_n;
	return x_sig_pred;
}

MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig_aug, double Dt) {
	MatrixXd Xsig_pred = MatrixXd(n_x, Xsig_aug.cols());

	// pass each augmented sigma point through the process model
	for(int i = 0; i < Xsig_aug.cols(); i++) {
		Xsig_pred.col(i) = ProcessModelPrediction(Xsig_aug.col(i), Dt);
	}

	return Xsig_pred;
}

void UKF::PredictMeanAndCovariance(MatrixXd Xsig_pred, VectorXd* x_pred_out, MatrixXd* P_pred_out) {
	VectorXd x_pred = VectorXd(n_x);  // predicted state vector
	MatrixXd P_pred = MatrixXd(n_x, n_x);  // predicted covariance matrix

	// create predicted state mean
	x_pred.fill(0);
	for(int i = 0; i < Xsig_pred.cols(); i++) {
		x_pred += weights(i) * Xsig_pred.col(i);
	}

	// predict state covariance matrix
	P_pred.fill(0);
	for(int i = 0; i < Xsig_pred.cols(); i++) {
		// state difference
		VectorXd x_diff = Xsig_pred.col(i) - x_pred;

		// angle normalization
		while(x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
		while(x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;

		P_pred += weights(i) * x_diff * x_diff.transpose();
	}

	*x_pred_out = x_pred;
	*P_pred_out = P_pred;
}

void UKF::Predict(double Dt, MatrixXd* Xsig_pred_out, VectorXd* x_pred_out, MatrixXd* P_pred_out) {
	// Generate (Augmented) Sigma point representation of previous state
	// cout << "Generating augmented sigma points..." << endl;
	MatrixXd Xsig_aug = GenerateAugmentedSigmaPoints();

	// Predict Sigma points (for time k+1) from augmented sigma points (for time k)
	// cout << "Predicting sigma points..." << endl;
	MatrixXd Xsig_pred = PredictSigmaPoints(Xsig_aug, Dt);
	*Xsig_pred_out = Xsig_pred;

	// Predict mean and covariance from predicted sigma points
	// cout << "Predicting mean and covariance..." << endl;
	VectorXd x_pred = VectorXd(n_x);
	MatrixXd P_pred = MatrixXd(n_x, n_x);
    PredictMeanAndCovariance(Xsig_pred, &x_pred, &P_pred);
	*x_pred_out = x_pred;
	*P_pred_out = P_pred;
}

// UPDATE LIDAR

void UKF::UpdateLidar(MeasurementPackage measurement_package) {
	/* TODO:
		* Use lidar data to update the belief about the object's position.
		* Modify the state vector, x_, and covariance, P_.
		* Calculate the lidar NIS.
	*/

	// Predict measurement


	// Update state

}

// UPDATE RADAR

VectorXd UKF::RadarMeasurementModelPrediction(VectorXd x_sig_pred) {
	float x = x_sig_pred(0);
	float y = x_sig_pred(1);
	float vel = x_sig_pred(2);
	float yaw = x_sig_pred(3);

	float c1 = sqrt(x * x + y * y);
	VectorXd z_pred = VectorXd(n_z);
	z_pred <<   c1,
				atan2(y, x),
				(x * vel * cos(yaw) + y * vel * sin(yaw)) / c1;
	return z_pred;
}

MatrixXd UKF::PredictRadarSigmaPoints(MatrixXd Xsig_pred) {
	MatrixXd Zsig_pred = MatrixXd(n_z, Xsig_pred.cols());

	for(int i = 0; i < Xsig_pred.cols(); i++) {
		Zsig_pred.col(i) = RadarMeasurementModelPrediction(Xsig_pred.col(i));
	}

	return Zsig_pred;
}

void UKF::PredictRadarMeanAndCovariance(MatrixXd Zsig_pred, VectorXd* z_pred_out, MatrixXd* S_pred_out) {
	VectorXd z_pred = VectorXd(n_z);  // measurement mean
	MatrixXd S_pred = MatrixXd(n_z, n_z);  // measurement covariance

	// calculate predicted measurement mean

	z_pred.fill(0);
	for(int i = 0; i < Zsig_pred.cols(); i++) {
		z_pred += weights(i) * Zsig_pred.col(i);
	}

	// calculate predicted measurement covariance

	// create measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R.fill(0);
	R(0, 0) = std_rad_r * std_rad_r;
	R(1, 1) = std_rad_phi * std_rad_phi;
	R(2, 2) = std_rad_rd * std_rad_rd;

	S_pred.fill(0);
	for(int i = 0; i < Zsig_pred.cols(); i++) {
		VectorXd z_diff = Zsig_pred.col(i) - z_pred;

		//normalize angle
		while(z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
		while(z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;

		S_pred += weights(i) * z_diff * z_diff.transpose();
	}
	S_pred += R;

	*z_pred_out = z_pred;
	*S_pred_out = S_pred;
}

void UKF::UpdateRadar(MatrixXd Xsig_pred, VectorXd x_pred, MatrixXd P_pred, VectorXd z) {
	/* TODO:
		* Calculate the radar NIS.
	*/

	// Predict measurement sigma points
	MatrixXd Zsig_pred = PredictRadarSigmaPoints(Xsig_pred);

	// Predict measurement mean and covariance
	VectorXd z_pred = VectorXd(n_z);
	MatrixXd S_pred = MatrixXd(n_z, n_z);
	PredictRadarMeanAndCovariance(Zsig_pred, &z_pred, &S_pred);

	// Update state

	// calculate cross correlation matrix
	MatrixXd Tc = MatrixXd(n_x, n_z);
	Tc.fill(0);
	for(int i = 0; i < Zsig_pred.cols(); i++) {
		VectorXd x_diff = Xsig_pred.col(i) - x_pred;
		// Normalize angle
		while(x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
		while(x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;

		VectorXd z_diff = Zsig_pred.col(i) - z_pred;
		// Normalize angle
		while(z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
		while(z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;

		Tc += weights(i) * x_diff * z_diff.transpose();
	}

	// calculate Kalman gain K
	MatrixXd K = Tc * S_pred.inverse();

	VectorXd z_diff = z - z_pred;
	// normalize angle
	while(z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
	while(z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;

	x = x_pred + K * z_diff;
	P = P_pred - K * S_pred * K.transpose();
}
