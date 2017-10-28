#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

	// initially set to false, set to true in first call of ProcessMeasurement
	bool is_initialized;

	// if this is false, laser measurements will be ignored (except for init)
	bool use_laser;
	// if this is false, radar measurements will be ignored (except for init)
	bool use_radar;

	// state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
	VectorXd x;
	// state covariance matrix
	MatrixXd P;
	// augmented sigma points matrix
	MatrixXd Xsig_aug;
	// predicted sigma points matrix
	MatrixXd Xsig_pred;

	// time when the state is true, in us
	long long time_us;

	// Process noise standard deviation longitudinal acceleration in m/s^2
	double std_a;
	// Process noise standard deviation yaw acceleration in rad/s^2
	double std_yawdd;

	// Laser measurement noise standard deviation position1 in m
	double std_las_px;
	// Laser measurement noise standard deviation position2 in m
	double std_las_py;

	// Radar measurement noise standard deviation radius in m
	double std_rad_r;
	// Radar measurement noise standard deviation angle in rad
	double std_rad_phi;
	// Radar measurement noise standard deviation radius change in m/s
	double std_rad_rd;

	// Weights of sigma points
	VectorXd weights;

	// State dimension
	int n_x;
	// Augmented state dimension
	int n_aug;

	// Sigma point spreading parameter
	double lambda;


	UKF();
	virtual ~UKF();

	/* ProcessMeasurement
	 * @param meas_package The latest measurement data of either radar or laser
	*/
	void ProcessMeasurement(MeasurementPackage meas_package);

	// void GenerateSigmaPoints(MatrixXd* Xsig_out);
	void GenerateAugmentedSigmaPoints(MatrixXd* Xsig_out);

	/* Passes the augmented sigma point vector Xsig_aug through the process model
	 * @param {VectorXd} Xsig_aug - augmented sigma points to be passed through the process model
	 * @param {float} Dt - elapsed time since last measurement
	 * return {VectorXd} - predicted sigma points calculated from process model
	*/
	VectorXd ProcessModelPrediction(VectorXd Xsig_aug, float Dt);

	/* Gives predicted sigma points (through process model) based on augmented sigma points
	 * @param {MatrixXd} Xsig_aug - augmented sigma points
	 * @param {float} Dt - elapsed time since last measurement
	 * return {MatrixXd} - predicted sigma points
	*/
	MatrixXd PredictSigmaPoints(MatrixXd Xsig_aug, float Dt);

	/* Predicts sigma points, the state, and the state covariance matrix
	 * @param Dt - Time between k and k+1 in seconds
	*/
	void Prediction(double Dt);

	/* Updates the state and the state covariance matrix using a laser measurement
	 * @param meas_package The measurement at k+1
	*/
	void UpdateLidar(MeasurementPackage meas_package);

	/* Updates the state and the state covariance matrix using a radar measurement
	 * @param meas_package The measurement at k+1
	*/
	void UpdateRadar(MeasurementPackage meas_package);

};

#endif /* UKF_H */
