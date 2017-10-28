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
	// Laser measurement dimension
	int n_z_las;
	// Radar measurement dimension
	int n_z_rad;
	// Measurement dimension (will change based on measurement received)
	int n_z;

	// Sigma point spreading parameter
	double lambda;


	UKF();
	virtual ~UKF();

	/* ProcessMeasurement
	 * @param meas_package The latest measurement data of either radar or laser
	*/
	void ProcessMeasurement(MeasurementPackage meas_package);

	// Generate augmented sigma point representation of current state (x and P)
	MatrixXd GenerateAugmentedSigmaPoints();

	/* Uses the process model to create predicted sigma points from augmented sigma points
	 * @param {VectorXd} x_sig_aug - augmented sigma points to be passed through the process model
	 * @param {float} Dt - elapsed time since last measurement
	 * return {VectorXd} - predicted sigma points calculated from process model
	*/
	VectorXd ProcessModelPrediction(VectorXd x_sig_aug, float Dt);

	/* Gives predicted sigma points (through process model) based on augmented sigma points
	 * @param {MatrixXd} Xsig_aug - augmented sigma points
	 * @param {float} Dt - elapsed time since last measurement
	 * return {MatrixXd} - predicted sigma points
	*/
	MatrixXd PredictSigmaPoints(MatrixXd Xsig_aug, float Dt);

	/* Fills out predicted state and covariance (x_pred and P_pred) from the predicted sigma points (Xsig_pred)
	 * @param {MatrixXd} Xsig_pred - predicted sigma points
	 * @param {VectorXd*} x_pred - predicted state mean to be filled out
	 * @param {MatrixXd*} P_pred - predicted covariance to be filled out
	*/
	void PredictMeanAndCovariance(MatrixXd Xsig_pred, VectorXd* x_pred_out, MatrixXd* P_pred_out);

	/* Predicts sigma points, the state, and the state covariance matrix
	 * @param Dt - Time between k and k+1 in seconds
	*/
	void Prediction(double Dt);

	/* Updates the state and the state covariance matrix using a laser measurement
	 * @param meas_package The measurement at k+1
	*/
	void UpdateLidar(MeasurementPackage meas_package);

	/* Maps a predicted sigma points vector to Radar measurement space
	 * @param {VectorXd} Xsig_pred - predicted sigma points to be passed through the Radar measurement model
	 * return {VectorXd} - Radar sigma points (sigma points in Radar measurement space)
	*/
	VectorXd RadarMeasurementModelPrediction(VectorXd Xsig_pred);

	/* Gives predicted sigma points mapped to Radar measurement space using Radar measurement model
	 * @param {MatrixXd} Xsig_pred - predicted sigma points
	 * return {MatrixXd} - Radar sigma points (sigma points in Radar measurement space)
	*/
	MatrixXd PredictRadarSigmaPoints(MatrixXd Xsig_pred);

	/* Fills out predicted measurement mean and covariance (z_pred and S_pred) from the predicted measurement sigma points (Zsig_pred)
	 * @param {MatrixXd} Zsig_pred - predicted measurement sigma points
	 * @param {VectorXd*} z_pred - predicted measurement mean to be filled out
	 * @param {MatrixXd*} S_pred - predicted measurement covariance to be filled out
	*/
	void PredictRadarMeanAndCovariance(MatrixXd Zsig_pred, VectorXd* z_pred_out, MatrixXd* S_pred_out);

	/* Updates the state and the state covariance matrix using a radar measurement
	 * @param meas_package The measurement at k+1
	*/
	void UpdateRadar(MeasurementPackage meas_package);

};

#endif /* UKF_H */
