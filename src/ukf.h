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
	long long last_time;

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
	int lambda;


	UKF();
	virtual ~UKF();

	/* ProcessMeasurement
	 * @param measurement_package - latest measurement data (either radar or laser)
	*/
	void ProcessMeasurement(MeasurementPackage measurement_package);

	// Generate augmented sigma point representation of current state (x and P)
	MatrixXd GenerateAugmentedSigmaPoints();

	/* Uses the process model to create predicted sigma points from augmented sigma points
	 * @param {VectorXd} x_sig_aug - augmented sigma points to be passed through the process model
	 * @param {double} Dt - elapsed time since last measurement
	 * return {VectorXd} - predicted sigma points calculated from process model
	*/
	VectorXd ProcessModelPrediction(VectorXd x_sig_aug, double Dt);

	/* Gives predicted sigma points (through process model) based on augmented sigma points
	 * @param {MatrixXd} Xsig_aug - augmented sigma points
	 * @param {double} Dt - elapsed time since last measurement
	 * return {MatrixXd} - predicted sigma points
	*/
	MatrixXd PredictSigmaPoints(MatrixXd Xsig_aug, double Dt);

	/* Predicts state mean and covariance from predicted sigma points (Xsig_pred)
	 * @param {MatrixXd} Xsig_pred - predicted sigma points
	 * @param {VectorXd*} x_pred_out - predicted state mean to be filled out
	 * @param {MatrixXd*} P_pred_out - predicted state covariance to be filled out
	*/
	void PredictMeanAndCovariance(MatrixXd Xsig_pred, VectorXd* x_pred_out, MatrixXd* P_pred_out);

	/* Predicts sigma points, state mean, and state covariance matrix
	 * @param Dt - Time between k and k+1 in seconds
	 * @param {MatrixXd*} Xsig_pred_out - predicted sigma points to be filled out
	 * @param {VectorXd*} x_pred_out - predicted state mean to be filled out
	 * @param {MatrixXd*} P_pred_out - predicted state covariance to be filled out
	*/
	void Predict(double Dt, MatrixXd* Xsig_pred_out, VectorXd* x_pred_out, MatrixXd* P_pred_out);

	/* Maps a predicted sigma points vector to Lidar measurement space
	 * @param {VectorXd} Xsig_pred - predicted sigma points to be passed through the Lidar measurement model
	 * return {VectorXd} - Lidar sigma points (sigma points in Lidar measurement space)
	*/
	VectorXd LidarMeasurementModelPrediction(VectorXd Xsig_pred);

	/* Gives predicted sigma points mapped to Lidar measurement space using Lidar measurement model
	 * @param {MatrixXd} Xsig_pred - predicted sigma points
	 * return {MatrixXd} - Lidar sigma points (sigma points in Lidar measurement space)
	*/
	MatrixXd PredictLidarSigmaPoints(MatrixXd Xsig_pred);

	/* Fills out predicted measurement mean and covariance (z_pred and S_pred) from the predicted measurement sigma points (Zsig_pred)
	 * @param {MatrixXd} Zsig_pred - predicted measurement sigma points
	 * @param {VectorXd*} z_pred - predicted measurement mean to be filled out
	 * @param {MatrixXd*} S_pred - predicted measurement covariance to be filled out
	*/
	void PredictLidarMeanAndCovariance(MatrixXd Zsig_pred, VectorXd* z_pred_out, MatrixXd* S_pred_out);

	/* Updates the state mean and covariance (x and P) using a radar measurement and predicted values for the elapsed time
	 * @param {MatrixXd} Xsig_pred - predicted state sigma points for the elapsed time
	 * @param {VectorXd} x_pred - predicted mean
	 * @param {MatrixXd} P_pred - predicted covariance
	 * @param {VectorXd} z - new radar measurement
	*/
	void UpdateLidar(MatrixXd Xsig_pred, VectorXd x_pred, MatrixXd P_pred, VectorXd z);

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

	/* Updates the state mean and covariance (x and P) using a radar measurement and predicted values for the elapsed time
	 * @param {MatrixXd} Xsig_pred - predicted state sigma points for the elapsed time
	 * @param {VectorXd} x_pred - predicted mean
	 * @param {MatrixXd} P_pred - predicted covariance
	 * @param {VectorXd} z - new radar measurement
	*/
	void UpdateRadar(MatrixXd Xsig_pred, VectorXd x_pred, MatrixXd P_pred, VectorXd z);

};

#endif /* UKF_H */
