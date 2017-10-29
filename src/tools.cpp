#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.empty()){
		cout << "No estimations provided. Unable to calculate RMSE." << endl;
		return rmse;
	}
	if(estimations.size() != ground_truth.size()){
		cout << "Estimations and ground truth samples not equal in size. Unable to calculate RMSE." << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(int i = 0; i < estimations.size(); ++i){
		VectorXd residual = estimations[i] - ground_truth[i];
		VectorXd residual_squared = residual.array() * residual.array();
		rmse += residual_squared;
	}

	rmse /= estimations.size();  // calculate the mean
	return rmse.array().sqrt();  // calculate the squared root
}
