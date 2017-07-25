#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // Check validity of ground truth data RWM
  if(estimations.size() != ground_truth.size()
  		|| estimations.size() == 0){
  	std::cout<< "Invalid: Ground truth data estimation incorrect" << endl;
  	return rmse;
  }

  // Accumulate squared residuals RWM
  for(unsigned int i=0; i < estimations.size(); i++){

  	VectorXd residual = estimations[i] - ground_truth[i];

  	// coefficient-wise multiplication RWM
  	residual = residual.array()*residual.array();
  	rmse += residual;
  }

  // Calculate mean RWM
  rmse = rmse/estimations.size();

  // Calculate the squared root RWM
  rmse = rmse.array().sqrt();
  //std::cout<<"RMSE = "<<rmse<<endl;
  return rmse;
}