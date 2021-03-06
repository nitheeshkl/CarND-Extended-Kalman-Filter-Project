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
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  for(unsigned int i=0; i<estimations.size(); ++i) {
      VectorXd diff = estimations[i] - ground_truth[i];
      diff = diff.array() * diff.array();
      rmse += diff;
  }

  rmse = rmse/estimations.size();
  rmse = sqrt(rmse.array());

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c1 = px*px + py*py;
  // fix divide by zero
  if (c1 < 0.00001) {
      px += 0.001;
      py += 0.001;
      c1 = px*px + py*py;
  }
  float c2 = sqrt(c1);
  float c3 = c1*c2;

  Hj << px/c2, py/c2, 0, 0,
        -py/c1, px/c1, 0, 0,
        (py*(vx*py - vy*px))/c3, (px*(vy*px - vx*py))/c3, px/c2, py/c2;

  return Hj;
}

double Tools::NormalizeAngle(double& angle) {
  return atan2(sin(angle), cos(angle));
}
