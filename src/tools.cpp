#include <iostream>
#include "tools.h"

using namespace std;
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
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    float prod1 = (px*px) + (py*py);
    float prod2 = sqrt(prod1);
    float prod3 = (prod1*prod2);

    //check division by zero
    if(fabs(prod1) < 0.0001){
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        Hj << 0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0;
        return Hj;
    }

    //compute the Jacobian matrix
    // avoid mulitple calcucations
    float vx_py = vx*py;
    float vy_px = vy*px;

    Hj << (px/prod2), (py/prod2), 0, 0,
          -(py/prod1), (px/prod1), 0, 0,
            py*(vx_py - vy_px)/prod3, px*(vy_px - vx_py)/prod3, px/prod2, py/prod2;

    return Hj;
}

/**
 * get values between -pi and pi
 */
float Tools::wrapMinMax(float x, float min, float max)
{
    if (x < min){
        x += 2*max;
    }
    else if (x > max){
        x -= 2*max;
    }
    return x;
}
