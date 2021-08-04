#include "kalman_filter.h"
#include <math.h>       /* -PI<atan2<+PI , M_PI (aka M_PI) , cos, sin */

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) 
  {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
  }

void KalmanFilter::Predict()
  {
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
  }

void KalmanFilter::Update(const VectorXd &z) 
  {
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
  
    //new estimate ( hence state update) 
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
  }

void KalmanFilter::UpdateEKF(const VectorXd &z)
  {
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
  
    float c1 = (px*px)+(py*py);
  
    float Rho = sqrt(c1);
    float Phi = M_PI/2;
    float Rho_dot=0;
  

    if (fabs(px) < 0.0001)
      {
      if (py < 0)
        {
         Phi = - M_PI/2 ; 
        }
      }
                
    else { Phi= atan2(py,px); }
  
  
    if (fabs(Rho) > 0.0001)
      { Rho_dot = ((px*vx) + (py*vy))/Rho; }
    else 
      { Rho_dot=0;}
      
    VectorXd z_pred(3);
    z_pred << Rho,Phi,Rho_dot;
    
    VectorXd y = z - z_pred;
  
    while ( y(1) < -M_PI ) 
      {
        y(1)=y(1)+(2*M_PI);
      }
      
    while ( y(1) >M_PI ) 
      {
        y(1)=y(1)-(2*M_PI);
      }
  
    MatrixXd Hj_ = H_;
    MatrixXd Hjt = Hj_.transpose();
    MatrixXd S = Hj_ * P_ * Hjt + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Hjt;
    MatrixXd K = PHt * Si;
  
    // new estimate (hence state update) 
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
  }
  


