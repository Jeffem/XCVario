/*
*) Refactor the code to remove reduentent part. 
*) Compiled for Linux with C++14 standard
Copyright (c) 2021 Balamurugan Kandan.
MIT License; See LICENSE.md for complete details
Author: 2021 Balamurugan Kandan
*/

/*
Updated to be a class, use Eigen, and compile as an Arduino library.
Added methods to get gyro and accel bias. Added initialization to
estimated angles rather than assuming IMU is level.

Copyright (c) 2016 - 2019 Regents of the University of Minnesota and Bolder Flight Systems Inc.
MIT License; See LICENSE.md for complete details
Author: Brian Taylor
*/

/*
Addapted from earlier version
Copyright 2011 Regents of the University of Minnesota. All rights reserved.
Original Author: Adhika Lie
*/

/* Commentaire GFM
 * Il faudra revoir la gestion des thread avec les contrôles Mutex
 * Il faudra tester la nouvelle fonction Matrice inverse
 */

#include "ekfNavINS.h"
#include "stdlib.h"
#include <math.h>
#include <utility>
   using namespace std;


   double determinant(Matrix Gl,int m);
   Matrix inverse(Matrix Gl,int m);
   double  *ptr;//to access the Cofactor matrix


   double determinant(Matrix Gl,int m)

    {Matrix new_mat;double  Cofactor[6][6];
     double pivot[6][6];double det=0;
     if (m==2)
      {det=0;
       det=(Gl(0,0)*Gl(1,1))-(Gl(0,1)*Gl(1,0));
       return det;
      }//end if(m==2)
     else //if m!=2
      {for (int row=0;row<=(m-1);row++)
       {for (int col=0;col<=(m-1);col++)
        { pivot[row][col]=Gl(row,col);
          if(row==0)
           {
            int r=0;
            for(int i=row+1;i<=(m-1);i++)
             {int s=0;
              if(col==0)
               {for(int j=col+1;j<=(m-1);j++)
                 {new_mat(r,s)=Gl(i,j);
                  s++;
                 }
               }
              else
               {for(int j=0;j<=(m-1);j++)
                 {if(j==col)
                   break;
                  else
                   {new_mat(r,s)=Gl(i,j);
                    s++;
                   }
                 }
                }
               r++;
              }
             }
          if(row!=0)
           {int r=0;
            for (int i=0;i<=(m-1);i++)
             {
              if (i!=row)
               {int s=0;
               if(col==0)
                {for(int j=col+1;j<=(m-1);j++)
                  {new_mat(r,s)=Gl(i,j);
                   s++;
                  }
                }
               if(col!=0)
                {for(int j=0;j<=(m-1);j++)
                  {
                   if(j!=col)
                    {new_mat(r,s)=Gl(i,j);
                     s++;
                    }
                   }
                 }
                r++;
               }
              }
             }
         Cofactor[row][col]=pow(-1,row+col+2)*determinant(new_mat,m-1);
         ptr=&Cofactor[row][0];
          }

      double det=0;
      for(int coll=0;coll<=(m-1);coll++)
       {
        det=det+Cofactor[row][coll]*pivot[row][coll];
       }
     }
    }// end of else
   return det;
   }//end determinant function

   Matrix inverse(Matrix Gl,int m)
   {double adj[6][6];Matrix inv;double cof[6][6];//double Cofactor[6][6];
    double det=determinant(Gl,m);
      int index=0;
      for (int row=0;row<=(m-1);row++)
       {for (int col=0;col<=(m-1);col++)
         {
          cof[row][col]=*(ptr+index);
          //transpose of cofactor matrix for adjoint matrix
          adj[row][col]=cof[col][row];
          index++;
         }
        }
      for (int row=0;row<=(m-1);row++)//now using adjoint values for inverse
       {for (int col=0;col<=(m-1);col++)
         {
          inv(row,col)=adj[row][col]*1/det;
         }
       }
      return inv;
   }//end inverse function

void ekfNavINS::ekf_init(uint64_t time, double vn,double ve,double vd,double lat,double lon,double alt,float p,float q,float r,float ax,float ay,float az,float hx,float hy, float hz) {
  // grab initial gyro values for biases
  gbx = p;
  gby = q;
  gbz = r;
  std::tie(theta,phi,psi) = getPitchRollYaw(ax, ay, az, hx, hy, hz);
  // euler to quaternion
  quat = toQuaternion(phi, theta, psi);
  // Assemble the matrices
  // ... gravity
  grav(2,0) = G;
  // ... H
  //H.block(0,0,5,5) = Eigen::Matrix<float,5,5>::Identity();
  H(0,0)=1;H(1,1)=1;H(2,2)=1;H(3,3)=1;H(4,4)=1;
  // ... Rw
  //Rw.block(0,0,3,3) = powf(SIG_W_A,2.0f) * Eigen::Matrix<float,3,3>::Identity();
  Rw(0,0)=SIG_W_A*SIG_W_A;Rw(1,1)=SIG_W_A*SIG_W_A;Rw(2,2)=SIG_W_A*SIG_W_A;
  //Rw.block(3,3,3,3) = powf(SIG_W_G,2.0f) * Eigen::Matrix<float,3,3>::Identity();
  Rw(3,3)=SIG_W_G*SIG_W_G;Rw(4,4)=SIG_W_G*SIG_W_G;Rw(5,5)=SIG_W_G*SIG_W_G;
  //Rw.block(6,6,3,3) = 2.0f * powf(SIG_A_D,2.0f) / TAU_A*Eigen::Matrix<float,3,3>::Identity();
  Rw(6,6)=SIG_A_D*SIG_A_D;Rw(7,7)=SIG_A_D*SIG_A_D;Rw(8,8)=SIG_A_D*SIG_A_D;
  //Rw.block(9,9,3,3) = 2.0f * powf(SIG_G_D,2.0f) / TAU_G*Eigen::Matrix<float,3,3>::Identity();
  Rw(9,9)=SIG_G_D*SIG_G_D;Rw(10,10)=SIG_G_D*SIG_G_D;Rw(11,11)=SIG_G_D*SIG_G_D;

  // ... P
  //P.block(0,0,3,3) = powf(P_P_INIT,2.0f) * Eigen::Matrix<float,3,3>::Identity();
  P(0,0)=P_P_INIT*P_P_INIT;P(1,1)=P_P_INIT*P_P_INIT;P(2,2)=P_P_INIT*P_P_INIT;
  //P.block(3,3,3,3) = powf(P_V_INIT,2.0f) * Eigen::Matrix<float,3,3>::Identity();
  P(3,3)=P_V_INIT*P_V_INIT;P(4,4)=P_V_INIT*P_V_INIT;P(5,5)=P_V_INIT*P_V_INIT;
  //P.block(6,6,2,2) = powf(P_A_INIT,2.0f) * Eigen::Matrix<float,2,2>::Identity();
  P(4,4)=P_A_INIT*P_A_INIT;P(7,7)=P_A_INIT*P_A_INIT;
  //P(8,8) = powf(P_HDG_INIT,2.0f);
  P(8,8)=P_HDG_INIT*P_HDG_INIT;
  //P.block(9,9,3,3) = powf(P_AB_INIT,2.0f) * Eigen::Matrix<float,3,3>::Identity();
  P(9,9)=P_AB_INIT*P_AB_INIT;P(10,10)=P_AB_INIT*P_AB_INIT;P(11,11)=P_AB_INIT*P_AB_INIT;
  //P.block(12,12,3,3) = powf(P_GB_INIT,2.0f) * Eigen::Matrix<float,3,3>::Identity();
  P(12,12)=P_GB_INIT*P_GB_INIT;P(13,13)=P_GB_INIT*P_GB_INIT;P(14,14)=P_GB_INIT*P_GB_INIT;
  // ... R
  //R.block(0,0,2,2) = powf(SIG_GPS_P_NE,2.0f) * Eigen::Matrix<float,2,2>::Identity();
  R(0,0)=SIG_GPS_P_NE*SIG_GPS_P_NE;R(1,1)=SIG_GPS_P_NE*SIG_GPS_P_NE;
  //R(2,2) = powf(SIG_GPS_P_D,2.0f);
  R(2,2)=SIG_GPS_P_D*SIG_GPS_P_D;
  //R.block(3,3,2,2) = powf(SIG_GPS_V_NE,2.0f) * Eigen::Matrix<float,2,2>::Identity();
  R(3,3) = SIG_GPS_V_NE*SIG_GPS_V_NE;R(4,4) = SIG_GPS_V_NE*SIG_GPS_V_NE;
  //R(5,5) = powf(SIG_GPS_V_D,2.0f);
  R(5,5) = SIG_GPS_V_D*SIG_GPS_V_D;
  // .. then initialize states with GPS Data
  lat_ins = lat;
  lon_ins = lon;
  alt_ins = alt;
  vn_ins = vn;
  ve_ins = ve;
  vd_ins = vd;
  // specific force
  f_b(0,0) = ax;
  f_b(1,0) = ay;
  f_b(2,0) = az;
  /* initialize the time */
  _tprev = time;
}

void ekfNavINS::ekf_update( uint64_t time/*, unsigned long TOW*/, double vn,double ve,double vd,double lat,double lon,double alt,
                          float p,float q,float r,float ax,float ay,float az,float hx,float hy, float hz ) {
  if (!initialized_) {
    ekf_init(time, vn, ve, vd, lat, lon, alt, p, q, r, ax, ay, az, hx, hy, hz);
    // initialized flag
    initialized_ = true;
  } else {
    // get the change in time
    float _dt = ((float)(time - _tprev)) / 1e6;
    // Update Gyro and Accelerometer biases
    updateBias(ax, ay, az, p, q, r);
    // Update INS values
    updateINS();
    // Attitude Update
    dq(0) = 1.0f;
    dq(1) = 0.5f*om_ib(0,0)*_dt;
    dq(2) = 0.5f*om_ib(1,0)*_dt;
    dq(3) = 0.5f*om_ib(2,0)*_dt;
    quat = qmult(quat,dq);
    double normeq = quat.normFro();
    quat=quat/normeq;
    // Avoid quaternion flips sign
    if (quat(0) < 0) {
      quat = -1.0f*quat;
    }
    // AHRS Transformations
    C_N2B = quat2dcm(quat);
    C_B2N = transpose(C_N2B);
    // obtain euler angles from quaternion
    std::tie(phi, theta, psi) = toEulerAngles(quat);
    // Velocity Update
    dx = C_B2N*f_b + grav;
    vn_ins += _dt*dx(0,0);
    ve_ins += _dt*dx(1,0);
    vd_ins += _dt*dx(2,0);
    // Position Update
    dxd = llarate(V_ins,lla_ins);
    lat_ins += _dt*dxd(0,0);
    lon_ins += _dt*dxd(1,0);
    alt_ins += _dt*dxd(2,0);
    // Jacobian update
    updateJacobianMatrix();
    // Update process noise and covariance time
    updateProcessNoiseCovarianceTime(_dt);
    // Gps measurement update
    //if ((TOW - previousTOW) > 0) {
    if ((time - _tprev) > 0) {
      //previousTOW = TOW;
      lla_gps(0,0) = lat;
      lla_gps(1,0) = lon;
      lla_gps(2,0) = alt;
      V_gps(0,0) = vn;
      V_gps(1,0) = ve;
      V_gps(2,0) = vd;
      // Update INS values
      updateINS();
      // Create measurement Y
      updateCalculatedVsPredicted();
      // Kalman gain
      K = inverse(transpose(P*H)*(transpose(H*P*H) + R),6);
      // Covariance update
      P =transpose( (I15-K*H)*P*(I15-K*H)) + transpose(K*R*K);
      // State update
      x = K*y;
      // Update the results
      update15statesAfterKF();
      _tprev = time;
    }
    // Get the new Specific forces and Rotation Rate
    updateBias(ax, ay, az, p, q, r);
  }
}

void ekfNavINS::ekf_update(uint64_t time) {
//  std::shared_lock lock(shMutex);
  ekf_update(time, /*0,*/ gpsVel.vN, gpsVel.vE, gpsVel.vD,
                      gpsCoor.lat, gpsCoor.lon, gpsCoor.alt,
                      imuDat.gyroX, imuDat.gyroY, imuDat.gyroZ,
                      imuDat.accX, imuDat.accY, imuDat.accZ,
                      imuDat.hX, imuDat.hY, imuDat.hZ);
}

void ekfNavINS::imuUpdateEKF(uint64_t time, imuData imu) {
  {
//    std::unique_lock lock(shMutex);
    imuDat = imu;
  }
  ekf_update(time);
}

void ekfNavINS::gpsCoordinateUpdateEKF(gpsCoordinate coor) {
//  std::unique_lock lock(shMutex);
  gpsCoor = coor;
}

void ekfNavINS::gpsVelocityUpdateEKF(gpsVelocity vel) {
//  std::unique_lock lock(shMutex);
  gpsVel = vel;
}

void ekfNavINS::updateINS() {
  // Update lat, lng, alt, velocity INS values to matrix
  lla_ins(0,0) = lat_ins;
  lla_ins(1,0) = lon_ins;
  lla_ins(2,0) = alt_ins;
  V_ins(0,0) = vn_ins;
  V_ins(1,0) = ve_ins;
  V_ins(2,0) = vd_ins;
}

std::tuple<float,float,float> ekfNavINS::getPitchRollYaw(float ax, float ay, float az, float hx, float hy, float hz) {
  // initial attitude and heading
  theta = asinf(ax/G);
  phi = -asinf(ay/(G*cosf(theta)));
  // magnetic heading correction due to roll and pitch angle
  Bxc = hx*cosf(theta) + (hy*sinf(phi) + hz*cosf(phi))*sinf(theta);
  Byc = hy*cosf(phi) - hz*sinf(phi);
  // finding initial heading
  psi = -atan2f(Byc,Bxc);
  return (std::make_tuple(theta,phi,psi));
}

void ekfNavINS::updateCalculatedVsPredicted() {
      // Position, converted to NED
      pos_ecef_ins = lla2ecef(lla_ins);
      pos_ecef_gps = lla2ecef(lla_gps);
      pos_ned_gps = ecef2ned(pos_ecef_gps - pos_ecef_ins, lla_ins);
      // Update the difference between calculated and predicted
      y(0,0) = (float)(pos_ned_gps(0,0));
      y(1,0) = (float)(pos_ned_gps(1,0));
      y(2,0) = (float)(pos_ned_gps(2,0));
      y(3,0) = (float)(V_gps(0,0) - V_ins(0,0));
      y(4,0) = (float)(V_gps(1,0) - V_ins(1,0));
      y(5,0) = (float)(V_gps(2,0) - V_ins(2,0));
}

void ekfNavINS::update15statesAfterKF() {
	  Matrix temp = Matrix::zeros(3,1);
	  temp(0,0)=x(0,0);temp(1,0)=x(1,0);temp(2,0)=x(2,0);
      //estmimated_ins = llarate ((x.block(0,0,3,1)).cast<double>(), lat_ins, alt_ins);
	  estmimated_ins = llarate (temp, lat_ins, alt_ins);lat_ins += estmimated_ins(0,0);
      lon_ins += estmimated_ins(1,0);
      alt_ins += estmimated_ins(2,0);
      vn_ins = vn_ins + x(3,0);
      ve_ins = ve_ins + x(4,0);
      vd_ins = vd_ins + x(5,0);
      // Attitude correction
      dq(0,0) = 1.0f;
      dq(1,0) = x(6,0);
      dq(2,0) = x(7,0);
      dq(3,0) = x(8,0);
      quat = qmult(quat,dq);
      double normeq = quat.normFro();
      quat=quat/normeq;
      // obtain euler angles from quaternion
      std::tie(phi, theta, psi) = toEulerAngles(quat);
      abx = abx + x(9,0);
      aby = aby + x(10,0);
      abz = abz + x(11,0);
      gbx = gbx + x(12,0);
      gby = gby + x(13,0);
      gbz = gbz + x(14,0);
}

void ekfNavINS::updateBias(float ax,float ay,float az,float p,float q, float r) {
  f_b(0,0) = ax - abx;
  f_b(1,0) = ay - aby;
  f_b(2,0) = az - abz;
  om_ib(0,0) = p - gbx;
  om_ib(1,0) = q - gby;
  om_ib(2,0) = r - gbz;
}

void ekfNavINS::updateProcessNoiseCovarianceTime(float _dt) {
  PHI = I15+Fs*_dt;
  // Process Noise
  Gs=Matrix::zeros(15,12);
  //Gs.block(3,0,3,3) = -C_B2N;
  Gs(3,0)=-C_B2N(0,0);Gs(3,1)=-C_B2N(0,1);Gs(3,2)=-C_B2N(0,2);
  Gs(4,0)=-C_B2N(1,0);Gs(4,1)=-C_B2N(1,1);Gs(4,2)=-C_B2N(1,2);
  Gs(5,0)=-C_B2N(2,0);Gs(5,1)=-C_B2N(2,1);Gs(5,2)=-C_B2N(2,2);
  //Gs.block(6,3,3,3) = -0.5f*I3;
  Gs(6,3)=-0.5f;Gs(7,4)=-0.5f;Gs(8,5)=-0.5f;
  //Gs.block(9,6,6,6) = I6;
  Gs(9,6)=1;Gs(10,7)=1;Gs(11,8)=1;Gs(12,9)=1;Gs(13,10)=1;Gs(14,11)=1;
  // Discrete Process Noise
  Q = transpose(PHI*_dt*Gs*Rw*Gs);
  Q = 0.5f*(Q+transpose(Q));
  // Covariance Time Update
  P = transpose(PHI*P*PHI)+Q;
  P = 0.5f*(P+transpose(P));
}

void ekfNavINS::updateJacobianMatrix() {
    // Jacobian
  Fs = Matrix::zeros(15,15);
  // ... pos2gs
  //Fs.block(0,3,3,3) = I3;
  Fs(0,0)=1;Fs(1,1)=1;Fs(2,2)=1;
  // ... gs2pos
  Fs(5,2) = -2.0f*G/EARTH_RADIUS;
  // ... gs2att
  //Fs.block(3,6,3,3) = -2.0f*C_B2N*sk(f_b);
  Fs(3,6)=-2.0f*C_B2N(0,0)*sk(f_b)(0,0);Fs(3,7)=-2.0f*C_B2N(0,1)*sk(f_b)(0,1);Fs(3,8)=-2.0f*C_B2N(0,2)*sk(f_b)(0,2);
  Fs(4,6)=-2.0f*C_B2N(1,0)*sk(f_b)(1,0);Fs(4,7)=-2.0f*C_B2N(1,1)*sk(f_b)(1,1);Fs(4,8)=-2.0f*C_B2N(1,2)*sk(f_b)(1,2);
  Fs(5,6)=-2.0f*C_B2N(2,0)*sk(f_b)(2,0);Fs(5,7)=-2.0f*C_B2N(2,1)*sk(f_b)(2,1);Fs(5,8)=-2.0f*C_B2N(2,2)*sk(f_b)(2,2);
  // ... gs2acc
  //Fs.block(3,9,3,3) = -C_B2N;
  Fs(3,9)=-C_B2N(0,0);Fs(3,10)=-C_B2N(0,1);Fs(3,11)=-C_B2N(0,2);
  Fs(4,9)=-C_B2N(1,0);Fs(4,10)=-C_B2N(1,1);Fs(4,11)=-C_B2N(1,2);
  Fs(5,9)=-C_B2N(2,0);Fs(5,10)=-C_B2N(2,1);Fs(5,11)=-C_B2N(2,2);
  // ... att2att
  //Fs.block(6,6,3,3) = -sk(om_ib);
  Fs(6,6)=-sk(om_ib)(0,0);Fs(6,7)=-sk(om_ib)(0,1);Fs(6,8)=-sk(om_ib)(0,2);
  Fs(7,6)=-sk(om_ib)(1,0);Fs(7,7)=-sk(om_ib)(1,1);Fs(7,8)=-sk(om_ib)(1,2);
  Fs(8,6)=-sk(om_ib)(2,0);Fs(8,7)=-sk(om_ib)(2,1);Fs(8,8)=-sk(om_ib)(2,2);

  // ... att2gyr
  //Fs.block(6,12,3,3) = -0.5f*I3;
  Fs(6,12)=-0.5f;Fs(7,13)=-0.5f;Fs(7,13)=-0.5f;
  // ... Accel Markov Bias
  //Fs.block(9,9,3,3) = -1.0f/TAU_A*I3;
  Fs(9,9)=-1.0f/TAU_A;Fs(10,10)=-1.0f/TAU_A;Fs(11,11)=-1.0f/TAU_A;
  //Fs.block(12,12,3,3) = -1.0f/TAU_G*I3;
  Fs(12,2)=-1.0f/TAU_G;Fs(13,13)=-1.0f/TAU_G;Fs(14,14)=-1.0f/TAU_G;
}

// This function gives a skew symmetric matrix from a given vector w
Matrix ekfNavINS::sk(Matrix w) {
  Matrix C(3,3);
  C(0,0) = 0.0f;    C(0,1) = -w(2,0); C(0,2) = w(1,0);
  C(1,0) = w(2,0);  C(1,1) = 0.0f;    C(1,2) = -w(0,0);
  C(2,0) = -w(1,0); C(2,1) = w(0,0);  C(2,2) = 0.0f;
  return C;
}
/*constexpr std::pair<double, double> ekfNavINS::earthradius(double lat) {
  double denom = fabs(1.0 - (ECC2 * pow(sin(lat),2.0)));
  double Rew = EARTH_RADIUS / sqrt(denom);
  double Rns = EARTH_RADIUS * (1.0-ECC2) / (denom*sqrt(denom));
  return (std::make_pair(Rew, Rns));
}
*/
// This function calculates the rate of change of latitude, longitude, and altitude.
Matrix ekfNavINS::llarate(Matrix V,Matrix lla) {
  //double Rew, Rns;
  Matrix lla_dot;
  //std::tie(Rew, Rns) = earthradius(lla(0,0));
  lla_dot(0,0) = V(0,0)/(Rns + lla(2,0));
  lla_dot(1,0) = V(1,0)/((Rew + lla(2,0))*cos(lla(0,0)));
  lla_dot(2,0) = -V(2,0);
  return lla_dot;
}

// This function calculates the rate of change of latitude, longitude, and altitude.
Matrix ekfNavINS::llarate(Matrix V, double lat, double alt) {
  Matrix lla;
  lla(0,0) = lat;
  lla(1,0) = 0.0; // Not used
  lla(2,0) = alt;
  return llarate(V, lla);
}

// This function calculates the ECEF Coordinate given the Latitude, Longitude and Altitude.
Matrix ekfNavINS::lla2ecef(Matrix lla) {
  //double Rew;
  Matrix ecef;
  //std::tie(Rew, std::ignore) = earthradius(lla(0,0));
  ecef(0,0) = (Rew + lla(2,0)) * cos(lla(0,0)) * cos(lla(1,0));
  ecef(1,0) = (Rew + lla(2,0)) * cos(lla(0,0)) * sin(lla(1,0));
  ecef(2,0) = (Rew * (1.0 - ECC2) + lla(2,0)) * sin(lla(0,0));
  return ecef;
}

// This function converts a vector in ecef to ned coordinate centered at pos_ref.
Matrix ekfNavINS::ecef2ned(Matrix ecef,Matrix pos_ref) {
  Matrix ned;
  ned(1,0)=-sin(pos_ref(1,0))*ecef(0,0) + cos(pos_ref(1,0))*ecef(1,0);
  ned(0,0)=-sin(pos_ref(0,0))*cos(pos_ref(1,0))*ecef(0,0)-sin(pos_ref(0,0))*sin(pos_ref(1,0))*ecef(1,0)+cos(pos_ref(0,0))*ecef(2,0);
  ned(2,0)=-cos(pos_ref(0,0))*cos(pos_ref(1,0))*ecef(0,0)-cos(pos_ref(0,0))*sin(pos_ref(1,0))*ecef(1,0)-sin(pos_ref(0,0))*ecef(2,0);
  return ned;
}

// quaternion to dcm
Matrix ekfNavINS::quat2dcm(Matrix q) {
  Matrix C_N2B;
  C_N2B(0,0) = 2.0f*powf(q(0,0),2.0f)-1.0f + 2.0f*powf(q(1,0),2.0f);
  C_N2B(1,1) = 2.0f*powf(q(0,0),2.0f)-1.0f + 2.0f*powf(q(2,0),2.0f);
  C_N2B(2,2) = 2.0f*powf(q(0,0),2.0f)-1.0f + 2.0f*powf(q(3,0),2.0f);

  C_N2B(0,1) = 2.0f*q(1,0)*q(2,0) + 2.0f*q(0,0)*q(3,0);
  C_N2B(0,2) = 2.0f*q(1,0)*q(3,0) - 2.0f*q(0,0)*q(2,0);

  C_N2B(1,0) = 2.0f*q(1,0)*q(2,0) - 2.0f*q(0,0)*q(3,0);
  C_N2B(1,2) = 2.0f*q(2,0)*q(3,0) + 2.0f*q(0,0)*q(1,0);

  C_N2B(2,0) = 2.0f*q(1,0)*q(3,0) + 2.0f*q(0,0)*q(2,0);
  C_N2B(2,1) = 2.0f*q(2,0)*q(3,0) - 2.0f*q(0,0)*q(1,0);
  return C_N2B;
}

// quaternion multiplication
Matrix ekfNavINS::qmult(Matrix p, Matrix q) {
  Matrix r;
  r(0,0) = p(0,0)*q(0,0) - (p(1,0)*q(1,0) + p(2,0)*q(2,0) + p(3,0)*q(3,0));
  r(1,0) = p(0,0)*q(1,0) + q(0,0)*p(1,0) + p(2,0)*q(3,0) - p(3,0)*q(2,0);
  r(2,0) = p(0,0)*q(2,0) + q(0,0)*p(2,0) + p(3,0)*q(1,0) - p(1,0)*q(3,0);
  r(3,0) = p(0,0)*q(3,0) + q(0,0)*p(3,0) + p(1,0)*q(2,0) - p(2,0)*q(1,0);
  return r;
}

// bound angle between -180 and 180
float ekfNavINS::constrainAngle180(float dta) {
  if(dta >  M_PI) dta -= (M_PI*2.0f);
  if(dta < -M_PI) dta += (M_PI*2.0f);
  return dta;
}

// bound angle between 0 and 360
float ekfNavINS::constrainAngle360(float dta){
  dta = fmod(dta,2.0f*M_PI);
  if (dta < 0)
    dta += 2.0f*M_PI;
  return dta;
}

Matrix ekfNavINS::toQuaternion(float yaw, float pitch, float roll) {
    float cy = cosf(yaw * 0.5f);
    float sy = sinf(yaw * 0.5f);
    float cp = cosf(pitch * 0.5f);
    float sp = sinf(pitch * 0.5f);
    float cr = cosf(roll * 0.5f);
    float sr = sinf(roll * 0.5f);
    Matrix q;
    q(0) = cr * cp * cy + sr * sp * sy; // w
    q(1) = cr * cp * sy - sr * sp * cy; // x
    q(2) = cr * sp * cy + sr * cp * sy; // y
    q(3) = sr * cp * cy - cr * sp * sy; // z
    return q;
}

std::tuple<float, float, float> ekfNavINS::toEulerAngles(Matrix quat) {
    float roll, pitch, yaw;
    // roll (x-axis rotation)
    float sinr_cosp = 2.0f * (quat(0,0)*quat(1,0)+quat(2,0)*quat(3,0));
    float cosr_cosp = 1.0f - 2.0f * (quat(1,0)*quat(1,0)+quat(2,0)*quat(2,0));
    roll = atan2f(sinr_cosp, cosr_cosp);
    // pitch (y-axis rotation)
    double sinp = 2.0f * (quat(0,0)*quat(2,0) - quat(1,0)*quat(3,0));
    //angles.pitch = asinf(-2.0f*(quat(1,0)*quat(3,0)-quat(0,0)*quat(2,0)));
    if (std::abs(sinp) >= 1)
        pitch = std::copysign(M_PI / 2.0f, sinp); // use 90 degrees if out of range
    else
        pitch = asinf(sinp);
    // yaw (z-axis rotation)
    float siny_cosp = 2.0f * (quat(1,0)*quat(2,0)+quat(0,0)*quat(3,0));
    float cosy_cosp = 1.0f - 2.0f * (quat(2,0)*quat(2,0)+quat(3,0)*quat(3,0));
    yaw = atan2f(siny_cosp, cosy_cosp);
    return std::make_tuple(roll, pitch, yaw);
}

