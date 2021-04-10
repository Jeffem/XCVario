/*
*) Refactor the code to remove reduentent part and improve the readabilty. 
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

#pragma once

#include <stdint.h>
#include <math.h>
#include <tuple>
#include <mutex>
#include <shared_mutex>
#include "../components/Linear-Algebra/src/Matrix.hpp"
#include <utility>
   using namespace std;

constexpr float SIG_W_A = 0.05f;
// Std dev of gyro output noise (rad/s)
constexpr float SIG_W_G = 0.00175f;
// Std dev of Accelerometer Markov Bias
constexpr float SIG_A_D = 0.01f;
// Correlation time or time constant
constexpr float TAU_A = 100.0f;
// Std dev of correlated gyro bias
constexpr float SIG_G_D = 0.00025;
// Correlati1on time or time constant
constexpr float TAU_G = 50.0f;
// GPS measurement noise std dev (m)
constexpr float SIG_GPS_P_NE = 3.0f;
constexpr float SIG_GPS_P_D = 6.0f;
// GPS measurement noise std dev (m/s)
constexpr float SIG_GPS_V_NE = 0.5f;
constexpr float SIG_GPS_V_D = 1.0f;
// Initial set of covariance
constexpr float P_P_INIT = 10.0f;
constexpr float P_V_INIT = 1.0f;
constexpr float P_A_INIT = 0.34906f;
constexpr float P_HDG_INIT = 3.14159f;
constexpr float P_AB_INIT = 0.9810f;
constexpr float P_GB_INIT = 0.01745f;
// acceleration due to gravity
constexpr float G = 9.807f;
// major eccentricity squared
constexpr double ECC2 = 0.0066943799901;
// earth semi-major axis radius (m)
constexpr double EARTH_RADIUS = 6378137.0;

class gpsCoordinate {
    public:
        double lat;
        double lon;
        double alt;
};

class gpsVelocity {
    public:
        double vN;
        double vE;
        double vD;
};

class imuData {
    public:
        float gyroX;
        float gyroY;
        float gyroZ;
        float accX;
        float accY;
        float accZ;
        float hX;
        float hY;
        float hZ;
};

class ekfNavINS {
  public:
    // ekf_update
    void ekf_update( uint64_t time/*, unsigned long TOW*/,   /* Time, Time of the week from GPS */
                    double vn, double ve, double vd,    /* Velocity North, Velocity East, Velocity Down */
                    double lat, double lon, double alt, /* GPS latitude, GPS longitude, GPS/Barometer altitude */
                    float p, float q, float r,          /* Gyro P, Q and R  */
                    float ax, float ay, float az,       /* Accelarometer X, Y and Z */
                    float hx, float hy, float hz        /* Magnetometer HX, HY, HZ */ );
    // returns whether the INS has been initialized
    bool initialized()          { return initialized_; }
    // returns the pitch angle, rad
    float getPitch_rad()        { return theta; }
    // returns the roll angle, rad
    float getRoll_rad()         { return phi; }
    // returns the heading angle, rad
    float getHeadingConstrainAngle180_rad()      { return constrainAngle180(psi); }
    float getHeading_rad()      { return psi; }
    // returns the INS latitude, rad
    double getLatitude_rad()    { return lat_ins; }
    // returns the INS longitude, rad
    double getLongitude_rad()   { return lon_ins; }
    // returns the INS altitude, m
    double getAltitude_m()      { return alt_ins; }
    // returns the INS north velocity, m/s
    double getVelNorth_ms()     { return vn_ins; }
    // returns the INS east velocity, m/s
    double getVelEast_ms()      { return ve_ins; }
    // returns the INS down velocity, m/s
    double getVelDown_ms()      { return vd_ins; }
    // returns the INS ground track, rad
    float getGroundTrack_rad()  { return atan2f((float)ve_ins,(float)vn_ins); }
    // returns the gyro bias estimate in the x direction, rad/s
    float getGyroBiasX_rads()   { return gbx; }
    // returns the gyro bias estimate in the y direction, rad/s
    float getGyroBiasY_rads()   { return gby; }
    // returns the gyro bias estimate in the z direction, rad/s
    float getGyroBiasZ_rads()   { return gbz; }
    // returns the accel bias estimate in the x direction, m/s/s
    float getAccelBiasX_mss()   { return abx; }
    // returns the accel bias estimate in the y direction, m/s/s
    float getAccelBiasY_mss()   { return aby; }
    // returns the accel bias estimate in the z direction, m/s/s
    float getAccelBiasZ_mss()   { return abz; }
    // return pitch, roll and yaw
    std::tuple<float,float,float> getPitchRollYaw(float ax, float ay, float az, float hx, float hy, float hz);
    void imuUpdateEKF(uint64_t time, imuData imu);
    void gpsCoordinateUpdateEKF(gpsCoordinate coor);
    void gpsVelocityUpdateEKF(gpsVelocity vel);

  private:
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// member variables /////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    gpsCoordinate gpsCoor;
    gpsVelocity   gpsVel;
    imuData       imuDat;
//    mutable std::shared_mutex shMutex;
    // initialized
    bool initialized_ = false;
    // timing
    uint64_t _tprev;
    //float _dt;
    unsigned long previousTOW;
    // estimated attitude
    float phi, theta, psi;
    // estimated NED velocity
    double vn_ins, ve_ins, vd_ins;
    // estimated location
    double lat_ins, lon_ins, alt_ins;
    // magnetic heading corrected for roll and pitch angle
    float Bxc, Byc;
    // accelerometer bias
    float abx = 0.0, aby = 0.0, abz = 0.0;
    // gyro bias
    float gbx = 0.0, gby = 0.0, gbz = 0.0;
    // earth radius at location
    double Re, Rn, denom;
    Matrix I15 = Matrix::identity(15);
    Matrix I6 = Matrix::identity(6);
    Matrix I3 = Matrix::identity(3);
    // State matrix
    Matrix Fs = Matrix::identity(15);
    // State transition matrix
    Matrix PHI = Matrix::zeros(15,15);
    // Covariance matrix
    Matrix P = Matrix::zeros(15,15);
    // For process noise transformation
    Matrix Gs = Matrix::zeros(15,12);
    Matrix Rw = Matrix::zeros(12,12);
    // Process noise matrix
    Matrix Q = Matrix::zeros(15,15);
    // Gravity model
    Matrix grav = Matrix::zeros(3,1);
    // Rotation rate
    Matrix om_ib = Matrix::zeros(3,1);
    // Specific force
    Matrix f_b = Matrix::zeros(3,1);
    // DCM
    Matrix C_N2B = Matrix::zeros(3,3);
    // DCM transpose
    Matrix C_B2N = Matrix::zeros(3,3);
    // Temporary to get dxdt
    Matrix dx = Matrix::zeros(3,1);
    Matrix dxd = Matrix::zeros(3,1);
    // Estimated INS
    Matrix estmimated_ins = Matrix::zeros(3,1);
    // NED velocity INS
    Matrix V_ins = Matrix::zeros(3,1);
    // LLA INS
    Matrix lla_ins = Matrix::zeros(3,1);
    // NED velocity GPS
    Matrix V_gps = Matrix::zeros(3,1);
    // LLA GPS
    Matrix lla_gps = Matrix::zeros(3,1);
    // Position ECEF INS
    Matrix pos_ecef_ins = Matrix::zeros(3,1);
    // Position NED INS
    Matrix pos_ned_ins = Matrix::zeros(3,1);
    // Position ECEF GPS
    Matrix pos_ecef_gps = Matrix::zeros(3,1);
    // Position NED GPS
    Matrix pos_ned_gps = Matrix::zeros(3,1);
    // Quat
    Matrix quat = Matrix::zeros(4,1);
    // dquat
    Matrix dq = Matrix::zeros(4,1);
    // difference between GPS and INS
    Matrix y = Matrix::zeros(6,1);
    // GPS measurement noise
    Matrix R = Matrix::zeros(6,6);
    Matrix x = Matrix::zeros(15,1);
    // Kalman Gain
    Matrix K = Matrix::zeros(15,6);
    Matrix H = Matrix::zeros(6,15);
    // skew symmetric
    Matrix sk(Matrix w);

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// member functions /////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // ekf_init
    void ekf_init(uint64_t time, 
                 double vn,double ve,double vd, 
                 double lat,double lon,double alt,
                 float p,float q,float r,
                 float ax,float ay,float az,
                 float hx,float hy, float hz);
    // lla rate
    Matrix llarate(Matrix V, Matrix lla);
    Matrix llarate(Matrix V, double lat, double alt);
    // lla to ecef
    Matrix lla2ecef(Matrix lla);
    // ecef to ned
    Matrix ecef2ned(Matrix ecef, Matrix pos_ref);
    // quaternion to dcm
    Matrix quat2dcm(Matrix q);
    // quaternion multiplication
    Matrix qmult(Matrix p, Matrix q);
    // maps angle to +/- 180
    float constrainAngle180(float dta);
    // maps angle to 0-360
    float constrainAngle360(float dta);
    // Returns Radius - East West and Radius - North South
    // la paire génère des erreurs, remplacées par constante pour lat=45°
    // voir si on peut revenir à une fonction avec un tuple au lieu d'une paire
    double Rew = 6421050;
    double Rns = 6367381;
    //constexpr std::pair<double, double> earthradius(double lat);

    // Yaw, Pitch, Roll to Quarternion
    Matrix toQuaternion(float yaw, float pitch, float roll);
    // Quarternion to Yaw, Pitch, Roll
    std::tuple<float, float, float> toEulerAngles(Matrix quat);
    // Update Jacobian matrix
    void updateJacobianMatrix();
    // Update Process Noise and Covariance Time
    void updateProcessNoiseCovarianceTime(float _dt);
    // Update Gyro and Accelerometer Bias
    void updateBias(float ax,float ay,float az,float p,float q, float r);
    // Update 15 states after KF state update
    void update15statesAfterKF();
    // Update difference between predicted and calculated GPS and IMU values
    void updateCalculatedVsPredicted();
    void ekf_update(uint64_t time);
    void updateINS();
};
