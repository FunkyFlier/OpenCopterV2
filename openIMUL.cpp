// Implementation of Madgwick's IMU and AHRS algorithms.
// See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/
//This code is provided under the GNU General Public Licence
//12/19/2012
#include "openIMUL.h"
//#include <Streaming.h>



openIMU::openIMU(float *gyroX, float *gyroY, float *gyroZ, float *accX, float *accY, 
float *accZ, float *scAccX, float *scAccY, float *scAccZ, float *magX, float *magY, float *magZ, float *xGPS , float *yGPS, float *zBaro,float *xGPSVel ,float *yGPSVel,float *zBaroVel,float *G_Dt){

  //constructor for the 10DOF system
  gx = gyroX;
  gy = gyroY;
  gz = gyroZ;

  ax = accX;
  ay = accY;
  az = accZ;

  mx = magX;
  my = magY;
  mz = magZ;

  sax = scAccX;
  say = scAccY;
  saz = scAccZ;

  XPosGPS = xGPS;
  YPosGPS = yGPS;
  ZPosBaro = zBaro;
  XVelGPS = xGPSVel;
  YVelGPS = yGPSVel;
  ZVelBaro = zBaroVel;
  dt = G_Dt;

  q0.val = 1;
  q1.val = 0;
  q2.val = 0;
  q3.val = 0;

  //Kalman stuff
  XEst.val = 0;
  velX.val = 0;

  YEst.val = 0;
  velY.val = 0;

  velZ.val = 0;
  ZEst.val = 0;


  kPosGPS = 0.20;
  kVelGPS = 0.36;
  kAccGPS = 0.00003;

  kPosBaro = 0.09;
  kVelBaro = 0.09;
  kAccBaro = 1e-6;

  for (int i = 0; i < LAG_SIZE; i++){
    XEstHist[i] = 0;
    YEstHist[i] = 0;
    XVelHist[i] = 0;
    YVelHist[i] = 0;
  }
  for (int i = 0; i < LAG_SIZE_BARO; i++){
    ZEstHist[i] = 0;
    ZVelHist[i] = 0;
  }


}

void openIMU::GetGravOffset(void){
  initialAccMagnitude.val = -1.0 * sqrt(*sax * *sax + *say * *say + *saz * *saz) ;
  //Serial<<*sax<<","<<*say<<","<<*saz<<"\r\n";
  //Serial.println(initialAccMagnitude);
}
void openIMU::GetInertial(void){

  inertialX.val = R11 * *sax + R21 * *say + R31 * *saz;// - inertialXOffSet.val;
  inertialY.val = R12 * *sax + R22 * *say + R32 * *saz;// - inertialYOffSet.val;
  inertialZ.val = R13 * *sax + R23 * *say + R33 * *saz - initialAccMagnitude.val;// - inertialZOffSet.val;
}


void openIMU::UpdateLagIndex(void){
  /*
  currentEstIndex++;
   
   lagIndex++;
   if (currentEstIndex == LAG_SIZE){
   currentEstIndex = 0;
   }
   if (lagIndex == LAG_SIZE){
   lagIndex = 0;
   
   }*/
   

  currentEstIndex++;
  if (currentEstIndex == (LAG_SIZE)){
    currentEstIndex = 0;
  }
  lagIndex = currentEstIndex - 53;
  /*
  if (lagAmt > 0){
    lagIndex = currentEstIndex - 53;
  }
  else{
    lagIndex = currentEstIndex;
  }*/
  
  
  if (lagIndex < 0){
    lagIndex = LAG_SIZE + lagIndex;
  }
  
  
  
  currentEstIndex_z++;
  if (currentEstIndex_z == (LAG_SIZE_BARO)){
    currentEstIndex_z = 0;
  }
  lagIndex_z = currentEstIndex_z - 23;

  if (lagIndex_z < 0){
    lagIndex_z = LAG_SIZE_BARO + lagIndex_z;
  }
}



void openIMU::InitialQuat(){
  //calculate the ypr from sensors convert to quaternion and rotation matrix
  radPitch.val = atan2(*sax,sqrt(*say * *say + *saz * *saz));
  radRoll.val = atan2(-1.0 * *say, -1.0 * *saz);
  //Serial.println(ToDeg(radPitch.val));
  //Serial.println(ToDeg(radRoll.val));

  R11 = cos(radPitch.val);
  R13 = -1.0 * sin(radPitch.val);
  R21 = sin(radRoll.val)*sin(radPitch.val);
  R22 = cos(radRoll.val);
  R23 = cos(radPitch.val)*sin(radRoll.val);

  bx = *mx * R11 + *mz * R13;
  by = *mx * R21 + *my * R22 + *mz * R23;
  radYaw.val = atan2(-1.0 * by, bx) - DECLINATION;
  //Serial.println(ToDeg(radYaw.val));
  //Serial.println(DECLINATION);
  //Serial.println(ToDeg(DECLINATION));

  //Serial<<*ax<<","<<*ay<<","<<*az<<"\r\n";
  //Serial<<*sax<<","<<*say<<","<<*saz<<"\r\n";

  q0.val = cos(radYaw.val/2.0)*cos(radPitch.val/2.0)*cos(radRoll.val/2.0) + sin(radYaw.val/2.0)*sin(radPitch.val/2.0)*sin(radRoll.val/2.0); 
  q1.val = cos(radYaw.val/2.0)*cos(radPitch.val/2.0)*sin(radRoll.val/2.0) - sin(radYaw.val/2.0)*sin(radPitch.val/2.0)*cos(radRoll.val/2.0); 
  q2.val = cos(radYaw.val/2.0)*sin(radPitch.val/2.0)*cos(radRoll.val/2.0) + sin(radYaw.val/2.0)*cos(radPitch.val/2.0)*sin(radRoll.val/2.0); 
  q3.val = sin(radYaw.val/2.0)*cos(radPitch.val/2.0)*cos(radRoll.val/2.0) - cos(radYaw.val/2.0)*sin(radPitch.val/2.0)*sin(radRoll.val/2.0);
  magnitude.val = sqrt(q0.val *  q0.val + q1.val *  q1.val + q2.val *  q2.val + q3.val *  q3.val); 
  q0.val = q0.val / magnitude.val;
  q1.val = q1.val / magnitude.val;
  q2.val = q2.val / magnitude.val;
  q3.val = q3.val / magnitude.val;


  q0q0 = q0.val*q0.val;
  q1q1 = q1.val*q1.val;
  q2q2 = q2.val*q2.val;
  q3q3 = q3.val*q3.val;

  q0q1 = q0.val*q1.val;
  q0q2 = q0.val*q2.val;
  q0q3 = q0.val*q3.val;

  q1q2 = q1.val*q2.val;
  q1q3 = q1.val*q3.val;

  q2q3 = q2.val*q3.val;
  //generate rotation matrix
  R11 = 2*(q0q0-0.5+q1q1);
  R12 = 2*(q1q2+q0q3);
  R13 = 2*(q1q3-q0q2);
  R21 = 2*(q1q2-q0q3);
  R22 = 2*(q0q0-0.5+q2q2);
  R23 = 2*(q2q3+q0q1);
  R31 = 2*(q1q3+q0q2);
  R32 = 2*(q2q3-q0q1);
  R33 = 2*(q0q0-0.5+q3q3);  

  //rotate by declination 
  COS_DEC = cos(DECLINATION);
  SIN_DEC = sin(DECLINATION);
  R11 = R11*COS_DEC - R12*SIN_DEC;
  R12 = R12*COS_DEC + R11*SIN_DEC;

  R21 = R21*COS_DEC - R22*SIN_DEC;
  R22 = R22*COS_DEC + R21*SIN_DEC;

  R31 = R31*COS_DEC - R32*SIN_DEC;
  R32 = R32*COS_DEC + R31*SIN_DEC;


  //Serial<<_FLOAT(R11,7)<<","<<_FLOAT(R12,7)<<","<<_FLOAT(R13,7)<<"\r\n"<<_FLOAT(R21,7)<<","<<_FLOAT(R22,7)<<","<<_FLOAT(R23,7)<<"\r\n"<<_FLOAT(R31,7)<<","<<_FLOAT(R32,7)<<","<<_FLOAT(R33,7)<<"\r\n";
  //Serial<<"+++\r\n";
  rotError = (R11*R21 + R12*R22 + R13*R23) * 0.5;
  //Serial<<_FLOAT(rotError,7)<<"\r\n";
  xOrtho[0] = R11 - R21 * rotError;
  xOrtho[1] = R12 - R22 * rotError;
  xOrtho[2] = R13 - R23 * rotError;

  //Serial<<_FLOAT(xOrtho[0],7)<<","<<_FLOAT(xOrtho[1],7)<<","<<_FLOAT(xOrtho[2],7)<<"\r\n";

  yOrtho[0] = R21 - R11 * rotError;
  yOrtho[1] = R22 - R12 * rotError;
  yOrtho[2] = R23 - R13 * rotError;
  //Serial<<_FLOAT(yOrtho[0],7)<<","<<_FLOAT(yOrtho[1],7)<<","<<_FLOAT(yOrtho[2],7)<<"\r\n";

  zOrtho[0] = xOrtho[1] * yOrtho[2] - xOrtho[2] * yOrtho[1];
  zOrtho[1] = xOrtho[2] * yOrtho[0] - xOrtho[0] * yOrtho[2];
  zOrtho[2] = xOrtho[0] * yOrtho[1] - xOrtho[1] * yOrtho[0];
  //Serial<<_FLOAT(zOrtho[0],7)<<","<<_FLOAT(zOrtho[1],7)<<","<<_FLOAT(zOrtho[2],7)<<"\r\n";

  //normScale = 0.5 * (3 - xOrtho[0] * xOrtho[0] + xOrtho[1] * xOrtho[1] + xOrtho[2] * xOrtho[2] );
  //normScale = InvSqrt(xOrtho[0] * xOrtho[0] + xOrtho[1] * xOrtho[1] + xOrtho[2] * xOrtho[2]);
  normScale = 1.0/sqrt(xOrtho[0] * xOrtho[0] + xOrtho[1] * xOrtho[1] + xOrtho[2] * xOrtho[2]);
  //Serial<<_FLOAT(normScale,7)<<"\r\n";
  R11 = xOrtho[0] * normScale;
  R12 = xOrtho[1] * normScale;
  R13 = xOrtho[2] * normScale;

  //normScale = 0.5 * (3 - yOrtho[0] * yOrtho[0] + yOrtho[1] * yOrtho[1] + yOrtho[2] * yOrtho[2] );
  //normScale = InvSqrt(yOrtho[0] * yOrtho[0] + yOrtho[1] * yOrtho[1] + yOrtho[2] * yOrtho[2]);
  normScale = 1.0/sqrt(yOrtho[0] * yOrtho[0] + yOrtho[1] * yOrtho[1] + yOrtho[2] * yOrtho[2]);
  //Serial<<_FLOAT(normScale,7)<<"\r\n";
  R21 = yOrtho[0] * normScale;
  R22 = yOrtho[1] * normScale;
  R23 = yOrtho[2] * normScale;

  //normScale = 0.5 * (3 - zOrtho[0] * zOrtho[0] + zOrtho[1] * zOrtho[1] + zOrtho[2] * zOrtho[2] );
  //normScale = InvSqrt(zOrtho[0] * zOrtho[0] + zOrtho[1] * zOrtho[1] + zOrtho[2] * zOrtho[2]);
  normScale = 1.0/sqrt(zOrtho[0] * zOrtho[0] + zOrtho[1] * zOrtho[1] + zOrtho[2] * zOrtho[2]);
  //Serial<<_FLOAT(normScale,7)<<"\r\n";
  R31 = zOrtho[0] * normScale;
  R32 = zOrtho[1] * normScale;
  R33 = zOrtho[2] * normScale;

  //Serial<<_FLOAT(R11,7)<<","<<_FLOAT(R12,7)<<","<<_FLOAT(R13,7)<<"\r\n"<<_FLOAT(R21,7)<<","<<_FLOAT(R22,7)<<","<<_FLOAT(R23,7)<<"\r\n"<<_FLOAT(R31,7)<<","<<_FLOAT(R32,7)<<","<<_FLOAT(R33,7)<<"\r\n";
  //Serial<<"---\r\n";
  recipNorm = InvSqrt(*mx * *mx + *my * *my + *mz * *mz);
  *mx *= recipNorm;
  *my *= recipNorm;
  *mz *= recipNorm;

  hx = R11 * *mx + R21 * *my + R31 * *mz;
  hy = R12 * *mx + R22 * *my + R32 * *mz;
  hz = R13 * *mx + R23 * *my + R33 * *mz;

  bx = sqrt(hx * hx + hy * hy);
  bz = hz;


}

void openIMU::Predict(void){



  biasedX = (*sax - accelBiasX.val);
  biasedY = (*say - accelBiasY.val);
  biasedZ = (*saz - accelBiasZ.val);
  inertialXBiased.val = R11 * biasedX + R21 * biasedY + R31 * biasedZ;//  - inertialXOffSet.val;
  inertialYBiased.val = R12 * biasedX + R22 * biasedY + R32 * biasedZ;// - inertialYOffSet.val;
  inertialZBiased.val = R13 * biasedX + R23 * biasedY + R33 * biasedZ - initialAccMagnitude.val;// - inertialZOffSet.val; 



  velX.val = velX.val + inertialXBiased.val * *dt;
  velY.val = velY.val + inertialYBiased.val * *dt;
  velZ.val = velZ.val + inertialZBiased.val * *dt;


  XEst.val = XEst.val + velX.val * *dt;
  YEst.val = YEst.val + velY.val * *dt;
  ZEst.val = ZEst.val + velZ.val * *dt;



  XEstHist[currentEstIndex] = XEst.val;
  YEstHist[currentEstIndex] = YEst.val;
  
  
  XVelHist[currentEstIndex] = velX.val;
  YVelHist[currentEstIndex] = velY.val;
  
  
  ZEstHist[currentEstIndex_z] = ZEst.val;
  ZVelHist[currentEstIndex_z] = velZ.val;
  lagZVel.val = -1.0 * ZVelHist[lagIndex_z];
  //lagEstX.val = XVelHist[lagIndex];
  //lagEstY.val = YVelHist[lagIndex];


  ZEstUp.val = -1.0 * ZEst.val;
  velZUp.val = -1.0 * velZ.val;


}
void openIMU::CorrectGPS(void){
  xPosError.val = XEstHist[lagIndex] - *XPosGPS;
  yPosError.val = YEstHist[lagIndex] - *YPosGPS;

  xVelError.val = XVelHist[lagIndex] - *XVelGPS;
  yVelError.val = YVelHist[lagIndex] - *YVelGPS;

  XEst.val = XEst.val - kPosGPS * xPosError.val;
  YEst.val = YEst.val - kPosGPS * yPosError.val;

  velX.val = velX.val - kVelGPS * xVelError.val;
  velY.val = velY.val - kVelGPS * yVelError.val;




  //Serial<<_FLOAT(accelBiasX.val,5)<<","<<_FLOAT(accelBiasY.val,5)<<","<<_FLOAT(accelBiasZ.val,5)<<"\r\n";
  accelBiasXEF = R11 * accelBiasX.val + R21 * accelBiasY.val + R31 * accelBiasZ.val;
  accelBiasYEF = R12 * accelBiasX.val + R22 * accelBiasY.val + R32 * accelBiasZ.val;
  accelBiasZEF = R13 * accelBiasX.val + R23 * accelBiasY.val + R33 * accelBiasZ.val;

  //Serial<<_FLOAT(R11,4)<<","<<_FLOAT(R12,4)<<","<<_FLOAT(R13,4)<<"\r\n"<<_FLOAT(R21,4)<<","<<_FLOAT(R22,4)<<","<<_FLOAT(R23,4)<<"\r\n"<<_FLOAT(R31,4)<<","<<_FLOAT(R32,4)<<","<<_FLOAT(R33,4)<<"\r\n";
  //Serial<<_FLOAT(accelBiasXEF,5)<<","<<_FLOAT(accelBiasYEF,5)<<","<<_FLOAT(accelBiasZEF,5)<<","<<_FLOAT(xVelError.val,5)<<","<< _FLOAT(yVelError.val,5)<<"\r\n";
  accelBiasXEF = accelBiasXEF + kAccGPS * xVelError.val;
  accelBiasYEF = accelBiasYEF + kAccGPS * yVelError.val;
  //Serial<<_FLOAT(accelBiasXEF,5)<<","<<_FLOAT(accelBiasYEF,5)<<","<<_FLOAT(accelBiasZEF,5)<<"\r\n";

  accelBiasX.val = R11*accelBiasXEF + R12*accelBiasYEF + R13*accelBiasZEF;
  accelBiasY.val = R21*accelBiasXEF + R22*accelBiasYEF + R23*accelBiasZEF;
  accelBiasZ.val = R31*accelBiasXEF + R32*accelBiasYEF + R33*accelBiasZEF;
  //Serial<<_FLOAT(accelBiasX.val,5)<<","<<_FLOAT(accelBiasY.val,5)<<","<<_FLOAT(accelBiasZ.val,5)<<"\r\n";
  //Serial<<"*******\r\n";

}

void openIMU::CorrectAlt(void){
  zPosError.val = ZEstHist[lagIndex_z] + *ZPosBaro;
  zVelError.val = ZVelHist[lagIndex_z] + *ZVelBaro;
  
  ZEst.val = ZEst.val - kPosBaro * zPosError.val;
  velZ.val = velZ.val - kVelBaro * zVelError.val;

  accelBiasXEF = R11*accelBiasX.val + R21*accelBiasY.val + R31*accelBiasZ.val;
  accelBiasYEF = R12*accelBiasX.val + R22*accelBiasY.val + R32*accelBiasZ.val;
  accelBiasZEF = R13*accelBiasX.val + R23*accelBiasY.val + R33*accelBiasZ.val;


  accelBiasZEF = accelBiasZEF + kAccBaro * zVelError.val;

  accelBiasX.val = R11*accelBiasXEF + R12*accelBiasYEF + R13*accelBiasZEF;
  accelBiasY.val = R21*accelBiasXEF + R22*accelBiasYEF + R23*accelBiasZEF;
  accelBiasZ.val = R31*accelBiasXEF + R32*accelBiasYEF + R33*accelBiasZEF;

  ZEstUp.val = -1.0 * ZEst.val;
  velZUp.val = -1.0 * velZ.val;
}
/*
void openIMU::AHRSStart(void){
 
 }
 void openIMU::AHRSEnd(){
 
 }*/
void openIMU::IntegrateGyro(){
  magnitude.val =  sqrt(*ax * *ax + *ay * *ay + *az * *az);
  magnitudeDifference.val = fabs(initialAccMagnitude.val +  magnitude.val);
  if (magnitudeDifference.val > FEEDBACK_LIMIT){
    skipFeedBack = true;
  }
  if (kiAcc > 0){

    *gx = *gx + integralFBX;
    *gy = *gy + integralFBY;
    *gz = *gz + integralFBZ;  
  }
  dtby2 = *dt * 0.5;
  q0.val += -1 * dtby2*(*gx * q1.val + *gy * q2.val + *gz * q3.val);
  q1.val +=      dtby2*(*gx * q0.val - *gy * q3.val + *gz * q2.val);
  q2.val +=      dtby2*(*gx * q3.val + *gy * q0.val - *gz * q1.val);
  q3.val +=      dtby2*(*gy * q1.val - *gx * q2.val + *gz * q0.val);


  //normalize the quaternion
  recipNorm = InvSqrt(q0.val * q0.val + q1.val * q1.val + q2.val * q2.val + q3.val * q3.val);
  q0.val *= recipNorm;
  q1.val *= recipNorm;
  q2.val *= recipNorm;
  q3.val *= recipNorm;

  q0q0 = q0.val*q0.val;
  q1q1 = q1.val*q1.val;
  q2q2 = q2.val*q2.val;
  q3q3 = q3.val*q3.val;

  q0q1 = q0.val*q1.val;
  q0q2 = q0.val*q2.val;
  q0q3 = q0.val*q3.val;

  q1q2 = q1.val*q2.val;
  q1q3 = q1.val*q3.val;

  q2q3 = q2.val*q3.val;
  //generate rotation matrix
  R11 = 2*(q0q0-0.5+q1q1);
  R12 = 2*(q1q2+q0q3);
  R13 = 2*(q1q3-q0q2);
  R21 = 2*(q1q2-q0q3);
  R22 = 2*(q0q0-0.5+q2q2);
  R23 = 2*(q2q3+q0q1);
  R31 = 2*(q1q3+q0q2);
  R32 = 2*(q2q3-q0q1);
  R33 = 2*(q0q0-0.5+q3q3);  
  //COS_DEC = cos(DECLINATION);
  //SIN_DEC = sin(DECLINATION);
  //rotate by declination 
  R11 = R11*COS_DEC - R12*SIN_DEC;
  R12 = R12*COS_DEC + R11*SIN_DEC;

  R21 = R21*COS_DEC - R22*SIN_DEC;
  R22 = R22*COS_DEC + R21*SIN_DEC;

  R31 = R31*COS_DEC - R32*SIN_DEC;
  R32 = R32*COS_DEC + R31*SIN_DEC;



}

void openIMU::AHRSupdate() {
  //normalize the sensor readings
  magnitude.val =  sqrt(*ax * *ax + *ay * *ay + *az * *az);
  magnitudeDifference.val = fabs(initialAccMagnitude.val +  magnitude.val);
  //Serial2<<initialAccMagnitude.val<<","<<magnitude.val<<","<<magnitudeDifference.val<<"\r\n";
  //if (magnitudeDifference.val < FEEDBACK_LIMIT && fabs(inertialX.val) < 0.3 && fabs(inertialY.val) < 0.3){
  //if (magnitudeDifference.val < FEEDBACK_LIMIT && sqrt(inertialX.val * inertialX.val + inertialY.val * inertialY.val) < 0.3 ){
    //if (magnitudeDifference.val < FEEDBACK_LIMIT && sqrt(*gx * *gx + *gy * *gy + *gz * *gz) < 0.018 ){
  if (magnitudeDifference.val < FEEDBACK_LIMIT ){
   //initialAccMagnitude.val = initialAccMagnitude.val * 0.9 - magnitude.val * 0.1;
    //initialAccMagnitude.val = initialAccMagnitude.val * 0.25 - sqrt(*sax * *sax + *say * *say + *saz * *saz) * 0.75;

    feedBack = 2;
    recipNorm = 1/magnitude.val;
    *ax *= recipNorm;
    *ay *= recipNorm;
    *az *= recipNorm;

    recipNorm = InvSqrt(*mx * *mx + *my * *my + *mz * *mz);
    *mx *= recipNorm;
    *my *= recipNorm;
    *mz *= recipNorm;

    hx = R11 * *mx + R21 * *my + R31 * *mz;
    hy = R12 * *mx + R22 * *my + R32 * *mz;
    hz = R13 * *mx + R23 * *my + R33 * *mz;

    /*bx_ = sqrt(hx * hx + hy * hy);
    bz_ = hz;

    bx = bx * 0.95 + bx_ * 0.05;
    bz = bz * 0.95 + bz_ * 0.05;*/
    bx = sqrt(hx * hx + hy * hy);
    bz = hz;


    wx = R11*bx + R13*bz;
    wy = R21*bx + R23*bz;
    wz = R31*bx + R33*bz;

    vx = R13;
    vy = R23;
    vz = R33;

    exm = (*my * wz - *mz * wy);
    eym = (*mz * wx - *mx * wz);
    ezm = (*mx * wy - *my * wx);

    exa = (*ay * vz - *az * vy);
    eya = (*az * vx - *ax * vz);
    eza = (*ax * vy - *ay * vx);

    kiDTAcc = kiAcc * *dt;
    kiDTMag = kiMag * *dt;
    if (kiAcc > 0){
      integralFBX += exa * kiDTAcc+ exm * kiDTMag;
      integralFBY += eya * kiDTAcc+ eym * kiDTMag;
      integralFBZ += eza * kiDTAcc+ ezm * kiDTMag;
      *gx = *gx + integralFBX;
      *gy = *gy + integralFBY;
      *gz = *gz + integralFBZ;  
    }
    else{
      integralFBX = 0;
      integralFBY = 0;
      integralFBZ = 0;  
    }
    *gx += exa * kpAcc + exm * kpMag;
    *gy += eya * kpAcc + eym * kpMag;
    *gz += eza * kpAcc + ezm * kpMag;
  }
  else{
    feedBack = 0;
  }

  dtby2 = *dt * 0.5;
  //Serial2<<*gx<<","<<*gy<<","<<*gz<<","<<_FLOAT(*dt,5)<<","<<dtby2<<"\r\n";
  q0.val += -1.0 * dtby2*(*gx * q1.val + *gy * q2.val + *gz * q3.val);
  q1.val +=      dtby2*(*gx * q0.val - *gy * q3.val + *gz * q2.val);
  q2.val +=      dtby2*(*gx * q3.val + *gy * q0.val - *gz * q1.val);
  q3.val +=      dtby2*(*gy * q1.val - *gx * q2.val + *gz * q0.val);


  //normalize the quaternion
  recipNorm = 1/sqrt(q0.val * q0.val + q1.val * q1.val + q2.val * q2.val + q3.val * q3.val);
  q0.val *= recipNorm;
  q1.val *= recipNorm;
  q2.val *= recipNorm;
  q3.val *= recipNorm;
  /*
  q0q0 = q0.val*q0.val;
   q1q1 = q1.val*q1.val;
   q2q2 = q2.val*q2.val;
   q3q3 = q3.val*q3.val;
   
   q0q1 = q0.val*q1.val;
   q0q2 = q0.val*q2.val;
   q0q3 = q0.val*q3.val;
   
   q1q2 = q1.val*q2.val;
   q1q3 = q1.val*q3.val;
   
   q2q3 = q2.val*q3.val;
   //generate rotation matrix
   R11 = 2*(q0q0-0.5+q1q1);
   R12 = 2.0*(q1q2+q0q3);
   R13 = 2.0*(q1q3-q0q2);
   R21 = 2.0*(q1q2-q0q3);
   R22 = 2.0*(q0q0-0.5+q2q2);
   R23 = 2.0*(q2q3+q0q1);
   R31 = 2.0*(q1q3+q0q2);
   R32 = 2.0*(q2q3-q0q1);
   R33 = 2.0*(q0q0-0.5+q3q3);  
   
   //rotate by declination 
   R11 = R11*COS_DEC - R12*SIN_DEC;
   R12 = R12*COS_DEC + R11*SIN_DEC;
   
   R21 = R21*COS_DEC - R22*SIN_DEC;
   R22 = R22*COS_DEC + R21*SIN_DEC;
   
   R31 = R31*COS_DEC - R32*SIN_DEC;
   R32 = R32*COS_DEC + R31*SIN_DEC;
   */
}

void openIMU::GenerateRotationMatrix(void){
  q0q0 = q0.val*q0.val;
  q1q1 = q1.val*q1.val;
  q2q2 = q2.val*q2.val;
  q3q3 = q3.val*q3.val;

  q0q1 = q0.val*q1.val;
  q0q2 = q0.val*q2.val;
  q0q3 = q0.val*q3.val;

  q1q2 = q1.val*q2.val;
  q1q3 = q1.val*q3.val;

  q2q3 = q2.val*q3.val;
  //generate rotation matrix
  R11 = 2*(q0q0-0.5+q1q1);
  R12 = 2.0*(q1q2+q0q3);
  R13 = 2.0*(q1q3-q0q2);
  R21 = 2.0*(q1q2-q0q3);
  R22 = 2.0*(q0q0-0.5+q2q2);
  R23 = 2.0*(q2q3+q0q1);
  R31 = 2.0*(q1q3+q0q2);
  R32 = 2.0*(q2q3-q0q1);
  R33 = 2.0*(q0q0-0.5+q3q3);  
  COS_DEC = cos(DECLINATION);
  SIN_DEC = sin(DECLINATION);
  //rotate by declination 
  R11 = R11*COS_DEC - R12*SIN_DEC;
  R12 = R12*COS_DEC + R11*SIN_DEC;

  R21 = R21*COS_DEC - R22*SIN_DEC;
  R22 = R22*COS_DEC + R21*SIN_DEC;

  R31 = R31*COS_DEC - R32*SIN_DEC;
  R32 = R32*COS_DEC + R31*SIN_DEC;


  rotError = (R11*R21 + R12*R22 + R13*R23) * 0.5;
  //rotError2 = (R11 * R31 + R12 * R32 + R13 * R33);
  if (rotError != 0.0 ){
    xOrtho[0] = R11 - R21 * rotError;
    xOrtho[1] = R12 - R22 * rotError;
    xOrtho[2] = R13 - R23 * rotError;

    yOrtho[0] = R21 - R11 * rotError;
    yOrtho[1] = R22 - R12 * rotError;
    yOrtho[2] = R23 - R13 * rotError;

    zOrtho[0] = xOrtho[1] * yOrtho[2] - xOrtho[2] * yOrtho[1];
    zOrtho[1] = xOrtho[2] * yOrtho[0] - xOrtho[0] * yOrtho[2];
    zOrtho[2] = xOrtho[0] * yOrtho[1] - xOrtho[1] * yOrtho[0];

    normScale = 1.0/sqrt(xOrtho[0] * xOrtho[0] + xOrtho[1] * xOrtho[1] + xOrtho[2] * xOrtho[2]);

    R11 = xOrtho[0] * normScale;
    R12 = xOrtho[1] * normScale;
    R13 = xOrtho[2] * normScale;

    normScale = 1.0/sqrt(yOrtho[0] * yOrtho[0] + yOrtho[1] * yOrtho[1] + yOrtho[2] * yOrtho[2]);

    R21 = yOrtho[0] * normScale;
    R22 = yOrtho[1] * normScale;
    R23 = yOrtho[2] * normScale;

    normScale = 1.0/sqrt(zOrtho[0] * zOrtho[0] + zOrtho[1] * zOrtho[1] + zOrtho[2] * zOrtho[2]);

    R31 = zOrtho[0] * normScale;
    R32 = zOrtho[1] * normScale;
    R33 = zOrtho[2] * normScale;
  }
  else{
    normScale = 1.0/sqrt(R11 * R11 + R12 * R12 + R13 * R13);

    R11 = R11 * normScale;
    R12 = R12 * normScale;
    R13 = R13 * normScale;

    normScale = 1.0/sqrt(R21 *R21 + R22 * R22 + R23 * R23);

    R21 = R21 * normScale;
    R22 = R22 * normScale;
    R23 = R23 * normScale;

    normScale = 1.0/sqrt(R31 * R31 + R32 * R32 + R33 * R33);

    R31 = R31 * normScale;
    R32 = R32 * normScale;
    R33 = R33 * normScale;
  }
/*normScale = 1.0/sqrt(R11 * R11 + R12 * R12 + R13 * R13);

    R11 = R11 * normScale;
    R12 = R12 * normScale;
    R13 = R13 * normScale;

    normScale = 1.0/sqrt(R21 *R21 + R22 * R22 + R23 * R23);

    R21 = R21 * normScale;
    R22 = R22 * normScale;
    R23 = R23 * normScale;

    normScale = 1.0/sqrt(R31 * R31 + R32 * R32 + R33 * R33);

    R31 = R31 * normScale;
    R32 = R32 * normScale;
    R33 = R33 * normScale;*/
  /*
  R11 = 1;
   R12 = 0;
   R13 = 0;
   R21 = 0;
   R22 = 1;
   R23 = 0;
   R31 = 0;
   R32 = 0;
   R33 = 1;  */
}

void openIMU::GetEuler(void){
  rawRoll.val = ToDeg(FastAtan2(2 * (q0.val * q1.val + q2.val * q3.val),1 - 2.0 * (q1.val * q1.val + q2.val * q2.val)));
  roll.val=  rawRoll.val - rollOffset.val;

  rawPitch.val = ToDeg(asin(2.0 * (q0.val * q2.val - q3.val * q1.val)));
  pitch.val =  rawPitch.val - pitchOffset.val;

  yaw.val = ToDeg(FastAtan2(2.0 * (q0.val * q3.val + q1.val * q2.val) , 1 - 2.0* (q2.val * q2.val + q3.val * q3.val)));

  if (yaw.val < 0){
    yaw.val +=360;
  }

}
void openIMU::GetPitch(void){
  rawPitch.val = ToDeg(asin(2.0 * (q0.val * q2.val - q3.val * q1.val)));
  pitch.val =  rawPitch.val - pitchOffset.val;
}

void openIMU::GetRoll(void){
  rawRoll.val = ToDeg(FastAtan2(2.0 * (q0.val * q1.val + q2.val * q3.val),1 - 2.0 * (q1.val * q1.val + q2.val * q2.val)));
  roll.val=  rawRoll.val - rollOffset.val;
}

void openIMU::GetYaw(void){
  yaw.val = ToDeg(FastAtan2(2.0 * (q0.val * q3.val + q1.val * q2.val) , 1 - 2.0* (q2.val * q2.val + q3.val * q3.val))) ;
  if (yaw.val < 0){
    yaw.val +=360;
  }
}













































