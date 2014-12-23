/*
void SonarInit(){

  DDRB |= (1<<PB5);
  TCCR1A = (1<<WGM11)|(1<<COM1A1);
  TCCR1B = (1<<WGM13)|(1<<WGM12)|(1<<CS11)|(1<<CS10);
  ICR1 = PERIOD_TRIG;  
  OCR1A = 2; 

  DDRB &= ~(1<<PB6);
  PORTK |= (1<<PB6);
  PCMSK0 |= 1<<PCINT6;
  PCICR |= 1<<0;

}*/



ISR(PCINT0_vect){

  if (((PINB & 1<<PB6)>>PB6) == 1){
    start = micros();

  }
  else{
    width = (micros() - start);
    if (width <= 17400){
      newPing = true;
      baroCorrect = false;
    }
    else{
      baroCorrect = true;
    }

  }

}


void GetInitialQuat(){
  float inertialSumX,inertialSumY,inertialSumZ;
  imu.InitialQuat();



  imu.kpAcc = kp_waypoint_position.val;
  imu.kiAcc = ki_waypoint_position.val;
  imu.kpMag = kd_waypoint_position.val;
  imu.kiMag = fc_waypoint_position.val;
  //imu.FEEDBACK_LIMIT = 0.5;
  imu.DECLINATION = ToRad(fc_cross_track.val);
  imu.COS_DEC = cos(imu.DECLINATION);
  imu.SIN_DEC = sin(imu.DECLINATION);
  //Serial<<imu.kpAcc<<","<<imu.kiAcc<<","<<imu.kpMag<<","<<imu.kiMag<<","<<imu.FEEDBACK_LIMIT<<","<<imu.DECLINATION<<"\r\n";
  imu.GenerateRotationMatrix();
  imu.GetGravOffset();
  imuTimer = micros();
  for (int i = 0; i < 500; i++){
    while(micros() - imuTimer < 13333){
    }

    imuDT = (micros() - imuTimer ) * 0.000001;
    imuTimer = micros();
    GetMag();
    GetAcc();
    GetGyro();

    accToFilterX = -1.0 * filtAccX.val;//if the value from the smoothing filter is sent it will not work when the algorithm normalizes the vector
    accToFilterY = -1.0 * filtAccY.val;
    accToFilterZ = -1.0 * filtAccZ.val;
    imu.AHRSupdate();
    imu.GenerateRotationMatrix();
    imu.GetEuler();
    //Serial<<i<<","<<imu.pitch.val<<","<<imu.roll.val<<","<<imu.yaw.val<<"\r\n";
  }

  imu.initialAccMagnitude.val = 0;
  inertialSumZ = 0;
  for(int i = 0; i < 200; i ++){
    GetAcc();
    imu.GetInertial();

    inertialSumZ += imu.inertialZ.val;
    delayMicroseconds(2500);
  }

  imu.initialAccMagnitude.val = inertialSumZ / 200.0;
  //imu.initialAccMagnitude.val = - 9.81;

  imu.accelBiasX.val = 0;
  imu.accelBiasY.val = 0;
  imu.accelBiasZ.val = 0;

  imu.currentEstIndex = (uint8_t)kp_waypoint_velocity.val;
  imu.lagIndex = 0;

}
void CalibrateSensors(){

  generalPurposeTimer = millis();
  while(1){
    if ( millis() - generalPurposeTimer >= 10){
      generalPurposeTimer = millis();
      /*
      SPI.setDataMode(SPI_MODE3);
       AccSSLow();
       SPI.transfer(DATAX0 | READ | MULTI);
       accX.buffer[0] = SPI.transfer(0x00);
       accX.buffer[1] = SPI.transfer(0x00);
       accY.buffer[0] = SPI.transfer(0x00);
       accY.buffer[1] = SPI.transfer(0x00);
       accZ.buffer[0] = SPI.transfer(0x00);
       accZ.buffer[1] = SPI.transfer(0x00);
       AccSSHigh();  
       accY.val *= -1;
       accZ.val *= -1;*/


      AccSSLow();
      SPI.transfer(OUT_X_L_A | READ | MULTI);
      accX.buffer[0] = SPI.transfer(0x00);
      accX.buffer[1] = SPI.transfer(0x00);
      accY.buffer[0] = SPI.transfer(0x00);
      accY.buffer[1] = SPI.transfer(0x00);
      accZ.buffer[0] = SPI.transfer(0x00);
      accZ.buffer[1] = SPI.transfer(0x00);
      accX.val = accX.val>>4;
      accY.val = accY.val>>4;
      accZ.val = accZ.val>>4;
      AccSSHigh();
      tempX = accX.val *  0.7071067 + accY.val * 0.7071067;
      tempY = accX.val * -0.7071067 + accY.val * 0.7071067;
      accX.val = tempX;
      accY.val = tempY;

      I2c.read(MAG_ADDRESS,LSM303_OUT_X_H,6);
      magX.buffer[1] = I2c.receive();//X
      magX.buffer[0] = I2c.receive();
      magZ.buffer[1] = I2c.receive();//Z
      magZ.buffer[0] = I2c.receive();
      magY.buffer[1] = I2c.receive();//Y
      magY.buffer[0] = I2c.receive();
      //magY.val *= -1;
      //magZ.val *= -1;

      PollPressure();
      if (newBaro == true){
        newBaro = false;
      } 
      if (sendCalibrationData == true){
        SendCalData();
      }
    }
    Radio();
  }
}

void SendCalData(){
  int16_u temp;

  txSum = 0;
  txDoubleSum = 0;
  radioPrint->write(0xAA);
  switch(cmdNum){
  case 0:
    radioPrint->write(7);


    radioPrint->write((uint8_t)0x00);
    txSum += 0;
    txDoubleSum += txSum;

    radioPrint->write(magX.buffer[0]);
    txSum += magX.buffer[0];
    txDoubleSum += txSum;

    radioPrint->write(magX.buffer[1]);
    txSum += magX.buffer[1];
    txDoubleSum += txSum;

    radioPrint->write(magY.buffer[0]);
    txSum += magY.buffer[0];
    txDoubleSum += txSum;

    radioPrint->write(magY.buffer[1]);
    txSum += magY.buffer[1];
    txDoubleSum += txSum;

    radioPrint->write(magZ.buffer[0]);
    txSum += magZ.buffer[0];
    txDoubleSum += txSum;

    radioPrint->write(magZ.buffer[1]);
    txSum += magZ.buffer[1];
    txDoubleSum += txSum;

    radioPrint->write(txSum);
    radioPrint->write(txDoubleSum);

    break;
  case 1:
    radioPrint->write(7);


    radioPrint->write(1);
    txSum += 1;
    txDoubleSum += txSum;

    radioPrint->write(accX.buffer[0]);
    txSum += accX.buffer[0];
    txDoubleSum += txSum;

    radioPrint->write(accX.buffer[1]);
    txSum += accX.buffer[1];
    txDoubleSum += txSum;

    radioPrint->write(accY.buffer[0]);
    txSum += accY.buffer[0];
    txDoubleSum += txSum;

    radioPrint->write(accY.buffer[1]);
    txSum += accY.buffer[1];
    txDoubleSum += txSum;

    radioPrint->write(accZ.buffer[0]);
    txSum += accZ.buffer[0];
    txDoubleSum += txSum;

    radioPrint->write(accZ.buffer[1]);
    txSum += accZ.buffer[1];
    txDoubleSum += txSum;

    radioPrint->write(txSum);
    radioPrint->write(txDoubleSum);
    break;
  case 2:
    radioPrint->write(17);

    radioPrint->write(2);
    txSum += 2;
    txDoubleSum += txSum;

    temp.val = rawRCVal[THRO];
    radioPrint->write(temp.buffer[0]);
    txSum += temp.buffer[0];
    txDoubleSum += txSum;
    radioPrint->write(temp.buffer[1]);
    txSum += temp.buffer[1];
    txDoubleSum += txSum;

    temp.val = rawRCVal[AILE];
    radioPrint->write(temp.buffer[0]);
    txSum += temp.buffer[0];
    txDoubleSum += txSum;
    radioPrint->write(temp.buffer[1]);
    txSum += temp.buffer[1];
    txDoubleSum += txSum;

    temp.val = rawRCVal[ELEV];
    radioPrint->write(temp.buffer[0]);
    txSum += temp.buffer[0];
    txDoubleSum += txSum;
    radioPrint->write(temp.buffer[1]);
    txSum += temp.buffer[1];
    txDoubleSum += txSum;

    temp.val = rawRCVal[RUDD];
    radioPrint->write(temp.buffer[0]);
    txSum += temp.buffer[0];
    txDoubleSum += txSum;
    radioPrint->write(temp.buffer[1]);
    txSum += temp.buffer[1];
    txDoubleSum += txSum;

    temp.val = rawRCVal[GEAR];
    radioPrint->write(temp.buffer[0]);
    txSum += temp.buffer[0];
    txDoubleSum += txSum;
    radioPrint->write(temp.buffer[1]);
    txSum += temp.buffer[1];
    txDoubleSum += txSum;

    temp.val = rawRCVal[AUX1];
    radioPrint->write(temp.buffer[0]);
    txSum += temp.buffer[0];
    txDoubleSum += txSum;
    radioPrint->write(temp.buffer[1]);
    txSum += temp.buffer[1];
    txDoubleSum += txSum;

    temp.val = rawRCVal[AUX2];
    radioPrint->write(temp.buffer[0]);
    txSum += temp.buffer[0];
    txDoubleSum += txSum;
    radioPrint->write(temp.buffer[1]);
    txSum += temp.buffer[1];
    txDoubleSum += txSum;

    temp.val = rawRCVal[AUX3];
    radioPrint->write(temp.buffer[0]);
    txSum += temp.buffer[0];
    txDoubleSum += txSum;
    radioPrint->write(temp.buffer[1]);
    txSum += temp.buffer[1];
    txDoubleSum += txSum;

    radioPrint->write(txSum);
    radioPrint->write(txDoubleSum);
    break;
  }
}

void StartUpAHRSRun(){
  if ( micros() - imuTimer >= 13333){  
    imuDT = (micros() - imuTimer) * 0.000001;
    imuTimer = micros();
    GetAcc();
    GetMag();
    GetGyro();
    PollPressure();
    if (newBaro == true){
      newBaro = false;
    } 
    imu.AHRSupdate();
  } 
}



void GPSStart(){
  uint8_t LEDState;
  gps.init();


  generalPurposeTimer = millis();
  while ((millis() - generalPurposeTimer < 1000) && (gps.newData == false)){
    gps.Monitor();
    if (gps.newData == true){
      GPSDetected = true;
    }
  }
  //Serial<<"gps det: "<<GPSDetected<<"\r\n";
  //to do add feed back with leds
  if (GPSDetected == true){
    gpsFailSafe = false;
    while (gps.data.vars.gpsFix != 3){

      gps.Monitor();
      if (millis() - generalPurposeTimer > 500){
        generalPurposeTimer = millis();
        LEDState++;
        if (LEDState == 4){
          LEDState = 0;
        }
      }
      switch (LEDState){
      case 0:
        digitalWrite(13,HIGH);
        digitalWrite(RED,LOW);
        digitalWrite(YELLOW,LOW);
        digitalWrite(GREEN,LOW);
        break;
      case 1:
        digitalWrite(13,LOW);
        digitalWrite(RED,HIGH);
        digitalWrite(YELLOW,LOW);
        digitalWrite(GREEN,LOW);
        break;
      case 2:
        digitalWrite(13,LOW);
        digitalWrite(RED,LOW);
        digitalWrite(YELLOW,HIGH);
        digitalWrite(GREEN,LOW);
        break;
      case 3:
        digitalWrite(13,LOW);
        digitalWrite(RED,LOW);
        digitalWrite(YELLOW,LOW);
        digitalWrite(GREEN,HIGH);
        break;
      }
    }
    while(gps.data.vars.hAcc > 500000){
      gps.Monitor();
      if (millis() - generalPurposeTimer > 500){
        generalPurposeTimer = millis();
        LEDState++;
        if (LEDState == 4){
          LEDState = 0;
        }
      }
      switch (LEDState){
      case 0:
        digitalWrite(13,HIGH);
        digitalWrite(RED,LOW);
        digitalWrite(YELLOW,LOW);
        digitalWrite(GREEN,LOW);
        break;
      case 1:
        digitalWrite(13,HIGH);
        digitalWrite(RED,HIGH);
        digitalWrite(YELLOW,LOW);
        digitalWrite(GREEN,LOW);
        break;
      case 2:
        digitalWrite(13,HIGH);
        digitalWrite(RED,LOW);
        digitalWrite(YELLOW,HIGH);
        digitalWrite(GREEN,LOW);
        break;
      case 3:
        digitalWrite(13,HIGH);
        digitalWrite(RED,LOW);
        digitalWrite(YELLOW,LOW);
        digitalWrite(GREEN,HIGH);
        break;
      }
    }

    gps.newData = false;
    while (gps.newData == false){
      gps.Monitor();
    }
    homeBase.lat.val = gps.data.vars.lat;
    homeBase.lon.val = gps.data.vars.lon;
  }  



}
void GetAltitude(float *press,float *pressInit, float *alti){

  float pressureRatio =  *press /  *pressInit;
  *alti = (1.0f - pow(pressureRatio, 0.190295f)) * 44330.0f;
}

void PollPressure(void){
  /*  if (millis() - baroPollTimer > POLL_RATE){
   switch (pressureState){
   case 0://read ut
   StartUT();
   pressureState = 1;
   baroTimer = millis();
   break;
   case 1://wait for ready signal
   if (millis() - baroTimer > 5){
   pressureState = 2;
   ut = ReadUT();
   StartUP();
   baroTimer = millis();
   }
   
   break;
   case 2://read up
   if (millis() - baroTimer > CONV_TIME){
   up = ReadUP();
   temperature = Temperature(ut);
   pressure.val = Pressure(up);
   pressureState = 0;
   newBaro = true;
   baroPollTimer = millis();
   }
   break;
   
   }
   }*/


  if (millis() - baroPollTimer >= BARO_CONV_TIME){
    switch(baroState){
    case 0://start temp conv
      BaroSSLow();
      SPI.transfer(CONVERT_D2_OSR4096);
      BaroSSHigh();
      baroState = 1;
      baroDelayTimer = millis();
      break;
    case 1:
      if (millis() - baroDelayTimer >= 10){
        BaroSSLow();
        SPI.transfer(ADC_READ);
        D2.buffer[2] = SPI.transfer(0x00);
        D2.buffer[1] = SPI.transfer(0x00);
        D2.buffer[0] = SPI.transfer(0x00);
        BaroSSHigh();
        baroState = 2;
      }
      break;
    case 2:
      BaroSSLow();
      SPI.transfer(CONVERT_D1_OSR4096);
      BaroSSHigh();
      baroState = 3;
      baroDelayTimer = millis();
      break;
    case 3:
      if (millis() - baroDelayTimer >= 10){
        BaroSSLow();
        SPI.transfer(ADC_READ);
        D1.buffer[2] = SPI.transfer(0x00);
        D1.buffer[1] = SPI.transfer(0x00);
        D1.buffer[0] = SPI.transfer(0x00);
        BaroSSHigh();
        baroState = 0;
        baroPollTimer = millis();
        GetBaro();
        newBaro = true;
      }
      break;
    }
  }  
}
void GetBaro(){


  dT = D2.val-(((uint32_t)C5.val)<<8);
  TEMP = (dT * C6.val)/8388608;
  OFF = C2.val * 65536.0 + (C4.val * dT) / 128;
  SENS = C1.val * 32768.0 + (C3.val * dT) / 256;

  if (TEMP < 0) {
    // second order temperature compensation when under 20 degrees C
    float T2 = (dT*dT) / 0x80000000;
    float Aux = TEMP*TEMP;
    float OFF2 = 2.5*Aux;
    float SENS2 = 1.25*Aux;
    TEMP = TEMP - T2;
    OFF = OFF - OFF2;
    SENS = SENS - SENS2;
  }

  P = (D1.val*SENS/2097152 - OFF)/32768;
  temperature = TEMP + 2000;
  pressure.val = P;

}

void BaroInit(void){
  /* pressureState = 0;
   newBaro = false;
   I2c.read(BMP085_ADDRESS,0xAA,22);
   msb = I2c.receive();
   lsb = I2c.receive();
   ac1 = (msb << 8) | lsb;
   
   msb = I2c.receive();
   lsb = I2c.receive();
   ac2 = (msb << 8) | lsb;
   
   msb = I2c.receive();
   lsb = I2c.receive();
   ac3 = (msb << 8) | lsb;
   
   msb = I2c.receive();
   lsb = I2c.receive();
   ac4 = (msb << 8) | lsb;
   
   msb = I2c.receive();
   lsb = I2c.receive();
   ac5 = (msb << 8) | lsb;
   
   msb = I2c.receive();
   lsb = I2c.receive();
   ac6 = (msb << 8) | lsb;
   
   msb = I2c.receive();
   lsb = I2c.receive();
   b1 = (msb << 8) | lsb;
   
   msb = I2c.receive();
   lsb = I2c.receive();
   b2 = (msb << 8) | lsb;
   
   msb = I2c.receive();
   lsb = I2c.receive();
   mb = (msb << 8) | lsb;
   
   msb = I2c.receive();
   lsb = I2c.receive();
   mc = (msb << 8) | lsb;
   
   msb = I2c.receive();
   lsb = I2c.receive();
   md = (msb << 8) | lsb;*/

  BaroSSLow();
  SPI.transfer(MS5611_RESET);
  BaroSSHigh();
  delay(5);
  BaroSSLow();
  SPI.transfer(MS5611_PROM_C1);
  C1.buffer[1] = SPI.transfer(0x00);
  C1.buffer[0] = SPI.transfer(0x00);
  BaroSSHigh();
  delay(1);
  BaroSSLow();
  SPI.transfer(MS5611_PROM_C2);
  C2.buffer[1] = SPI.transfer(0x00);
  C2.buffer[0] = SPI.transfer(0x00);
  BaroSSHigh();
  delay(1);
  BaroSSLow();
  SPI.transfer(MS5611_PROM_C3);
  C3.buffer[1] = SPI.transfer(0x00);
  C3.buffer[0] = SPI.transfer(0x00);
  BaroSSHigh();
  delay(1);
  BaroSSLow();
  SPI.transfer(MS5611_PROM_C4);
  C4.buffer[1] = SPI.transfer(0x00);
  C4.buffer[0] = SPI.transfer(0x00);
  BaroSSHigh();
  delay(1);
  BaroSSLow();
  SPI.transfer(MS5611_PROM_C5);
  C5.buffer[1] = SPI.transfer(0x00);
  C5.buffer[0] = SPI.transfer(0x00);
  BaroSSHigh();
  delay(1);
  BaroSSLow();
  SPI.transfer(MS5611_PROM_C6);
  C6.buffer[1] = SPI.transfer(0x00);
  C6.buffer[0] = SPI.transfer(0x00);
  BaroSSHigh();
  //this is to get the ground pressure for relative altitude
  //lower pressure than this means positive altitude
  //higher pressure than this means negative altitude
  baroCount = 0;
  baroSum = 0;
  while (baroCount < 10){//use a while instead of a for loop because the for loop runs too fast
    PollPressure();
    if (newBaro == true){
      newBaro = false;
      baroCount++;
      baroSum += pressure.val;
    }    
  }
  pressureInitial = baroSum / 10;   


  initialTemp.val = temperature;

  //use the line below for altitdue above sea level
  //pressureInitial = 101325;

}


void MagInit(){
  //continous conversion 220Hz
  I2c.write((uint8_t)MAG_ADDRESS,(uint8_t)LSM303_CRA_REG,(uint8_t)0x18);
  I2c.write((uint8_t)MAG_ADDRESS,(uint8_t)LSM303_CRB_REG,(uint8_t)0x60);
  I2c.write((uint8_t)MAG_ADDRESS,(uint8_t)LSM303_MR_REG,(uint8_t)0x00);
  I2c.read(MAG_ADDRESS,LSM303_OUT_X_H,6);
  magX.buffer[1] = I2c.receive();//X
  magX.buffer[0] = I2c.receive();
  magZ.buffer[1] = I2c.receive();//Z
  magZ.buffer[0] = I2c.receive();
  magY.buffer[1] = I2c.receive();//Y
  magY.buffer[0] = I2c.receive();
  //magY.val *= -1;
  //magZ.val *= -1;
  shiftedMagX  = magX.val - magOffSetX;
  shiftedMagY  = magY.val - magOffSetY;
  shiftedMagZ  = magZ.val - magOffSetZ;
  scaledMagX = magWInv00 * shiftedMagX + magWInv01 * shiftedMagY + magWInv02 * shiftedMagZ;
  scaledMagY = magWInv10 * shiftedMagX + magWInv11 * shiftedMagY + magWInv12 * shiftedMagZ;
  scaledMagZ = magWInv20 * shiftedMagX + magWInv21 * shiftedMagY + magWInv22 * shiftedMagZ;
  calibMagX.val = scaledMagX;
  calibMagY.val = scaledMagY;
  calibMagZ.val = scaledMagZ;
  for (uint8_t i = 0; i < 100; i++){
    GetMag();
    delay(5);
  }

}

void AccInit(){
  /*
  SPI.setDataMode(SPI_MODE3);
   
   AccSSLow();
   SPI.transfer(WRITE | SINGLE | BW_RATE);
   SPI.transfer(0x0C);
   AccSSHigh();
   
   AccSSLow();
   SPI.transfer(WRITE | SINGLE | POWER_CTL);
   SPI.transfer(0x08);//start measurment
   AccSSHigh();
   
   AccSSLow();
   SPI.transfer(WRITE | SINGLE | DATA_FORMAT);
   SPI.transfer(0x08);//full resolution + / - 16g
   AccSSHigh();*/
  AccSSLow();
  SPI.transfer(CTRL_REG1_A | WRITE | SINGLE);
  SPI.transfer(0x77);//400Hz all axes enabled
  AccSSHigh();

  AccSSLow();
  SPI.transfer(CTRL_REG2_A | WRITE | SINGLE);
  SPI.transfer(0x00);//high pass filter not used
  AccSSHigh();

  AccSSLow();
  SPI.transfer(CTRL_REG3_A | WRITE | SINGLE);
  SPI.transfer(0x00);//not using interrupts for polling sensor
  AccSSHigh();

  AccSSLow();
  SPI.transfer(CTRL_REG4_A | WRITE | SINGLE);
  SPI.transfer(0x18);//little endian
  AccSSHigh();

  AccSSLow();
  SPI.transfer(CTRL_REG5_A | WRITE | SINGLE);
  SPI.transfer(0x00);//not using interrupts for polling sensor
  AccSSHigh();

  AccSSLow();
  SPI.transfer(CTRL_REG6_A | WRITE | SINGLE);
  SPI.transfer(0x00);//not using interrupts for polling sensor
  AccSSHigh();


  GetAcc();

  if (accX.val > 0){
    scaledAccX.val = accX.val * accXScalePos;
  }
  else{
    scaledAccX.val = accX.val * accXScaleNeg;
  }
  if (accY.val > 0){
    scaledAccY.val = accY.val  * accYScalePos;
  }
  else{
    scaledAccY.val = accY.val  * accYScaleNeg;
  }
  if (accZ.val > 0){
    scaledAccZ.val = accZ.val * accZScalePos;
  }
  else{
    scaledAccZ.val = accZ.val * accZScaleNeg;
  }
  filtAccX.val = scaledAccX.val;
  filtAccY.val = scaledAccY.val;
  filtAccZ.val = scaledAccZ.val; 


  for (uint16_t i = 0;i < 100; i++){
    GetAcc();
    delayMicroseconds(2500);
  }

}

void GyroInit(){
  SPI.setDataMode(SPI_MODE0);
  GyroSSLow();
  SPI.transfer(L3G_CTRL_REG2 | WRITE | SINGLE);
  SPI.transfer(0x00); //high pass filter disabled
  GyroSSHigh();

  GyroSSLow();
  SPI.transfer(L3G_CTRL_REG3 | WRITE | SINGLE);
  SPI.transfer(0x00); //not using interrupts
  GyroSSHigh();

  GyroSSLow();
  SPI.transfer(L3G_CTRL_REG4 | WRITE | SINGLE);
  SPI.transfer(0x20); //2000dps scale
  GyroSSHigh();

  GyroSSLow();
  SPI.transfer(L3G_CTRL_REG5 | WRITE | SINGLE);
  SPI.transfer(0x02); //out select to use the second LPF
  //not using HPF or interrupts
  GyroSSHigh();

  GyroSSLow();
  SPI.transfer(L3G_CTRL_REG1 | WRITE | SINGLE);
  SPI.transfer(0x8F);
  GyroSSHigh();
  //this section takes an average of 500 samples to calculate the offset
  //if this step is skipped the IMU will still work, but this simple step gives better results
  gyroOffsetX = 0;
  gyroOffsetY = 0;
  gyroOffsetZ = 0;
  gyroSumX = 0;
  gyroSumY = 0;
  gyroSumZ = 0;
  for (uint16_t j = 0; j < 500; j ++){
    GyroSSLow();
    SPI.transfer(L3G_OUT_X_L  | READ | MULTI);
    gyroX.buffer[0] = SPI.transfer(0x00);
    gyroX.buffer[1] = SPI.transfer(0x00);
    gyroY.buffer[0] = SPI.transfer(0x00);
    gyroY.buffer[1] = SPI.transfer(0x00);
    gyroZ.buffer[0] = SPI.transfer(0x00);
    gyroZ.buffer[1] = SPI.transfer(0x00);

    GyroSSHigh();
    delay(3);
  }
  for (uint16_t j = 0; j < 50; j ++){
    GyroSSLow();
    SPI.transfer(L3G_OUT_X_L  | READ | MULTI);
    gyroX.buffer[0] = SPI.transfer(0x00);
    gyroX.buffer[1] = SPI.transfer(0x00);
    gyroY.buffer[0] = SPI.transfer(0x00);
    gyroY.buffer[1] = SPI.transfer(0x00);
    gyroZ.buffer[0] = SPI.transfer(0x00);
    gyroZ.buffer[1] = SPI.transfer(0x00);

    GyroSSHigh();
    gyroSumX += gyroX.val;
    gyroSumY += (gyroY.val);
    gyroSumZ += (gyroZ.val);
    //Serial<<gyroX.val<<","<<gyroY.val<<","<<gyroZ.val<<"\r\n";
    
    delay(3);
  }
  gyroOffsetX = gyroSumX / 50;
  gyroOffsetY = gyroSumY / 50;
  gyroOffsetZ = gyroSumZ / 50;
  //Serial<<"++++++++++\r\n";
  //Serial<<gyroOffsetX<<","<<gyroOffsetY<<","<<gyroOffsetZ<<"\r\n";
  GetGyro();

}

void GetMag(){
  I2c.read(MAG_ADDRESS,LSM303_OUT_X_H,6);
  magX.buffer[1] = I2c.receive();//X
  magX.buffer[0] = I2c.receive();
  magZ.buffer[1] = I2c.receive();//Z
  magZ.buffer[0] = I2c.receive();
  magY.buffer[1] = I2c.receive();//Y
  magY.buffer[0] = I2c.receive();
  //magY.val *= -1;
  //magZ.val *= -1;
  //sign swap due to different compass
  //deltaTemp.val = temperature - calibTempMag.val;

  /* magX.val -= deltaTemp.val * xSlopeMag.val;
   magY.val += deltaTemp.val * ySlopeMag.val;
   magZ.val += deltaTemp.val * zSlopeMag.val;*/

  shiftedMagX  = magX.val - magOffSetX;
  shiftedMagY  = magY.val - magOffSetY;
  shiftedMagZ  = magZ.val - magOffSetZ;
  scaledMagX = magWInv00 * shiftedMagX + magWInv01 * shiftedMagY + magWInv02 * shiftedMagZ;
  scaledMagY = magWInv10 * shiftedMagX + magWInv11 * shiftedMagY + magWInv12 * shiftedMagZ;
  scaledMagZ = magWInv20 * shiftedMagX + magWInv21 * shiftedMagY + magWInv22 * shiftedMagZ;
  calibMagX.val = calibMagX.val * 0.75 + scaledMagX * 0.25;
  calibMagY.val = calibMagY.val * 0.75 + scaledMagY * 0.25;
  calibMagZ.val = calibMagZ.val * 0.75 + scaledMagZ * 0.25;
  magToFiltX = calibMagX.val;
  magToFiltY = calibMagY.val;
  magToFiltZ = calibMagZ.val;
}
void UpdateOffset(){
   gyroOffsetX = 0;
  gyroOffsetY = 0;
  gyroOffsetZ = 0;
  gyroSumX = 0;
  gyroSumY = 0;
  gyroSumZ = 0;
 
  for (uint16_t j = 0; j < 50; j ++){
    GyroSSLow();
    SPI.transfer(L3G_OUT_X_L  | READ | MULTI);
    gyroX.buffer[0] = SPI.transfer(0x00);
    gyroX.buffer[1] = SPI.transfer(0x00);
    gyroY.buffer[0] = SPI.transfer(0x00);
    gyroY.buffer[1] = SPI.transfer(0x00);
    gyroZ.buffer[0] = SPI.transfer(0x00);
    gyroZ.buffer[1] = SPI.transfer(0x00);

    GyroSSHigh();
    gyroSumX += gyroX.val;
    gyroSumY += (gyroY.val);
    gyroSumZ += (gyroZ.val);
    //Serial<<gyroX.val<<","<<gyroY.val<<","<<gyroZ.val<<"\r\n";
    
    delay(3);
  }
  gyroOffsetX = gyroSumX / 50;
  gyroOffsetY = gyroSumY / 50;
  gyroOffsetZ = gyroSumZ / 50;
}
void GetGyro(){

  SPI.setDataMode(SPI_MODE0);
  GyroSSLow();
  SPI.transfer(L3G_OUT_X_L  | READ | MULTI);
  gyroX.buffer[0] = SPI.transfer(0x00);
  gyroX.buffer[1] = SPI.transfer(0x00);
  gyroY.buffer[0] = SPI.transfer(0x00);
  gyroY.buffer[1] = SPI.transfer(0x00);
  gyroZ.buffer[0] = SPI.transfer(0x00);
  gyroZ.buffer[1] = SPI.transfer(0x00);
  GyroSSHigh();


  //deltaTemp.val = temperature - initialTemp.val;
  /* deltaTemp.val = 0;
   
   
   gyroX.val -= ( (deltaTemp.val * xSlopeGyro.val) + gyroOffsetX); 
   gyroY.val -= ( (deltaTemp.val * ySlopeGyro.val) + gyroOffsetY);
   gyroZ.val -= ( (deltaTemp.val * zSlopeGyro.val) + gyroOffsetZ);*/


  gyroX.val -= gyroOffsetX; 
  gyroY.val -= gyroOffsetY;
  gyroZ.val -= gyroOffsetZ;
  tempX = gyroX.val *  0.7071067 + gyroY.val * 0.7071067;
  tempY = gyroX.val * -0.7071067 + gyroY.val * 0.7071067;
  gyroX.val = tempX;
  gyroY.val = tempY;
  //Serial<<gyroX.val<<","<<gyroY.val<<","<<gyroZ.val<<"\r\n";
  degreeGyroX.val = gyroX.val * 0.07;
  degreeGyroY.val = gyroY.val * 0.07;
  degreeGyroZ.val = gyroZ.val * 0.07;

  radianGyroX = ToRad(degreeGyroX.val);
  radianGyroY = ToRad(degreeGyroY.val);
  radianGyroZ = ToRad(degreeGyroZ.val);


}
void GetAcc(){
  AccSSLow();
  SPI.transfer(OUT_X_L_A | READ | MULTI);
  accX.buffer[0] = SPI.transfer(0x00);
  accX.buffer[1] = SPI.transfer(0x00);
  accY.buffer[0] = SPI.transfer(0x00);
  accY.buffer[1] = SPI.transfer(0x00);
  accZ.buffer[0] = SPI.transfer(0x00);
  accZ.buffer[1] = SPI.transfer(0x00);

  AccSSHigh();
  accX.val = accX.val>>4;
  accY.val = accY.val>>4;
  accZ.val = accZ.val>>4;

  tempX = accX.val *  0.7071067 + accY.val * 0.7071067;
  tempY = accX.val * -0.7071067 + accY.val * 0.7071067;
  accX.val = tempX;
  accY.val = tempY;

  //shiftedAccX.val = sqrt((float)accX.val * (float)accX.val + (float)accY.val * (float)accY.val + (float)accZ.val * (float)accZ.val);
  //deltaTemp.val = temperature - calibTempAcc.val;
  //deltaTemp.val =0;

  /*  accXScalePos = 9.8 / kp_pitch_rate.val ;
   accYScalePos = 9.8 / ki_pitch_rate.val ;
   accZScalePos = 9.8 / kd_pitch_rate.val ;
   
   accXScaleNeg = 9.8 / fc_pitch_rate.val ;
   accYScaleNeg = 9.8 / kp_roll_rate.val ;
   accZScaleNeg = 9.8 / ki_roll_rate.val ;*/


  /*
  accX.val -= (xAccOffset.val + (deltaTemp.val) * xSlopeAcc.val);
   accY.val -= (yAccOffset.val + (deltaTemp.val) * ySlopeAcc.val);
   accZ.val -= (deltaTemp.val) * zSlopeAcc.val;*/
  //accX.val -= (deltaTemp.val * xSlopeAcc.val);
  //accY.val -= (deltaTemp.val * ySlopeAcc.val);
  //accZ.val -= (deltaTemp.val * zSlopeAcc.val);

  if (accX.val > 0){
    scaledAccX.val = accX.val * accXScalePos;
    //scaledAccX.val = accX.val * 0.019232445;
  }
  else{
    scaledAccX.val = accX.val * accXScaleNeg;
    //scaledAccX.val = accX.val * 0.018987648;
  }
  if (accY.val > 0){
    scaledAccY.val = accY.val * accYScalePos;
    //scaledAccY.val = accY.val * 0.01898699;
  }
  else{
    scaledAccY.val = accY.val * accYScaleNeg;
    //scaledAccY.val = accY.val * 0.019236658;
  }
  if (accZ.val > 0){
    scaledAccZ.val = accZ.val * accZScalePos;
    //scaledAccZ.val = accZ.val * 0.019563049;
  }
  else{
    scaledAccZ.val = accZ.val * accZScaleNeg;
    //scaledAccZ.val = accZ.val * 0.018522824;
  } 

  /*shiftedAccX.val  = accX.val - -3.915546;
   shiftedAccY.val  = accY.val - -5.206482;
   shiftedAccZ.val  = accZ.val - -13.800706;
   scaledAccX.val = 0.019069 * shiftedAccX.val + 0.000018 * shiftedAccY.val + -0.000029 * shiftedAccZ.val;
   scaledAccY.val = 0.000018 * shiftedAccX.val + 0.019152 * shiftedAccY.val + 0.000017 * shiftedAccZ.val;
   scaledAccZ.val = -0.000029 * shiftedAccX.val + 0.000017 * shiftedAccY.val + 0.018951 * shiftedAccZ.val; */
  /*  shiftedAccX.val  = accX.val - -3.28;
   shiftedAccY.val  = accY.val - -4.79;
   shiftedAccZ.val  = accZ.val - -13.76;
   scaledAccX.val = shiftedAccX.val * 0.019086;
   scaledAccY.val = shiftedAccY.val * 0.019138;
   scaledAccZ.val = shiftedAccZ.val * 0.018958;*/


  //Filter(&scaledAccX,&scaledAccY,&scaledAccZ,&filtAccX.val,&filtAccY.val,&filtAccZ.val);
  filtAccX.val = filtAccX.val * 0.9 + scaledAccX.val * 0.1;
  filtAccY.val = filtAccY.val * 0.9 + scaledAccY.val * 0.1;
  filtAccZ.val = filtAccZ.val * 0.9 + scaledAccZ.val * 0.1;

}

void WaitForTempStab(){
  boolean stabTemp = false;
  boolean tog;
  uint8_t tempState = 0;
  while (stabTemp == false){
    if (millis() - ledTimer >= 1000){
      ledTimer = millis();
      tog = ~tog;
      digitalWrite(RED,tog);
      digitalWrite(YELLOW,tog);
      digitalWrite(GREEN,tog);
      digitalWrite(13,tog);
    }
    PollPressure();
    if (newBaro == true){
      newBaro = false;
      switch(tempState){
      case 0://set final temperature
        initialTemp.val  = temperature;
        tempState = 1;
        generalPurposeTimer = millis();
        digitalWrite(RED,LOW);

        break;
      case 1://wait 
        if (abs(temperature - initialTemp.val ) > 200){
          generalPurposeTimer = millis();//reset timer if temp has changed by more than a degree
          initialTemp.val = temperature;
        }
        if (millis() - generalPurposeTimer > 3000){
          //if (millis() - generalPurposeTimer > 1){
          tempState = 2;
        }
        digitalWrite(YELLOW,LOW);

        break;
      case 2:
        if (abs(temperature - initialTemp.val ) <= 200){
          initialTemp.val  = temperature;
          stabTemp = true;
        }
        else{
          tempState = 0;
        }
        digitalWrite(GREEN,LOW);
        break;
      }
    }
  }
  baroCount = 0;
  baroSum = 0;
  while (baroCount < 10){//use a while instead of a for loop because the for loop runs too fast
    PollPressure();
    if (newBaro == true){
      newBaro = false;
      baroCount++;
      baroSum += pressure.val;
    }    
  }
  pressureInitial = baroSum / 10;    
}






























