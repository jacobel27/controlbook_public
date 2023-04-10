#include "arduino_setup.h"
#include "time_utilities.h"
#include "sensor_utilities.h"
#include "motor_utilities.h"
#include "reference_utilities.h"
#include "ctrlPID.h"

//=============================================================================
// declare global structures
//=============================================================================

// instantiate structures
Reference reference;

// instantiate classes as global variables
TimeUtilities timing;
SensorUtilities sensors;
MotorUtilities rotors;
ReferenceUtilities signal_generator;
ReferenceUtilities signal2;
CtrlPID controller;

//=============================================================================
// arduino setup function (runs once at start of simulation)
//=============================================================================
void setup()
{
  // start serial communication
  Serial.begin(9600);

  // initialize all classes
  timing.init();  // initialize current time and sample rate
  sensors.init();  // initialize sensors
  controller.init();  // initialize controller
  signal_generator.init(15*3.14/180, 0.01, 0.0);
  signal2.init(15*3.14/180, 0.015, 0.0);

  // initialize buttons on breakout board and watchdog timers
  initialize_buttons();
}

//=============================================================================
// arduino loop function (loops forever)
//=============================================================================
void loop()
{
  timing.update();  // update current time and sample rate
  sensors.update();  // update sensors
  //float psi_ref = signal_generator.square_signal(timing.current);
  //reference.psi = signal_generator.joystick_yaw();
  //reference.theta = signal_generator.joystick_pitch();
  reference.psi = signal_generator.square_signal(timing.current);
  reference.theta = signal2.square_signal(timing.current);

  controller.update(reference, sensors, rotors, timing.Ts);  // update controller

  // zero encoders if zero button pushed
  zeroButton.update();  
  if ( zeroButton.pressed() ) sensors.zero(); 
   
  // calibrate rotors if calibrate button pushed
  calButton.update();
  if ( calButton.pressed() ) rotors.init();
  
  // reset watchdog timer
  wdt_reset();
  digitalWrite(LED_RX, LOW);
}
