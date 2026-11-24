/*
  QT Py RP2040 + BNO085 (Wire1 / STEMMA QT, I2C 0x4A) -> raw IMU stream over USB serial.

  Output: one ASCII line per sensor event, nothing is paired or overwritten on the board.
    A,<seq>,<t_us>,<x>,<y>,<z>          raw accelerometer   [ADC counts]
    G,<seq>,<t_us>,<x>,<y>,<z>,<temp>   raw gyroscope       [ADC counts]
    C,<seq>,<t_us>,<x>,<y>,<z>          calibrated acc      [m/s^2]  (sensor's own processing)
    U,<seq>,<t_us>,<x>,<y>,<z>          uncalibrated gyro   [rad/s]  (sensor's own processing)
    [Info] ... / [ERR] ...              status messages
  seq  = per-report 8-bit sequence number from the sensor; a jump of >1 means a lost report.
  t_us = A/G: sample time from the sensor's own clock (uint32, wraps after ~71.6 min).
         C/U: board clock (micros()) when the event was received. The library's own event
              time is unusable over I2C. Pair C with A and U with G by sequence number.
  C/U have to be enabled anyway: the BNO08x only produces raw reports while the matching
  processed report is running. They serve as a scale reference for the raw counts.

  Changes vs. the original sketch (firmware/backup/sensor_read_out_code_original.ino.txt):
  - Events are handled in our own sh2 callback. Adafruit_BNO08x::getSensorEvent() keeps only the
    last event of each sh2_service() call, which silently dropped samples whenever the sensor
    delivered several reports at once.
  - No 5 s startup delay and no 5 s deadline for opening the port: the board waits until the
    host opens the port, and closing/reopening the port works without replugging.

  LED: orange/black = waiting for host, green/white = streaming,
       red+purple = sensor not found, red+yellow = setup failed, red+blue = sensor stopped responding
*/
#include <Adafruit_BNO08x.h>
#include <Adafruit_NeoPixel.h>
#include <Wire.h>

#define REPORT_INTERVAL_US 5000  // 200 Hz for all four reports

#define FIRMWARE_VERSION "imu_stream 2.0"

// COLOR CODES
#define ORANGE 0xffa000
#define BLACK 0x000000
#define PURPLE 0x800080
#define RED 0xff0000
#define GREEN 0x00ff00
#define BLUE 0x0000ff
#define YELLOW 0xffff00
#define WHITE 0xffffff

#define STATE_WAIT_HOST 1
#define STATE_SENSOR_INIT 2
#define STATE_SETUP 3
#define STATE_RUN 4
#define STATE_ERR -1

#define LED_BLINK_TIME 500 //millisec
#define MAX_TIME_SENSOR_INIT 2000 //millisec
#define MAX_TIME_SENSOR_READOUT 5000 //millisec
#define MAX_SENSOR_CONNECTION_RETRIES 5

#define SENSOR_INIT_ERROR 2
#define SENSOR_SETUP_ERROR 3
#define SENSOR_UNRESPONSIVE_ERROR 4

#define BNO08X_RESET -1

int state = STATE_WAIT_HOST;
int sensor_connection_retries = 0;
unsigned long time_sensor_init_start;
unsigned long time_last_event;
bool host_connected = false;

unsigned long time_led_color_change;
uint32_t led_color_1 = ORANGE;
uint32_t led_color_2 = BLACK;
char current_color = 0;

Adafruit_BNO08x bno08x(BNO08X_RESET);
Adafruit_NeoPixel led(1, PIN_NEOPIXEL, NEO_GRB + NEO_KHZ800);

void setup() {
  // LED setup
  #if defined(NEOPIXEL_POWER)
  pinMode(NEOPIXEL_POWER, OUTPUT);
  digitalWrite(NEOPIXEL_POWER, HIGH);
  #endif

  Wire1.setClock(400000);

  led.begin();
  led.setBrightness(50);

  Serial.begin(115200);  // baud rate is ignored for USB CDC
}

void led_update() {
  if ((millis() - time_led_color_change) > LED_BLINK_TIME) {
    led.fill(current_color == 0 ? led_color_2 : led_color_1);
    led.show();
    current_color = !current_color;
    time_led_color_change = millis();
  }
}

void throw_error(int reason) {
  state = STATE_ERR;
  led_color_1 = RED;
  Serial.print("[ERR] ");
  switch (reason) {
    case SENSOR_INIT_ERROR:
      Serial.println("Sensor connection failed!");
      led_color_2 = PURPLE;
      break;
    case SENSOR_SETUP_ERROR:
      Serial.println("Failed sensor setup!");
      led_color_2 = YELLOW;
      break;
    case SENSOR_UNRESPONSIVE_ERROR:
      Serial.println("No data from sensor is being received!");
      led_color_2 = BLUE;
      break;
  }
}

void begin_sensor_init() {
  time_sensor_init_start = millis();
  state = STATE_SENSOR_INIT;
}

void state_wait_host() {
  led_color_1 = ORANGE;
  led_color_2 = BLACK;
  if (Serial) {
    host_connected = true;
    Serial.println("[Info] Serial connection established.");
    Serial.println("[Info] Firmware: " FIRMWARE_VERSION);
    begin_sensor_init();
  }
}

void state_sensor_init() {
  if (!bno08x.begin_I2C(0x4A, &Wire1)) {
    if ((millis() - time_sensor_init_start) > MAX_TIME_SENSOR_INIT) {
      throw_error(SENSOR_INIT_ERROR);
    }
  } else {
    Serial.println("[Info] Sensor connection established.");
    state = STATE_SETUP;
  }
}

// Events are queued by the sh2 callback (which runs inside the library's I2C handling)
// and printed afterwards from the main loop.
struct Event {
  uint8_t id, seq;
  uint32_t t_us;
  int16_t x, y, z, temp;
  float fx, fy, fz;
};
#define QUEUE_LEN 64  // power of two
Event queue[QUEUE_LEN];
volatile uint32_t q_head = 0, q_tail = 0;
uint32_t q_overflows = 0;

// Called by sh2_service() once per sensor event, possibly several times per call.
static void on_sensor_event(void *cookie, sh2_SensorEvent_t *event) {
  sh2_SensorValue_t v;
  if (sh2_decodeSensorEvent(&v, event) != SH2_OK) {
    return;
  }
  time_last_event = millis();
  if (q_head - q_tail >= QUEUE_LEN) {
    q_overflows++;
    return;
  }
  Event &e = queue[q_head % QUEUE_LEN];
  e.id = v.sensorId;
  e.seq = v.sequence;
  switch (v.sensorId) {
    case SH2_RAW_ACCELEROMETER:
      e.t_us = v.un.rawAccelerometer.timestamp;
      e.x = v.un.rawAccelerometer.x; e.y = v.un.rawAccelerometer.y; e.z = v.un.rawAccelerometer.z;
      break;
    case SH2_RAW_GYROSCOPE:
      e.t_us = v.un.rawGyroscope.timestamp;
      e.x = v.un.rawGyroscope.x; e.y = v.un.rawGyroscope.y; e.z = v.un.rawGyroscope.z;
      e.temp = v.un.rawGyroscope.temperature;
      break;
    case SH2_ACCELEROMETER:
      e.t_us = micros();
      e.fx = v.un.accelerometer.x; e.fy = v.un.accelerometer.y; e.fz = v.un.accelerometer.z;
      break;
    case SH2_GYROSCOPE_UNCALIBRATED:
      e.t_us = micros();
      e.fx = v.un.gyroscopeUncal.x; e.fy = v.un.gyroscopeUncal.y; e.fz = v.un.gyroscopeUncal.z;
      break;
    default:
      return;
  }
  q_head++;
}

void print_queued_events() {
  char line[96];
  while (q_tail != q_head) {
    const Event &e = queue[q_tail % QUEUE_LEN];
    int n = 0;
    if (e.id == SH2_RAW_ACCELEROMETER) {
      n = snprintf(line, sizeof(line), "A,%u,%lu,%d,%d,%d\n", e.seq, (unsigned long)e.t_us,
                   e.x, e.y, e.z);
    } else if (e.id == SH2_RAW_GYROSCOPE) {
      n = snprintf(line, sizeof(line), "G,%u,%lu,%d,%d,%d,%d\n", e.seq, (unsigned long)e.t_us,
                   e.x, e.y, e.z, e.temp);
    } else if (e.id == SH2_ACCELEROMETER) {
      n = snprintf(line, sizeof(line), "C,%u,%lu,%.5f,%.5f,%.5f\n", e.seq, (unsigned long)e.t_us,
                   e.fx, e.fy, e.fz);
    } else if (e.id == SH2_GYROSCOPE_UNCALIBRATED) {
      n = snprintf(line, sizeof(line), "U,%u,%lu,%.6f,%.6f,%.6f\n", e.seq, (unsigned long)e.t_us,
                   e.fx, e.fy, e.fz);
    }
    q_tail++;
    if (n > 0 && host_connected) {
      Serial.write(line, n);
    }
  }
  if (q_overflows && host_connected) {
    Serial.print("[Info] Event queue overflow, lost events: ");
    Serial.println(q_overflows);
    q_overflows = 0;
  }
}

void state_setup() {
  // begin_I2C() registers the library's own callback; replace it with ours.
  sh2_setSensorCallback(on_sensor_event, NULL);

  // The processed reports must be enabled too, otherwise the raw reports stay silent.
  if (!bno08x.enableReport(SH2_ACCELEROMETER, REPORT_INTERVAL_US) ||
      !bno08x.enableReport(SH2_GYROSCOPE_UNCALIBRATED, REPORT_INTERVAL_US) ||
      !bno08x.enableReport(SH2_RAW_ACCELEROMETER, REPORT_INTERVAL_US) ||
      !bno08x.enableReport(SH2_RAW_GYROSCOPE, REPORT_INTERVAL_US)) {
    throw_error(SENSOR_SETUP_ERROR);
    return;
  }
  Serial.println("[Info] Sensor setup successful.");
  Serial.println("[Info] Sending data...");
  time_last_event = millis();
  state = STATE_RUN;
}

void state_run() {
  led_color_1 = GREEN;
  led_color_2 = WHITE;

  // Track the host: stop writing while the port is closed, resume when it is reopened.
  bool connected = Serial;
  if (connected && !host_connected) {
    Serial.println("[Info] Serial connection established.");
    Serial.println("[Info] Firmware: " FIRMWARE_VERSION);
    Serial.println("[Info] Sending data...");
  }
  host_connected = connected;
  if (!connected) {
    led_color_1 = ORANGE;
  }

  if (bno08x.wasReset()) {
    if (host_connected) {
      Serial.println("[Info] Sensor reset, re-enabling reports.");
    }
    state = STATE_SETUP;
    return;
  }

  sh2_service();  // calls on_sensor_event() for every pending event
  print_queued_events();

  if ((millis() - time_last_event) > MAX_TIME_SENSOR_READOUT) {
    if (sensor_connection_retries < MAX_SENSOR_CONNECTION_RETRIES) {
      sensor_connection_retries += 1;
      if (host_connected) {
        Serial.println("[Info] No sensor data, reconnecting.");
      }
      begin_sensor_init();
    } else {
      throw_error(SENSOR_UNRESPONSIVE_ERROR);
    }
  }
}

void loop() {
  while (1) {
    led_update();
    switch (state) {
      case STATE_WAIT_HOST:
        state_wait_host();
        break;
      case STATE_SENSOR_INIT:
        state_sensor_init();
        break;
      case STATE_SETUP:
        state_setup();
        break;
      case STATE_RUN:
        state_run();
        break;
      default:
        break;  // error state: only the LED blinks
    }
  }
}
