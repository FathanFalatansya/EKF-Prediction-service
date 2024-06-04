#include <iostream>
#include <string>
#include <mosquitto.h>
#include <cmath>

#define MQTT_HOST "172.22.155.65"
#define MQTT_PORT 1883
#define MQTT_TOPIC "radar_data"

struct mosquitto *mosq = NULL;

// Definisi ukuran vektor status dan vektor pengukuran
const int state_size = 4;  // (x, y, V_a, chi_a)
const int measurement_size = 2;  // (x, y)

// Langkah waktu
const double dt = 1.0;  // Asumsikan langkah waktu 1 detik

// Matriks noise proses dan noise pengukuran
double Q[state_size][state_size] = {{0.1, 0, 0, 0},
                                    {0, 0.1, 0, 0},
                                    {0, 0, 0.1, 0},
                                    {0, 0, 0, 0.1}};
double R[measurement_size][measurement_size] = {{0.1, 0},
                                                {0, 0.1}};

// Vektor status
double x[state_size] = {0, 0, 10, 0};  // Inisialisasi awal (x, y, V_a, chi_a)

// Matriks kovariansi status
double P[state_size][state_size] = {{1, 0, 0, 0},
                                    {0, 1, 0, 0},
                                    {0, 0, 1, 0},
                                    {0, 0, 0, 1}};

// Matriks model transisi status
double F[state_size][state_size];

// Matriks model pengukuran
double H[measurement_size][state_size] = {{1, 0, 0, 0},
                                           {0, 1, 0, 0}};

// Vektor pengukuran
double z[measurement_size];

// Inisialisasi Extended Kalman Filter
void initializeEKF(double x[], double P[][state_size]) {
    // Tidak ada operasi spesifik yang diperlukan di sini
}

// Fungsi untuk menghitung turunan model dinamis
void f(const double x[], double at, double an, double x_pred[]) {
    x_pred[0] = x[0] + x[2] * cos(x[3]) * dt;
    x_pred[1] = x[1] + x[2] * sin(x[3]) * dt;
    x_pred[2] = x[2] + at * dt;
    x_pred[3] = x[3] + (an / x[2]) * dt;
}

// Fungsi untuk menghitung Jacobian dari model dinamis
void calculateJacobianF(const double x[], double an, double F[][state_size]) {
    F[0][0] = 1;
    F[0][1] = 0;
    F[0][2] = cos(x[3]) * dt;
    F[0][3] = -x[2] * sin(x[3]) * dt;
    F[1][0] = 0;
    F[1][1] = 1;
    F[1][2] = sin(x[3]) * dt;
    F[1][3] = x[2] * cos(x[3]) * dt;
    F[2][0] = 0;
    F[2][1] = 0;
    F[2][2] = 1;
    F[2][3] = 0;
    F[3][0] = 0;
    F[3][1] = 0;
    F[3][2] = -an / (x[2] * x[2]) * dt;
    F[3][3] = 1;
}

// Fungsi prediksi
void predict(double x[], double P[][state_size], double at, double an) {
    // Prediksi status
    double x_pred[state_size];
    f(x, at, an, x_pred);

    // Hitung Jacobian dari model transisi status
    double F[state_size][state_size];
    calculateJacobianF(x, an, F);

    // Prediksi kovariansi status
    double P_temp[state_size][state_size];
    for (int i = 0; i < state_size; ++i) {
        for (int j = 0; j < state_size; ++j) {
            P_temp[i][j] = 0;
            for (int k = 0; k < state_size; ++k) {
                P_temp[i][j] += F[i][k] * P[k][j];
            }
        }
    }

    for (int i = 0; i < state_size; ++i) {
        for (int j = 0; j < state_size; ++j) {
            P[i][j] = 0;
            for (int k = 0; k < state_size; ++k) {
                P[i][j] += P_temp[i][k] * F[j][k];
            }
        }
    }

    for (int i = 0; i < state_size; ++i) {
        for (int j = 0; j < state_size; ++j) {
            P[i][j] += Q[i][j];
        }
    }
}

// Fungsi update
void update(double x[], double P[][state_size], const double z[]) {
    // Prediksi pengukuran
    double z_pred[measurement_size];
    for (int i = 0; i < measurement_size; ++i) {
        z_pred[i] = 0;
        for (int j = 0; j < state_size; ++j) {
            z_pred[i] += H[i][j] * x[j];
        }
    }

    // Residual pengukuran
    double y[measurement_size];
    for (int i = 0; i < measurement_size; ++i) {
        y[i] = z[i] - z_pred[i];
    }

    // Kovariansi residual
    double S[measurement_size][measurement_size];
    for (int i = 0; i < measurement_size; ++i) {
        for (int j = 0; j < measurement_size; ++j) {
            S[i][j] = 0;
            for (int k = 0; k < state_size; ++k) {
                S[i][j] += H[i][k] * P[k][j];
            }
        }
    }

    for (int i = 0; i < measurement_size; ++i) {
        for (int j = 0; j < measurement_size; ++j) {
            S[i][j] += R[i][j];
        }
    }

    // Gain Kalman
    double K[state_size][measurement_size];
    for (int i = 0; i < state_size; ++i) {
        for (int j = 0; j < measurement_size; ++j) {
            K[i][j] = 0;
            for (int k = 0; k < measurement_size; ++k) {
                K[i][j] += P[i][k] * H[k][j];
            }
        }
    }

    double S_inv[measurement_size][measurement_size];
    double determinant = S[0][0] * S[1][1] - S[0][1] * S[1][0];
    S_inv[0][0] = S[1][1] / determinant;
    S_inv[0][1] = -S[0][1] / determinant;
    S_inv[1][0] = -S[1][0] / determinant;
    S_inv[1][1] = S[0][0] / determinant;

    double temp[state_size][measurement_size];
    for (int i = 0; i < state_size; ++i) {
        for (int j = 0; j < measurement_size; ++j) {
            temp[i][j] = 0;
            for (int k = 0; k < measurement_size; ++k) {
                temp[i][j] += P[i][k] * H[k][j];
            }
        }
    }

    double K_temp[state_size][measurement_size];
    for (int i = 0; i < state_size; ++i) {
        for (int j = 0; j < measurement_size; ++j) {
            K_temp[i][j] = 0;
            for (int k = 0; k < measurement_size; ++k) {
                K_temp[i][j] += temp[i][k] * S_inv[k][j];
            }
        }
    }

    // Update status
    for (int i = 0; i < state_size; ++i) {
        for (int j = 0; j < measurement_size; ++j) {
            x[i] += K_temp[i][j] * y[j];
        }
    }

    // Update kovariansi status
    double P_temp[state_size][state_size];
    for (int i = 0; i < state_size; ++i) {
        for (int j = 0; j < state_size; ++j) {
            P_temp[i][j] = 0;
            for (int k = 0; k < measurement_size; ++k) {
                P_temp[i][j] += K_temp[i][k] * H[k][j];
            }
        }
    }

    for (int i = 0; i < state_size; ++i) {
        for (int j = 0; j < state_size; ++j) {
            P[i][j] -= P_temp[i][j];
        }
    }
}

// Function to handle the received message
void on_message(struct mosquitto *mosq, void *userdata, const struct mosquitto_message *message) {
    if (message->payloadlen) {
        std::string payload = std::string(static_cast<char*>(message->payload), message->payloadlen);
        std::cout << "Received message: " << payload << std::endl;

        // Split the message based on comma
        size_t firstCommaPos = payload.find(',');
        if (firstCommaPos != std::string::npos) {
            size_t secondCommaPos = payload.find(',', firstCommaPos + 1);
            size_t thirdCommaPos = payload.find(',', secondCommaPos + 1);
            if (secondCommaPos != std::string::npos && thirdCommaPos != std::string::npos) {
                double range = std::stod(payload.substr(0, firstCommaPos));
                double bearing = std::stod(payload.substr(firstCommaPos + 1, secondCommaPos - firstCommaPos - 1));
                double speed = std::stod(payload.substr(secondCommaPos + 1, thirdCommaPos - secondCommaPos - 1));
                double acceleration = std::stod(payload.substr(thirdCommaPos + 1));

                // Calculate tangential and normal acceleration
                double at = acceleration * cos(bearing);
                double an = acceleration * sin(bearing);

                // Process the received data with EKF
                // Input tangential acceleration and normal acceleration
                predict(x, P, at, an);

                // Output estimated status
                std::cout << "Estimated status: ";
                for (int i = 0; i < state_size; ++i) {
                    std::cout << x[i] << " ";
                }
                std::cout << std::endl;
            } else {
                std::cerr << "Invalid message format" << std::endl;
            }
        } else {
            std::cerr << "Invalid message format" << std::endl;
        }
    } else {
        std::cerr << "Received empty message" << std::endl;
    }
}

int main() {
    mosquitto_lib_init();

    mosq = mosquitto_new(NULL, true, NULL);
    if (!mosq) {
        std::cerr << "Error: Out of memory." << std::endl;
        return 1;
    }

    // Set callback to handle messages
    mosquitto_message_callback_set(mosq, on_message);

    int rc = mosquitto_connect(mosq, MQTT_HOST, MQTT_PORT, 60);
    if (rc) {
        std::cerr << "Unable to connect to MQTT broker: " << mosquitto_strerror(rc) << std::endl;
        return 1;
    }

    rc = mosquitto_subscribe(mosq, NULL, MQTT_TOPIC, 0);
    if (rc) {
        std::cerr << "Error subscribing to topic: " << mosquitto_strerror(rc) << std::endl;
        return 1;
    }

    std::cout << "Subscribed to topic: " << MQTT_TOPIC << std::endl;

    rc = mosquitto_loop_forever(mosq, -1, 1);
    if (rc != MOSQ_ERR_SUCCESS) {
        std::cerr << "Error in the loop: " << mosquitto_strerror(rc) << std::endl;
    }

    mosquitto_destroy(mosq);
    mosquitto_lib_cleanup();

    return 0;
}

