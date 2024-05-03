#include <cmath>
#include <cstdio>

struct XYZ {
    double x;
    double y;
    double z;
};

// Use these values for D65 whitepoint.
// published values vary slightly, and will affect the results.
XYZ d65WhitePoint = { 95.0489, 100.0, 108.8840};


// This function should implement the CIE dE 2000 colour difference
double deltaE_CIE2000(const XYZ& a, const XYZ& b);


struct TestVector {
    int id;
    XYZ a;
    XYZ b;
    double dE;
};

TestVector tests[] = {
    {  1, {       18.0,  18.421875,   98.71875 }, {       17.5,  18.421875, 103.328125 }, 2.0428693215228901 },
    {  2, {   18.09375,  18.421875,  94.953125 }, {       17.5,  18.421875, 103.328125 }, 2.8704685611441016 },
    {  3, {   18.03125,  18.421875,  90.171875 }, {       17.5,  18.421875, 103.328125 }, 3.44367359586365 },
    {  4, {      17.25,  18.421875, 105.765625 }, {       17.5,  18.421875, 103.328125 }, 0.98953039312009172 },
    {  5, {   17.28125,  18.421875, 106.609375 }, {       17.5,  18.421875, 103.328125 }, 1.0029552436525526 },
    {  6, {   17.34375,  18.421875, 107.765625 }, {       17.5,  18.421875, 103.328125 }, 0.97417324220258006 },
    {  7, {       17.5,  18.421875,    20.0625 }, {  17.328125,  18.421875,  19.015625 }, 2.3224208445484056 },
    {  8, {  17.328125,  18.421875,  19.015625 }, {       17.5,  18.421875,    20.0625 }, 2.3224208445484056 },
    {  9, {   17.96875,  18.421875,    20.0625 }, {  17.046875,  18.421875,  20.046875 }, 7.2394063666101118 },
    { 10, {   17.96875,  18.421875,    20.0625 }, {  17.046875,  18.421875,  20.046875 }, 7.2394063666101118 },
    { 11, {   17.96875,  18.421875,    20.0625 }, {  17.046875,  18.421875,  20.046875 }, 7.2394063666101118 },
    { 12, {   17.96875,  18.421875,    20.0625 }, {  17.046875,  18.421875,  20.046875 }, 7.2394063666101118 },
    { 13, {       17.5,  18.421875,  18.765625 }, {       17.5,  18.421875,   21.40625 }, 4.8166832023426167 },
    { 14, {       17.5,  18.421875,  18.765625 }, {       17.5,  18.421875,   21.40625 }, 4.8166832023426167 },
    { 15, {       17.5,  18.421875,  18.765625 }, {       17.5,  18.421875,   21.40625 }, 4.8166832023426167 },
    { 16, {   17.96875,  18.421875,    20.0625 }, {       17.5,  18.421875,   21.40625 }, 4.3294648525052608 },
    { 17, {   17.96875,  18.421875,    20.0625 }, {     51.875,  45.171875,   68.59375 }, 27.160687903094356 },
    { 18, {   17.96875,  18.421875,    20.0625 }, {    26.5625,      29.25,  15.203125 }, 22.874523554449674 },
    { 19, {   17.96875,  18.421875,    20.0625 }, {  17.296875,   23.90625,   27.96875 }, 31.784604209825009 },
    { 20, {   17.96875,  18.421875,    20.0625 }, {  30.671875,   25.96875,  19.421875 }, 19.476261849382428 },
    { 21, {   17.96875,  18.421875,    20.0625 }, {   18.09375,  18.421875,      19.75 }, 0.99603164263017707 },
    { 22, {   17.96875,  18.421875,    20.0625 }, {     18.125,  18.421875,    20.0625 }, 1.0414923919236823 },
    { 23, {   17.96875,  18.421875,    20.0625 }, {  17.859375,  18.421875,      19.75 }, 0.95305102475035275 },
    { 24, {   17.96875,  18.421875,    20.0625 }, {  18.109375,  18.421875,     19.875 }, 0.99326866414194748 },
    { 25, {  19.453125,   28.40625,      11.75 }, {  19.609375,  28.640625,  10.734375 }, 1.2804838792912299 },
    { 26, {   22.53125,   31.59375,  39.046875 }, {  22.640625,     31.375,  37.328125 }, 1.2730197313816478 },
    { 27, {    29.0625,  29.578125,   36.28125 }, {   28.84375,  29.734375,     36.125 }, 1.8229382999799022 },
    { 28, {    4.15625,   8.546875,   8.140625 }, {   4.421875,   8.515625,   8.765625 }, 1.841290832234854 },
    { 29, {    4.96875,    3.71875,     19.875 }, {   4.671875,     3.8125,  18.046875 }, 2.075123140549433 },
    { 30, {  15.640625,       9.25,    5.09375 }, {  15.953125,    9.15625,   4.453125 }, 1.4155796660121873 },
    { 31, {    73.1875,  78.046875,       83.0 }, {     74.125,    78.8125,  85.765625 }, 1.458222602968178 },
    { 32, {    74.1875,    78.3125,    86.5625 }, {   69.34375,   73.40625,  80.890625 }, 1.5635548124666301 },
    { 33, {   0.703125,       0.75,   0.984375 }, {   0.609375,    0.65625,   0.859375 }, 0.6553561450639207 },
    { 34, {    0.21875,   0.234375,   0.328125 }, {    0.09375,    0.09375,   0.140625 }, 1.024061208946461}
};


int main() {
    const int numTests = sizeof(tests) / sizeof(tests[0]);
    const double tolerance = 0.01;
    
    int numFailures = 0;
    for (int i = 0; i < numTests; i++) {
        const TestVector& test = tests[i];
        double diff = deltaE_CIE2000(test.a, test.b);
        if (std::abs(diff - test.dE) > tolerance) {
            std::fprintf(stderr, "Failed: %d, expected %f - got %f (delta: %f)\n", test.id, test.dE, diff, std::abs(diff - test.dE));
            numFailures++;
        }
    }

    if (numFailures == 0) {
        std::fprintf(stderr, "Passed\n");
    }
}

/*
 Put your implementation of CIE dE2000 here...
*/

/*
The formula for CIE dE2000 was taken from:
https://en.wikipedia.org/wiki/Color_difference

The conversion from XYZ to LAB was taken from:
https://en.wikipedia.org/wiki/CIELAB_color_space

Variables are named according to the wikipedia article
*/

// Struct to store the L*a*b* colour space
// This colour space should be used in the CIE dE 2000 calculation
struct Lab {
    double L;
    double a;
    double b;
};

// Convert XYZ colour space to L*a*b*
Lab XYZtoLAB(const XYZ& xyz) {
    // Normalise the XYZ values based on the reference white point
    double X = xyz.x / d65WhitePoint.x;
    double Y = xyz.y / d65WhitePoint.y;
    double Z = xyz.z / d65WhitePoint.z;

    // f lambda function to apply to variables
    auto f = [](double t) {
        const double delta = 6.0 / 29.0;
        if (t > pow(delta, 3))
            return cbrt(t);
        else
            return t / (3 * delta * delta) + 4.0 / 29.0;
    };


    // Calculate the L*a*b* values
    Lab lab;
    lab.L = 116 * f(Y) - 16;
    lab.a = 500 * (f(X) - f(Y));
    lab.b = 200 * (f(Y) - f(Z));

    return lab;
}

// utility function to calculate the atan2 in degrees 0-360
// edge case for atan2(0, 0) results in 0
double degree_atan2(double b, double a) {
    double radians = atan2(b, a);
    if (b == 0 || a == 0) {
        return 0;
    }
    double degrees = radians * 180 / M_PI;
    if (degrees < 0) {
        degrees = 360 + degrees;
    }
    return degrees;
}

// utility functions in degrees
double degree_cos(double angle) {
    return cos(angle * M_PI / 180);
}

double degree_sin(double angle) {
    return sin(angle * M_PI / 180);
}

// Function to calculate the delta h value, inputs are in degrees
double calc_delta_h_prime(double h1p, double h2p) {
    double diff = h2p - h1p;
    double absDiff = abs(diff);
    if (absDiff > 180) {
        if (diff > 0) {
            return diff - 360;
        } else {
            return diff + 360;
        }
    }
    return diff;
}

// function to calculate H bar prime
double calc_H_bar_prime(double h1p, double h2p, double C_1_prime, double C_2_prime) {
    double absDiff = abs(h1p - h2p);
    double result = 0;
    if (absDiff > 180){
        if (h1p + h2p < 360)
            result = (h1p + h2p + 360) / 2;
        else
            result = (h1p + h2p - 360) / 2;
    } else
        result = (h1p + h2p) / 2;
    
    // if one angle is indeterminate, then the other angle is the average of the two
    if (C_1_prime == 0 || C_2_prime == 0)
        result = result * 2;
    
    return result;
}


// Function to calculate the CIE dE 2000 colour difference
// Note: k_L, k_C, and k_H are set to 1
double deltaE_CIE2000(const XYZ& a, const XYZ& b) {
    // Convert the XYZ values to L*a*b*
    Lab labA = XYZtoLAB(a);
    Lab labB = XYZtoLAB(b);

    // Define all double variables
    double L1_star, a1_star, b1_star; // L*a*b* values for colour a
    double L2_star, a2_star, b2_star; // L*a*b* values for colour b
    double L_bar, C_bar; // Average values
    double a_1_prime, a_2_prime, C_1_prime, C_2_prime, h_1_prime, h_2_prime; // Prime values
    double C_bar_prime, H_bar_prime, delta_h_prime, delta_C_prime, delta_L_prime, delta_H_prime; // Prime values

    // Assign the L*a*b* values for colour a
    L1_star = labA.L;
    a1_star = labA.a;
    b1_star = labA.b;

    // Assign the L*a*b* values for colour b
    L2_star = labB.L;
    a2_star = labB.a;
    b2_star = labB.b;


    // Calculate the average values
    L_bar = (L1_star + L2_star) / 2;
    C_bar = (sqrt(pow(a1_star, 2) + pow(b1_star, 2)) + sqrt(pow(a2_star, 2) + pow(b2_star, 2))) / 2;

    // calculate the prime values
    delta_L_prime = L2_star - L1_star;

    a_1_prime = a1_star + (a1_star / 2) * (1 - sqrt(pow(C_bar, 7) / (pow(C_bar, 7) + pow(25, 7))));
    a_2_prime = a2_star + (a2_star / 2) * (1 - sqrt(pow(C_bar, 7) / (pow(C_bar, 7) + pow(25, 7))));

    C_1_prime = sqrt(pow(a_1_prime, 2) + pow(b1_star, 2));
    C_2_prime = sqrt(pow(a_2_prime, 2) + pow(b2_star, 2));

    C_bar_prime = (C_1_prime + C_2_prime) / 2;
    delta_C_prime = C_2_prime - C_1_prime;

    h_1_prime = degree_atan2(b1_star, a_1_prime);
    h_2_prime = degree_atan2(b2_star, a_2_prime);

    delta_h_prime = calc_delta_h_prime(h_1_prime, h_2_prime);
    delta_H_prime = 2 * sqrt(C_1_prime * C_2_prime) * degree_sin(delta_h_prime / 2);
    H_bar_prime = calc_H_bar_prime(h_1_prime, h_2_prime, C_1_prime, C_2_prime);
    
    // calculate final values
    double T = 1 - 0.17 * degree_cos(H_bar_prime - 30) + 0.24 * degree_cos(2 * H_bar_prime) + 0.32 * degree_cos(3 * H_bar_prime + 6) - 0.20 * degree_cos(4 * H_bar_prime - 63);

    double S_L = 1 + (0.015 * pow(L_bar - 50, 2)) / sqrt(20 + pow(L_bar - 50, 2));
    double S_C = 1 + 0.045 * C_bar_prime;
    double S_H = 1 + 0.015 * C_bar_prime * T;

    double R_T = -2 * sqrt(pow(C_bar_prime, 7) / (pow(C_bar_prime, 7) + pow(25, 7))) * degree_sin(60 * exp(-pow((H_bar_prime - 275) / 25, 2)));

    double result = sqrt(pow(delta_L_prime / S_L, 2) + pow(delta_C_prime / S_C, 2) + pow(delta_H_prime / S_H, 2) + R_T * delta_C_prime * delta_H_prime / (S_C * S_H));
    return result;
}
