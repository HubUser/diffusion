#include <iostream>
#include <string>

extern "C" {
    void get_diffusion(float energy, float ratio, float kp, float* coefA, float* coefB);
}
    void radbelInputTet2DstDoubleQ(float, float, std::string, std::string, float, float, int, float*, float*);

void get_diffusion(float energy, float ratio, float kp, float* coefA, float* coefB) {
    const float PI = 3.14159;

    float gam = 1 + energy / 511.;
    std::string gamS = std::to_string(ratio);
    std::string gamG = std::to_string(energy);
    radbelInputTet2DstDoubleQ(gam, ratio, gamS, gamG, PI * 45. / 180., 2.5, kp, coefA, coefB); 

    std::cout << "Finish" << std::endl;
}
