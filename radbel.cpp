#include <vector>
#include <future>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "ctpl.h"
#include "coefs.h"
#include "FC.h"

extern "C" {
    void localdtet2dst(float*, float*, float*, int*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, int*);
    float alpha(float* a, float* b);
}

Coefs parallel_exec(int id, const Args a);
void printResults(Args a, Coefs c);
void radbelInputTet2DstDoubleQ(float gam, float s, std::string nameS, std::string nameG, float lM, float alphaL, int kp, float* coefA, float* coefB) {

    std::string nA = std::to_string(alphaL);
    std::string nameB = std::to_string(kp);

    std::string fileHeader = "alp n0 n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 n13 n14 n15 n16 n17 n18 n19 n20 n21 n22 n23 n24 n25 n26 n27 n28 n29 n30 n 31 n32 n33 n34 n35 n36 n37 n38 n39 n40 n41 n42 n43 n44 n45 n46 n47 n48 n49 n50 n51 n52 n53 n54 n55 n56 n57 n58 n59 n60 n-1 n-2 n-3 n -4 n-5 n-6 n-7 n-8 n-9 n-10 n-11 n-12 n-13 n-14 n-15 n-16 n-17 n-18 n-19 n-20 n-21 n-22 n-23 n-24 n-25 n-26 n-27 n-28 n-29 n-30 n-31 n-32 n-33 n-34 n-35 n-36 n-37 n-38 n-39 n-40 n-41 n-42 n-43 n-44 n-45 n-46 n-47 n-48 n-49 n-50 n-51 n-52 n-53 n-54 n-55 n-56 n-57 n-58 n-59 n-60 sum sum5 sum10 sum20 sum40";
    std::vector<std::string> labels = { "_ee_mag(", "_aa_mag(", "_aa_mgs(", "_ee_mgs(" };
    for (auto label: labels) {
        std::ostringstream buffer;
        buffer << nameG << label << nameS << ")_" << nA << "_Kp" << nameB << ".txt";
        std::string fileName = buffer.str();
        std::ofstream fileUnit(fileName, std::ofstream::out|std::ofstream::trunc);
        fileUnit << fileHeader << std::endl;
        fileUnit.close();
    }

    int numN = 61;
    const int nalp = 100;
    ctpl::thread_pool p(2);
    std::vector<std::future<Coefs>> results(nalp);
    for (int ia = 0; ia < nalp; ia++) {
        Args a(gam, s, nameS, nameG, nameB, nA, lM, alphaL, kp, coefA, coefB, ia, numN, nalp);
        results[ia] = p.push(parallel_exec, a);
    }
    for (int ia = 0; ia < nalp; ia++) {
        Args a(gam, s, nameS, nameG, nameB, nA, lM, alphaL, kp, coefA, coefB, ia, numN, nalp);
        printResults(a, results[ia].get());
    }
}

Coefs parallel_exec(int id, const Args a) {
    const float alphaMax = 90.0;
    float coefDaaLoc[2 * a.numN + 5], coefDeeB[2 * a.numN + 5];
    float coefDaaSB[2 * a.numN + 5], coefDeeSB[2 * a.numN + 5];
    float* coefDaaB = new float[2 * a.numN + 5];
    /* float coefDaaB[2 * a.numN + 5], coefDaaSB[2 * a.numN + 5], coefDeeSB[2 * a.numN + 5]; */
    float coefDeaLoc[2 * a.numN + 5], coefDeeLoc[2 * a.numN + 5];
    const float dl = 0.01;
    const float lM0 = 3.14159;
    float gam = a.gam;
    int numN = a.numN;

    const float valueTm01=11.5372, a1Tm1=14.3157, a2Tm1=-8.18009, a3Tm1=1.38532, a4Tm1=0.0, a5Tm1=0.0;
    const float valueTm02=66.0, a1Tm2=1.0, a2Tm2=0.0, a3Tm2=0.0, a4Tm2=0.0, a5Tm2=0.0;

    const float valueTw01=5.68019, a1Tw1=4.62479, a2Tw1=3.04756, a3Tw1=-5.05697, a4Tw1=1.83305, a5Tw1=-2.02501e-001;
    const float valueTw02=5.68019, a1Tw2=4.62479, a2Tw2=3.04756, a3Tw2=-5.05697, a4Tw2=1.83305, a5Tw2=-2.02501e-001;

    const float L = 5.;
    const float norm = 184. / pow(L / 4.5, 3.);
    float alpEQ = (5.0 + (alphaMax - 5.0) * a.ia / (a.nalp - 1.0) ) * 3.14159/180.;
    float Tl = 1.3 - 0.56 * sin(alpEQ);

    coefDeeB[0] = alpEQ * 180.0 / 3.14159;
    coefDaaB[0] = alpEQ * 180.0 / 3.14159;
    coefDaaSB[0] = alpEQ * 180.0 / 3.14159;
    coefDeeSB[0] = alpEQ * 180.0 / 3.14159;
    for (int i = 1; i < 2*a.numN+5; i++) {
        coefDeeSB[i] = 0;
        coefDeeB[i] = 0;
        coefDaaB[i] = 0;
        coefDaaSB[i] = 0;
    }
    float lam = 0.0;
    float alp = alpha(&lam, &alpEQ);
    while (alp > 0) {
        float b = sqrt(1 + 3 * pow(sin(lam), 2))/pow(cos(lam), 6);
        float ne0 = pow(cos(lam), (-2 * a.alphaL));
        float ldeg = lam * 180 / 3.14159 / 10.0;

        float valueTm1=(valueTm01+a1Tm1*ldeg+a2Tm1*pow(ldeg, 2)+a3Tm1*pow(ldeg, 3)+a4Tm1 * pow(ldeg, 4)+a5Tm1*pow(ldeg, 5))*3.14159/180.0;
        float valueTm2=(valueTm02+a1Tm2*ldeg+a2Tm2*pow(ldeg, 2)+a3Tm2*pow(ldeg, 3)+a4Tm2 *pow(ldeg, 4)+a5Tm2*pow(ldeg, 5))*3.14159/180.0;
        float valueTw1=sqrt(2.0)*(valueTw01+a1Tw1*ldeg+a2Tw1*pow(ldeg, 2)+a3Tw1*pow(ldeg, 3)+a4Tw1*pow(ldeg, 4)+a5Tw1*pow(ldeg, 5))*3.14159/180.0;
        float valueTw2=sqrt(2.0)*(valueTw02+a1Tw2*ldeg+a2Tw2*pow(ldeg, 2)+a3Tw2*pow(ldeg, 3)+a4Tw2*pow(ldeg, 4)+a5Tw2*pow(ldeg, 5))*3.14159/180.0;

        float Bw = 0;
        for (int j1 = 0; j1 < 4; j1++) {
            for (int i1 = 0; i1 < 4; i1++) {
                Bw += a.coefB[j1 * 4 + i1]*pow(a.kp, i1)*pow(ldeg*10, j1);
            }
        }
        Bw *= 1000.0;

        if (a.kp > 50 && ldeg > 3.0 && Bw < 3.0) Bw = 3.0;
        if (Bw < 1.0) Bw = 1.0;

        float factorQ = 0.0;
        if (a.kp < 65) {
            for (int j1 = 0; j1 < 4; j1++) {
                for (int i1 = 0; i1 < 4; i1++) {
                    factorQ += a.coefA[j1 * 4 + i1] * pow(a.kp, i1) * pow(ldeg * 10, j1);
                }
            }
            factorQ = pow(10.0, factorQ);
            if (factorQ < 0.005) factorQ = 0.0;
        }

        float valueWm = 0.35; //  wave frequency

        float valueTm[2];
        float valueTw[2];
        if (lam < a.lM) {
            valueTm[0] = valueTm1;
            valueTm[1] = valueTm2;
            valueTw[0] = valueTw1;
            valueTw[1] = valueTw2;

            float temp = a.s * sqrt(ne0) / b;
            float temp1 = 1.0;
            float temp2 = 1.571;
            int temp3 = 1;
            localdtet2dst(&gam, &temp, &alp, &numN, coefDaaLoc, coefDeaLoc, coefDeeLoc, valueTm, valueTw, &temp1, &factorQ, &temp2, &b, &valueWm, &temp3);
            /* for (int ii = 0; ii < 2 * numN + 5; ii++) { */
            /*     std::cout << coefDaaLoc[ii] << "  "; */
            /* } */
            /* std::cout << std::endl; */
            /* std::cout.flush(); */
            /* int *x = 0; */
            /* *x = 0; */

            for (int i = 1; i < 2 * a.numN + 5; i++) {
                coefDaaB[i] += coefDaaLoc[i]*cos(alp)/pow(cos(alpEQ), 2)*pow(cos(lam), 7)*dl/Tl/norm*Bw*Bw/100.0/100.0;
                coefDeeB[i] += coefDeeLoc[i]/cos(alp)*sqrt(1+3*pow(sin(lam), 3))*cos(lam)*dl/Tl/norm*Bw*Bw/100.0/100.0;
            }
        }

        if (lam < lM0) {
            valueTm[0] = 0.0;
            valueTm[1] = 0.0;
            valueTw[0] = 0.523;
            valueTw[1] = 1.0;

            float temp = a.s * sqrt(ne0) / b;
            float temp1 = 1.0;
            float temp2 = 1.571;
            int temp3 = 1;
            float temp4 = 0.0;
            localdtet2dst(&gam, &temp, &alp, &numN, coefDaaLoc, coefDeaLoc, coefDeeLoc, valueTm, valueTw, &temp1, &temp4, &temp2, &b, &valueWm, &temp3);

            for (int i = 1; i < 2 * a.numN + 5; i++) {
                coefDaaSB[i] += coefDaaLoc[i]*cos(alp)/pow(cos(alpEQ), 2)*pow(cos(lam), 7)*dl/Tl/norm*Bw*Bw/100.0/100.0;
                coefDeeSB[i] += coefDeeLoc[i]/cos(alp)*sqrt(1+3 * pow(sin(lam), 3))*cos(lam)*dl/Tl/norm*Bw*Bw/100.0/100.0;
            }
        }

        lam = lam + dl;
        alp = alpha(&lam, &alpEQ);
    }
    /* if (id ==0) { */
    /*     std::cout << a.print() << std::endl; */
    /*     for (int ii = 0; ii < 2 * numN + 5; ii++) { */
    /*         std::cout << coefDaaB[ii] << "  "; */
    /*     } */
    /*     std::cout << std::endl; */
    /*     std::cout.flush(); */
    /*     int *x = 0; */
    /*     *x = 0; */
    /* } */
    Coefs res = Coefs(coefDeeB, coefDaaB, coefDaaSB, coefDeeSB, a.numN);
    delete[] coefDaaB;
    return res;
}

void printResults(Args a, Coefs c) {
    std::vector<float*> coefs = { c.coefDeeB, c.coefDaaB, c.coefDaaSB, c.coefDeeSB };
    std::vector<std::string> labels = { "_ee_mag(", "_aa_mag(", "_aa_mgs(", "_ee_mgs(" };
    for (int l = 0; l < labels.size(); l++) {
        std::ostringstream buffer;
        buffer << a.nameG << labels[l] << a.nameS << ")_" << a.nA << "_Kp" << a.nameB << ".txt";
        std::string fileName = buffer.str();
        std::ofstream fileUnit(fileName, std::ofstream::out|std::ofstream::app);
        fileUnit << std::scientific;
        for (int c = 0; c < 2 * a.numN + 5; c++) {
            fileUnit << coefs[l][c] << "  ";
        }
        fileUnit << std::endl;
        fileUnit.close();
    }
    std::cout << "Progress: " << a.ia << "/" << a.nalp << std::endl;
}
