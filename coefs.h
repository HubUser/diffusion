struct Args {
    Args(float gam, float s, std::string nameS, std::string nameG, std::string nameB, std::string nA, float lM, float alphaL, 
            int kp, float* coefA, float* coefB, int ia, int numN, int nalp) {
        this->gam = gam;
        this->s = s;
        this->nameS = nameS;
        this->nameG = nameG;
        this->nameB = nameB;
        this->nA = nA;
        this->lM = lM;
        this->alphaL = alphaL;
        this->kp = kp;
        this->coefA = coefA;
        this->coefB = coefB;
        this->ia = ia;
        this->numN = numN;
        this->nalp = nalp;
    }
    float gam;
    float s;
    std::string nameS;
    std::string nameG;
    std::string nameB;
    std::string nA;
    float lM;
    float alphaL;
    int kp;
    float* coefA;
    float* coefB;
    int ia;
    int numN;
    int nalp;
    std::string print() const {
        std::ostringstream buffer;
        buffer << gam << "  " << s << "  " << lM << "  " << ia << "  " << numN << "  " << nalp << std::endl;
        return buffer.str();
    }
};

struct Coefs {
    Coefs(float* coefDeeB, float* coefDaaB, float* coefDaaSB, float* coefDeeSB, const int& numN_) 
    : numN(numN_) {
        this->coefDeeB = new float[2 * numN + 5];
        this->coefDaaB = new float[2 * numN + 5];
        this->coefDaaSB = new float[2 * numN + 5];
        this->coefDeeSB = new float[2 * numN + 5];
        for (int ii = 0; ii < 2 * numN + 5; ii++) {
            this->coefDeeB[ii] = coefDeeB[ii];
            this->coefDaaB[ii] = coefDaaB[ii];
            this->coefDaaSB[ii] = coefDaaSB[ii];
            this->coefDeeSB[ii] = coefDeeSB[ii];
        }
    }
    Coefs(const Coefs& c) : numN(c.numN) {
        this->coefDeeB = new float[2 * numN + 5];
        this->coefDaaB = new float[2 * numN + 5];
        this->coefDaaSB = new float[2 * numN + 5];
        this->coefDeeSB = new float[2 * numN + 5];
        for (int ii = 0; ii < 2 * numN + 5; ii++) {
            this->coefDeeB[ii] = c.coefDeeB[ii];
            this->coefDaaB[ii] = c.coefDaaB[ii];
            this->coefDaaSB[ii] = c.coefDaaSB[ii];
            this->coefDeeSB[ii] = c.coefDeeSB[ii];
        }
    }
    ~Coefs() {
        delete[] coefDeeB;
        delete[] coefDaaB;
        delete[] coefDaaSB;
        delete[] coefDeeSB;
    }
    float* coefDeeB;
    float* coefDaaB;
    float* coefDaaSB;
    float* coefDeeSB;
    const int numN;
};
