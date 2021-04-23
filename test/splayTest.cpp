#include "config.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <sys/time.h>
#include <chrono>
using namespace std;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::system_clock;

default_random_engine e;
uniform_int_distribution<int> random4(0, 3);
uniform_real_distribution<double> random_double(0, 1.0);

void initRandom() {
    unsigned seed = (unsigned)time(NULL);
    e = std::default_random_engine(seed);
    srand(seed);
}

double randomDouble() {
    return random_double(e);
}
int randomInt(int n) {
    return rand() % n;
}

const double beta = 10.0;
const double V = 1.0;
const double K = 1.0;

double U = 1.0;
double mu = 0.5 * U;

double D = 1.0;
double myGamma = -1.0;

int mJ = 12, mL = 12, mF = 3;

vector<double> randomVector(int length) {
    vector<double> res(length);
    for (int i = 0; i < length; ++i) {
        res[i] = randomDouble();
    }
    return res;
}

void testConfiguration(bool efficiencyTest = false, int minimumN = 500, int dataCount = 10000) {
    vector<double> aJ = randomVector(mJ + 1);
    // vector<double> aJ(mJ + 1, 0.0);
    vector<double> aL = randomVector(mL + 1);
    // vector<double> aL = {1.0, 1.0, 1.0};
    vector<double> aF = randomVector(mF + 1);
    // SplayWeightConfiguration c(aJ, aL, aF);
    SplayWeightConfiguration c(aJ, aL, aF);
    double dWeightJL, dWeightF;
    double totalError = 0.0;
    // int dataCount = 10000;
    double averageN = 0;

    auto millisecBeforeSimulation = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    for (int i = 0; i < dataCount; ++i) {
        bool isInsert = randomInt(2);
        if (isInsert) {
            bool spin = randomInt(2);
            // bool spin = true;
            double tau = randomDouble() * beta;
            // cout << "insert spin = " << spin << ", tau = " << tau << endl;
            dWeightJL = c.insertDeltaWeightJL(spin, tau);
            dWeightF = c.insertDeltaWeightF();
            // cout << "after calculating weights" << endl;
            // c.weight += c.insertDeltaWeight(spin, tau);
            // cout << "finished weight calculation" << endl;
            c.insertSpin(spin, tau);
            // cout << "dWeightJL = " << dWeightJL << endl;
            // cout << "dWeightF = " << dWeightF << endl;
            c.updateWeight(dWeightJL, dWeightF);
        } else if (c.n > minimumN) {
            int idx = randomInt(c.n);
            // cout << "remove spin = " << c.c[idx].first << ", tau = " << c.c[idx].second << endl;
            dWeightJL = c.removeDeltaWeightJL(c.c[idx].first, c.c[idx].second);
            dWeightF = c.removeDeltaWeightF();
            // c.weight += c.removeDeltaWeight(c.c[idx].first, c.c[idx].second);
            c.updateWeight(dWeightJL, dWeightF);
            c.removeSpin(idx);
        }
        // cout << c.getWeight() << ' ' <<  c.bruteForceMBH() << endl;
        // printf("%.12lf\n", c.getWeight() - c.bruteForceMBH());
        if (!efficiencyTest) {
            double weight = c.getWeight(), mbhWeight = c.bruteForceMBH();
            totalError += (weight - mbhWeight) / mbhWeight;
        }
        // cout << "weight = " << weight << endl;
        // cout << "mbh weight = " << mbhWeight << endl;
        // cout << "error = " << ((weight - mbhWeight) / mbhWeight) << endl;
        // c.printPF();

        averageN += c.n;
    }
    auto millisecAfterSimulation = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    averageN /= dataCount;
    // cout << "average n = " << averageN << endl;
    if (efficiencyTest) {
        cout << "Test the efficiency for a simulation of " << dataCount << " steps." << endl;
        cout << (millisecAfterSimulation - millisecBeforeSimulation) * 1000000.0 / dataCount << " milliseconds for 10^6 steps, with average n = " << averageN << ", (mL, mJ, mF) = (" << mL << ", " << mJ << ", " << mF << ")." << endl;
    }
    if (!efficiencyTest) {
        double averageError = totalError / dataCount;
        cout << "average error for m^2 configuration = " << averageError << endl;
    }
    // c.printA();
}

int main() {
    initRandom();
    testConfiguration(false, 30, 10000);
    testConfiguration(true, 500, 100000);
    return 0;
}