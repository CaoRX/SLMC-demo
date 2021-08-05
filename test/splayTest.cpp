#include "config.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <sys/time.h>
#include <chrono>
#include <fstream>
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
const double range = 2.0;
double generateDataPoint() {
    return (randomDouble() - 0.5) * range;
}
// const double V = 1.0;
// const double K = 1.0;

// double U = 1.0;
// double mu = 0.5 * U;

// double D = 1.0;
// double myGamma = -1.0;

int mJ = 12, mL = 12, mF = 3;

vector<double> randomVector(int length) {
    vector<double> res(length);
    for (int i = 0; i < length; ++i) {
        res[i] = randomDouble();
    }
    return res;
}

void testConfigurationCorrectness(int dataCount = 10000, int minimumN = 30) {
    // run a test of dataCount steps, n cannot be decreased if lower than minimumN
    // test the correctness by comparing with O(n^2) weight calculation
    vector<double> aJ = randomVector(mJ + 1);
    vector<double> aL = randomVector(mL + 1);
    vector<double> aF = randomVector(mF + 1);
    bool efficiencyTest = false;
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
            double tau = generateDataPoint();
            dWeightJL = c.insertDeltaWeightJL(spin, tau);
            dWeightF = c.insertDeltaWeightF();
            c.insertSpin(spin, tau);
            c.updateWeight(dWeightJL, dWeightF);
        } else if (c.n > minimumN) {
            int idx = randomInt(c.n);
            dWeightJL = c.removeDeltaWeightJL(c.c[idx].first, c.c[idx].second);
            dWeightF = c.removeDeltaWeightF();
            c.updateWeight(dWeightJL, dWeightF);
            c.removeSpin(idx);
        }
        if (!efficiencyTest) {
            double weight = c.getWeight(), mbhWeight = c.bruteForceMBH();
            totalError += fabs((weight - mbhWeight) / mbhWeight);
        }

        averageN += c.n;
    }
    auto millisecAfterSimulation = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    averageN /= dataCount;
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

void testEfficiency(int dataCount = 100000, int minimumN = 500, int maximumN = 1000, int label = -1, int loop = 1) {
    // run a test of dataCount steps, n cannot be decreased if lower than minimumN
    // test the efficiency(and also correctness) by comparing with O(n) update
    if (maximumN < minimumN) {
        maximumN = minimumN * 2;
    }
    if (label == -1) {
        label = minimumN;
    }
    string fileName = "n" + std::to_string(label) + ".txt";
    cout << "output has been written to " << fileName << endl;
    fstream fout(fileName, std::ios::out);
    fout << "efficiency test: (mJ, mL, mF) = (" << mJ << ", " << mL << ", " << mF << "), step = " << dataCount << ", minimumN = " << minimumN << std::endl;
    // printf("efficiency test: (mJ, mL, mF) = (%d, %d, %d), step = %d, minimumN = %d\n", mJ, mL, mF, dataCount, minimumN);

    double effTime = 0, simpleTime = 0;
    double totalAverageN = 0;
    double totalErrorSum = 0.0;
    double maximumError = 0.0, minimumV = -1.0;

    for (int lp = 0; lp < loop; ++lp) {
        vector<double> aJ = randomVector(mJ + 1);
        vector<double> aL = randomVector(mL + 1);
        vector<double> aF = randomVector(mF + 1);

        SplayWeightConfiguration c(aJ, aL, aF);
        SimpleConfiguration cS(aJ, aL, aF);
        double dWeightJL, dWeightF;
        double dWeightJLS, dWeightFS;
        double totalError = 0.0;
        // int dataCount = 10000;
        double averageN = 0;

        vector<pair<bool, double> > vertices;
        vector<int> removeIdx;
        vector<bool> opts; // false for remove, true for insert

        vector<double> weightEff, weightSimple;
        int currN = 0;
        for (int i = 0; i < dataCount; ++i) {
            bool opt;
            if (currN <= minimumN) {
                opt = true;
            } else if (currN >= maximumN) {
                opt = false;
            } else {
                opt = randomInt(2);
            }
            if (opt) {
                opts.push_back(opt);
                vertices.push_back({randomInt(2), generateDataPoint()});
                removeIdx.push_back(-1);
                ++currN;
            } else {
                opts.push_back(opt);
                vertices.push_back({0, 0});
                removeIdx.push_back(randomInt(currN));
                --currN;
            }
        }

        auto millisecBeforeSimulation = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < dataCount; ++i) {
            bool isInsert = opts[i];
            if (isInsert) {
                bool spin = vertices[i].first;
                double tau = vertices[i].second;

                dWeightJL = c.insertDeltaWeightJL(spin, tau);
                dWeightF = c.insertDeltaWeightF();
                c.insertSpin(spin, tau);
                c.updateWeight(dWeightJL, dWeightF);
            } else {
                int idx = removeIdx[i];
                dWeightJL = c.removeDeltaWeightJL(c.c[idx].first, c.c[idx].second);
                dWeightF = c.removeDeltaWeightF();
                c.updateWeight(dWeightJL, dWeightF);
                c.removeSpin(idx);
            }
            double weight = c.getWeight();
            weightEff.push_back(weight);
            averageN += c.n;
        }
        auto millisecAfterSimulation = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        averageN /= dataCount;
        effTime += (millisecAfterSimulation - millisecBeforeSimulation);

        averageN = 0;

        millisecBeforeSimulation = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < dataCount; ++i) {
            bool isInsert = opts[i];
            if (isInsert) {
                bool spin = vertices[i].first;
                double tau = vertices[i].second;

                dWeightJLS = cS.insertDeltaWeightJL(spin, tau);
                dWeightFS = cS.insertDeltaWeightF();
                cS.insertSpin(spin, tau);
                cS.updateWeight(dWeightJLS, dWeightFS);
            } else {
                int idx = removeIdx[i];
                dWeightJLS = cS.removeDeltaWeightJL(cS.c[idx].first, cS.c[idx].second);
                dWeightFS = cS.removeDeltaWeightF();
                cS.updateWeight(dWeightJLS, dWeightFS);
                cS.removeSpin(idx);
            }
            double weight = cS.getWeight();
            weightSimple.push_back(weight);
            averageN += cS.n;
        }
        averageN /= dataCount;
        totalAverageN += averageN;
        millisecAfterSimulation = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

        simpleTime += millisecAfterSimulation - millisecBeforeSimulation;

        totalError = 0.0;

        for (int i = 0; i < dataCount; ++i) {
            double error = fabs(weightEff[i] - weightSimple[i]) / fabs(weightSimple[i]);
            totalError += error;
            maximumError = max(maximumError, error);
            minimumV = (minimumV < 0) ? fabs(weightSimple[i]) : min(fabs(weightSimple[i]), minimumV);
            // if (minimumN == 6) {
            //     fout << error << endl;
            // }
            // totalError += fabs(weightEff[i] - weightSimple[i]) / fabs(weightSimple[i]);
            // if 
        }
        totalError /= dataCount;
        totalErrorSum += totalError;
    }

    totalAverageN /= loop;
    simpleTime /= loop;
    effTime /= loop;
    totalErrorSum /= loop;
    
    fout << "total time for efficient update on " << dataCount << " steps with average n = " << totalAverageN << " is " << effTime << "ms." << endl;
    fout << "total time for simple update on " << dataCount << " steps with average n = " << totalAverageN << " is " << simpleTime << "ms." << endl;
    fout << "average weight error between efficient update and simple update is " << totalErrorSum << endl;

    // fout << "minimum absolute weight = " << minimumV << ", maximum error = " << maximumError << endl;

    fout.close();
}

int main(int argc, char **argv) {
    initRandom();
    vector<int> averageNs = {10, 15, 30, 50, 100, 300, 500, 1000, 3000};
    // vector<int> averageNs = {10, 15, 30};
    int loop = 1;
    if (argc > 1) {
        loop = atoi(argv[1]);
    }
    // testConfigurationCorrectness(1000, 30);
    for (int averageN: averageNs) {
        int low = int(averageN - sqrt(averageN));
        int high = int(averageN + sqrt(averageN));
        testEfficiency(100000, low, high, averageN, loop);
    }
    return 0;
}
