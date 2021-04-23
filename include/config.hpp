#pragma once 

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

#include "splay.hpp"

using namespace std;

class SplayWeightConfiguration {
public: 
    vector<pair<bool, double> > c;

    int mJ, mL, mF;
    int mMax;
    vector<double> aJ, aL, aF;
    // we need to maintain:
    // 1. f(n) (calculated once only)
    // 2. \sum_{ij} L(x(|c[i].second - c[j].second|))
    // 3. \sum_{ij} J(x(|c[i].second - c[j].second|))
    vector<double> pf, spf;
    int n;

    vector<double> fVal;
    vector<vector<double> > cFac;
    int nCounter;

    double weightJL, weightF;
    Splay *momentSplay[2];

    double getSpin(bool s) {
        return s ? 1.0 : (-1.0);
    }

    void prepareCFac() {
        cFac = vector<vector<double> >(mMax + 1);
        for (int i = 0; i <= mMax; ++i) {
            cFac[i] = vector<double>(i + 1);
            cFac[i][0] = 1;
            for (int j = 1; j < i; ++j) {
                cFac[i][j] = cFac[i - 1][j] + cFac[i - 1][j - 1];
            }
            cFac[i][i] = 1;
        }
    }

    void initConfiguration() {
        n = 0; c.clear();
        weightJL = 0.0; weightF = fFunc(0); // must run after randomInit
    }

    void initPF() {
        mMax = max(mJ, mL); 
        pf = vector<double>(mMax + 1, 0.0);
        spf = vector<double>(mMax + 1, 0.0);
    }

    void initF() {
        fVal.clear();
        fVal.push_back(fFunc(0)); nCounter = 0;
    }

    SplayWeightConfiguration(vector<double> &_aJ, vector<double> &_aL, vector<double> &_aF) {
        // init aJ, aL, aF; mJ, mL, mF;

        aJ = _aJ; aL = _aL; aF = _aF;
        mJ = int(aJ.size()) - 1;
        mL = int(aL.size()) - 1;
        mF = int(aF.size()) - 1;

        initConfiguration();
        initPF();
        initF();
        prepareCFac();

        momentSplay[0] = new Splay(mMax);
        momentSplay[1] = new Splay(mMax);
    }

    double polyFunc(double x, int order, vector<double> &a) {
        double base = 1.0;
        double res = 0.0;
        for (int i = 0; i <= order; ++i) {
            res += a[i] * base;
            base *= x;
        }
        return res;
    }
    double JFunc(double dt) {
        return polyFunc(dt, mJ, aJ);
    }
    double LFunc(double dt) {
        return polyFunc(dt, mL, aL);
    }
    double fFunc(int n) {
        return polyFunc(n, mF, aF);
    }

    double deltaF(int n) {
        // fVal[n + 1] - f[n]
        while (nCounter < n + 1) {
            ++nCounter; fVal.push_back(fFunc(nCounter));
        }
        return fVal[n + 1] - fVal[n];
    }


    double bruteForceMBH() {
        if (n == 0) {
            return fFunc(n);
        }
        // -beta H^{eff}
        double res = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                res += (JFunc(fabs(c[i].second - c[j].second)) * getSpin(c[i].first) * getSpin(c[j].first) + LFunc(fabs(c[i].second - c[j].second)));
            }
        }
        res /= n;
        res += fFunc(n);
        return res;
    }

    // our task: to maintain a structure for insert / remove (t, s)
    // so for insert / remove, we can update weight by O(m)

    void addSpinWeight(bool spin, double tau) {
        double base = 1.0;
        for (int i = 0; i <= mMax; ++i) {
            pf[i] += base;
            spf[i] += (spin ? 1 : -1) * base;
            base *= tau;
        }
    }

    void removeSpinWeight(int idx) {
        double fac = c[idx].first ? 1 : -1;
        double tau = c[idx].second;

        double base = 1.0;
        for (int i = 0; i <= mMax; ++i) {
            pf[i] -= base;
            spf[i] -= fac * base;
            base *= tau;
        }
    }

    void removeSpinWeight(bool spin, double tau) {
        double fac = spin ? 1 : -1;

        double base = 1.0;
        for (int i = 0; i <= mMax; ++i) {
            pf[i] -= base;
            spf[i] -= fac * base;
            base *= tau;
        }
    }

    void insertSpin(bool spin, double tau) {
        addSpinWeight(spin, tau);
        c.push_back({spin, tau});
        // double s = getSpin(spin);
        // momentJSplay->insert(getSpin(spin) * tau);
        // cout << "after normal insert" << endl;
        momentSplay[spin]->insert(tau);

        ++n;
    }

    void removeSpin(int idx) {
        removeSpinWeight(idx);
        swap(c[idx], c.back());
        momentSplay[c.back().first]->remove(c.back().second);
        c.pop_back();
        --n;
    }

    double insertDeltaWeightJL(bool spin, double tau) {
        double res = 0.0;
        vector<double> tbases(mMax + 1);
        double base = 1.0;
        for (int i = 0; i <= mMax; ++i) {
            tbases[i] = base; base *= tau;
        }

        double fac = (spin ? 1 : -1);
        double tFac, tiFac;

        vector<double> moment[2];
        vector<double> momentJ(mMax + 1), momentL(mMax + 1);
        for (int i = 0; i < 2; ++i) {
            moment[i] = momentSplay[i]->queryLower(tau);
        }
        for (int i = 0; i <= mMax; ++i) {
            momentJ[i] = moment[1][i] - moment[0][i];
            momentL[i] = moment[1][i] + moment[0][i];
        }
        // cout << "momentJ: ";
        // for (int i = 0; i < momentJ.size(); ++i) {
        //     cout << momentJ[i] << ' ';
        // }
        // cout << endl;

        // cout << "momentL: ";
        // for (int i = 0; i < momentL.size(); ++i) {
        //     cout << momentL[i] << ' ';
        // }
        // cout << endl;
        // moment for J: moment[spin] - momont[spin ^ 1]
        // cout << "after makeing moment" << endl;

        for (int i = 0; i <= mJ; ++i) {
            tFac = (i & 1) ? (-1) : 1;
            for (int j = 0; j + i <= mJ; ++j) {
                tiFac = (j & 1) ? (-1) : 1;
                res += 2 * tFac * aJ[i + j] * cFac[i + j][i] * spf[j] * tbases[i];
                if ((i + j) & 1) {
                    res += 2 * (2 * tiFac) * aJ[i + j] * cFac[i + j][i] * momentJ[j] * tbases[i];
                }
            }
        }
        res *= fac;

        // cout << "res = " << res << endl;

        for (int i = 0; i <= mL; ++i) {
            tFac = (i & 1) ? (-1) : 1;
            // consider the order tau^i in L(tau)
            for (int j = 0; j + i <= mL; ++j) {
                // consider the order-(i + j) term
                tiFac = (j & 1) ? (-1) : 1;
                res += 2 * tFac * aL[i + j] * cFac[i + j][i] * pf[j] * tbases[i];
                if ((i + j) & 1) {
                    res += 2 * (2 * tiFac) * aL[i + j] * cFac[i + j][i] * momentL[j] * tbases[i];
                }
            }
        }
        // cout << "mL = " << mL << endl;
        // cout << "res = " << res << endl;

        // res *= 2;
        res += (aJ[0] + aL[0]);

        // cout << "dWeight = " << res << endl;

        return res;

        // for D weightJL
    }

    double insertDeltaWeightF() {
        return deltaF(n);
    }
    double removeDeltaWeightJL(bool spin, double tau) {
        removeSpinWeight(spin, tau);
        --n; 

        double res = insertDeltaWeightJL(spin, tau);

        ++n;
        addSpinWeight(spin, tau);
        
        return -res;
    }
    double removeDeltaWeightF() {
        return -deltaF(n - 1);
    }

    void updateWeight(double dWeightJL, double dWeightF) {
        weightJL += dWeightJL; weightF += dWeightF;
    }
    double getWeight() {
        return (n == 0) ? weightF : (weightJL / n + weightF);
    }

    void printPF() {
        for (int i = 0; i <= mMax; ++i) {
            cout << pf[i] << ' ';
        }
        cout << endl;
        for (int i = 0; i <= mMax; ++i) {
            cout << spf[i] << ' ';
        }
        cout << endl;
    }
    void printA() {
        cout << "aJ: ";
        for (int i = 0; i < aJ.size(); ++i) {
            cout << aJ[i] << ' ';
        }
        cout << endl;
        cout << "aL: ";
        for (int i = 0; i < aL.size(); ++i) {
            cout << aL[i] << ' ';
        }
        cout << endl;
    }
    // Config toConfig() {
    //     return Config(c);
    // }

    ~SplayWeightConfiguration() {
        if (momentSplay[0]) {
            delete momentSplay[0];
        }
        if (momentSplay[1]) {
            delete momentSplay[1];
        }
    }
};