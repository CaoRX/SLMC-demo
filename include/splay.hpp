#pragma once
#include <iostream>
#include <vector>

using namespace std;
#define LL long long

class Node {
public:
    Node *ch[2];
    Node *f;
    double val;
    int s;
    vector<double> moments;
    int momentN;
    Node(double _val, int _momentN): val(_val), momentN(_momentN) {
        ch[0] = ch[1] = nullptr;
        f = nullptr;
        s = 1; 
        moments = vector<double>(momentN + 1, 0.0);
        moments[0] = 1.0;
        for (int i = 1; i <= momentN; ++i) {
            moments[i] = moments[i - 1] * val;
        }
    }
    void update() {
        s = 1 + chs(0) + chs(1);
        double moment = 1.0;
        for (int i = 0; i <= momentN; ++i) {
            moments[i] = moment + chMoment(0, i) + chMoment(1, i);
            moment *= val;
        }
    }
    void addMoment(double val) {
        double moment = 1.0;
        for (int i = 0; i <= momentN; ++i) {
            moments[i] += moment;
            moment *= val;
        }
    }
    int chs(int chd) {
        if (ch[chd]) {
            return ch[chd]->s;
        } else {
            return 0;
        }
    }
    double chMoment(int chd, int idx) {
        if (ch[chd]) {
            return ch[chd]->moments[idx];
        } else {
            return 0;
        }
    }
    int chi(Node *c) {
        return (c == ch[1]);
    }
};

class Splay {
public:
    Node *rt;
    int momentN;
    Splay(int _momentN): momentN(_momentN) {
        rt = nullptr;
    }
    void rot(Node *x) {
        Node *f = x->f;
        int chd = (f->chi(x)); // x == f->ch[chd]
        if (f->f) {
            f->f->ch[f->f->chi(f)] = x;
        }
        x->f = f->f; f->f = x;
        f->ch[chd] = x->ch[!chd]; 
        if (x->ch[!chd]) {
            x->ch[!chd]->f = f;
        }
        x->ch[!chd] = f;
        x->s = f->s; f->update(); x->update();
    }

    void splay(Node *x, Node *top = nullptr) {
        // rotate x until x->f == top
        // when top = nullptr: 

        Node *f, *ff;

        while (x->f && x->f != top) {
            f = x->f; ff = f->f;
            if (ff == top) {
                rot(x);
            } else if (f->chi(x) ^ ff->chi(f)) {
                rot(x); rot(x);
            } else {
                rot(f); rot(x);
            }
        }

        if (top == nullptr) {
            rt = x;
        }
    }
    void insert(double val) {
        Node *cur = rt;
        Node **loc = &rt;
        Node *f = nullptr;
        while (cur) {
            cur->s += 1;
            cur->addMoment(val);
            if (val > cur->val) {
                loc = &cur->ch[1];
            } else {
                loc = &cur->ch[0];
            }
            f = cur;
            cur = *loc;
        }
        *loc = new Node(val, momentN);
        (*loc)->f = f;
        splay(*loc);
    }
    Node *find(double val) {
        Node *cur = rt;
        while (cur && cur->val != val) {
            if (val < cur->val) {
                cur = cur->ch[0];
            } else {
                cur = cur->ch[1];
            }
        }
        return cur;
    }
    Node *find(Node *subtree, double val) {
        if (!subtree) {
            return nullptr;
        }
        if (subtree->val < val) {
            return find(subtree->ch[1], val);
        } else if (subtree->val > val) {
            return find(subtree->ch[0], val);
        } else {
            Node *lch = find(subtree->ch[0], val);
            return lch ? lch : subtree;
        }
    }
    void remove(double val) {
        Node *rd = find(val);
        splay(rd);
        rt = combine(rd->ch[0], rd->ch[1]);
        delete rd;
    }
    Node *largest(Node *x) {
        Node *cur;
        for (cur = x; cur->ch[1]; cur = cur->ch[1]) ;
        return cur;
    }
    Node *combine(Node *l, Node *r) {
        if (l) {
            l->f = nullptr;
        }
        if (r) {
            r->f = nullptr;
        }
        if (!l) {
            return r;
        }
        if (!r) {
            return l;
        }
        Node *nrt = largest(l);
        splay(nrt);
        nrt->ch[1] = r;
        r->f = nrt;
        nrt->update();
        return nrt;
    }

    pair<Node *, bool> prevX(Node *u, double x) {
        // get the largest key < x
        // if found, bool = true, otherwise false
        if (u == nullptr) {
            return make_pair(nullptr, false);
        }
        if (u->val >= x) {
            return prevX(u->ch[0], x);
        }
        pair<Node *, bool> rPrev = prevX(u->ch[1], x);
        if (rPrev.second) {
            return rPrev;
        } else {
            return make_pair(u, true);
        }
    }

    pair<Node *, bool> nextX(Node *u, double x) {
        // get the smallest key > x
        // if found, bool = true, otherwise false
        if (u == nullptr) {
            return make_pair(nullptr, false);
        }
        if (u->val <= x) {
            return nextX(u->ch[1], x);
        }
        pair<Node *, bool> lNext = nextX(u->ch[0], x);
        if (lNext.second) {
            return lNext;
        } else {
            return make_pair(u, true);
        }
    }

    vector<double> queryLower(double x) {
        if (rt == nullptr) {
            return vector<double>(momentN + 1, 0.0);
        }
        pair<Node *, bool> nxt = nextX(rt, x);
        if (!nxt.second) {
            return rt->moments;
        }
        Node *nxtp = nxt.first;
        splay(nxtp);
        if (nxtp->ch[0]) {
            return nxtp->ch[0]->moments;
        } else {
            return vector<double>(momentN + 1, 0.0);
        }
    }
    void clear() {
        while (rt) {
            rt = combine(rt->ch[0], rt->ch[1]);
        }
    }
};