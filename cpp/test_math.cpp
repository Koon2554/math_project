#include <bits/stdc++.h>
#include "test.hpp"
using namespace std;
#define ld long double
const ld EPS = 1e-11;

class RandomGenerator {
private:
    mt19937 engine;
public:
    RandomGenerator() : engine(random_device{}()) {}

    ld getRandomInt(ld min, ld max) {
        uniform_int_distribution<int> dist(min, max);
        return dist(engine);
    }
};

vector<ld> generatePolynomial(const vector<ld>& roots) {
    ld n = roots.size();
    vector<ld> poly = {1, -roots[0]};

    for (size_t i = 1; i < roots.size(); ++i) {
        vector<ld> next_poly(poly.size() + 1, 0.0);
        ld root = roots[i];
        
        for (size_t j = 0; j < poly.size(); ++j) {
            next_poly[j] += poly[j];
            next_poly[j + 1] -= poly[j] * root;
        }
        poly = next_poly;
    }
    return poly;
}

int test() {
    RandomGenerator r_1;
    ld cnt = r_1.getRandomInt(1, 5);
    vector<ld> root_res(cnt);
    for (ld i = 0; i < cnt; i++) {
        root_res[i] = r_1.getRandomInt(-10, 10);
    }

    // ld cnt = 3;
    vector<ld> p;

    sort(root_res.begin(), root_res.end());

    // vector<ld> root_res = {-4, 1, 4};
    cout << '\n';
    p = generatePolynomial(root_res);

    root_res.erase(unique(root_res.begin(), root_res.end()), root_res.end());
    for (auto i : root_res) cout << i << ' ';
    cout << '\n';
    for (auto i : p) cout << i << ' ';
    cout << '\n';

    vector<ld> ans = find_root(cnt, p);
    for (auto i : ans) cout << i << ' ';
    cout << '\n';
    if ((ld)ans.size() != (ld)root_res.size()) return 0;
    // bool c = true;
    for (ld i = 0; i < (ld)ans.size(); i++) {
        if (abs(ans[i] - root_res[i]) > EPS) {
            cout << '\n';
            return 0;
        }

        cout << ans[i] - root_res[i] << ' ';
    }
    cout << '\n';
    return 1;
}

int main() {
    ld cnt = 0, max = 5;
    for (ld i = 0; i < max; i++) {
        // cnt += test
        cout << test() << '\n';
    }

    // cout << cnt;
    return 0;
}