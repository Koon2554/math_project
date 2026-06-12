#include <bits/stdc++.h>
using namespace std;

using ld = long double;

const ld EPS = 1e-12;

struct Poly {
    vector<ld> a; // highest degree -> constant

    Poly() {}

    Poly(vector<ld> v) : a(v) {
        normalize();
    }

    void normalize() {
        while (!a.empty() && fabsl(a[0]) < EPS)
            a.erase(a.begin());

        if (a.empty())
            a.push_back(0);
    }

    int deg() const {
        return (int)a.size() - 1;
    }

    ld eval(ld x) const {
        ld res = 0;
        for (ld c : a)
            res = res * x + c;
        return res;
    }
};

Poly derivative(const Poly& P) {
    vector<ld> res;

    int n = P.deg();

    for (int i = 0; i < n; i++) {
        res.push_back(
            (n - i) * P.a[i]
        );
    }

    return Poly(res);
}

Poly remainder_poly(Poly A, Poly B) {

    while (A.deg() >= B.deg() &&
           !(A.a.size() == 1 &&
             fabsl(A.a[0]) < EPS))
    {
        ld factor = A.a[0] / B.a[0];

        int diff = A.deg() - B.deg();

        vector<ld> temp = B.a;

        for (auto &x : temp)
            x *= factor;

        temp.insert(temp.end(), diff, 0);

        int sz = max(
            (int)A.a.size(),
            (int)temp.size()
        );

        vector<ld> AA(sz, 0);
        vector<ld> TT(sz, 0);

        int shiftA = sz - A.a.size();
        int shiftT = sz - temp.size();

        for (int i = 0; i < A.a.size(); i++)
            AA[i + shiftA] = A.a[i];

        for (int i = 0; i < temp.size(); i++)
            TT[i + shiftT] = temp[i];

        for (int i = 0; i < sz; i++)
            AA[i] -= TT[i];

        A = Poly(AA);
    }

    return A;
}

struct Sturm {

    vector<Poly> seq;

    Sturm(const Poly& P) {

        seq.push_back(P);
        seq.push_back(derivative(P));

        while (true) {

            Poly R =
                remainder_poly(
                    seq[seq.size()-2],
                    seq[seq.size()-1]
                );

            bool zero = true;

            for (auto x : R.a)
                if (fabsl(x) > EPS)
                    zero = false;

            if (zero)
                break;

            for (auto &x : R.a)
                x = -x;

            seq.push_back(R);
        }
    }

    int variations(ld x) const {

        int ans = 0;

        int prev = 0;

        for (auto &P : seq) {

            ld val = P.eval(x);

            if (fabsl(val) < EPS)
                continue;

            int cur =
                (val > 0 ? 1 : -1);

            if (prev && cur != prev)
                ans++;

            prev = cur;
        }

        return ans;
    }

    int roots_in(ld l, ld r) const {
        return variations(l) -
               variations(r);
    }
};

vector<ld> roots;

void isolate(
    const Sturm& S,
    const Poly& P,
    ld l,
    ld r
) {

    int cnt = S.roots_in(l, r);

    if (cnt == 0)
        return;

    if (cnt == 1) {

        ld L = l;
        ld R = r;

        ld fL = P.eval(L);
        ld fR = P.eval(R);

        if (fabsl(fL) < EPS) {
            roots.push_back(L);
            return;
        }

        if (fabsl(fR) < EPS) {
            roots.push_back(R);
            return;
        }

        for (int it = 0; it < 200; it++) {

            ld mid =
                (L + R) / 2.0L;

            ld fM =
                P.eval(mid);

            if ((fL > 0) ==
                (fM > 0))
            {
                L = mid;
                fL = fM;
            }
            else {
                R = mid;
            }
        }

        roots.push_back(
            (L + R) / 2.0L
        );

        return;
    }

    ld mid = (l + r) / 2.0L;

    isolate(S, P, l, mid);
    isolate(S, P, mid, r);
}

int main() {

    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;

    vector<ld> coef(n + 1);

    for (auto &x : coef)
        cin >> x;

    Poly P(coef);

    ld M = 1;

    for (int i = 1; i <= n; i++) {
        M = max(
            M,
            fabsl(coef[i] / coef[0])
        );
    }

    M += 1;

    Sturm S(P);

    isolate(
        S,
        P,
        -M,
        M
    );

    sort(
        roots.begin(),
        roots.end()
    );

    cout << roots.size() << '\n';

    cout << fixed
         << setprecision(12);

    for (auto x : roots)
        cout << x << ' ';

    cout << '\n';
}