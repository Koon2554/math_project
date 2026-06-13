#include <bits/stdc++.h>
using namespace std;

using ld = long double;
const ld EPS = 1e-11; // ปรับค่า EPS ให้เหมาะสมกับการคำนวณซ้ำซ้อน

struct Poly {
    vector<ld> a;

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

    bool is_zero() const {
        return a.size() == 1 && fabsl(a[0]) < EPS;
    }

    ld eval(ld x) const {
        ld res = 0;
        for (ld c : a)
            res = res * x + c;
        return res;
    }
};

Poly derivative(const Poly& P) {
    if (P.deg() == 0) return Poly({0}); // ดักจับการดิฟค่าคงที่
    vector<ld> res;
    int n = P.deg();
    for (int i = 0; i < n; i++) {
        res.push_back((n - i) * P.a[i]);
    }
    return Poly(res);
}

Poly remainder_poly(Poly A, Poly B) {
    if (B.is_zero()) return Poly({0}); // ป้องกัน Division by Zero

    // บังคับให้ Degree ลดลงเสมอในแต่ละรอบ ป้องกัน Infinite Loop
    while (A.deg() >= B.deg() && !A.is_zero()) {
        ld factor = A.a[0] / B.a[0];
        vector<ld> AA = A.a;
        
        for (int i = 0; i <= B.deg(); i++) {
            AA[i] -= factor * B.a[i];
        }
        
        // ลบพจน์นำหน้าทิ้งอย่างเด็ดขาดเพื่อหลีกเลี่ยงปัญหาเศษทศนิยมตกค้าง
        AA.erase(AA.begin());
        A = Poly(AA);
    }

    return A;
}

struct Sturm {
    vector<Poly> seq;

    Sturm(const Poly& P) {
        if (P.is_zero() || P.deg() == 0) return;

        seq.push_back(P);
        seq.push_back(derivative(P));

        while (true) {
            Poly R = remainder_poly(seq[seq.size()-2], seq[seq.size()-1]);

            if (R.is_zero())
                break;

            for (auto &x : R.a)
                x = -x;

            seq.push_back(R);
        }
    }

    int variations(ld x) const {
        int ans = 0;
        int prev = 0;

        for (const auto &P : seq) {
            ld val = P.eval(x);
            if (fabsl(val) < EPS) continue;

            int cur = (val > 0 ? 1 : -1);
            if (prev && cur != prev) ans++;
            prev = cur;
        }
        return ans;
    }

    int roots_in(ld l, ld r) const {
        return variations(l) - variations(r);
    }
};

void isolate(const Sturm& S, const Poly& P, ld l, ld r, vector<ld>& roots) {
    int cnt = S.roots_in(l, r);

    if (cnt == 0) return;

    if (cnt == 1) {
        ld L = l, R = r;
        ld fL = P.eval(L), fR = P.eval(R);

        // กรณีเจอรากพอดีที่ขอบเขต
        if (fabsl(fL) < EPS) { roots.push_back(L); return; }
        if (fabsl(fR) < EPS) { roots.push_back(R); return; }

        // Bisection Method
        for (int it = 0; it < 100; it++) {
            ld mid = (L + R) / 2.0L;
            
            // ใช้กฎเครื่องหมายสลับหากเป็นรากเดี่ยว (ทำงานเร็ว O(N))
            if ((fL > 0) != (fR > 0)) {
                ld fM = P.eval(mid);
                if ((fL > 0) == (fM > 0)) { L = mid; fL = fM; }
                else { R = mid; }
            } 
            // หากเครื่องหมายเหมือนกัน (Multiple root ที่สัมผัสแกน X) 
            // ให้ใช้ Sturm ในการไกด์ Bisection แทน
            else {
                if (S.roots_in(L, mid) == 1) { R = mid; } 
                else { L = mid; }
            }
        }
        roots.push_back((L + R) / 2.0L);
        return;
    }

    ld mid = (l + r) / 2.0L;
    isolate(S, P, l, mid, roots);
    isolate(S, P, mid, r, roots);
}

vector<ld> find_root(int n, vector<ld> coef) {
    vector<ld> roots; // เปลี่ยนเป็น Local Variable ป้องกันข้อมูลเก่าปนเปื้อน
    Poly P(coef);
    
    if (P.is_zero() || P.deg() == 0) return roots;

    ld M = 1;
    // ใช้ P.a ที่ทำ normalize แล้ว ป้องกัน coef[0] เป็น 0
    for (size_t i = 1; i < P.a.size(); i++) {
        M = max(M, fabsl(P.a[i] / P.a[0]));
    }
    M += 1;

    Sturm S(P);
    isolate(S, P, -M, M, roots);

    sort(roots.begin(), roots.end());
    
    // กำจัดรากที่ซ้ำกันเกินไปจากการเหลื่อมล้ำของทศนิยม (Deduplication)
    vector<ld> unique_roots;
    for (ld r : roots) {
        if (unique_roots.empty() || fabsl(r - unique_roots.back()) > 1e-7) {
            unique_roots.push_back(r);
        }
    }

    return unique_roots;
}