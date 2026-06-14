#include <bits/stdc++.h>
#include "test.hpp"
using namespace std;

using ld = long double;
const ld EPS = 1e-14;
const ld EPS_ROOT = 1e-10;

struct Poly {
    vector<ld> a;

    Poly() {}

    Poly(vector<ld> v) : a(v) {
        normalize();
    }

    // ลบสัมประสิทธิ์นำที่มีค่าเข้าใกล้ 0 ทิ้ง
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

    // Horner's method สำหรับหาค่า P(x)
    ld eval(ld x) const {
        ld res = 0;
        for (ld c : a)
            res = res * x + c;
        return res;
    }
    
    // Horner's method สำหรับหาค่า P'(x) โดยตรง
    ld eval_derivative(ld x) const {
        if (deg() == 0) return 0;
        ld res = 0;
        int n = deg();
        for (int i = 0; i < n; i++) {
            res = res * x + (n - i) * a[i];
        }
        return res;
    }
};

// หาอนุพันธ์ของพหุนาม
Poly derivative(const Poly& P) {
    if (P.deg() == 0) return Poly({0});
    vector<ld> res;
    int n = P.deg();
    for (int i = 0; i < n; i++) {
        res.push_back((n - i) * P.a[i]);
    }
    return Poly(res);
}

// หาเศษเหลือจากการหารพหุนาม (A mod B)
Poly remainder_poly(Poly A, Poly B) {
    if (B.is_zero()) return Poly({0});

    while (A.deg() >= B.deg() && !A.is_zero()) {
        ld factor = A.a[0] / B.a[0];
        vector<ld> AA = A.a;
        
        for (int i = 0; i <= B.deg(); i++) {
            AA[i] -= factor * B.a[i];
        }
        
        AA.erase(AA.begin());
        A = Poly(AA);
    }

    return A;
}

// Sturm Sequence
struct Sturm {
    vector<Poly> seq;

    Sturm() {}

    Sturm(const Poly& P) {
        if (P.is_zero() || P.deg() == 0) return;

        seq.push_back(P);
        seq.push_back(derivative(P));

        while (true) {
            Poly R = remainder_poly(seq[seq.size()-2], seq[seq.size()-1]);

            if (R.is_zero())
                break;

            // กลับเครื่องหมายเศษ R -> -R
            for (auto &x : R.a)
                x = -x;

            seq.push_back(R);
        }
    }

    // นับจำนวนการเปลี่ยนเครื่องหมาย (Sign variations)
    int variations(ld x) const {
        if (seq.empty()) return 0;
        int ans = 0;
        int prev = 0;

        for (const auto &P : seq) {
            ld val = P.eval(x);
            if (fabsl(val) < EPS) continue;

            int cur = (val > 0 ? 1 : -1);
            if (prev != 0 && cur != prev) ans++;
            prev = cur;
        }
        return ans;
    }

    // จำนวนรากในช่วง [l, r]
    int roots_in(ld l, ld r) const {
        return variations(l) - variations(r);
    }
};

// ค้นหารากแบบ Divide and Conquer (Bisection)
void isolate(const Sturm& S, const Poly& P, ld l, ld r, vector<ld>& roots) {
    int cnt = S.roots_in(l, r);

    if (cnt == 0) return;

    // ถ้ารากมีตัวเดียวในช่วงนี้ ใช้ Bisection ในการหาจุดบรรจบ
    if (cnt == 1) {
        ld L = l, R = r;
        ld fL = P.eval(L), fR = P.eval(R);

        // ตรวจสอบรากที่ขอบเขต
        if (fabsl(fL) < EPS_ROOT) { 
            roots.push_back(L); 
            return; 
        }
        if (fabsl(fR) < EPS_ROOT) { 
            roots.push_back(R); 
            return; 
        }

        // Bisection Method
        for (int it = 0; it < 150; it++) {
            ld mid = L + (R - L) / 2.0L;
            ld fM = P.eval(mid);
            
            if (fabsl(fM) < EPS_ROOT) {
                roots.push_back(mid);
                return;
            }
            
            // กรณี simple root (เครื่องหมายต่าง)
            if ((fL > 0) != (fR > 0)) {
                if ((fL > 0) == (fM > 0)) { 
                    L = mid; fL = fM; 
                } else { 
                    R = mid; fR = fM; 
                }
            } 
            // กรณี multiple root (เครื่องหมายเหมือน) ต้องเช็คจำนวนรากย่อย
            else {
                if (S.roots_in(L, mid) == 1) { 
                    R = mid; 
                } else { 
                    L = mid; 
                }
            }
            
            // ตรวจสอบ convergence
            if (R - L < 1e-15) break;
        }
        
        roots.push_back(L + (R - L) / 2.0L);
        return;
    }

    // ถ้ามีมากกว่า 1 ราก ให้แบ่งครึ่งช่วงหาต่อไป (Divide)
    ld mid = l + (r - l) / 2.0L;
    isolate(S, P, l, mid, roots);
    isolate(S, P, mid, r, roots);
}

// การปรับปรุงรากให้แม่นยำยิ่งขึ้นโดยใช้ Newton-Raphson
ld newton_raphson_refine(const Poly& P, ld x0, int max_iter = 50) {
    ld x = x0;
    ld tol = 1e-15;
    
    for (int i = 0; i < max_iter; i++) {
        ld fx = P.eval(x);
        ld dfx = P.eval_derivative(x);
        
        // ป้องกันการหารด้วยศูนย์
        if (fabsl(dfx) < EPS) break;
        
        ld dx = fx / dfx;
        x = x - dx;
        
        if (fabsl(dx) < tol) break;
    }
    
    return x;
}

// ฟังก์ชันหลักในการค้นหาราก
vector<ld> find_root(int n, vector<ld> coef) {
    vector<ld> roots;
    Poly P(coef);
    
    if (P.is_zero() || P.deg() == 0) return roots;

    // คำนวณขอบเขต Cauchy bound ที่ถูกต้อง
    ld max_ratio = 0;
    for (size_t i = 1; i < P.a.size(); i++) {
        max_ratio = max(max_ratio, fabsl(P.a[i] / P.a[0]));
    }
    ld M = 1.0L + max_ratio;

    // สร้าง Sturm Sequence
    Sturm S(P);
    
    // ค้นหารากโดย Divide and Conquer
    isolate(S, P, -M, M, roots);

    // จัดเรียงผลลัพธ์
    sort(roots.begin(), roots.end());
    
    // // ปรับปรุงความแม่นยำด้วย Newton-Raphson
    // for (auto& r : roots) {
    //     r = newton_raphson_refine(P, r, 100);
    // }
    
    // Deduplication: ลบรากที่ซ้ำกันออก
    // vector<ld> unique_roots;
    // for (ld r : roots) {
    //     if (unique_roots.empty() || fabsl(r - unique_roots.back()) > 1e-11) {
    //         unique_roots.push_back(r);
    //     }
    // }

    return roots;
}

// ===== Test =====
int temp() {
    // ปรับแต่ง I/O stream ให้เร็วขึ้น
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    cout << fixed << setprecision(10);
    
    // Test 1: x^3 - 2x + 1 = 0
    {
        cout << "=== Test 1: x^3 - 2x + 1 = 0 ===\n";
        vector<ld> coef = {1, 0, -2, 1};
        auto roots = find_root(3, coef);
        
        cout << "Found " << roots.size() << " roots:\n";
        Poly P(coef);
        for (ld x : roots) cout << "x = " << x << ", P(x) = " << P.eval(x) << "\n";
        cout << "\n";
    }
    
    // Test 2: x^2 - 2 = 0
    {
        cout << "=== Test 2: x^2 - 2 = 0 ===\n";
        vector<ld> coef = {1, 0, -2};
        auto roots = find_root(2, coef);
        
        cout << "Found " << roots.size() << " roots:\n";
        Poly P(coef);
        for (ld x : roots) cout << "x = " << x << ", P(x) = " << P.eval(x) << "\n";
        cout << "\n";
    }
    
    // Test 3: x^4 - 5x^2 + 4 = 0
    {
        cout << "=== Test 3: x^4 - 5x^2 + 4 = 0 ===\n";
        vector<ld> coef = {1, 0, -5, 0, 4};
        auto roots = find_root(4, coef);
        
        cout << "Found " << roots.size() << " roots:\n";
        Poly P(coef);
        for (ld x : roots) cout << "x = " << x << ", P(x) = " << P.eval(x) << "\n";
        cout << "\n";
    }
    
    // Test 4: (x-1)^3 = 0 (triple root)
    {
        cout << "=== Test 4: (x-1)^3 = 0 (triple root) ===\n";
        vector<ld> coef = {1, -3, 3, -1};
        auto roots = find_root(3, coef);
        
        cout << "Found " << roots.size() << " roots:\n";
        Poly P(coef);
        for (ld x : roots) cout << "x = " << x << ", P(x) = " << P.eval(x) << "\n";
        cout << "\n";
    }
    
    // Test 5: x^5 - 5x^3 + 4x = 0
    {
        cout << "=== Test 5: x^5 - 5x^3 + 4x = 0 ===\n";
        vector<ld> coef = {1, 0, -5, 0, 4, 0};
        auto roots = find_root(5, coef);
        
        cout << "Found " << roots.size() << " roots:\n";
        Poly P(coef);
        for (ld x : roots) cout << "x = " << x << ", P(x) = " << P.eval(x) << "\n";
        cout << "\n";
    }
    
    return 0;
}