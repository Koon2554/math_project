import sympy as sp
from mpmath import mp

mp.dps = 100

def make_poly(coefs, x):
    expr = sum(
        sp.Rational(c) * x ** (len(coefs) - 1 - i)
        for i, c in enumerate(coefs)
    )
    return sp.Poly(expr, x)

def poly_value(coefs, x):
    res = mp.mpf('0')
    for c in coefs:
        res = res * x + mp.mpf(str(c))
    return res

def E(theta, coefs):
    t = mp.tan(theta)
    n = len(coefs) - 1
    s = mp.mpf('0')
    
    for i in range(n):
        s += ((-1) ** i) * mp.mpf(str(coefs[i])) * t ** (n - i)
        
    s += ((-1) ** n) * mp.mpf(str(coefs[n]))
    
    power = (n - 1) // 2
    return ((-1) ** power) * s

def build_sturm(P):
    seq = [P]
    D = P.diff()
    if D.is_zero: return seq
    seq.append(D)
    
    A = P
    B = D
    while not B.is_zero:
        _, R = sp.div(A, B)
        if R.is_zero: break
        R = -R
        seq.append(R)
        A = B
        B = R
    return seq

def count_sign_variations(sturm_seq, x_val):
    x_val = sp.Float(str(x_val), 150)
    prev = None
    count = 0
    for P in sturm_seq:
        val = sp.N(P.eval(x_val), 150)
        if val == 0: continue
        sign = 1 if val > 0 else -1
        if prev is not None and sign != prev:
            count += 1
        prev = sign
    return count

def roots_in_interval(sturm_seq, L, R):
    return count_sign_variations(sturm_seq, L) - count_sign_variations(sturm_seq, R)

def refine_root_via_E(L, R, coefs):
    theta_L = mp.atan(-R)
    theta_R = mp.atan(-L)
    
    fL = E(theta_L, coefs)
    fR = E(theta_R, coefs)
    
    eps = mp.mpf('1e-105')
    
    if abs(fL) <= eps: return -mp.tan(theta_L)
    if abs(fR) <= eps: return -mp.tan(theta_R)
    
    if fL * fR > 0:
        x_L, x_R = L, R
        fx_L = poly_value(coefs, x_L)
        for _ in range(400):
            x_M = (x_L + x_R) / 2
            fx_M = poly_value(coefs, x_M)
            if abs(fx_M) <= eps: return x_M
            if fx_L * fx_M <= 0:
                x_R = x_M
            else:
                x_L = x_M
                fx_L = fx_M
        return (x_L + x_R) / 2

    for _ in range(400):
        theta_M = (theta_L + theta_R) / 2
        fM = E(theta_M, coefs)
        
        if abs(fM) <= eps:
            return -mp.tan(theta_M)
            
        if fL * fM <= 0:
            theta_R = theta_M
            fR = fM
        else:
            theta_L = theta_M
            fL = fM
            
    theta = (theta_L + theta_R) / 2
    return -mp.tan(theta)

def isolate(sturm_seq, coefs, L, R):
    cnt = roots_in_interval(sturm_seq, L, R)
    if cnt == 0: return []
    if cnt == 1:
        root = refine_root_via_E(L, R, coefs)
        return [root]
        
    M = (L + R) / 2
    
    if M == L or M == R:
        return [refine_root_via_E(L, R, coefs)]
        
    left = isolate(sturm_seq, coefs, L, M)
    right = isolate(sturm_seq, coefs, M, R)
    return left + right

def find_roots_lills_method(coefs):
    x = sp.Symbol('x')
    P = make_poly(coefs, x)
    
    sqf_factors = P.sqf_list()[1]
    all_roots = []
    
    for factor, multiplicity in sqf_factors:
        if factor.degree() == 0: continue
        sturm_seq = build_sturm(factor)
        
        coeffs_factor = factor.all_coeffs()
        a0 = abs(coeffs_factor[0])
        mx = max(abs(c) / a0 for c in coeffs_factor[1:]) if len(coeffs_factor) > 1 else 0
        
        M = mp.mpf(str(1 + mx))
        
        unique_roots = isolate(
            sturm_seq,
            list(coeffs_factor),
            -M,
            M
        )
        
        for root in unique_roots:
            for _ in range(multiplicity):
                all_roots.append(root)
                
    all_roots.sort()
    return all_roots
