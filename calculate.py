from fractions import Fraction

ZERO_TOL = Fraction(0)        
MIN_WIDTH = Fraction(1, 10**7)
 
def remove_leading_zeros(c):
    c = list(c)
    while len(c) > 1 and c[0] == 0:
        c.pop(0)
    return c


def is_zero_poly(p):
    return all(c == 0 for c in p)


def poly_eval(poly, x):
    val = Fraction(0)
    for c in poly:
        val = val * x + c
    return val


def poly_deriv(poly):
    n = len(poly) - 1
    if n <= 0:
        return [Fraction(0)]
    return [poly[i] * (n - i) for i in range(n)]


def poly_div_remainder(num, den):
    num = remove_leading_zeros(num)
    den = remove_leading_zeros(den)
    if len(num) < len(den):
        return num

    num = list(num)
    len_n, len_d = len(num), len(den)
    for i in range(len_n - len_d + 1):
        t = num[i] / den[0]
        for j in range(len_d):
            num[i + j] -= t * den[j]

    rem = num[len_n - len_d + 1:]
    return remove_leading_zeros(rem) if rem else [Fraction(0)]


def poly_div_exact(f, d):
    f = remove_leading_zeros(list(f))
    d = remove_leading_zeros(list(d))
    quotient = []
    while len(f) >= len(d):
        lead = f[0] / d[0]
        quotient.append(lead)
        for j in range(len(d)):
            f[j] -= lead * d[j]
        f.pop(0)
    return quotient if quotient else [Fraction(0)]


def poly_gcd(f, g):
    f = remove_leading_zeros(f)
    g = remove_leading_zeros(g)
    while not is_zero_poly(g):
        f, g = g, poly_div_remainder(f, g)
    return f


def make_squarefree(coeffs):
    f = remove_leading_zeros(coeffs)
    if len(f) <= 1:
        return f
    fp = poly_deriv(f)
    if is_zero_poly(fp):
        return f
    g = poly_gcd(f, fp)
    if len(g) <= 1:
        return f
    return remove_leading_zeros(poly_div_exact(f, g))



def cauchy_root_bound(coeffs):
    coeffs = remove_leading_zeros(coeffs)
    a0 = abs(coeffs[0])
    if len(coeffs) == 1:
        return Fraction(1)
    max_ratio = max(abs(c) for c in coeffs[1:]) / a0
    return Fraction(1) + max_ratio



def build_sturm_sequence(P):
    P = remove_leading_zeros(P)
    seq = [P, poly_deriv(P)]
    while True:
        rem = poly_div_remainder(seq[-2], seq[-1])
        rem = [-c for c in rem]
        if is_zero_poly(rem):
            break
        seq.append(rem)
        if len(seq) > len(P) + 5:   
            break
    return seq


def count_sign_variations(seq, x):
    count = 0
    last_sign = 0
    for p in seq:
        v = poly_eval(p, x)
        if v != 0:
            sign = 1 if v > 0 else -1
            if last_sign != 0 and sign != last_sign:
                count += 1
            last_sign = sign
    return count


def sturm_root_count(seq, a, b):
    nudge = Fraction(1, 10**6)
    if poly_eval(seq[0], a) == 0:
        a += nudge
    if poly_eval(seq[0], b) == 0:
        b -= nudge
    return count_sign_variations(seq, a) - count_sign_variations(seq, b)



def find_roots_in_interval(seq, L, R, min_width, known_k=-1):
    k = sturm_root_count(seq, L, R) if known_k == -1 else known_k
    if k == 0:
        return []
    if (R - L) < min_width:
        return [(L + R) / 2]

    mid = (L + R) / 2
    k_left = sturm_root_count(seq, L, mid)
    left = find_roots_in_interval(seq, L, mid, min_width, k_left)
    right = find_roots_in_interval(seq, mid, R, min_width, k - k_left)
    return left + right

def find_polynomial_roots(coeffs):
    coeffs = [Fraction(c) for c in coeffs]
    coeffs = remove_leading_zeros(coeffs)
    if len(coeffs) == 1:
        return []

    coeffs_sf = make_squarefree(coeffs)   
    M = cauchy_root_bound(coeffs)
 
    L = -(M + 1) - Fraction(1, 3)
    R = (M + 1) + Fraction(2, 7)

    seq = build_sturm_sequence(coeffs_sf) 
    roots = find_roots_in_interval(seq, L, R, MIN_WIDTH)
    return sorted(round(float(r)) for r in roots)
