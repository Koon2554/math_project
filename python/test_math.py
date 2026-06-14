import random
import sympy as sp
import calculate as cal

x = sp.Symbol('x')

for _ in range(10):
    roots = [random.randint(-10,10)
             for _ in range(random.randint(1,10))]

    poly = 1
    for r in roots:
        poly *= (x-r)

    poly = sp.expand(poly)

    coeffs = sp.Poly(poly, x).all_coeffs()

    print("roots =", sorted(roots))
    print("coeff =", coeffs)
    res = cal.find_roots_lills_method(coeffs)
    print("res   =", res)
    if sorted(roots) == res:
        print("we can it!!")
    else :
        print("OMG")
    # print("sympy =", sp.nroots(poly))
    print()