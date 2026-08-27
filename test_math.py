import random
import time
import sympy as sp
import calculate

x = sp.Symbol('x')

TEST = 1000

correct = 0
total_time = 0.0

correct_degree = [0] * 11
total_degree = [0] * 11


for i in range(1, 11, 1):
    for _ in range(TEST):
        roots = [random.randint(-100, 100)
                 for _ in range(i)]

        poly = 1
        for r in roots:
            poly *= (x - r)

        poly = sp.expand(poly)
        coeffs = sp.Poly(poly, x).all_coeffs()

        expected = sorted(set(roots))

        degree = len(coeffs) - 1
        total_degree[degree] += 1

        start = time.perf_counter()
        result = calculate.find_polynomial_roots(coeffs)
        elapsed = time.perf_counter() - start
        total_time += elapsed

        if expected == result:
            correct += 1
            correct_degree[degree] += 1
        else:
            print("=" * 50)
            print("Wrong Answer")
            print("Degree   :", degree)
            print("Coeff    :", coeffs)
            print("Expected :", expected)
            print("Found    :", result)
            print("=" * 50)

print("\n========== SUMMARY ==========")
print(f"Correct : {correct}/{TEST*10}")
print(f"Accuracy: {correct/(TEST*10)*100:.2f}%")
print(f"Average Time: {total_time/(TEST*10):.6f} sec\n")

print("Accuracy by Degree")
print("-------------------------------")
print("Degree | Correct | Total | Accuracy")
print("-------------------------------")

for d in range(1, 11):
    if total_degree[d] == 0:
        continue

    acc = 100 * correct_degree[d] / total_degree[d]
    print(f"{d:6d} | {correct_degree[d]:7d} | {total_degree[d]:5d} | {acc:7.2f}%")
