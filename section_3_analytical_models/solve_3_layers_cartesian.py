from sympy import *

D1, D2 = symbols("D1, D2")
S1, S2 = symbols("S1, S2")
e, L = symbols("e, L")
c_0, c_1 = symbols("c_0, c_1")
a1, a2, a3, b1, b2, b3 = symbols("a1, a2, a3, b1, b2, b3")
x_0 = 0
x = Symbol("x")

u1 = a1 * x + b1
u2 = a2 * x + b2
u3 = a3 * x + b3


list_of_equations = [
    u1.subs(x, e) / S1 - u2.subs(x, e) / S2,
    u2.subs(x, e + L) / S2 - u3.subs(x, e + L) / S1,
    D1 * diff(u1, x).subs(x, e) - D2 * diff(u2, x).subs(x, e),
    D2 * diff(u2, x).subs(x, e + L) - D1 * diff(u3, x).subs(x, e + L),
    u1.subs(x, x_0) - c_0,
    u3.subs(x, L + (2 * e)) - c_1,
]


res = solve(list_of_equations, a1, a2, a3, b1, b2, b3)

print(f"a1 = {latex(res[a1])}")
print(f"a2 = {latex(res[a2])}")
print(f"a3 = {latex(res[a3])}")
print(f"b1 = {latex(res[b1])}")
print(f"b2 = {latex(res[b2])}")
print(f"b3 = {latex(res[b3])}")


def compute_flux(e, D1, D2, S1, S2, L, c_0, c_1):
    new_a = res[a2]
    new_a = new_a.subs("e", e)
    new_a = new_a.subs("D1", D1)
    new_a = new_a.subs("S1", S1)
    new_a = new_a.subs("c_0", c_0)
    new_a = new_a.subs("c_1", c_1)
    new_a = new_a.subs("L", L)
    new_a = new_a.subs("D2", D2)
    new_a = new_a.subs("S2", S2)
    flux = -D2 * new_a

    return flux


def compute_PRF(e, D1, D2, S1, S2, L, pressure):
    flux_with_barrier = compute_flux(e, D1, D2, S1, S2, L, pressure**0.5 * S1, 0)
    print(f"flux with barrier = {latex(flux_with_barrier)}")

    flux_wo_barrier = compute_flux(0, D2, D2, S2, S2, L, pressure**0.5 * S2, 0)
    print(f"flux without barrier = {latex(flux_wo_barrier)}")

    return flux_wo_barrier / flux_with_barrier


PRF_expression = simplify(
    compute_PRF(e=e, D1=D1, S1=S1, D2=D2, S2=S2, L=L, pressure=Symbol("P"))
)

if __name__ == "__main__":
    print(f"PRF = {latex(PRF_expression)}")
