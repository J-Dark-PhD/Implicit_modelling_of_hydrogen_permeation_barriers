from sympy import *


def grad(u):
    return diff(u, r)


D1, D2 = symbols("D1, D2")
S1, S2 = symbols("S1, S2")
e, L = symbols("e L")
r, r_0 = symbols("r r_0")
u_0, u_L = symbols("c_0, c_L")
a1, a2, a3, b1, b2, b3 = symbols("a1, a2, a3, b1, b2, b3")

u1 = a1 * log(r) + b1
u2 = a2 * log(r) + b2
u3 = a3 * log(r) + b3

r_interface_1 = r_0 + e
r_interface_2 = r_0 + e + L
list_of_equations = [
    u1.subs(r, r_interface_1) / S1 - u2.subs(r, r_interface_1) / S2,
    u2.subs(r, r_interface_2) / S2 - u3.subs(r, r_interface_2) / S1,
    D1 * grad(u1).subs(r, r_interface_1) - D2 * grad(u2).subs(r, r_interface_1),
    D2 * grad(u2).subs(r, r_interface_2) - D1 * grad(u3).subs(r, r_interface_2),
    u1.subs(r, r_0) - u_0,
    u3.subs(r, r_0 + L + 2 * e) - u_L,
]

# for some reason takes a bit of time...
res = solve(list_of_equations, a1, a2, a3, b1, b2, b3)


def compute_flux(e, r_0, D1, D2, S1, S2, L, c_0):
    new_a = res[a2]
    new_a = new_a.subs("e", e)
    new_a = new_a.subs("r_0", r_0)
    new_a = new_a.subs("D1", D1)
    new_a = new_a.subs("S1", S1)
    new_a = new_a.subs("c_0", c_0)
    new_a = new_a.subs("c_L", 0)
    new_a = new_a.subs("L", L)
    new_a = new_a.subs("D2", D2)
    new_a = new_a.subs("S2", S2)
    flux = -D2 * new_a
    return flux


def compute_PRF(e, r_0, D1, D2, S1, S2, L, pressure):
    flux_with_barrier = compute_flux(e, r_0 - e, D1, D2, S1, S2, L, pressure**0.5 * S1)
    flux_wo_barrier = compute_flux(0, r_0, D2, D2, S2, S2, L, pressure**0.5 * S2)

    return flux_wo_barrier / flux_with_barrier


PRF_expression = simplify(
    compute_PRF(e=e, r_0=r_0, D1=D1, S1=S1, D2=D2, S2=S2, L=L, pressure=Symbol("P"))
)

if __name__ == "__main__":
    print(PRF_expression)
    print(latex(PRF_expression))
