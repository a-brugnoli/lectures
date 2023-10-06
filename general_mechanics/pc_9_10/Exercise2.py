from sympy import symbols, Eq, solve, diff, log, simplify

r, p, D, a, nu, A1, C1, A2, B2, C2 = symbols('r p D a nu A1 C1 A2 B2 C2')

w_1_A1C1 = - p /(64*D)*r**4 + A1*r**2/4 + C1
w_2_A2B2C2 = - p*r**4/(64*D) + p*a**2/(2*D)*(r**2*log(r) - r**2) + A2*r**2/4 + B2*log(r) + C2

# Support at r = a for w1
w_1_A1C1_at_a = w_1_A1C1.subs(r, a)

support_w1_at_a = Eq(w_1_A1C1_at_a, 0)

C1_fun_A1 = solve(support_w1_at_a, C1)[0]
print(f"C1 in terms of A1:  {C1_fun_A1}")
w_1_A1 = w_1_A1C1.subs(C1, C1_fun_A1)

print(f"w1 in terms of A1: {w_1_A1}")

# Support at r = a for w2
w_2_A2B2C2_at_a = w_2_A2B2C2.subs(r, a)

support_w2_at_a = Eq(w_2_A2B2C2_at_a, 0)

C2_fun_AB2 = solve(support_w2_at_a, C2)[0]
print(f"C2 in terms of A2 and B2:  {C2_fun_AB2}")
w_2_A2B2 = w_2_A2B2C2.subs(C2, C2_fun_AB2)

print(f"w2 in terms of A2 and B2: {w_2_A2B2}")

# Free condition at r = 2a
dw_2_A2B2_dr = diff(w_2_A2B2, r)
ddw_2_A2B2_ddr = diff(dw_2_A2B2_dr, r)

Mrr_2 = - D*(ddw_2_A2B2_ddr + nu/r*dw_2_A2B2_dr)

Mrr_2_at_2a = Mrr_2.subs(r, 2*a)

free_w2_at_2a = Eq(Mrr_2_at_2a, 0)

A2_fun_B2 = solve(free_w2_at_2a, A2)[0]
print(f"A2 in terms of B2:  {A2_fun_B2}")

# Continuity rotation in r = a
dw_1_A1_dr = diff(w_1_A1, r)

dw_1_A1_dr_at_a = dw_1_A1_dr.subs(r, a)
print(f"dw1_dr_at_a: {dw_1_A1_dr_at_a}")

w_2_B2 = w_2_A2B2.subs(A2, A2_fun_B2)
dw_2_B2_dr = diff(w_2_B2, r)
dw_2_B2_dr_at_a = dw_2_B2_dr.subs(r, a)
print(f"dw2_dr_at_a: {dw_2_B2_dr_at_a}")

continuity_rot_at_a = Eq(dw_1_A1_dr_at_a, dw_2_B2_dr_at_a)

A1_fun_B2 = solve(continuity_rot_at_a, A1)[0]
print(f"A1 in terms of B2:  {A1_fun_B2}")

# Continuity moment in r = a
w_1_B2 = w_1_A1.subs(A1, A1_fun_B2)

dw_1_B2_dr = diff(w_1_B2, r)
ddw_1_B2_ddr = diff(dw_1_B2_dr, r)
ddw_1_B2_ddr_at_a = ddw_1_B2_ddr.subs(r, a)

ddw_2_B2_ddr = diff(dw_2_B2_dr, r)
ddw_2_B2_ddr_at_a = ddw_2_B2_ddr.subs(r, a)

continuity_mom_at_a = Eq(ddw_1_B2_ddr_at_a, ddw_2_B2_ddr_at_a)

B2_sol = solve(continuity_mom_at_a, B2)[0]

w_1_sol = simplify(w_1_B2.subs(B2, B2_sol))
w_2_sol = simplify(w_2_B2.subs(B2, B2_sol))

print(f"Solution for w1:  {w_1_sol}")
print(f"Solution for w2:  {w_2_sol}")







