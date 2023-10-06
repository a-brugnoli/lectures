from sympy import symbols, Eq, solve, diff, log

r, p, D, a, nu, A1, C1, A2, B2, C2 = symbols('r p D a nu A1 C1 A2 B2 C2')


w_1_AC = - p /(64*D)*r**4 + A1*r**2/4 + C1
w_2_ABC = - p*r**4/(64*D) + p*a**2/(2*D)*(r**2*log(r) - r**2) + A2*r**2/4 + B2*log(r) + C2

# Support at r = a
w_1_AC_at_a = w_1_AC.subs(r, a)

support_w1_at_a = Eq(w_1_AC_at_a, 0)

C1_fun_A1 = solve(support_w1_at_a, C1)
print(f"C1 in terms of A1:  {C1_fun_A1}")
w_1_A = w_1_AC.subs(C1, C1_fun_A1)

w_2_ABC_at_a = w_2_ABC.subs(r, a)

support_w2_at_a = Eq(w_2_ABC_at_a, 0)

C2_fun_AB2 = solve(support_w2_at_a, C2)
print(f"C2 in terms of A2:  {C2_fun_AB2}")
w_2_AB = w_2_ABC.subs(C2, C2_fun_AB2)

# Free condition at r = 2a
dw_2_AB_dr = diff(w_2_AB, r)

ddw_AB_ddr = diff(dw_2_AB_dr, r)

Mrr_2 = - D*(ddw_AB_ddr + nu/r*dw_2_AB_dr)

Mrr_2_at_2a = Mrr_2.subs(r, 2*a)

free_w2_at_2a = Eq(Mrr_2_at_2a, 0)

A2_fun_B2 = solve(free_w2_at_2a, A2)
print(f"A2 in terms of B2:  {A2_fun_B2}")


