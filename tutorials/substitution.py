import sympy

# Define symbols (variables)
r, a = sympy.symbols('r a')

# Define the expression f(r)
expression = r**2 + 3*r - 5

# Substitute a for r in the expression
expression_in_terms_of_a = expression.subs(r, a)

# Print the expression in terms of a
print("Expression in terms of a:", expression_in_terms_of_a)
