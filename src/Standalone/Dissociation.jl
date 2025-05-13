# Symbolic Manipulation to Solve for pO₂
using Symbolics
using Nemo

# Define the symbolic variables
@variables pO₂ XO₂ K₁ α₁ XCO₂ K₂ α₂ β₁ β₂



# Original equation:
# pO₂ = X_O₂ * K₁ * (1 + α₁ * (X_CO₂ * (K₂ * (1 + α₂ * pO₂)) / (1 + β₂ * pO₂))) /
#                        (1 + β₁ * (X_CO₂ * (K₂ * (1 + α₂ * pO₂)) / (1 + β₂ * pO₂)))

# Clear the fractions
lhs = pO₂ * (1 + β₁ * XCO₂ * K₂ * (1 + α₂ * pO₂) / (1 + β₂ * pO₂))
rhs = XO₂ * K₁ * (1 + α₁ * XCO₂ * K₂ * (1 + α₂ * pO₂) / (1 + β₂ * pO₂))

# Expand the terms
expanded_lhs = expand(lhs)
expanded_rhs = expand(rhs)

# Move everything to one side
equation = expanded_lhs - expanded_rhs ~ 0

# Solve symbolically
solutions = symbolic_solve(equation, pO₂)

# Expand and collect the solutions
full_solution = expand(solutions)

# Display the solutions
println("Symbolic Solution for pO₂:")
@show full_solution[1]

"""
CO₂
"""

@variables pCO₂ XO₂ K₁ α₁ XCO₂ K₂ α₂ β₁ β₂


lhs = pCO₂ * (1 + β₂ * XO₂ * K₁ * (1 + α₁ * pCO₂) / (1 + β₁ * pCO₂))
rhs = XCO₂ * K₂ * (1 + α₂ * XO₂ * K₁ * (1 + α₁ * pCO₂) / (1 + β₁ * pCO₂))

# Expand the terms
expanded_lhs = expand(lhs)
expanded_rhs = expand(rhs)

# Move everything to one side
equation = expanded_lhs - expanded_rhs ~ 0

# Solve symbolically
solutions = symbolic_solve(equation, pCO₂)

# Expand and collect the solutions
full_solution = expand(solutions)

# Display the solutions
println("Symbolic Solution for pCO₂:")
@show full_solution[1]

"""Helper Functions"""




function pressure2conc(pO₂, pCO₂; CₛₐₜO₂=(9*0.0227), CₛₐₜCO₂=(86.11*0.0227), h₁=0.3836, h₂=1.819, K₁=14.99, K₂=194.4, α₁=0.03198, α₂=0.05591, β₁=0.008275, β₂=0.03255)
    # Perform calculations using inputs and parameters
    XO₂ = pO₂ * (1 + β₁ * pCO₂) / (K₁ * (1 + α₁ * pCO₂))
    XCO₂ = pCO₂ * (1 + β₂ * pO₂) / (K₂ * (1 + α₂ * pO₂))

    cO₂ = CₛₐₜO₂ * (XO₂)^(1/h₁)/(1 + (XO₂)^(1/h₁))
    cCO₂ = CₛₐₜCO₂ * (XCO₂)^(1/h₂)/(1 + (XCO₂)^(1/h₂))

    return cO₂, cCO₂
end

# Example usage:
a, b = pressure2conc(40, 45)
println("cO₂: ", a)
println("cCO₂: ", b)

function conc2pressure(cO₂, cCO₂; CₛₐₜO₂=(9*0.0227), CₛₐₜCO₂=(86.11*0.0227), h₁=0.3836, h₂=1.819, K₁=14.99, K₂=194.4, α₁=0.03198, α₂=0.05591, β₁=0.008275, β₂=0.03255)
    # Perform calculations using inputs and parameters


    XO₂ = ((cO₂ / CₛₐₜO₂) / (1 - (cO₂ / CₛₐₜO₂)))^h₁
    XCO₂ = ((cCO₂ / CₛₐₜCO₂) / (1 - (cCO₂ / CₛₐₜCO₂)))^h₂

    pO₂ = (-1 + √(1 + 2K₁*XO₂*β₂ + 2K₂*XCO₂*β₁ + (K₁^2)*(XO₂^2)*(β₂^2) - 2K₁*K₂*XCO₂*XO₂*α₁*α₂ + 4K₁*K₂*XCO₂*XO₂*α₁*β₂ + 4K₁*K₂*XCO₂*XO₂*α₂*β₁ - 2K₁*K₂*XCO₂*XO₂*β₁*β₂ + (K₂^2)*(XCO₂^2)*(β₁^2) + 2(K₁^2)*K₂*XCO₂*(XO₂^2)*α₁*α₂*β₂ + 2K₁*(K₂^2)*(XCO₂^2)*XO₂*α₁*α₂*β₁ + (K₁^2)*(K₂^2)*(XCO₂^2)*(XO₂^2)*(α₁^2)*(α₂^2)) + K₁*XO₂*β₂ - K₂*XCO₂*β₁ + K₁*K₂*XCO₂*XO₂*α₁*α₂) / (2β₂ + 2K₂*XCO₂*α₂*β₁)

    pCO₂ = (-1 + √(1 + 2K₁*XO₂*β₂ + 2K₂*XCO₂*β₁ + (K₁^2)*(XO₂^2)*(β₂^2) - 2K₁*K₂*XCO₂*XO₂*α₁*α₂ + 4K₁*K₂*XCO₂*XO₂*α₁*β₂ + 4K₁*K₂*XCO₂*XO₂*α₂*β₁ - 2K₁*K₂*XCO₂*XO₂*β₁*β₂ + (K₂^2)*(XCO₂^2)*(β₁^2) + 2(K₁^2)*K₂*XCO₂*(XO₂^2)*α₁*α₂*β₂ + 2K₁*(K₂^2)*(XCO₂^2)*XO₂*α₁*α₂*β₁ + (K₁^2)*(K₂^2)*(XCO₂^2)*(XO₂^2)*(α₁^2)*(α₂^2)) - K₁*XO₂*β₂ + K₂*XCO₂*β₁ + K₁*K₂*XCO₂*XO₂*α₁*α₂) / (2β₁ + 2K₁*XO₂*α₁*β₂)

    # cO₂ = CₛₐₜO₂ * (XO₂)^(1/h₁)/(1 + (XO₂)^(1/h₁))
    # cCO₂ = CₛₐₜCO₂ * (XCO₂)^(1/h₂)/(1 + (XCO₂)^(1/h₂))
    return pO₂, pCO₂
end

c, d = conc2pressure(a, b)