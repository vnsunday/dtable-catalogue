
# 
Table Distribution:
    Table:
        Row1(Label=X): x1, x2, ..., xn
        Row2(Lebel=P): p1, p2, ..., pn
    Where:
        p[i]=p(x[i])=P(X=x[i])

Distribution Function (Cumulative)(CMF)
    F(x) = P(X < x), x ∈ R (Real-Number)

Probability Density function (PDF)
    f(x) = F'(x). (The 1st derivative of F.)
    =>
    F(x) = integral{-∞,x}f(x)dx

# Special Numbers
Expectation: E(X)
    (a) X: a discrete Random variables, probability function p(x)
        E(X) = sum{i}x[i]p[i]
    (b) X: a continous random variables. PDF is f(x)
        E(X) = integral{-∞,∞}xf(x)dx

Variance of X: VX
    Let μ = EX
    VX = E((X - μ)²)
