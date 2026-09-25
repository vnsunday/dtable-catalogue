
# Sampling
Random samples


Definition:
    Random sampling size n from Ground-Set having Random variable X is a set variables X1,X2,...,Xn. Sastified
        (a) Indepence statistically
        (b) Having the same statistical distribution to X.

    Sastified X[i] are called random variables who independent and distributed identically.

Definition of Statistic:
    A statistic is a function Y = g(X1,X2,...,Xn).

    Note: a statistic is a measured function, and specified (not including unknown parameter).

# Confidence Interval
Confidence Interval: Interval Estimation.

Input: 
    x1,x2,...,xn follows a Distribution D[θ]
        θ: Unknown parameter
    1 - α: Confidence Level
Output:
    Find (θ1,θ2), an estimation range for θ.
    With a given confidence level = 1 - α.

Terminologies:

Approach:
    Pick (choose) a statistic (a function)
        G=G(x1,x2,...,xn,θ) ~ D2
        G follows a Distribution D2 (without unknown parameter θ).
            => D2 is totally specified.
        Note: To D2, θ is variable rather than Parameter.

    With (1 - α) Confidence Level; choose α1 and α2 such that:
        α1 + α2 = α.
        P(G < g[α1]) = α1
        P(G > g[1-α2]) = α2
        <=>
        P(g[α1] < G(x1,x2,...,xn,θ) < g[1 - α2]) = 1 - α1 - α2 = 1 - α
        Call θ1=g[α1]; θ2 = g[1 - α2]
        We have 
            p(θ1 < θ < θ2) = 1 - α

## Confidence Interval for Mean
Input: 
    X ~ N(μ,σ²)
    Observing samples (originally by X) (known values):
        x1, x2, ..., xn
    1 - α: Confidence Level
    Where
        μ: unknown
        σ = σ0: known
        1 - α: Known

Question:
    Interval estimation (μ1, μ2) for (1 - α) Confidence Level.
    
Solution:
    X¯ = (x1 + x2 + ... + xn) / n

    Choose the statistic.
        G = Z = 
            (X¯-μ)√n/σ0

    By the proved 
        Z ~ N(0,1)

    By the proved result.


