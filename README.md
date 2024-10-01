# ALibm 
Fall 2023
Sebastian (Alex) Schott


Full Report: https://docs.google.com/document/d/18aQoqCwM4V8EnPv3WZaG5QA2jsTiSpEXtvwtrYlNPJk/edit?usp=sharing


## Abstract 

The accuracy of mathematical computation in scientific and technological fields depends on the ability of math libraries to evaluate elementary functions (such as ln(x) or ex). Unfortunately, standard math libraries like math.h which are widely used fail to produce the correctly rounded results for all inputs. This rounding issue, while seemingly minor at each calculation, can accumulate to yield significant errors, especially in domains that rely heavily on precise floating point computations. Typically, standard math libraries use polynomial approximations to estimate the real value of each elementary function by minimizing the absolute error across all inputs. However, this technique struggles when the exact real number solution falls on the border of two floating-point values. In such cases, to achieve the correctly rounded result would require extreme precision, which the polynomial cannot handle. RLibm [1] is a math library which solves this problem to produce correctly rounded results for every input. Their method also uses polynomial approximation, but instead of approximating the real number value of the function, they approximate the rounded result. This allows for more flexibility in creating the polynomial, and allows the final polynomial to produce correctly rounded results across various floating point representations. However, there are some elementary functions which RLibm does not implement, such as ex-1 due to difficulty in range reduction. This project extends RLibm to more elementary functions. First I implemented the algorithm described in the original RLibm paper and used it to generate a polynomial approximation of log2(x) and expm1(x) which was correctly rounded on all float inputs over small domains. Using multiple polynomials, I was able to implement log2 in a correctly rounded 32 bit floating point library.


[1] An Approach to Generate Correctly Rounded Math Libraries for New Floating Point Variants.Jay P. Lim, Mridul Aanjaneya, John Gustafson, and Santosh Nagarakatte. 2021.48th ACM SIGPLAN Symposium on Principles of Programming Languages (POPL 2021).
