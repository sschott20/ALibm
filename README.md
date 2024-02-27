# ALibm 
Fall 2023
Sebastian (Alex) Schott


Full Report: https://docs.google.com/document/d/e/2PACX-1vRlNT80sajOJU-dzRMy9EdidS_jz_HgkvCDGyu_IAHdwPP8Ai1_ZGORA_tk_nn_74wZqJpdJWkGJekf/pub


## BACKGROUND
The world runs on real numbers. All scientific endeavors use math libraries to approximate elementary functions (e.g. ln(x) or ex). Unfortunately, the implementations in the most widely used math libraries (e.g. math.h) do not produce the correctly rounded result for all inputs. At each calculation, the error is small, but when compounded multiple times the error can be significant. 

Although correctly rounded math libraries do exist, they come with performance drawbacks which make them less enticing. In addition, given the importance of floating point calculations in matrix multiplication and machine learning, there are many modern floating point representations (e.g. Bfloat16, TensorFloat32 and Posits) which are commonly used instead of IEEE and do not have correctly rounded math libraries at all. 


## PROJECT

Existing math libraries use polynomial approximations which approximate the real value of f(x). With minimax, the polynomial is made to minimize the maximum error compared to the real value of f(x) across all inputs. However, when the correct real number solution lies directly between two floating point values, the polynomial must be extremely precise to avoid incorrect rounding. In practice, this is ignored, and polynomials will not have the correct rounding for that input. 
The solution posed by Prof. Lim and others is to build a polynomial approximation which approximates the correctly rounded result, which is not equal to the correct real number result. 

This approach allows more freedom in the polynomial generation, and can produce correctly rounded results for all floating point inputs across multiple representations, while also being faster than the standard library implementations.


[1] An Approach to Generate Correctly Rounded Math Libraries for New Floating Point Variants.Jay P. Lim, Mridul Aanjaneya, John Gustafson, and Santosh Nagarakatte. 2021.48th ACM SIGPLAN Symposium on Principles of Programming Languages (POPL 2021).
