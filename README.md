# ALibm - CPSC 290 Proposal
Fall 2023
Sebastian (Alex) Schott

Proposal for CPSC 290 research project with Prof. Jay Lim 

1 BACKGROUND
  The world runs on real numbers. All scientific endeavors use math libraries to approximate elementary functions (e.g. ln(x) or ex). Unfortunately, the implementations in the most widely used math libraries (e.g. math.h) do not produce the correctly rounded result for all inputs. At each calculation, the error is small, but when compounded multiple times the error can be significant. Although correctly rounded math libraries do exist, they come with performance drawbacks which make them less enticing. In addition, given the importance of floating point calculations in matrix multiplication and machine learning, several new floating point representations have been proposed (e.g. Bfloat16, TensorFloat32 and Posits) which did not have correctly rounded math libraries at all. 
Existing math libraries use polynomial approximations which approximate the real value of f(x). With minimax, the polynomial is made to minimize the maximum error compared to the real value of f(x) across all inputs. However, when the correct real number solution lies directly between two floating point values, the polynomial must be extremely precise to avoid incorrect rounding. In practice, this is ignored, and polynomials will not have the correct rounding for that input. 
The solution posed by Prof. Lim and others is to build a polynomial approximation which approximates the correctly rounded result, which is not equal to the correct real number result. This approach allows more freedom in the polynomial generation, and can produce correctly rounded results for all floating point inputs across multiple representations, while also being faster than the standard library implementations.

2 PROJECT
	My project will be to work with Prof. Lim and find further improvements towards making a fully functional correctly rounded math library. There are two natural continuations for this project. The first would be to add new approximations for functions that have not yet been done. Functions like tanh, ex-e-x, or log(x+1) have not yet been added to the library as they do not have as natural of a range reduction as other elementary functions. The second would be to explore compositions of functions, which currently do not produce the correctly rounded results, even when the individual components do. 
In addition, this CPSC 290 project will lay the groundwork for future projects and  ideally I will be able to continue working with Prof. Lim on similar projects in another semester of CPSC 290, or even as potential for a senior thesis. 

[1] An Approach to Generate Correctly Rounded Math Libraries for New Floating Point Variants.Jay P. Lim, Mridul Aanjaneya, John Gustafson, and Santosh Nagarakatte. 2021.48th ACM SIGPLAN Symposium on Principles of Programming Languages (POPL 2021).

