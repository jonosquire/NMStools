# NMStools

Mathematica package NMStools.m contains a number of tools that are useful for nonmodal stability calculations.

1) EigenNDSolve - This is a function to numerically calculate the eigenmodes of a 1-D linear differential equation with general boundary conditions. The syntax is very similar to the native Mathematica NDSolve function, with the symbolic input and boundary conditions automatically discretized internally using Chebyshev polynomials with the tau method.

2) plotEValues - Nice plotting of the spectrum found through EigenNDSolve.

3) pseudoModes - Given the spectrum and the desired inner product, finds the fastest growing linear structures [see Squire & Bhattacharjee, Astrophys. J., 797 67 (2014)]

I also include a tutorial notebook, ExampleNMS, which goes through an example based on the advection equation. If you end up using it for research, I just ask that you cite [Squire & Bhattacharjee, Astrophys. J., 797 67 (2014)].
