This package computes Galois representations appearing in the Jacobian of a given curve C, given the char.poly. of the Frobenius at p for p a prime of good reduction of C. See my preprint "Hensel-lifting torsion points on Jacobians and Galois representations" (arXiv 1808.03939) for more information.

It requires PARI/GP (developement version 2.12) and the gcc compiler.

At the moment, hyperelliptic curves, superelliptic curves, and smooth plane curves are supported natively. It is possible to work with general curves, but one must provide Riemann-Roch data.

A few examples are provided:
* ExHyper.gp: the representation is afforded by a piece of the 7-torsion of a hyperelliptic curve, namely the modular curve X1(13) (genus 2).
* ExSmooth.gp: the representation is that afforded by the whole 2-torsion of a smooth plane curve, namely the Klein quartic (genus 3).
* ExSuper.gp: the representation is that afforded by the whole 2-torsion of the superelliptic curve y^4=x^3+x+1 (genus 3).
* ExHyperAlg.gp, ExSmoothAlg.gp, ExSuperAlg.gp demonstrate the syntax when using points on C whire are defined over number fields (this is slow, so avoid this if sufficiently many rational points are available). 

These examples are intended to demonstrate the functionalities of this package.
To run them, first compile the package by typing "make all", and then start GP and type 
read("ExX.gp")
where X is "Hyper" or "Smooth" or... etc.

This package now contains new experimental code to compute Galois representations attached to modular forms. In order to use this functionality, install the main package as above, then go to the MakMod/ subdirectory, type "make all" again, and try the examples Ex16.gp, ExD13.gp, ExD17.gp, or ExD19.gp. See my preprint "Moduli-friendly Eisenstein series over the p-adics and the computation of modular Galois representations" (arXiv 2004.14683) for a description the algorithm.

This package is provided in the hope it will be useful, but comes without any guarantee whatsoever.

Nicolas Mascot, May 6, 2020

