This package requires the branch "aurel-matmod-ZpXQ" of the Git version of PARI/GP, and the gcc compiler.

It computes Galois representations appearing in the Jacobian of a given curve C, given the char.poly. of the Frobenius at p for p a prime of good reduction of C.

At the moment, only the cases of hyperelliptic curves and of smooth plane curves are implemented.

Two examples are provided:
* ExHyper.gp: the representation is afforded by a piece of the 7-torsion of a hyperelliptic curve, namely the modular curve X1(13) (genus 2).
* ExSmooth.gp: the representation is that afforded by the whole 2-torsion of a smooth plane curve, namely the Klein quartic (genus 3).

These examples are intended to demonstrate the functionalities of this package.

To run them, first compile the package by typing "make all", and then start GP (version = "aurel-matmod-ZpXQ" git branch) and type 
read("ExX.gp")
where X is either "Hyper" or "Smooth".

This package is provided in the hope it will be useful but comes without any guarantee whatsoever.

Nicolas Mascot, September 1, 2018
