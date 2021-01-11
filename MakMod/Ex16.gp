read("install.gp");

S=mfinit([16,2,Mod(5,16)],0); \\ New cuspidal space of level 16, weight 2, nebentypus Conrey index 5
f=mfeigenbasis(S)[1]; \\ The unique eigenform (dim S = 1)
X = mfgalrep(f,5,[[3,1]],[30,100],10,3,0,1);
\\ Compute Gal rep attached to f mo the prime l|5 such that a3 = 1 mod l
\\ [30,100]: work p-adically, where 30 <= p <= 100 is chosen by the algorithm
\\ 10: choose p-adic accuracy so as to be able to identify rationals of height up to 10
\\ 3: use q-adic precision O(q^3) (technical, governs the number of candidate polynomials to describe the representation)
\\ 0: Do not compute extra data for Tp action (not needed in this example, and slows things down)
\\ 1: Use only 1 elliptic curve to contruct modular Jacobian (default; more curves not needed in this example, would slow things down)
X[1]


