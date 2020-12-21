read("install.gp");

\\ TODO this does not work: Frob_p does not isolate the rep space for any p.
\\ TODO use Hecke

S = mfinit([5,6],0); \\ Level 7, weight 8
f = mfeigenbasis(S)[1]; \\ Has a companion mod 13
print(Ser(mfcoefs(f,5),'q)); \\ Display q-exp
X = mfgalrep(f,13,[],300,50,3);
\\ Compute Gal rep attached to f mod 13
\\ []: no need to specify which prime above 13 is considered
\\ 300: work p-adically, where p<=300 is chosen by the algorithm
\\ 50: choose p-adic accuracy so as to be able to identify rationals of height up to 50
\\ 3: use q-adic precision O(q^3) (technical, governs the number of candidate polynomials to describe the representation)
X[1]
