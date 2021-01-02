read("install.gp");

default(threadsize,200M);
default(threadsizemax,2G);
default(parisize,2G);
default(parisizemax,100G);

S = mfinit([5,6],0); \\ Level 5, weight 6
f = mfeigenbasis(S)[1]; \\ Supersing mod 13
print(Ser(mfcoefs(f,5),'q)); \\ Display q-exp
X = mfgalrep(f,13,[],1000,300,6,1);
\\ Compute Gal rep attached to f mod 13
\\ []: no need to specify which prime above 13 is considered
\\ 1000: work p-adically, where p<=1000 is chosen by the algorithm
\\ 300: choose p-adic accuracy so as to be able to identify rationals of height up to 300
\\ 6: use q-adic precision O(q^6) (technical, governs the number of candidate polynomials to describe the representation)
\\ 1: Isolate rep by Tp instead of Frob_p
X[1]
