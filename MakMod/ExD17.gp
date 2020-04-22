read("install.gp");

f = mfDelta();
X = mfgalrep(f,17,[],1000,150,6);
\\ Compute Gal rep attached to Delta mod 17
\\ []: no need to specify which prime above 17 is considered
\\ 1000: work p-adically, where p<=1000 is chosen by the algorithm
\\ 150: choose p-adic accuracy so as to be able to identify rationals of height up to 150
\\ 6: use q-adic precision O(q^6) (technical, governs the number of candidate polynomials to describe the representation)
X[1]


