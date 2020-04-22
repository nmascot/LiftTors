read("install.gp");

f = mfDelta();
X = mfgalrep(f,19,[],1000,250,6);
\\ Compute Gal rep attached to Delta mod 19
\\ []: no need to specify which prime above 19 is considered
\\ 1000: work p-adically, where p<=1000 is chosen by the algorithm
\\ 250: choose p-adic accuracy so as to be able to identify rationals of height up to 250
\\ 6: use q-adic precision O(q^6) (technical, governs the number of candidate polynomials to describe the representation)
X[1]


