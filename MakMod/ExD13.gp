read("install.gp");

f = mfDelta();
X = mfgalrep(f,13,[],200,50,3);
\\ Compute Gal rep attached to Delta mod 13
\\ []: no need to specify which prime above 13 is considered
\\ 200: work p-adically, where p<=200 is chosen by the algorithm
\\ 50: choose p-adic accuracy so as to be able to identify rationals of height up to 50
\\ 3: use q-adic precision O(q^3) (technical, governs the number of candidate polynomials to describe the representation)
X[1]
