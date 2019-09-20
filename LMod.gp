LMod(N,H,p)= 
{ /* Local factor Lp of the modular curve X_H(N) */
	my(x='x,y='y,t='t,L=1,L1,G,GChi,chi,S,d,o,Z,z,Tp,ep);
	G = znstar(N,1); \\ (Z/NZ)*
	GChi = apply(c->znchar([G,c])[2],chargalois(G)); \\ Galois orbits of chars, old and new
	if(H==0,
		GChi = [GChi[1]]
	,
		H1=if(H==1,[-1],concat([-1],H)); \\ Want these elements in Ker chi
	);
	for(i=1,#GChi,
		chi = GChi[i];
		if(H,
			for(j=1,#H1,
				if(chareval(G,chi,H1[j]),
					print("Dropping chi=",chi, " because chi(",H1[j],")!=1");
					next(2)
				)
			)
		); \\ Drop chi if h not in Ker chi
		S = mfinit([N,2,[G,chi]],1);
		d = mfdim(S);
		print("Dropping chi=",chi, " because dim=0");
		if(d==0,next);
		o = charorder(G,chi);
		print("Char ",chi," , order ",o," , dim ",d);
		Z = polcyclo(o,t);
		z = Mod(t,Z);
		ep = chareval(G,chi,p,[z,o]);
		Tp = mfheckemat(S,p);
		L1 = polresultant(charpoly(Tp,y),x^2-y*x+p*ep,y);
		L1 = polresultant(Z,liftpol(L1),t);
		print(L1);
		L *= L1;
	);
	L;
}
