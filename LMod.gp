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
					if(default(debug),print("Dropping chi=",chi, " because chi(",H1[j],")!=1"));
					next(2)
				)
			)
		); \\ Drop chi if h not in Ker chi
		S = mfinit([N,2,[G,chi]],1);
		d = mfdim(S);
		if(default(debug),print("Dropping chi=",chi, " because dim=0"));
		if(d==0,next);
		o = charorder(G,chi);
		if(default(debug),print("Char ",chi," , order ",o," , dim ",d));
		Z = polcyclo(o,t);
		z = Mod(t,Z);
		ep = chareval(G,chi,p,[z,o]);
		Tp = mfheckemat(S,p);
		L1 = polresultant(charpoly(Tp,y),x^2-y*x+p*ep,y);
		L1 = polresultant(Z,liftpol(L1),t);
		if(default(debug),print(L1));
		L *= L1;
	);
	L;
}

LMod_multi(N,H,plist)=
{ /* Local factor Lp of the modular curve X_H(N) */
  my(x='x,y='y,t='t,np=#plist,L,sort,sort1,P,L1,G,GChi,chi,S,d,o,Z,z,Tp,ep);
	L = vector(np,i,1);
	sort = vecsort(plist,,5);
	sort1 = vector(np);
	for(i=1,np,sort1[sort[i]]=i);
	P = vector(np,i,plist[sort[i]]);
	\\ TODO sort by decr p
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
          if(default(debug),print("Dropping chi=",chi, " because chi(",H1[j],")!=1"));
          next(2)
        )
      )
    ); \\ Drop chi if h not in Ker chi
    S = mfinit([N,2,[G,chi]],1);
    d = mfdim(S);
    if(default(debug),print("Dropping chi=",chi, " because dim=0"));
    if(d==0,next);
    o = charorder(G,chi);
    if(default(debug),print("Char ",chi," , order ",o," , dim ",d));
    Z = polcyclo(o,t);
    z = Mod(t,Z);
    L1 = parapply(p->polresultant(Z,liftpol(polresultant(charpoly(mfheckemat(S,p),y),x^2-y*x+p*chareval(G,chi,p,[z,o]),y)),t),P);
    L = parvector(np,j,L[j]*L1[j]);
  );
  vector(np,i,L[sort1[i]]);
}

