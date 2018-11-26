\\ The first 4 functions establish a correspondence Fl^d <-> {0,1,...,l^d-1}
\\ and give the action of a matrix and vector addition under this correspondence

c2i(c,l)= \\ Conversion coordinates -> index in an Fl-space of finite dim
{
	my(i);
	i=subst(Pol(c),'x,l);
	if(i,i,l^#c);
}

i2c(i,l,d)= \\ Conversion index -> corrdinates in an Fl-space of dim d
{
	my(j=i,c=vector(d)~);
	for(n=1,d,
		c[d+1-n] = j%l;
		j = j\l
	);
  c;
}

ActOni(M,i,l)= \\ Action of the matrix M on Fl^d, in terms of indices
{
	my(c,d);
	d = #M;
	c = i2c(i,l,d);
	c = lift(Mod(M*c,l));
	c2i(c,l);
}

Chordi(i1,i2,l,d)= \\ u,v -> -(u+v) in  terms of indices
{
	my(c1,c2,c3);
	c1 = i2c(i1,l,d);
	c2 = i2c(i2,l,d);
	c3 = apply(x->x%l,-(c1+c2));
	c2i(c3,l);
}

ProjPol(Z0,l,d,op,T,pe)=
\\ Given conjugate algebraic numbers indexed by A^n,
\\ gathers them along vector lines in a symmetric way
\\ --> conjugate algebraic numbers indexed by P^(n-1).
\\ Useful to go from a linear representation to its projectivisation
{
	my(Z,PZ,done,j,v,z,PF,PA,t);
	Z=Mod(Mod(Z0,pe),T);
	PZ=List();
	PV=List();
	done=vector(l^d-1);
	for(i=1,l^d-1,
		if(done[i],next);
		v=i2c(i,l,d);
		listput(PV,v);
		z=op;
		for(k=1,l-1,
		 j=c2i(lift(Mod(k*v,l)),l);
		 done[j]=1;
		 z=if(op,z*Z[j],z+Z[j]);
		);
		listput(PZ,z);
	);
	PZ=Vec(PZ);
	PF=factorback(apply(z->'x-z,PZ));
	PA=bestappr(liftpol(PF));
	t=variable(T);
	if(poldegree(PA,t)==0,PA=subst(PA,t,0));
	[PA,PF,liftall(PZ),Vec(PV)];
}

TorsSpaceFrobGen(J,l,B,matFrob)=
\\ Given the basis B of a subspace T of J[l] stable under Frob,
\\ and the matrix of Frob w.r.t. B,
\\ Finds a minimal generating set WB of T as an Fl[Frob]-module
\\ and returns it, along as the coordinates of the generators on B
{
	my(F,Q,NF,WB,cB,n,c);
	[F,Q] = matfrobenius(Mod(matFrob,l),2); \\ F = rational canonical form of matFrob, Q = transition matrix
	Q = lift(Q^(-1));
	NF = apply(poldegree,matfrobenius(F,1)); \\ List of degrees of elementary divisors
	WB = List(); cWB = List();
	n = 1;
	for(i=1,#NF,
  	c=Q[,n];
  	listput(cWB,c);
  	c=centerlift(Mod(c,l));
  	listput(WB,PicLC(J,c,B));
  	n+=NF[i]
	);
	WB = Vec(WB); \\ Generating set of T under Frob
	cWB = Vec(cWB); \\ Coordinates of these generators on B
	[WB,cWB];
}

TorsSpaceFrob(J,gens,cgens,l,matFrob)=
\\ Given a list gens of generators of a Galois-submodule T of J[l],
\\ as well as the matrix of Frob and the coords of these generators (w.r.t. the same basis),
\\ computes a system of representatives of the non-zero orbits of T under Frob
\\ and returns it along as a vector of indices giving the position of these representatives in T
{
	my(d,ld,V,done,ndone,todo,c,i,W,ImodF);
	d = #matFrob;
	ld = l^d;
	V = vector(ld);
	ImodF = List(); \\ Reps of {i}/Frob
  done = vector(ld); \\ 1 = got it, -1 = will get it next time, 0 = haven't got it. Each are stable under Frob.
  c = vector(d);
  i = c2i(c,l);
  V[i] = JgetW0(J);
  done[i] = 1;
  ndone = 1;
  for(n=1,#gens,
		W = gens[n];
    c = cgens[n];
    i = c2i(c,l);
		listput(ImodF,i);
		while(1,
    	V[i] = W;
    	done[i] = 1;
    	ndone += 1;
			i = ActOni(matFrob,i,l);
			if(done[i]==1,break);
			W = PicFrob(J,W)
		)
  );
  while(ndone<ld,
    todo = List();
    for(j=1,ld,
      for(k=j,ld,
        if(done[j]==1 && done[k]==1,
          i = Chordi(j,k,l,d);
          if(done[i] == 0,
            listput(todo,[j,k,i]);
						while(done[i] == 0,
            	done[i] = -1;
							i = ActOni(matFrob,i,l);
						)
          )
        )
      )
    );
    todo = Vec(todo);
		print("Computing ",#todo," new points");
    todo = parapply(t->[PicChord(J,V[t[1]],V[t[2]],0),t[3]],todo);
    for(n=1,#todo,
			W = todo[n][1];
			i = todo[n][2];
			listput(ImodF,i);
			while(1,
      	V[i] = W;
      	done[i] = 1;
      	ndone +=1;
				i = ActOni(matFrob,i,l);
				if(done[i] == 1,break);
				W = PicFrob(J,W)
			)
    )
  );
  [V,Vec(ImodF)];
}

TorsSpaceFrobEval(J,TI,U,l,d,matFrob)=
\\ Given a submodule T of J[l] of Fl-dimension d specified by TI, 
\\ where TI is formed of a list ofrepresentatives of orbits of T under Frob
\\ as output by the previous function,
\\ as well as the matrix of Frob on T
\\ and evaluation data U obtained by PicEvalInit(),
\\ computes the evaluation of every nonzero point of T
{
  my(J=J,T=TI[1],ImodF=TI[2],Z,ZmodF,i,z);
  Z = vector(l^d-1,i,[]);
  ZmodF = parapply(i->PicEval(J,T[i],U),ImodF);
  for(n=1,#ImodF,
    i = ImodF[n];
    z = ZmodF[n];
    while(1,
      Z[i] = z;
      i = ActOni(matFrob,i,l);
      if(Z[i] != [],break);
      z = apply(x->Frob(x,JgetFrobMat(J),JgetT(J),Jgetpe(J)),z);
    )
  );
  Z;
}
