c2i(c,l)=
{
	my(i);
	i=subst(Pol(c),'x,l);
	if(i,i,l^#c);
}

i2c(i,l,d)=
{
	my(j=i,c=vector(d)~);
	for(n=1,d,
		c[d+1-n] = j%l;
		j = j\l
	);
  c;
}

ActOni(M,i,l)=
{
	my(c,d);
	d = #M;
	c = i2c(i,l,d);
	c = lift(Mod(M*c,l));
	c2i(c,l);
}

Chordi(i1,i2,l,d)=
{
	my(c1,c2,c3);
	c1 = i2c(i1,l,d);
	c2 = i2c(i2,l,d);
	c3 = apply(x->x%l,-(c1+c2));
	c2i(c3,l);
}

A2P1(A,l,op,T,pe)=
{
	my(P,vecop);
	\\vecop(v)=liftall(if(op,vecprod,vecsum)(Mod(v,T)*Mod(1,pe)));
	P = vector(l+1);
	for(s=0,l-1,
		P[s+1] = vecprod(vector(l-1,i,A[i+l*((s*i)%l)]))
	);
	P[l+1] = vecprod(vector(l-1,i,A[l*i]));
	P;
}

TorsSpace(J,B,l)=
{
	my(d=#B,ld,V,done,ndone=0,todo,c,i);
	ld = l^d;
	V = vector(ld);
	done = vector(ld);
	c = vector(d);
	i = c2i(c,l);
	V[i] = JgetW0(J);
	done[i] = 1;
	ndone = 1;
	for(n=1,#B,
		c = vector(d);
		c[n] = 1;
		i = c2i(c,l);
		V[i] = B[n];
		done[i] = 1;
		ndone += 1
	);
	while(ndone<ld,
		todo = List();
		for(j=1,ld,
			for(k=j,ld,
				if(done[j]==1 && done[k]==1,
					i = Chordi(j,k,l,d);
					if(done[i] == 0,
						listput(todo,[j,k,i]);
						done[i] = -1;
					)
				)
			)
		);
		todo = Vec(todo);
		print("Computing ",#todo," new points");
		todo = parapply(t->[PicChord(J,V[t[1]],V[t[2]],0),t[3]],todo);
		for(n=1,#todo,
			V[todo[n][2]] = todo[n][1];
			done[todo[n][2]] = 1;
			ndone +=1;
		)
	);
	V;
}

TorsSpaceFrob(J,gens,cgens,l,matFrob)=
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

