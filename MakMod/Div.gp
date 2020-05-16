DivAdd1(A,B,dimres,p,excess,flag)=
{ \\ Mult RR spaces A and B.
  \\ Takes dimres + excess products A[u]*B[v] with random u,v
  \\ If flag, also return the list of pairs [u,v]
  my(nA=#A,nB=#B,m,n,C,uv,r);
  m = #A[,1];
  n = dimres+excess;
  while(1,
		uv = vector(n,j,[1+random(nA),1+random(nB)]);
		\\C = matrix(m,n,i,j,A[i,uv[j][1]]*B[i,uv[j][2]]);
		C = matconcat(parvector(n,j,vector(m,i,A[i,uv[j][1]]*B[i,uv[j][2]])~));
    r = matindexrank(Mod(liftint(C),p))[2];
    if(#r == dimres,
      C = vecextract(C,r);
      return(if(flag,[C,vecextract(uv,r)],C))
    );
    print1("@",#r,"/",dimres," ");
  );
}

BalancedDiv(d,degs)=
{ /* Let degs = [a1,..,an]. Find balanced b1,..,bn such that sum ai*bi = d. */
  my(n=#degs,s=vecsum(degs),D,q);
  q = d\s;
  d -= q*s;
  D = vector(n,i,q);
	while(d,
  	for(i=1,n,
    	if(d>=degs[i],
      	d -= degs[i];
      	D[i] +=1
    	)
		)
	);
  D;
}

DivPerturb(D,degs)=
{
  my(n=#degs,D2,i);
  d = sum(j=1,n,D[j]*degs[j]);
  D2 = BalancedDiv(d-1,degs);
  i=n;
  while(degs[i]==1,
    if(D2[i]+1!=D[i],
      D2[i]+=1;
      return(D2)
    );
    i-=1;
  );
  error("I don't know how to perturb this divisor");
}

Divo2Div(Do,Orbs,tags,n)=
{
  my(D,nO,o,no);
  nO = #Orbs;
  D = vector(n);
  for(i=1,nO,
    o = Orbs[i];
    no = #o;
    for(j=1,no,
      D[GetCoef(tags,o[j])] = Do[i];
    )
  );
  D;
}

MRRsubspace(M4qexps,D,B,T,pe,p,e)=
{
  my(K,i=1,nD=vecsum(D),ncusps=#D);
	if(B==0,B=vector(ncusps));
  K = matrix(nD,#M4qexps[1]);
  for(s=1,ncusps,
    for(j=1,D[s],
      K[i,]=liftall(M4qexps[s][j+B[s],]);
      i++;
    )
  );
  matkerpadic(K,T,pe,p,e);
}
