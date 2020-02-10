DivAdd1(A,B,dimres,p,excess,flag)=
{ \\ Mult RR spaces A and B.
  \\ Takes dimres + excess products A[u]*B[v] with random u,v
  \\ If flag, also return the list of pairs [u,v]
  my(nA=#A,nB=#B,m,n,C,u,v,uv,a,b,r);
  m = #A[,1];
  n = dimres+excess;
  C = matrix(m,n);
  uv = vector(n);
  while(1,
    for(j=1,n,
      u = 1+random(nA);
      v = 1+random(nB);
      a = A[,u];
      b = B[,v];
      uv[j] = [u,v];
      for(i=1,m,
        C[i,j] = a[i]*b[i]
      )
    );
    r = matindexrank(Mod(C,p))[2];
    if(#r == dimres,
      C = vecextract(C,r);
      uv = vecextract(uv,r);
      return(if(flag,[C,uv],C))
    );
    print1("@",#r,"/",dimres," ");
  );
}

BalancedDiv(d,degs)=
{
  my(n=#degs,s=vecsum(degs),D,q);
  q = d\s;
  d -= q*s;
  D = vector(n,i,q);
  for(i=1,n,
    if(d>=degs[i],
      d -= degs[i];
      D[i] +=1
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

MRRsubspace(M4qexps,D,T,p,e)=
{
  my(K,i=1,nD=vecsum(D),ncusps=#D);
  K = matrix(nD,#M4qexps[1]);
  for(s=1,ncusps,
    for(j=1,D[s],
      K[i,]=liftall(M4qexps[s][j,]);
      i++;
    )
  );
  matkerpadic_safe(K,T,p,e);
}
