makfn(X,v)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,K,V2);
 V2=RReval(Z,3*d0,df);
 K=matkerMod(matconcat([V2,v]),p,T,e)[,1];
 \\print(-centerlift(K));
 -(sum(i=0,18,K[i+1]*('x^i))+'y*sum(i=0,15,K[i+20]*('x^i)))/K[36];
}

makfneval(X,v,P)=subst(subst(makfn(X,v),'x,P[1]),'y,P[2]);

makfnratzeros(X,v)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,K,V2,a,b,zx);
 V2=RReval(Z,3*d0,df);
 K=matkerMod(matconcat([V2,v]),p,T,e)[,1];
 \\print(-centerlift(K));
 a=-sum(i=0,18,K[i+1]*('x^i))/K[36];
 b=-sum(i=0,15,K[i+20]*('x^i))/K[36];
 zx=polrootsmod(a^2-b^2*f,[p,T]);
 vector(#zx,i,[zx[i],-subst(a,'x,zx[i])/subst(b,'x,zx[i])]);
}

makWratzeros(X,W)=
{
 my(Z,Z1);
 Z=Set(makfnratzeros(X,W[,1]));
 for(j=2,#W,Z=setintersect(Z,Set(makfnratzeros(X,W[,j]))));
 Vec(Z);
}


/*matechelon(A,p)=
{ \\ Row echelon form, assuming it has the same shape mod p
 my(M=A,m=#A~,n=#A,r=0,C=vector(m),D=vector(n),d);
 for(k=1,n,
  for(j=1,m,
   if(M[j,k]!=0 && C[j]==0,
    d=-1/M[j,k];
    M[j,k]=-1;
    for(s=k+1,n,M[j,s]*=d);
    for(i=1,m,
     if(i!=j,
      d=M[i,k];
      for(s=k+1,n,M[i,s]+=d*M[j,s]);
      C[j]=k;
      D[k]=j;
     )
    );
    next(2);
   )
  );
  r+=1;
  D[k]=0;
 );
 M;
}*/

matechelon(A,p,pn)=
{
 my(M=A,m=#A~,n=#A,piv=vector(m),nonpiv=List());
 i=1;
 for(j=1,n,
  \\print("j=",j," i=",i);
  \\ Look for nonzero in col j under i
  for(k=i,m,
   if(Mod(M[k,j],p),
    \\print("Row ",k," nonzero");
    if(k>i, \\ Need to swap rows i and k ?
     \\print("Swap ",i," ",k);
     for(s=j,n,
      t=M[i,s];
      M[i,s]=M[k,s];
      M[k,s]=t
     );
     \\printp(M)
    );
    \\ Now M[i,j] nonzero
    piv[i]=j; \\ The pivot for line i is at col j
    \\ Scale line by 1/t
   \\ print("Scaling row");
    t=1/M[i,j];
    M[i,j]=Mod(1,pn);
    for(s=j+1,n,M[i,s]*=t);
    \\printp(M);
    \\ Use pivot to clear other rows
    \\print("Clearing other rows");
    for(l=1,m,
     if(l!=i && M[l,j]!=0, \\ Need work?
      \\print(l);
      t=M[l,j];
      M[l,j]=Mod(0,pn);
      for(s=j+1,n,M[l,s]-=(t*M[i,s]))
     )
    );
    \\printp(M);
    i+=1;
    next(2);
   )
  );
  \\ the j column is 0
  if(M[i..m,j]!=0,
   print("Matrix has bad reduction!");
  ,
   \\print("Whole col ",j," is 0");
   listput(nonpiv,j);
  )
 );
 [M,piv,Vec(nonpiv)];
}

matechelon(A,p,pn)=
{
 my(M=A,m=#A~,n=#A,piv=vector(m),nonpiv=List());
 i=1;
 for(j=1,n,
  \\print("j=",j," i=",i);
  \\ Look for nonzero in col j under i
  for(k=i,m,
   if(Mod(M[k,j],p),
    \\print("Row ",k," nonzero");
    if(k>i, \\ Need to swap rows i and k ?
     \\print("Swap ",i," ",k);
     for(s=j,n,
      t=M[i,s];
      M[i,s]=M[k,s];
      M[k,s]=t
     );
     \\printp(M)
    );
    \\ Now M[i,j] nonzero
    piv[i]=j; \\ The pivot for line i is at col j
    \\ Scale line by 1/t
   \\ print("Scaling row");
    t=1/M[i,j];
    M[i,j]=Mod(1,pn);
    for(s=j+1,n,M[i,s]*=t);
    \\printp(M);
    \\ Use pivot to clear other rows
    \\print("Clearing other rows");
    for(l=1,m,
     if(l!=i && M[l,j]!=0, \\ Need work?
      \\print(l);
      t=M[l,j];
      M[l,j]=Mod(0,pn);
      for(s=j+1,n,M[l,s]-=(t*M[i,s]))
     )
    );
    \\printp(M);
    i+=1;
    next(2);
   )
  );
  \\ the j column is 0
  if(M[i..m,j]!=0,
   print("Matrix has bad reduction!");
  ,
   \\print("Whole col ",j," is 0");
   listput(nonpiv,j);
  )
 );
 [M,piv,Vec(nonpiv)];
}


/*nchooser(n,r)=
{
 my(v,a,j);
 v=List(vector(n,i,i));
 a=vector(r);
 for(i=1,r,
  j=1+random(n-r+1);
  a[i]=v[j];
  listpop(v,j)
 );
 [a,Vec(v)];
}*/

PropZeros(A)=
{
 my([m,n]=matsize(A),z=0);
 for(i=1,m,
  for(j=1,n,
   if(A[i,j]==0,z+=1)
  )
 );
 z/(m*n)+0.;
}

/*FindMinor(A,r)=
{ \\ TODO use echelon instead
 my([m,n]=matsize(A),B=A,C,k,U=V=List(),u=List([1..m]),v=List([1..n]));
 while(m>r,
  k=1+random(m);
  C=matrix(m-1,n);
  for(j=1,n,
   for(i=1,m-1,
    C[i,j]=B[if(i<k,i,i+1),j]
   )
  );
  if(matrank(C)>=r,
   B=C;
   listput(U,u[k]);
   listpop(u,k);
   m-=1;
   \\print(m," ",r)
  )
 );
 while(n>r,
  k=1+random(n);
  C=matrix(m,n-1);
  for(j=1,n-1,
   for(i=1,m,
    C[i,j]=B[i,if(j<k,j,j+1)]
   )
  );
  if(matrank(C)>=r,
   B=C;
   listput(V,v[k]);
   listpop(v,k);
   n-=1;
   \\print(n," ",r)
  )
 );
 [[Vecsmall(u),Vecsmall(U)],[Vecsmall(v),Vecsmall(V)]];
}*/

vecxpress(A,v)=my(k);k=matker(matconcat([A,v]));-k[1..#A,1]/k[#A+1,1];
