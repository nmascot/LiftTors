TorsPlane(X,W1,W2,l)=
{
 my(A=vector(l-1),B=vector(l-1),V=List());
 A[1]=W1;
 B[1]=W2;
 for(i=2,l-1,
  A[i]=PicAdd(X,W1,A[i-1]);
  B[i]=PicAdd(X,W2,B[i-1])
 );
 for(i=0,l-1,
  for(j=0,l-1,
   if(i==0 && j==0,next);
   if(j==0,listput(V,A[i]);next);
   if(i==0,listput(V,B[j]);next);
   listput(V,PicAdd(X,A[i],B[j]))
  )
 );
 Vec(V);
}

Z2pol(Z)=factorback(apply(u->'x-u,Z));

mordroot(f,p)=
{
 my(x=variable(f),N,fa,l,v);
 if(issquarefree(Mod(f,p))==0,error("Not squarefree!"));
 N=lcm(factormod(f,p,1)[,1]);
 N=p^N-1;
 fa=factor(N);
 for(i=1,#fa~,
  [l,v]=fa[i,];
  while(v,
   if(Mod(x^(N/l)-1,p)%Mod(f,p),break);
   N/=l;
   v-=1
  )
 );
 N;
}
