RandVec(A,typ)= \\ Random vector in the col span of A
 sum(i=1,#A,random(typ)*A[,i]);

\\ Ca ca marche (mais pourquoi?)
Hsort(A,p)=
{
 my(B=List());
 \\print("Hsort in ",#A);
 for(j=1,#A,
  if(Mod(A[,j],p),listput(B,A[,j]))
 );
 \\print("Hsort out ",#B);
 matconcat(Vec(B));
}
matkerMod(A,p,T,n)=if(n==1,Mod(liftint(matker(A)),p),Mod(Hsort(matkerzq(liftall(A),T,p,n),p),p^n)*Mod(1,T));
matimageMod(A,p,T,n)=if(n==1,Mod(liftint(matimage(A)),p),Mod(Hsort(matimagezq(liftall(A),T,p,n),p),T)*Mod(1,p^n));
matdetMod(A,p,T,n)=if(n==1,Mod(liftint(matdet(A)),p),Mod(matdetzq(liftall(A),T,p,n),T)*Mod(1,p^n));
matinvMod(A,p,T,n)=if(n==1,Mod(liftint(A^(-1)),p),Mod(matinvzq(liftall(A),T,p,n),T)*Mod(1,p^n));

matF(A,p,T,n)=
{ \\ transpose of L matrix of transpose of A, assuming the rows of A are indep mod p
 my(B);
 B=matsupplement(Mod(A,p));
 B=matconcat([A,B[,#A+1..#A~]]);
 B=matinvzq(liftall(B),T,p,n);
 Mod(B[1..#A,],T)*Mod(1,p^n);
}

inline(mat2col(A)=
{
 my(m=#A~,n=#A,v);
 v=vector(m*n)~;
 for(i=1,m,
  for(j=1,n,
   v[(i-1)*n+j]=A[i,j]
  )
 );
 v;
});

col2mat(v,m,n)=
{
 my(A=matrix(m,n));
 for(i=1,m,
  for(j=1,n,
   A[i,j]=v[(i-1)*n+j]
  )
 );
 A;
}

VecsmallCompl(v,n)=Vecsmall(setminus(Set([1..n]),Set(v)));

FindMinorCompl(A)=
{
 my([m,n]=matsize(A),u,U,v,V);
 [u,v]=matindexrank(A);
 U=VecsmallCompl(u,m);
 V=VecsmallCompl(v,n);
 [[u,U],[v,V]];
}

inline(M2ABCD(M,uv)=
{
 my([uU,vV]=uv,[u,U]=uU,[v,V]=vV);
 [
  matrix(#u,#v,i,j,M[u[i],v[j]]),
  matrix(#u,#V,i,j,M[u[i],V[j]]),
  matrix(#U,#v,i,j,M[U[i],v[j]]),
  matrix(#U,#V,i,j,M[U[i],V[j]])
 ];
});
