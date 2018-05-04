install(ZpXQ_inv,GGGL);
Zqinv(a,p,T,e)=Mod(ZpXQ_inv(liftall(a),T,p,e),p^e)*Mod(1,T);

WeiRed(f,h)=
{
 my(F=f+(h/2)^2);
 2^poldegree(F)*subst(F,x,x/2);
}

ordJ(f,p,e)= \\ #Jac(y²=f(x))(F_{p^e})
{
 subst(polresultant(subst(hyperellcharpoly(Mod(f,p)),x,x*y),x^e-1),y,1);
}

RandPt(f,p,T,e)= \\ Random (x,y) on y²=f(x)
{
 my(typ,v=variable(f),x,y,z,w);
 typ=Mod(T,p^e);
 while(1,
  x=Mod(random(typ),T);
  y=subst(f,v,x);
  if(Mod(y,p)==0,next);
  fa=factor('x^2-Mod(y,p));
  if(#fa~==2,
   z=polcoeff(fa[1,1],0);
   z=padicappr('x^2-liftint(y),liftint(z)+O(p^e))[1];
   z=Mod(z,p^e); 
   return([x,if(random(2),-z,z)]);
  )
 );
}

RReval1(P,n,df)= \\ Row of values of x^i and x^i*y at P
{
 my([x,y]=P);
 concat(vector(n+1,i,x^(i-1)),y*vector(n-df/2+1,i,x^(i-1)));
}

RReval(Ps,n,df)= \\ Matrix of values of x^i and x^i*y and the points in Ps 
{
 matconcat(apply(P->RReval1(P,n,df),Ps)~);
}

PicInit(f,p,a,e)= \\ Over Witt(F_{p^a})/p^e
{
 my(t,T,Frob,df,d0,g,nZ,Z,typ,n,Zp,P,Q,lcyc,FrobCyc,W0,V,KV);
 t=varlower("t",variable(f));
 T=liftint(ffinit(p,a,t));
 Frob=FrobPoly(p,T,e);
 df=poldegree(f);
 g=df/2-1; \\ df=2*g+2
 d0=df;
 nZ=6*d0+1;
 Zp=Set();
 n=0;
 Z=List();
 FrobCyc=List();
 \\ Find >=nZ distinct points on the curve, stable under Frob as a set
 while(n<nZ,
  P=RandPt(f,p,T,e);
  Zp=setunion(Zp,Set([Mod(P,p)]));
  if(#Zp>n,
   Q=P;
   lcyc=1;
   while(1,
    listput(Z,Q);
    n+=1;
    Q=Frobapply(Q,T,Frob);
    if(Q==P,break);
    Zp=setunion(Zp,Set([Mod(Q,p)]));
    lcyc+=1
   );
   listput(FrobCyc,lcyc);
  );
 );
 Z=Vec(Z);
 FrobCyc=Vec(FrobCyc);
 nZ=#Z;
 W0=RReval(Z,d0,df);
 V=RReval(Z,3/2*d0,df);
 KV=matkerMod(V~,p,T,e)~;
 [f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob];
}

PicInput(X,P1,P2)= \\ W representing [P1+P2-OO]
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,A,B);
 A=RReval(Z,2*g+3,df);
 B=RReval([P1,P2],2*g+3,df);
 A*matkerMod(B,p,T,e);
}

PicRand(X)=
{
 \\ TODO not generic
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,P1,P2);
 P1=RandPt(f,p,T,e);
 P2=RandPt(f,p,T,e);
 PicInput(X,P1,P2);
}

PicRand2(X)=
{
 \\ TODO still not generic
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,E1,E2,W1,W2);
 E1=vector(d0,i,RandPt(f,p,T,e));
 E2=vector(d0,i,RandPt(f,p,T,e));
 W1=V*matkerMod(RReval(E1,3/2*d0,df),p,T,e);
 W2=V*matkerMod(RReval(E2,3/2*d0,df),p,T,e);
 PicChord(X,W1,W2);
}

PicFrob(X,W)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,W2,i=0,c);
 W2=matrix(#Z,#W);
 for(o=1,#FrobCyc,
  c=FrobCyc[o];
  for(j=1,c-1,
   W2[i+j+1,]=Frobapply(W[i+j,],T,Frob)
  );
  W2[i+1,]=Frobapply(W[i+c,],T,Frob);
  i+=c
 );
 W2;
}
 
DivAdd(WA,WB,d,p,T,e,excess=0)=
{
 my(typ=Mod(T,p^e),WAB,s,t,r,nZ);
 nZ=matsize(WA)[1];
 while(1,
  WAB=matrix(nZ,d+excess);
  for(n=1,d+excess,
   s=RandVec(WA,typ);
   t=RandVec(WB,typ);
   for(k=1,nZ,
    WAB[k,n]=s[k]*t[k]
   )
  );
  r=matrank(Mod(WAB,p));
  if(r==d,
   if(excess,WAB=matimageMod(WAB,p,T,e));
   return(WAB)
  );
  print1("@");
  print1([r,d]);
 );
}

DivSub(WA,WB,KV,d,p,T,e,nIGS=2)=
{
 my(typ=Mod(T,p^e),K,KB,K2,s,W,r);
 KB=matkerMod(WB~,p,T,e)~;
 while(1,
  K=KV;
  for(n=1,nIGS,
   s=RandVec(WA,typ);
   K2=KB;
   for(i=1,#s,K2[,i]*=s[i]);
   K=matconcat([K;K2])
  );
  r=#matker(Mod(K,p));
  if(r==d,return(matkerMod(K,p,T,e)));
  print1("#");
  print1([r,d])
 );
}

PicChord(X,WA,WB,flag=0)= \\ flag = randomize s + 2 * return s
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,typ=Mod(T,p^e),WAWB,WAB,s,sV,WC);
 WAWB=DivAdd(WA,WB,4*d0+1-g,p,T,e); \\ H0(6D0-A-B)
 WAB=DivSub(V,WAWB,KV,d0+1-g,p,T,e); \\ H0(3D0-A-B)
 if((flag>>2)%2,
  for(i=1,#Z,
   if(WAB[i,]==0,print("Line ",i," of WAB is zero"))
  )
 );
 s=if(flag%2,RandVec(WAB,Mod(T,p)),WAB[,1]); \\ [s] = -3D0 + A + B + C 
 sV=V;
 for(i=1,#s,sV[i,]*=s[i]);
 WC=DivSub(WAB,sV,KV,2*d0+1-g,p,T,e); \\ H0(3D0-C)
 if((flag>>1)%2,[WC,s],WC);
}

PicAdd(X,WA,WB)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X);
 PicChord(X,PicChord(X,WA,WB),W0);
}

PicSub(X,WA,WB)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X);
 PicChord(X,PicChord(X,WB,W0),WA);
}

PicMul(X,W,n,flag=0)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,W2);
 if(n==0,return(W0));
 if(n==1,return(W));
 if(n==-1,return(PicChord(X,W,W0)));
 if(n%2,
  W2=PicMul(X,W,(n+1)\2,flag); 
  W2=PicChord(X,W2,W2,flag%2); \\ =-(n+1)*W
  PicChord(X,W2,W,flag%2);
 ,
  W2=PicMul(X,W,-n\2,flag);
  PicChord(X,W2,W2,flag%2)
 );
}

PicFrobPoly(X,W,F)= \\ F(Frob).W
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,d,FW,res);
 d=poldegree(F);
 FW=W;
 res=PicMul(X,W,polcoeff(F,0));
 for(i=1,d,
  FW=PicFrob(X,FW);
  res=PicAdd(X,res,PicMul(X,FW,polcoeff(F,i)))
 );
 res;
}

PicEq(X,WA,WB)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,s,sWB,KsWB,K,K2);
 s=WA[,1];
 sWB=WB;
 for(i=1,#s,sWB[i,]*=s[i]);
 KsWB=matkerMod(sWB~,p,T,e)~;
 K=KV;
 for(j=1,#WA,
  K2=KsWB;
  s=WA[,j];
  for(i=1,#s,K2[,i]*=s[i]);
  K=matconcat([K,K2]~)
 );
 (#matkerMod(K,p,T,e))>0; \\ TODO correct ???
}

PicIsZero(X,W)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X);
 PicEq(X,W,W0);
}

PicRandTors(X,l,C=0)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,N,c,W,lW,v);
 N=ordJ(f,p,poldegree(T));
 v=valuation(N,l);
 if(v==0,error(Str("No rational ",l,"-torsion")));
 N/=l^v;
 if(C,
  c=hyperellcharpoly(Mod(f,p));
  c=subst(c,'x,variable(C));
  if(c%C,error("Incorrect characteristic polynomial"));
  c=c/C;
  if(c%C==0,error("Eigen multiplicity"))
 );   
 while(1,
  W=PicRand(X);
  W=PicMul(X,W,N);
  if(PicIsZero(X,W),
   print("Found an ",l,"-power torsion point but it is 0");
   next
  );
  if(C,
   W=PicFrobPoly(X,W,c);
   if(PicIsZero(X,W),
    print("Found an ",l,"-power torsion point but its projection on the eigenspace is 0");
    next
   )
  );
  lW=PicMul(X,W,l);
  v=1;
  while(PicIsZero(X,lW)==0,
   v+=1;
   W=lW;
   lW=PicMul(X,W,l)
  );
  print("Order was initially ",l,"^",v);
  return(PicChord(X,W,W0,1)) \\ TODO flag for randomization
 );
} 

PicEval(X,W)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,EqU,K,s,sV,U,U2);
 if(g!=2,error("This is ad hoc test code only for genus 2"));
 U=RReval(Z,4,df);
 EqU=matkerMod(U~,p,T,e)~;
 K=matkerMod(EqU*W,p,T,e);
 \\if(#K>1,error("Genericity 1 failed"));
 s=(W*K)[,1];
 sV=V;
 for(i=1,#Z,sV[i,]*=s[i]);
 U=DivSub(W,sV,KV,5,p,T,e);
 U2=matrix(#Z,3);
 for(i=1,#Z,
  U2[i,1]=Mod(1,p);
  U2[i,2]=Z[i][1];
  U2[i,3]=Z[i][2]
 );
 EqU=matkerMod(U2~,p,T,e)~;
 K=matkerMod(EqU*U,p,T,e);
 if(#K!=1,error("Genericity 2 failed"));
 s=K[,1];
 s=matkerMod(matconcat([U2,U*s]),p,T,e)[,1];
 s=liftint(s)*(1+O(p^e));
 [s[1],s[2]]/-s[3];
}
PicEval2(X,W)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,EqU,K,s,sV,U,U2);
 if(g!=2,error("This is ad hoc test code only for genus 2"));
 U=RReval(Z,4,df);
 EqU=matkerMod(U~,p,T,e)~;
 K=matkerMod(EqU*W,p,T,e);
 \\if(#K>1,error("Genericity 1 failed"));
 s=(W*K)[,1];
 sV=V;
 for(i=1,#Z,sV[i,]*=s[i]);
 U=DivSub(W,sV,KV,5,p,T,e);
 U2=matrix(#Z,3);
 for(i=1,#Z,
  U2[i,1]=Mod(1,p);
  U2[i,2]=Z[i][1];
  U2[i,3]=Z[i][1]^2;
 );
 EqU=matkerMod(U2~,p,T,e)~;
 K=matkerMod(EqU*U,p,T,e);
 if(#K!=1,error("Genericity 2 failed"));
 s=K[,1];
 s=matkerMod(matconcat([U2,U*s]),p,T,e)[,1];
 s=liftint(s)*(1+O(p^e));
 [s[1],s[2]]/-s[3];
}

PicChart(X,W)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,n1,n2,s,sV,U);
 n1=2*d0-g;
 n2=d0-g;
 s=W*matkerMod(W[1..n1,],p,T,e);
 if(#s!=1,error("Genericity 1 failed"));
 sV=matrix(#Z,#V);
 for(j=1,#V,
  for(P=n1+1,#Z,
   sV[P,j]=s[P,1]*V[P,j]
  )
 );
 U=DivSub(W,sV,KV,d0+1-g,p,T,e);
 s=U*matkerMod(U[n1+1..n1+n2,],p,T,e);
 if(#s!=1,error("Genericity 2 failed"));
 s[n1+n2+1..#Z,1];
}
