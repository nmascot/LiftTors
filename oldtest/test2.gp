read("LinAlg.gp");
read("Frob.gp");
read("MakHyper.gp");
f=x^6-x+1;
X=PicInit(f,1009,7,1);
[f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X;
\\V2=DivAdd(V,V,6*d0+1-g,p,T,1);
E=vector(d0,i,RandPt(f,p,T,1));
W=V*matker(RReval(E,3/2*d0,df));
\\WV=DivAdd(W,V,5*d0+1-g,p,T,1);
\\F=vector(#Z,i,Z[i][1]); \\ Function

/*vecxpress(A,v)=my(k);k=matker(matconcat([A,v]));-k[1..#A,1]/k[#A+1,1];

A=matconcat(vector(#W,i,vecxpress(V,W[,i])));
A2=matsupplement(A);
V=V*A2;
A=matconcat(vector(#WV,i,vecxpress(V2,WV[,i])));
A2=matsupplement(A);
V2=V2*A2;

M=matrix(6,6);
{
for(i=1,6,
 s=V[,i+#W];
 s=vecxpress(V2,s);
 M[,i]=s[1+#WV..#V2];
);
}

MF=matrix(6,6);
{
for(i=1,6,
 s=V[,i+#W];
 for(P=1,#Z,s[P]*=F[P]);
 s=vecxpress(V2,s);
 MF[,i]=s[1+#WV..#V2];
);
}

res1=matdet(MF)/matdet(M);
res2=prod(i=1,#E,E[i][1]);*/
