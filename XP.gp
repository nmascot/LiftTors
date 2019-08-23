read("install.gp");
read("Hyper2RR.gp");
read("GalRep.gp");
f = x^5+x^4; h = x^3+x+1; \\ X1(13)
l = 7; \\ The representation occurs in the 7-torsion of the Jacobian
p = 17; \\ We choose to get it 17-adically
e = 256; \\ Target p-adic accuracy is O(17^32)
chi = x^2-x-2; \\ Char.poly. of the Frobenius at p (another possible choice is x^2-2*x-1)
Lp = hyperellcharpoly(Mod([f,h],p));
P1 = [-1,1]; P2 = [0,0];

  C = Hyper2RR([f,h],P1,P2);
  C=concat(C,['y]);

[f,g,d0,L,LL,L1,L2,Bad]=C;
 L = RR_rescale(L,p);
  LL = RR_rescale(LL,p);
  L1 = RR_rescale(L1,p);
  L2 = RR_rescale(L2,p);
  Bad *= lcm(apply(S->lcm(apply(f->denominator(content(f)),S)),[L,L1,L2]));
{  if(chi,
    print("T = part of J[",l,"] where Frob_",p," acts by ",chi);
    d = poldegree(chi); \\ Dimension of representation
    a = mordroot(chi,l) \\ q = p^a
  ,
    print("T = all of J[",l,"]");
    d=2*g;
    a = mordroot(Lp,l)
  );}
  J=PicInit(f,g,d0,L,LL,Bad,p,a,e);

NJ = polresultant(Lp,x^a-1);
M = NJ/l^4;
e1=1;
J1=PicRed(J,e1);
{
while(1,
W1=PicRand(J1);
W1=PicMul(J1,W1,M,0);
if(PicIsZero(J1,W1)==0,break);
);
}
W2=PicLiftTors(J,W1,1,l);