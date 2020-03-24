#include<pari/pari.h>

GEN ZetaFromPointCount(GEN N, ulong p, ulong g)
{
	pari_sp av = avma;
	GEN Z,D,L,Pi;
	ulong i;
	Z = cgetg(g+2,t_SER);
	Z[1] = 0;
  setsigne(Z,1);
  setvarn(Z,0);
  setvalp(Z,1);
  for(i=1;i<=g;i++) gel(Z,i+1) = gdiv(stoi(N[i]),utoi(i));
  Z = gexp(Z,0);
  D = mkpoln(2,gen_m1,gen_1);
	D[1] = 0;
	setsigne(D,1);
  setvarn(D,0);
  Z = gmul(Z,D);
  D = mkpoln(2,gneg(utoi(p)),gen_1);
	D[1] = 0;
	setsigne(D,1);
  setvarn(D,0);
  Z = gmul(Z,D);
  L = cgetg(2*g+3,t_POL);
	L[1] = 0;
  setvarn(L,0);
  setsigne(L,1);
  gel(L,2*g+2) = gen_1;
  for(i=1;i<=g;i++) gel(L,2*g+2-i) = gel(Z,i+2);
  Pi = gen_1;
  for(i=1;i<=g;i++)
  {
    Pi = muliu(Pi,p);
    gel(L,g+2-i) = mulii(gel(L,g+2+i),Pi);
  }
  return gerepilecopy(av,L);
}

GEN PlaneZeta(GEN f, ulong p)
{
	pari_sp av = avma, av1;
	ulong d1,d2,df,g,a,pa,i,ix,m,n;
	GEN P,Pa,N,t,T,vars,Mf,x,fx,f00,L;

	vars = variables_vecsmall(f);
  d1 = poldegree(f,vars[1]);
  d2 = poldegree(f,vars[2]);
  if(d1>d2) df = d1;
  else df = d2;
  g = (df-1)*(df-2);
  g = g/2;
  t = varlower("t",vars[2]);
	Mf = RgXX_to_RgM(f,df+1);
	f00 = cgetg(df+3,t_POL);
	f00[1] = 0;
	setsigne(f00,1);
	setvarn(f00,0);
	for(i=0;i<=df;i++)
	{
		gel(f00,i+2) = gcoeff(Mf,i+1,df+1-i);
	}
	P = utoi(p);

	N = cgetg(g+1,t_VECSMALL); /* Point counts */
	x = cgetg(g+1,t_VEC);
	for(i=1;i<=g;i++) gel(x,i) = gen_0;
	Pa = gen_1;
	for(a=1;a<=g;a++) /* Count C(F_{p^a}) */
	{
		if(a==1) T = t;
		else T=liftint(ffinit(P,a,varn(t)));
		n = 0;
		Pa = mulii(Pa,P);
		pa = itou(Pa); /* pa = p^a = q, detects overflows */
		av1 = avma;
		for(ix=0;ix<pa;ix++) /* Loop over F_q */
		{
			avma = av1;
			m = ix;
			for(i=1;i<=a;i++)
			{
				gel(x,i) = utoi(m%p);
				m = m/p;
			}
			fx = poleval(f,gtopolyrev(x,varn(T)));
			n += lg(polrootsmod(fx,mkvec2(T,P)))-1;
		}
		if(itou(gel(f00,df+2))%p == 0) n++;
		n += lg(polrootsmod(f00,mkvec2(T,P)))-1;
		N[a] = n;
	}
	L = ZetaFromPointCount(N,p,g);
	return gerepileupto(av,L);
}

GEN SuperZeta(GEN f, ulong m, ulong p) /* y^m = f(x), assumes f sqfree and gcd(deg(f),m)=1 */
{
	pari_sp av = avma, av1;
	GEN N,P,Pa,x,T,t,fx,L;
	ulong d,g,a,pa,ix,i,n,mx,e1,e2;

	d = degree(f);
	g = (m-1)*(d-1);
	g = g/2;

	P = utoi(p);
	t = varlower("t",gvar(f));
  N = cgetg(g+1,t_VECSMALL); /* Point counts */
  x = cgetg(g+1,t_VEC); /* Value of x in Fq */
  for(i=1;i<=g;i++) gel(x,i) = gen_0;
  Pa = gen_1;

  for(a=1;a<=g;a++) /* Count C(F_{p^a}) */
  {
		Pa = mulii(Pa,P);
    pa = itou(Pa); /* pa = p^a = q, detects overflows */
		e1 = ugcd(pa-1,m);
		if(e1==1)
		{ /* y -> y^m is bijective on Fq so point counting is trivial */
			N[a] = pa+1;
			continue;
		}
		e2 = (pa-1)/e1; /* x in Fq* is an m-th power iff. x^e2==1, in which case it has e1 m-th roots */
    if(a==1) T = t;
    else T=liftint(ffinit(P,a,varn(t)));
    n = 1; /* Point at oo */
    av1 = avma;
    for(ix=0;ix<pa;ix++) /* Loop over Fq */
    {
      avma = av1;
      mx = ix;
      for(i=1;i<=a;i++)
      {
        gel(x,i) = utoi(mx%p);
        mx = mx/p;
      }
      fx = poleval(f,gmodulo(gmodulo(gtopolyrev(x,varn(T)),P),T));
			fx = liftall(fx);
			if(gequal0(fx))
			{
				n++;
			}
			else
			{
				fx = Fq_powu(fx,e2,T,P);
				if(gequal1(fx))
					n += e1;
			}
    }
    N[a] = n;
  }
	L = ZetaFromPointCount(N,p,g);
	return gerepileupto(av,L);
}
