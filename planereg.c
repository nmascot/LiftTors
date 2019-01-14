#include<pari/pari.h>

GEN PlaneZeta(GEN f, ulong p)
{
	pari_sp av = avma, av1;
	ulong d1,d2,df,g,a,pa,i,ix,m,n;
	GEN P,N,t,T,vars,Mf,x,fx,f00,Z,D,L,Pi;

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
	setsigne(f00,1);
	setvarn(f00,0);
	for(i=0;i<=df;i++)
	{
		gel(f00,i+2) = gcoeff(Mf,i+1,df+1-i);
	}
	P = utoi(p);

	N = cgetg(g+1,t_VECSMALL);
	x = cgetg(g+2,t_POL);
	setvarn(x,varn(t));
	for(i=1;i<=g;i++) gel(x,i+1) = gen_0;
	pa = 1;
	for(a=1;a<=g;a++)
	{
		if(a==1) T = t;
		else T=liftint(ffinit(P,a,varn(t)));
		n = 0;
		pa *= p;
		av1 = avma;
		for(ix=0;ix<pa;ix++)
		{
			avma = av1;
			m = ix;
			for(i=1;i<=a;i++)
			{
				gel(x,i+1) = utoi(m%p);
				m = m/p;
			}
			fx = poleval(f,x);
			n += lg(polrootsmod(fx,mkvec2(T,P)))-1;
		}
		if(itou(gel(f00,df+2))%p == 0) n++;
		n += lg(polrootsmod(f00,mkvec2(T,P)))-1;
		N[a] = n;
	}
	Z = cgetg(g+2,t_SER);
	setsigne(Z,1);
	setvarn(Z,0);
	setvalp(Z,1);
	for(i=1;i<=g;i++) gel(Z,i+1) = gdiv(stoi(N[i]),utoi(i));
	Z = gexp(Z,0);
	D = mkpoln(2,gen_m1,gen_1);
	setvarn(D,0);
	Z = gmul(Z,D);
	D = mkpoln(2,gneg(P),gen_1);
  setvarn(D,0);
  Z = gmul(Z,D);
	L = cgetg(2*g+3,t_POL);
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
