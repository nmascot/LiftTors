CycleDecomp(s)=
{
	my(n=#s,done=vector(n),j,cyc=List(),c);
	for(i=1,n,
		if(done[i]==0,
			j=i;
			c=List();
			until(j==i,
				listput(c,j);
				done[j]=1;
				j=s[j]
			);
			listput(cyc,Vecsmall(c));
		)
	);
	Vec(cyc);
}




