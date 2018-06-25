c2i(c,l)=
{
	my(i);
	i=subst(Pol(c),'x,l);
	if(i,i,l^#c);
}

i2c(i,l,d)=
{
	my(j=i,c=vector(d));
	for(n=1,d,
		c[d+1-n] = j%l;
		j = j\l
	);
  c;
}

TorsSpace(J,B,l)=
{
	my(d=#B,ld,V,done,ndone=0,todo,c,i);
	ld = l^d;
	V = vector(ld);
	done = vector(ld);
	c = vector(d);
	i = c2i(c,l);
	V[i] = JgetW0(J);
	done[i] = 1;
	ndone = 1;
	for(n=1,#B,
		c = vector(d);
		c[n] = 1;
		i = c2i(c,l);
		V[i] = B[n];
		done[i] = 1;
		ndone += 1
	);
	while(ndone<ld,
		todo = List();
		for(j=1,ld,
			for(k=j,ld,
				if(done[j]==1 && done[k]==1,
					cj = i2c(j,l,d);
					ck = i2c(k,l,d);
					c = apply(x->x%l,-(cj+ck));
					i = c2i(c,l);
					if(done[i] == 0,
						listput(todo,[j,k,i]);
						done[i] = -1;
					)
				)
			)
		);
		todo = Vec(todo);
		print("Computing ",#todo," new points");
		todo = parapply(t->[PicChord(J,V[t[1]],V[t[2]],0),t[3]],todo);
		for(n=1,#todo,
			V[todo[n][2]] = todo[n][1];
			done[todo[n][2]] = 1;
			ndone +=1;
		)
	);
	V;
}

