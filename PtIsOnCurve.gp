PtIsOnCurve(f,P)=
{
	my(vars=variables(f),u);
	u=subst(subst(f,vars[1],P[1]),vars[2],P[2]);
	if(#vars==3,
		u=subst(u,vars[3],if(#P==3,P[3],1))
	);
	u==0;
}

PtIsOnHyperCurve(F,P)=
{
	my(f,h=0);
	if(type(F)=="t_VEC",
		[f,h]=F
	,
		f=F
	);
	PtIsOnCurve('y^2+h*'y-f,P);
}

PtIsOnSuperCurve(f,m,P)=PtIsOnCurve('y^m-f,P);
	
