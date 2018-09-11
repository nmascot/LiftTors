i=Mod('w,'w^2+1);
bar(z)=trace(z+0*i)-z;

A=[[1+2*i,3],[-1-4*i,5],[1+4*i,7],[-7-10*i,11],[-1+4*i,13],[7+0*i,17],[1-14*i,19],[17-4*i,23],[-9-12*i,29],[1+0*i,31],[-25+28*i,37],[-5+0*i,41],[17+40*i,47],[23-20*i,53],[-39+22*i,59],[63+20*i,61],[65-22*i,67]];

char1(X)=my([a,p]=X,c);c=x^3-a*x^2+p*bar(a)*x-p^3;if(p%8>4,-subst(c,x,-x),c);
char2(X)=my([a,p]=X);subst(liftpol(char1([a,p])*char1([bar(a),p])),'w,0);

\\ Data for l=3
L_3_5 = polrecip(78125*x^14 + 31250*x^13 + 9375*x^12 - 2500*x^11 - 375*x^10 + 350*x^9 + 675*x^8 + 200*x^7 + 135*x^6 + 14*x^5 - 3*x^4 - 4*x^3 + 3*x^2 + 2*x + 1);
L_3_7 = polrecip(823543*x^14 + 470596*x^13 + 16807*x^12 - 57624*x^11 - 31213*x^10 - 9604*x^9 + 2163*x^8 + 2352*x^7 + 309*x^6 - 196*x^5 - 91*x^4 - 24*x^3 + x^2 + 4*x + 1);
L_3_11 = polrecip(19487171*x^14 - 7086244*x^13 + 3382071*x^12 - 585640*x^11 + 347391*x^10 - 73084*x^9 + 49731*x^8 - 11440*x^7 + 4521*x^6 - 604*x^5 + 261*x^4 - 40*x^3 + 21*x^2 - 4*x + 1);
