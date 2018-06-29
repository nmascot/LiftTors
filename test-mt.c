#include <pari/pari.h>

GEN bnf3(GEN n)
{
	/*pari_sp av=avma;*/
	printf("Into worker\n");
	GEN pol,K;
	pol = mkpoln(4,gen_1,gen_0,gen_0,n);
	pari_printf("%Ps\n",pol);
	K = bnfinit0(pol,0,NULL,100);
	return K;
	/*return gerepileupto(av,K);*/
}

GEN testmt(ulong t, GEN N0)
{
	pari_sp av=avma;
	GEN vN,done,worker,output;
	long i,workid,pending;
	struct pari_mt pt;
	static entree ep_worker={"_workernico",0,(void*)bnf3,0,"G",""};
	static int entree_added = 0;
	if(entree_added==0) {pari_add_function(&ep_worker); entree_added = 1;}
	vN = cgetg(t+1,t_VEC);
	output = cgetg(t+1,t_VEC);
	for(i=1;i<=t;i++)
	{
		gel(vN,i) = addis(N0,i);
	}
	pari_printf("%Ps\n",vN);
	pending = 0;
	worker = strtofunction("_workernico");
	printf("Appel parallele\n");
	mt_queue_start(&pt, worker);
	for (i=1; i<=t || pending; i++)
  {
		printf("i=%ld\n",i);
    /* submit job (input) */
		printf("sub\n");
    mt_queue_submit(&pt, i, i<=t?mkvec(gel(vN,i)):NULL);
    /* get result (output) */
		printf("get\n");
    done = mt_queue_get(&pt, &workid, &pending);
    if (done)
		{
			printf("stor %lu\n",workid);
			gel(output,workid) = done;
		}
  }
	mt_queue_end(&pt);

	return gerepilecopy(av,output);
}

