#include <R.h>
#include <Rdefines.h>

     void mixloglike(double *jdens, double *pmix, double *wt, int *nobs, int *ng,double *mll)
     {
       int i, j;
       double pmll;
       mll[0] = 0.0;
     
       for(j = 0; j < *nobs; j++)
       {
	       pmll=0.0;
	       for(i = 0; i < *ng; i++)
	       {
		       pmll += jdens[i * *nobs + j] * pmix[i * *nobs + j];
	       }
	       mll[0] += -wt[j]*log(pmll);
       }
     }

void jointlike(double *mdens,int *nobs, int *nvar, double *jll)
{
	int i,j;
        for(j = 0; j < *nobs; j++)
	{ 
		jll[j] = 1.0;
		for(i = 0; i < *nvar; i++)
			jll[j] = jll[j]*mdens[i* *nobs + j];
	}
}

void mixpost(double *jdens, double *pmix, double *wt, int *nobs, int *ng,double *post)
{
	int i, j;
	double pmll;
     
	for(j = 0; j < *nobs; j++)
	{
		pmll=0.0;
		for(i = 0; i < *ng; i++)
		{
			post[i * *nobs + j] = jdens[i * *nobs + j] * pmix[i * *nobs + j];
			pmll += post[i * *nobs + j];
		}
		for(i = 0; i < *ng; i++)
			post[i * *nobs + j] = post[i * *nobs + j] / pmll;
	}
}

