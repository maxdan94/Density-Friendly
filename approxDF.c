/*
Maximilien Danisch
Mai 2016
http://bit.ly/maxdan94
maximilien.danisch@telecom-paristech.fr

Info:
Feel free to use these lines as you wish. This program computes an approximation of the density-friendly decomposition.

To compile:
gcc approxDF.c -fopenmp -o approxDF -O9

To execute:
./approxDF nthreads iter net.txt rates.txt pavafit.txt cuts.txt

- nthreads is the number of threads to use
- iter is the number of iterations over all edges to perform
- net.txt should contain the graph (one edge on each line: 2 unsigned separated by a space)
- rates.txt will contain the density value for each node
- pavafit.txt will contain the profile given by the PAVA fit, that is the week approximation of the density-friendly ("size density density-upperbound" on each line)
- cuts.txt will contain the profile given by correct cuts, that is the strong approximation of the density-friendly ("size density density-upperbound" on each line)
Some information will be printed in the terminal.

*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define NLINKS 500000000//25690705119//maximum number of edges for memory allocation, will increase if needed

typedef struct {
	unsigned s;
	unsigned t;
  double a;//alpha value of edge (s,t)
} edge;

typedef struct {
	unsigned n;
	double r;//rate value of node n
} node;

typedef struct {
	unsigned n;//number of nodes
	unsigned long long e;//number of edges
	unsigned *map;//correspondance between old and new nodeID
	edge *edges;//list of all edges
	node *nodes;//value associated to each node
	node *nodes2;//value associated to each node
  double *ne;//ne[i]=number of edges from i to nodes before (used for pava)
	unsigned *cd;//cumulative degree
	unsigned *cuts;
  unsigned iter;//number of iterations
} optim;

//compute the maximum of three unsigned
inline unsigned max3(unsigned a,unsigned b,unsigned c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

optim* readedgelist(char* edgelist){
	unsigned long long e1=NLINKS;
	optim *opt=(optim*)malloc(sizeof(optim));
	FILE *file;
	opt->n=0;
	opt->e=0;
	file=fopen(edgelist,"r");
	opt->edges=(edge*)malloc(e1*sizeof(edge));
	while (fscanf(file,"%u %u", &(opt->edges[opt->e].s), &(opt->edges[opt->e].t))==2) {
		opt->n=max3(opt->n,opt->edges[opt->e].s,opt->edges[opt->e].t);
		if (opt->e++==e1) {
			e1+=NLINKS;
			opt->edges=(edge*)realloc(opt->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	opt->n++;
	opt->edges=(edge*)realloc(opt->edges,opt->e*sizeof(edge));
	return opt;
}

void relabel(optim *opt) {
	unsigned long long i,j;
	unsigned *newlabel;

	newlabel=(unsigned*)malloc(opt->n*sizeof(unsigned));
	for (i=0;i<opt->n;i++) {
		newlabel[i]=opt->n;
	}
	opt->map=(unsigned*)malloc(opt->n*sizeof(unsigned));
	j=0;
	for (i=0;i<opt->e;i++) {
		if (newlabel[opt->edges[i].s]==opt->n){
			newlabel[opt->edges[i].s]=j;
			opt->map[j++]=opt->edges[i].s;
		}
		if (newlabel[opt->edges[i].t]==opt->n){
			newlabel[opt->edges[i].t]=j;
			opt->map[j++]=opt->edges[i].t;
		}
		opt->edges[i].s=newlabel[opt->edges[i].s];
		opt->edges[i].t=newlabel[opt->edges[i].t];
	}
	opt->n=j;
	free(newlabel);
	opt->map=(unsigned*)realloc(opt->map,opt->n*sizeof(unsigned));
}

//Step 1: Frank-Wolf gradiant descent

//initialize the optim datastructure
void init(optim *opt){
	unsigned long long k;
  opt->iter=0;
	opt->nodes=(node*)malloc(opt->n*sizeof(node));
  for (k=0;k<opt->n;k++){
    opt->nodes[k].n=k;
    opt->nodes[k].r=0;
  }
  #pragma omp parallel for private(k)
	for (k=0;k<opt->e;k++){
    opt->edges[k].a=.5;
    #pragma omp atomic update
		opt->nodes[opt->edges[k].s].r+=.5;
		#pragma omp atomic update
		opt->nodes[opt->edges[k].t].r+=.5;
	}
}

//one pass over all edges
void onepass(optim *opt){
	unsigned i,j;
	unsigned long long k;
	double gamma;
	opt->iter++;
	gamma=2./(2.+opt->iter);
	#pragma omp parallel for private(i,j,k)
	for (k=0;k<opt->e;k++){//parfor
		i=opt->edges[k].s;
		j=opt->edges[k].t;
		if (opt->nodes[i].r<opt->nodes[j].r){
			#pragma omp atomic update
			opt->nodes[i].r+=gamma*(1.-opt->edges[k].a);//careful
			#pragma omp atomic update
			opt->nodes[j].r-=gamma*(1.-opt->edges[k].a);//careful
			opt->edges[k].a=(1.-gamma)*opt->edges[k].a+gamma;
		}
		else if (opt->nodes[i].r>opt->nodes[j].r){
			#pragma omp atomic update
			opt->nodes[i].r-=gamma*opt->edges[k].a;//careful
			#pragma omp atomic update
			opt->nodes[j].r+=gamma*opt->edges[k].a;//careful
			opt->edges[k].a*=(1.-gamma);
		}
	}
}

//to print in file the value for each stub: NOT USED
void print_alphas(optim* opt, char* alphas){
	FILE *file=fopen(alphas,"w");
	unsigned long long k;
	for (k=0;k<opt->e;k++){
		fprintf(file,"%u %u %e\n",opt->map[opt->edges[k].s],opt->map[opt->edges[k].t],opt->edges[k].a);
	}
	fclose(file);
}

void freeoptim(optim *opt){
	free(opt->map);
	free(opt->edges);
	free(opt->nodes);
  free(opt->ne);
	free(opt->cuts);
	free(opt);
}

//Step 2: Isotonic regression with PAVA

//used for quicksort 
static int compare_nodes(void const *a, void const *b){
	node const *pa = a;
	node const *pb = b;
	if ((*pa).r<=(*pb).r)
		return 1;
	return -1;
}
//used for quicksort 
static int compare_edges(void const *a, void const *b){
	edge const *pa = a;
	edge const *pb = b;
	if ((*pa).s>(*pb).s)
		return 1;
	if ((*pa).s<(*pb).s)
		return -1;
	if ((*pa).t<=(*pb).t)
		return 1;
	return -1;
}

//prepare the array to fit with PAVA
void prepava(optim *opt){
  unsigned u,v;
  unsigned long long i;
  unsigned *newlabel=(unsigned*)malloc(opt->n*sizeof(unsigned));
  qsort(opt->nodes,opt->n,sizeof(node),compare_nodes);

  for (i=0;i<opt->n;i++){
    newlabel[opt->nodes[i].n]=i;
  }
  for (i=0;i<opt->e;i++){
		u=newlabel[opt->edges[i].s];
		v=newlabel[opt->edges[i].t];
		if (u<v){
			opt->edges[i].s=u;
			opt->edges[i].t=v;
		}
		else {
			opt->edges[i].s=v;
			opt->edges[i].t=u;
			opt->edges[i].a=1-opt->edges[i].a;
		}
	}
  free(newlabel);
  qsort(opt->edges,opt->e,sizeof(edge),compare_edges);

	opt->cd=(unsigned*)calloc((opt->n+1),sizeof(unsigned));
	for (i=0;i<opt->e;i++){
		opt->cd[opt->edges[i].s+1]++;
	}
	for (i=0;i<opt->n;i++){
		opt->cd[i+1]+=opt->cd[i];
	}

  opt->ne=(double*)calloc(opt->n,sizeof(double));
  for (i=0;i<opt->e;i++){
    u=opt->edges[i].s;
    v=opt->edges[i].t;
    opt->ne[(u>v)?u:v]++;
  }
}


//to print in file the value for each node
void print_rates(optim* opt,char* rates){
	FILE *file=fopen(rates,"w");
	unsigned i;
	for (i=0;i<opt->n;i++){
		fprintf(file,"%u %e\n",opt->map[opt->nodes[i].n],opt->nodes[i].r);
	}
	fclose(file);
}



//fit data structure:
typedef struct {
	unsigned n;//total number of aggregated points
	unsigned *nag;//nag[i]=number of points aggregated in i
	double *val;//val[i]=value of the aggregated points
	double *ub;
} isoreg;


//Pool Adjacent Violators Algorithm. Values to fit in vect and size of vect.
isoreg *pava(double *vect,unsigned n){
	isoreg *fit=(isoreg*)malloc(sizeof(isoreg));
	unsigned *nag=(unsigned*)malloc(n*sizeof(unsigned));
	double *val=(double*)malloc(n*sizeof(double));
	unsigned i,j;

	nag[0]=1;
	val[0]=vect[0];
	j=0;
	for (i=1;i<n;i++){
		j+=1;
		val[j]=vect[i];
		nag[j]=1;
		while ((j>0) && (val[j]>val[j-1]*0.999999)){
			val[j-1]=(nag[j]*val[j]+nag[j-1]*val[j-1])/(nag[j]+nag[j-1]);
			nag[j-1]+=nag[j];
			j--;
		}
	}
	fit->n=j+1;
	fit->nag=nag;
	fit->val=val;
	return fit;
}


//printing the result in file output: "nag val ub" on each line
void print_fit(isoreg *fit,char* output){
	FILE *file=fopen(output,"w");
	unsigned i;
	for (i=0;i<fit->n;i++){
		fprintf(file,"%u %e %e\n",fit->nag[i],fit->val[i],fit->ub[i]);
	}
	fclose(file);
}

void freeisoreg(isoreg *fit){
	free(fit->ub);
	free(fit->nag);
	free(fit->val);
	free(fit);
}


//computing the upperbound on the density for the weak approx
void upperbounds(isoreg* fit, optim *opt){
	unsigned i,j;
	unsigned long long k;
	double *r=(double*)malloc(opt->n*sizeof(double));
	edge ed;

	fit->ub=(double*)malloc(fit->n*sizeof(double));
	for (k=0;k<opt->n;k++){
		r[k]=opt->nodes[k].r;
	}
	opt->cuts=(unsigned*)malloc(opt->n*sizeof(unsigned));
	i=0;
	for (k=0;k<fit->n;k++){
		fit->ub[k]=r[i];
		for (j=0;j<fit->nag[k];j++){
			opt->cuts[i++]=k;
		}
	}
	for (k=0;k<opt->e;k++){
		ed=opt->edges[k];
		if (opt->cuts[ed.s]!=opt->cuts[ed.t]){
			r[ed.t]+=ed.a;
			if (fit->ub[opt->cuts[ed.t]]<r[ed.t]){
				fit->ub[opt->cuts[ed.t]]=r[ed.t];//carefull
			}
		}
	}
	free(r);
}


//Step 3: Checking if the cuts given by PAVA are correct

isoreg* mkcut(isoreg* fit, optim *opt){
	unsigned i,k,u,ncuts=0;
	unsigned long long j1,j2;
	double *r=(double*)malloc(opt->n*sizeof(double));
	double *r2=(double*)malloc(opt->n*sizeof(double));
	double *r_tmp;
	double min, max;
	edge ed;
	unsigned long long *d=(unsigned long long*)calloc(opt->n,sizeof(unsigned long long));
	unsigned long long *d2=(unsigned long long*)malloc(opt->n*sizeof(unsigned long long));
	isoreg *fit2=(isoreg*)malloc(sizeof(isoreg));
	fit2->nag=(unsigned*)calloc(opt->n,sizeof(unsigned));
	fit2->val=(double*)calloc(opt->n,sizeof(double));
	fit2->n=0;


	for (k=0;k<opt->n;k++){
		r[k]=opt->nodes[k].r;
		r2[k]=r[k];
	}

	j1=0;
	j2=0;
	for (i=0;i<fit->n;i++){
		fit2->nag[ncuts]+=fit->nag[i];
		fit2->val[ncuts]+=fit->nag[i]*fit->val[i];
		j2+=fit->nag[i];
		for (u=j1;u<j2;u++){
			d2[u]=0;
			for (k=opt->cd[u]+d[u];k<opt->cd[u+1];k++){
				ed=opt->edges[k];
				if (ed.t>=j2){
					r2[u]-=ed.a;
					r2[ed.t]+=ed.a;
					d2[u]++;
				}
				else{
					break;
				}
			}
		}
		min=opt->nodes[0].r;
		for (k=j1;k<j2;k++){
			min=(min<r2[k])?min:r2[k];
		}
		max=0;
		for (k=j2;k<opt->n;k++){
			max=(max>r2[k])?max:r2[k];
		}

		if (max<min){
			for (k=j1;k<j2;k++){
				opt->cuts[k]=ncuts;
			}
			fit2->val[ncuts]/=fit2->nag[ncuts];
			for (u=j1;u<j2;u++){
				d[u]+=d2[u];
			}
			for (k=j1;k<opt->n;k++){
				r[k]=r2[k];
			}
			ncuts++;
			j1=j2;
		}
		else{
			for (k=j1;k<opt->n;k++){
				r2[k]=r[k];
			}
		}
	}
	fit2->n=ncuts;//+1;
	free(r);
	free(r2);
	free(d);
	free(d2);
	return fit2;
}

//Computing error: mult0,add0,multMAX,addMAX,multAVE,addAVE
double* error(isoreg* fit){
	unsigned k;
	double add,mult;
	double *err=(double*)calloc(6,sizeof(double));
	err[0]=fit->ub[0]/fit->val[0]-1.;
	err[1]=fit->ub[0]-fit->val[0];
	for (k=0;k<fit->n;k++){
		mult=fit->ub[k]/fit->val[k]-1.;
		add=fit->ub[k]-fit->val[k];
		err[2]=(err[2]>mult)?err[2]:mult;
		err[3]=(err[3]>add)?err[3]:add;
		err[4]+=mult;
		err[5]+=add;
	}
	err[4]/=(double)(fit->n);
	err[5]/=(double)(fit->n);
	return err;
}


int main(int argc,char** argv){
	optim* opt;
  isoreg *fit,*fit2;
	unsigned nsgs,i,k;
	unsigned nthreads=atoi(argv[1]);
	unsigned rep=atoi(argv[2]);
	double *err;
  char* edgelist=argv[3];
  char* rates=argv[4];
  char* pavafit=argv[5];
	char* cuts=argv[6];

  omp_set_num_threads(nthreads);

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;
	printf("- Reading edgelist from file %s\n",edgelist);
	opt=readedgelist(edgelist);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	printf("- Building the datastructure\n");
	t1=time(NULL);
	relabel(opt);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	printf("- Building the datastructure\n");
	printf("- Number of nodes = %u\n",opt->n);
	printf("- Number of edges = %llu\n",opt->e);
	printf("- Computing the locally densest decomposition\n");

  printf("- Step 1: Frank-Wolf gradiant descent (%u iterations)\n",rep);
	t1=time(NULL);
	init(opt);
	for (i=0;i<rep;i++){
		printf("%u\n",i);
		onepass(opt);
	}
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
  
  //print_alphas(opt,alphas);////////////////////////////

  printf("- Step 2: Isotonic regression with PAVA\n");
	t1=time(NULL);
  prepava(opt);
  fit=pava(opt->ne,opt->n);
	upperbounds(fit,opt);
	print_fit(fit,pavafit);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

	err=error(fit);
	printf("- error for the weak approximation (i.e. without checking whether the cuts are correct):\n");
	printf("- additive (densest,max,average) = %e %e %e\n",err[1],err[3],err[5]);
	printf("- multiplicative (densest,max,average) = %e %e %e\n",err[0],err[2],err[4]);
	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

  print_rates(opt,rates);

  printf("- Step 3: Checking if the %u cuts given by PAVA are correct\n",fit->n);
	t1=time(NULL);
	fit2=mkcut(fit,opt);
	freeisoreg(fit);
	upperbounds(fit2,opt);
  print_fit(fit2,cuts);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

	err=error(fit2);
	printf("- error for the strong approximation (i.e. taking only the correct cuts):\n");
	printf("- additive (densest,max,average) = %e %e %e\n",err[1],err[3],err[5]);
	printf("- multiplicative (densest,max,average) = %e %e %e\n",err[0],err[2],err[4]);
	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));
	return 0;
}

