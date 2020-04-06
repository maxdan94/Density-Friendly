/*
Maximilien Danisch
Mai 2016
http://bit.ly/maxdan94
maximilien.danisch@telecom-paristech.fr

Info:
Feel free to use these lines as you wish. This program computes the exact density-friendly decomposition.

To compile:
g++ exactDF.cpp -fopenmp -fpermissive -o exactDF -O3

To execute:
./exactDF ncpu iter net.txt rates.txt pavafit.txt cuts.txt exact.txt

- nthreads is the number of threads to use
- iter is the number of iterations over all edges to perform
- net.txt should contain the graph (one edge on each line: 2 unsigned separated by a space)
- rates.txt will contain the density value for each node
- pavafit.txt will contain the profile given by the PAVA fit, that is the week approximation of the density-friendly ("size density density-upperbound" on each line)
- cuts.txt will contain the profile given by correct cuts, that is the strong approximation of the density-friendly ("size density density-upperbound" on each line)
- exact.txt will contain the exact density-friendly ("size density" on each line (it is not necesarily in decreasing order of density))
Some information will be printed in the terminal.


*/

#include <boost/graph/edge_list.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <algorithm>
#include <utility>
#include <stdio.h> 
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define NLINKS 500000000 //maximum number of edges for memory allocation, will increase if needed

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
	unsigned e;//number of edges
	unsigned *map;//correspondance between old and new nodeID
	edge *edges;//list of all edges
	node *nodes;//value associated to each node
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
	unsigned e1=NLINKS;
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
	unsigned i,j;
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
  opt->iter=0;
	opt->nodes=(node*)malloc(opt->n*sizeof(node));
  for (unsigned k=0;k<opt->n;k++){
    opt->nodes[k].n=k;
    opt->nodes[k].r=0;
  }
	for (unsigned k=0;k<opt->e;k++){
    opt->edges[k].a=.5;
		opt->nodes[opt->edges[k].s].r+=.5;
		opt->nodes[opt->edges[k].t].r+=.5;
	}
}

//one pass over all edges
void onepass(optim *opt){
	unsigned i,j,k;
	double gamma;
	opt->iter++;
	gamma=2./(2.+opt->iter);
	#pragma omp parallel for private(i,j,k)
	for (k=0;k<opt->e;k++){//parfor
		i=opt->edges[k].s;
		j=opt->edges[k].t;
		if (opt->nodes[i].r<opt->nodes[j].r){
			#pragma omp atomic update
			opt->nodes[i].r+=gamma*(1-opt->edges[k].a);//carefull
			#pragma omp atomic update
			opt->nodes[j].r-=gamma*(1-opt->edges[k].a);
			opt->edges[k].a=(1.-gamma)*opt->edges[k].a+gamma;
		}
		else if (opt->nodes[i].r>opt->nodes[j].r){
			#pragma omp atomic update
			opt->nodes[i].r-=gamma*opt->edges[k].a;//carefull
			#pragma omp atomic update
			opt->nodes[j].r+=gamma*opt->edges[k].a;
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

//used for quicksort (greatest hit in CS before this algorithm)
static int compare_nodes(void const *a, void const *b){
	if (((node*)a)->r <= ((node*)b)->r)
		return 1;
	return -1;
}
//used for quicksort (greatest hit in CS before this algorithm)
static int compare_edges(void const *a, void const *b){
	edge *pa = (edge *)a;
	edge *pb = (edge *)b;
	if ((*pa).s>(*pb).s)
		return 1;
	if ((*pa).s<(*pb).s)
		return -1;
	if ((*pa).t<=(*pb).t)
		return 1;
	return -1;
}

void prepava(optim *opt){
  unsigned u,v;
  unsigned *newlabel=(unsigned*)malloc(opt->n*sizeof(unsigned));
  qsort(opt->nodes,opt->n,sizeof(node),compare_nodes);

  for (unsigned i=0;i<opt->n;i++){
    newlabel[opt->nodes[i].n]=i;
  }
  for (unsigned i=0;i<opt->e;i++){
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
  //free(newlabel);
  qsort(opt->edges,opt->e,sizeof(edge),compare_edges);

	opt->cd=(unsigned*)calloc((opt->n+1),sizeof(unsigned));
	for (unsigned i=0;i<opt->e;i++){
		opt->cd[opt->edges[i].s+1]++;
	}
	for (unsigned i=0;i<opt->n;i++){
		opt->cd[i+1]+=opt->cd[i];
	}

  opt->ne=(double*)calloc(opt->n,sizeof(double));
  for (unsigned i=0;i<opt->e;i++){
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
		while ((j>0) && (val[j]>val[j-1]-1e-10)){//do val[j]>val[j-1] to have a non-increasing monotonic regression.
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

//printing the result in file output: "nag val" on each line
void print_fit(isoreg *fit,char* output){
	FILE *file=fopen(output,"w");
	for (unsigned i=0;i<fit->n;i++){
		fprintf(file,"%u %.10le\n",fit->nag[i],fit->val[i]);
	}
	fclose(file);
}

void freeisoreg(isoreg *fit){
	free(fit->nag);
	free(fit->val);
	free(fit);
}

//Step 3: Checking if the cuts given by PAVA are correct

//printing the result in file output: "nag val" on each line
isoreg* mkcut(isoreg* fit, optim *opt){
	unsigned j1,j2,ncuts=0;
	double *r=(double*)malloc(opt->n*sizeof(double));
	double *r2=(double*)malloc(opt->n*sizeof(double));
	double min, max;
	edge ed;
	unsigned *d=(unsigned*)calloc(opt->n,sizeof(unsigned));
	unsigned *d2=(unsigned*)malloc(opt->n*sizeof(unsigned));
	isoreg *fit2=(isoreg*)malloc(sizeof(isoreg));
	fit2->nag=(unsigned*)calloc(opt->n,sizeof(unsigned));
	fit2->val=(double*)calloc(opt->n,sizeof(double));
	fit2->n=0;

	opt->cuts=(unsigned*)malloc(opt->n*sizeof(unsigned));

	for (unsigned k=0;k<opt->n;k++){
		r[k]=opt->nodes[k].r;
		r2[k]=r[k];
	}

	j1=0;
	j2=0;
	for (unsigned i=0;i<fit->n;i++){
		fit2->nag[ncuts]+=fit->nag[i];
		fit2->val[ncuts]+=fit->nag[i]*fit->val[i];
		j2+=fit->nag[i];
		for (unsigned u=j1;u<j2;u++){
			d2[u]=0;
			for (unsigned k=opt->cd[u]+d[u];k<opt->cd[u+1];k++){
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
		for (unsigned k=j1;k<j2;k++){
			min=(min<r2[k])?min:r2[k];
		}
		max=0;
		for (unsigned k=j2;k<opt->n;k++){
			max=(max>r2[k])?max:r2[k];
		}
		if (max<min){
			for (unsigned k=j1;k<j2;k++){
				opt->cuts[k]=ncuts;
			}
			fit2->val[ncuts]/=fit2->nag[ncuts];
			for (unsigned u=j1;u<j2;u++){
				d[u]+=d2[u];
			}
			for (unsigned k=j1;k<opt->n;k++){///////////
				r[k]=r2[k];
			}
			ncuts++;
			j1=j2;
		}
		else{
			for (unsigned k=j1;k<opt->n;k++){///////////
				r2[k]=r[k];
			}
		}
	}
	fit2->n=ncuts;//+1;
	return fit2;
}

//Step 4: Density-Friendly in each indepent subgraphs

typedef struct {
	unsigned s;
	unsigned t;
} edge2;

typedef struct {
	unsigned n;//number of nodes
	unsigned e;//number of edges
	unsigned l;//number of loops
	edge2 *edges;//list of edges
	unsigned *loops;//loops[i]=number of loops of node i
	unsigned *deg;//deg[i]=degree of node i
} subgraph;

void freesubgraph(subgraph* sg){
	free(sg->edges);
	free(sg->loops);
	free(sg->deg);
	free(sg);
}

subgraph **mkindep(optim* opt,unsigned *nsub){
	*nsub=opt->cuts[opt->n-1]+1;
	subgraph **sgs=(subgraph**)malloc((*nsub)*sizeof(subgraph*));
	unsigned *newlabel=(unsigned*)malloc(opt->n*sizeof(unsigned));
	unsigned j,p,u,v,u2,v2;
	for (unsigned k=0;k<*nsub;k++){
		sgs[k]=(subgraph*)malloc(sizeof(subgraph));
		sgs[k]->n=0;
		sgs[k]->e=0;
		sgs[k]->l=0;
	}
	j=0;
	for (unsigned k=0;k<opt->n;k++){
		sgs[opt->cuts[k]]->n++;
		newlabel[k]=j++;
		if ((k<opt->n-1) && (opt->cuts[k] != opt->cuts[k+1])){
			j=0;
		}
	}
	for (unsigned k=0;k<opt->e;k++){
		u=opt->edges[k].s;
		v=opt->edges[k].t;
		p=opt->cuts[v];
		if (opt->cuts[u]==p){
			sgs[p]->e++;
		}
	}
	for (unsigned k=0;k<*nsub;k++){
		sgs[k]->edges=(edge2*)malloc(sgs[k]->e*sizeof(edge2));
		sgs[k]->loops=(unsigned*)calloc(sgs[k]->n,sizeof(unsigned));
		sgs[k]->deg=(unsigned*)calloc(sgs[k]->n,sizeof(unsigned));
		sgs[k]->e=0;
	}
	for (unsigned k=0;k<opt->e;k++){
		u=opt->edges[k].s;
		v=opt->edges[k].t;
		u2=newlabel[u];
		v2=newlabel[v];
		p=opt->cuts[v];
		if (opt->cuts[u]==p){
			sgs[p]->edges[sgs[p]->e].s=u2;
			sgs[p]->edges[sgs[p]->e++].t=v2;
			sgs[p]->deg[u2]++;
			sgs[p]->deg[v2]++;
		}
		else{
			sgs[p]->loops[v2]++;
			sgs[p]->l++;
		}
	}
	return sgs;
}

subgraph *allocsubgraph(unsigned n,unsigned e){
	subgraph *sg=(subgraph*)malloc(sizeof(subgraph));
	sg->n=n;
	sg->e=0;
	sg->l=0;
	sg->deg=(unsigned*)calloc(n,sizeof(unsigned));
	sg->loops=(unsigned*)calloc(n,sizeof(unsigned));
	sg->edges=(edge2*)malloc(e*sizeof(edge2));
	return sg;
}

using namespace boost;
typedef adjacency_list_traits < vecS, vecS, directedS > Traits; //the associated types of the adjacency_list class
typedef adjacency_list < vecS, vecS, directedS,
	    property < vertex_name_t, unsigned,
	    property < vertex_index_t, long,
	    property < vertex_color_t, boost::default_color_type,
	    property < vertex_distance_t, double,///////////???
	    property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,
	    property < edge_capacity_t, double,
	    property < edge_residual_capacity_t, double,
	    property < edge_reverse_t, Traits::edge_descriptor > > > > Graph; //the associated properties of the ajacency_list

void AddEdge(Traits::vertex_descriptor &v1, Traits::vertex_descriptor &v2, property_map < Graph, edge_reverse_t >::type &rev, const double capacity, Graph &g) {
  Traits::edge_descriptor e1 = add_edge(v1, v2, g).first;
  Traits::edge_descriptor e2 = add_edge(v2, v1, g).first;
  put(edge_capacity, g, e1, capacity);
   rev[e1] = e2;
   rev[e2] = e1;
}

void ReReadGraph(Traits::vertex_descriptor &s, Traits::vertex_descriptor &t, subgraph *sg, Graph &g) {
	graph_traits<Graph>::vertex_descriptor u, v;
	property_map < Graph, edge_reverse_t >::type rev = get(edge_reverse, g);
	double alpha=(sg->e+sg->l)/((double)(sg->n))+1./(((double)sg->n)*((double)sg->n));
	double capacity1=2*alpha, capacity2, capacity3=1.;
	for (unsigned i=0; i<sg->n; i++) {
		u=add_vertex(g);
		capacity2=(sg->deg[i]+2.*sg->loops[i]);
		AddEdge(s,u, rev, capacity2, g);
		AddEdge(u,t, rev, capacity1, g);
	}
	for(unsigned i=0; i<sg->e; i++){
		u=vertex(sg->edges[i].s+2,g);
		v=vertex(sg->edges[i].t+2,g);
		AddEdge(u,v,rev, capacity3, g);
		AddEdge(v,u,rev, capacity3, g);
	}
}

void maxflow(subgraph* sg,FILE *file) {
	Traits::vertex_descriptor s, t;
	Graph g;
	property_map < Graph, vertex_color_t >::type col = get(vertex_color, g);
	s=add_vertex(g);
  t=add_vertex(g);
	ReReadGraph(s, t, sg, g);
	double flow = boykov_kolmogorov_max_flow(g ,s, t);
	//std::cout << "number of nodes: " << sg->n << std::endl;
	//std::cout << "number of edges: " << sg->e << std::endl;
	//std::cout << "number of loops: " << sg->l << std::endl;
	//std::cout << "flow: " << flow << std::endl;

	bool *isin=(bool *)calloc(sg->n,sizeof(bool));
	unsigned n_nodes=0;
	for(unsigned i=0;i<sg->n;i++){
		if (col[vertex(i+2,g)]==4){
			isin[i]=1;
			n_nodes++;
			//std::cout << i << " " << col[vertex(i+2,g)] << std::endl;
		}
	}
	g.clear();
	if (n_nodes>0){
		subgraph *sg1=allocsubgraph(n_nodes,sg->e);
		subgraph *sg2=allocsubgraph(sg->n-n_nodes,sg->e);
		unsigned *newlabel=(unsigned*)malloc(sg->n*sizeof(unsigned));
		unsigned j1=0,j2=0;
		for(unsigned i=0;i<sg->n;i++){
			if (isin[i]){
				sg1->loops[j1]=sg->loops[i];
				sg1->l+=sg->loops[i];
				newlabel[i]=j1++;
			}
			else{
				sg2->loops[j2]=sg->loops[i];
				sg2->l+=sg->loops[i];
				newlabel[i]=j2++;
			}
		}
		j1=0,j2=0;
		for(unsigned i=0;i<sg->e;i++){
			if (isin[sg->edges[i].s] && isin[sg->edges[i].t]){
				sg1->edges[j1].s=newlabel[sg->edges[i].s];
				sg1->edges[j1].t=newlabel[sg->edges[i].t];
				sg1->deg[sg1->edges[j1].s]++;
				sg1->deg[sg1->edges[j1].t]++;
				j1++;
			}
			else if (isin[sg->edges[i].s] && (isin[sg->edges[i].t]==0)){
				sg2->loops[newlabel[sg->edges[i].t]]++;
				sg2->l++;
			}
			else if ((isin[sg->edges[i].s]==0) && isin[sg->edges[i].t]){
				sg2->loops[newlabel[sg->edges[i].s]]++;
				sg2->l++;
			}
			else{
				sg2->edges[j2].s=newlabel[sg->edges[i].s];
				sg2->edges[j2].t=newlabel[sg->edges[i].t];
				sg2->deg[sg2->edges[j2].s]++;
				sg2->deg[sg2->edges[j2].t]++;
				j2++;
			}
		}
		sg1->e=j1;
		sg2->e=j2;
		free(isin);
		freesubgraph(sg);
		maxflow(sg1,file);
		maxflow(sg2,file);
	}
	else {
		//<< "n,e,d = "
		#pragma omp critical
		{
			fprintf(file,"%u %le\n",sg->n,((double)(sg->e+sg->l))/sg->n);
			//fprintf(file,"%u %u %u %le\n",sg->n,sg->e,sg->l,((double)(sg->e+sg->l))/sg->n);
			//std::cout << sg->n << " " << sg->e << " " << sg->l << " " << ((double)(sg->e+sg->l))/sg->n << std::endl;
		}
		free(isin);
		freesubgraph(sg);
	}
}


int main(int argc,char** argv){
	optim* opt;
  isoreg *fit,*fit2;
	subgraph** sgs;
	unsigned nsgs;
	unsigned nthreads=atoi(argv[1]);
	unsigned rep=atoi(argv[2]);
  char* edgelist=argv[3];
  char* rates=argv[4];
  char* pavafit=argv[5];
	char* cuts=argv[6];
	char* exact=argv[7];

  omp_set_num_threads(nthreads);

	time_t t1,t2,t3;
	t1=time(NULL);
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
	printf("- Number of edges = %u\n",opt->e);
	printf("- Computing the locally densest decomposition\n");

  printf("- Step 1: Frank-Wolf gradiant descent (%u iterations)\n",rep);
	t1=time(NULL);
	init(opt);
	for (unsigned i=0;i<rep;i++){
		//printf("%u\n",i);
		onepass(opt);
	}
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
  //print_alphas(opt,alphas);

  printf("- Step 2: Isotonic regression with PAVA\n");
	t1=time(NULL);
  prepava(opt);
  fit=pava(opt->ne,opt->n);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	
  print_fit(fit,pavafit);
  print_rates(opt,rates);

  printf("- Step 3: Checking if the %u cuts given by PAVA are correct\n",fit->n);
	t1=time(NULL);
	fit2=mkcut(fit,opt);
	freeisoreg(fit);
  print_fit(fit2,cuts);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

  printf("- Step 4: Density-Friendly in each %u indepent subgraphs\n",fit2->n-1);
	freeisoreg(fit2);
	t1=time(NULL);
	sgs=mkindep(opt,&nsgs);
	freeoptim(opt);
	FILE *file=fopen(exact,"w");
	unsigned k;
	#pragma omp parallel for schedule(dynamic) private(k)
	for(k=0;k<nsgs;k++){
		maxflow(sgs[k],file);
	}
	fclose(file);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

	return 0;
}
