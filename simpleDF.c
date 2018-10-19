/*
Maximilien Danisch
October 2015
http://bit.ly/maxdan94
maximilien.danisch@telecom-paristech.fr

Info:
Feel free to use these lines as you wish.
This program computes an approximation of the density-friendly decomposition and a densest subgraph (highest average degree) in very large graphs.

To compile:
"gcc simpleDF.c -lm -O9 -o simpleDF".

To execute:
"./simpleDF edgelist.txt iter densest.txt decomp.txt".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
"iter" is the number of iterations, e.g. "100".
"densest.txt" contains the found subgraph: "upperbound density size node1 node2 node3...". "upperbound" is an upperbound on the maximum density
"decomp.txt" contains on each line a node and it's density-score

If you also want the statistics (upperbound density and size of the found subgraph at each iteration) type:
"./ds edgelist.txt iter densest.txt decomp.txt stat.txt".
"stat.txt" will contain these statistics in the following format: "iteration density size upperbound".
It will be a bit slower, but still fast ;).

*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define NLINKS 500000000 //maximum number of edges for memory allocation, will increase if needed


typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	unsigned id;
	unsigned idr;
	unsigned long val;
} idval;

typedef struct {
	unsigned size;
	double density;
	double uppbound;//upper bound on the maximum density
} subgraph;

typedef struct {
	unsigned n; //number of nodes
	unsigned e; //number of edges
	unsigned *map;//correspondance between labels
	edge *edges;//list of all edges
	long unsigned *y; //value associated to each node
	idval *rank;//used to sort and manage the values
	long unsigned *ne;//ne[i]=number of edges from node ranked i to earlier ranked nodes
	unsigned iter;//number of iterations performed
} optim;


//compute the maximum of three unsigned
inline unsigned max3(unsigned a,unsigned b,unsigned c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

optim* readedgelist(char* edgelist){
	unsigned e1=NLINKS;
	optim *opt=malloc(sizeof(optim));
	FILE *file;

	opt->n=0;
	opt->e=0;
	file=fopen(edgelist,"r");
	opt->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%u %u", &(opt->edges[opt->e].s), &(opt->edges[opt->e].t))==2) {
		opt->n=max3(opt->n,opt->edges[opt->e].s,opt->edges[opt->e].t);
		if (opt->e++==e1) {
			e1+=NLINKS;
			opt->edges=realloc(opt->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	opt->n++;

	opt->edges=realloc(opt->edges,opt->e*sizeof(edge));

	return opt;
}

void relabel(optim *opt) {
	unsigned i,j;
	unsigned *newlabel;

	newlabel=malloc(opt->n*sizeof(unsigned));
	for (i=0;i<opt->n;i++) {
		newlabel[i]=opt->n;
	}

	opt->map=malloc(opt->n*sizeof(unsigned));
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
	opt->map=realloc(opt->map,opt->n*sizeof(unsigned));
	opt->y=calloc(opt->n,sizeof(long unsigned));
	opt->rank=malloc(opt->n*sizeof(idval));
	opt->ne=malloc(opt->n*sizeof(long unsigned));
	opt->iter=0;
}


//allocating memory for the subgraph
subgraph *allocsg(optim *opt){
	subgraph *g=malloc(sizeof(subgraph));
	g->uppbound=(double)opt->e;
	return g;
}

//one path over all edges
inline void onepass(optim *opt){
	unsigned i,j,k;
	opt->iter++;
	for (k=0;k<opt->e;k++){
		i=opt->edges[k].s;
		j=opt->edges[k].t;
		if (opt->y[i]<=opt->y[j]){
			opt->y[i]++;
		}
		else{
			opt->y[j]++;
		}
	}
}

//used for quicksort (greatest hit in CS before this algorithm)
static int compare (void const *a, void const *b){
	idval const *pa = a;
	idval const *pb = b;

	if ((*pa).val<=(*pb).val)
		return 1;
	return -1;
}

//make a swap along the ranking to find the densest subgraph
void mkdensest(optim *opt, subgraph *g){
	unsigned i, j, k;
	unsigned m;
	double tmp=0.;

	for (i=0;i<opt->n;i++){
		opt->rank[i].id=i;
		opt->rank[i].val=opt->y[i];
		opt->ne[i]=0;
	}

	qsort(opt->rank,opt->n,sizeof(idval),compare);

	for (i=0;i<opt->n;i++){
		opt->rank[opt->rank[i].id].idr=i;
	}

	for (k=0;k<opt->e;k++){
		i=opt->rank[opt->edges[k].s].idr;
		j=opt->rank[opt->edges[k].t].idr;
		opt->ne[(i>j)?i:j]++;
	}

	if (g->uppbound>opt->rank[0].val/((double)opt->iter)){
		g->uppbound=opt->rank[0].val/((double)opt->iter);
	}

	g->density=0.;
	g->size=0;
	for (i=1;i<opt->n;i++){
		opt->ne[i]+=opt->ne[i-1];
		tmp=((double)(opt->ne[i]))/((double)i+1.);
		if (tmp>g->density){
			g->density=tmp;
			g->size=i+1;
		}
	}
}



void freeoptim(optim *opt){
	free(opt->map);
	free(opt->edges);
	free(opt->y);
	free(opt->rank);
	free(opt->ne);
	free(opt);
}

int main(int argc,char** argv){
	optim* opt;
	subgraph* g;
	unsigned i,rep=atoi(argv[2]);
	FILE *file;

	printf("Reading edgelist from file %s\n",argv[1]);
	opt=readedgelist(argv[1]);
	printf("Building the datastructure\n");
	relabel(opt);
	printf("Number of nodes = %u\n",opt->n);
	printf("Number of edges = %u\n",opt->e);
	printf("Computing the densest subgraph\n");

	if (argc==6){
		file=fopen(argv[5],"w");
		setvbuf(file,NULL,_IONBF,0);
		g=allocsg(opt);
		for (i=0;i<rep;i++){
			onepass(opt);
			mkdensest(opt,g);
			printf("%u %le %u %le\n",i+1,g->density,g->size,g->uppbound);
			fprintf(file,"%u %le %u %le\n",i+1,g->density,g->size,g->uppbound);
		}
		fclose(file);
	}

	else{
		for (i=0;i<rep;i++){
			printf("%u\n",i+1);
			onepass(opt);
		}
		g=allocsg(opt);
		mkdensest(opt,g);
	}

	printf("density = %le\n",g->density);
	printf("size = %u\n",g->size);
	printf("density upper bound = %le\n",g->uppbound);

	file=fopen(argv[3],"w");
	fprintf(file,"%le %le %u",g->uppbound,g->density,g->size);
	for (i=0;i<g->size;i++){
		fprintf(file," %u",opt->map[opt->rank[i].id]);
	}
	fprintf(file,"\n");
	fclose(file);

	file=fopen(argv[4],"w");
	for (i=0;i<opt->n;i++){
		fprintf(file,"%u %le\n",opt->map[opt->rank[i].id],((double)(opt->rank[i].val))/((double)rep));
	}
	fprintf(file,"\n");
	fclose(file);


	freeoptim(opt);
	free(g);
	return 0;
}
