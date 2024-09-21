#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

typedef struct treenode treenode;

double answer;
int num_direct;

struct treenode
{
	//boundary box and centre
	double xmin, xmax;
	double centre;

	//the start and end element of the array which contains particles in this node
	int start;
	int end;

	// pointers to the child nodes
	treenode *leftchild;
	treenode *rightchild;

};

void swap(float *array, int Index1, int Index2){

	float dummy;

	dummy = array[Index2];
	array[Index2] = array[Index1];
	array[Index1] = dummy;

}

int partition(float *array, int lo, int hi){
	
	int pivotIndex = rand() % (hi-lo) + lo ; 
	
	float pivotValue = array[pivotIndex];

	swap(array, pivotIndex, hi);

	int storeIndex = lo;

	int ii;

	for(ii=lo;ii<hi;ii++){

		if(array[ii] <= pivotValue){
			
			swap(array, ii, storeIndex);

			storeIndex = storeIndex + 1;			

		}

	}

	swap(array, storeIndex, hi);

	return storeIndex;

}

void QuickSort(float *array,int lo, int hi){

	int p;

	if(lo<hi){
	
		p = partition(array,lo, hi);
		QuickSort(array,lo,p-1);
		QuickSort(array,p+1,hi);

	}

}


treenode* kdtree(float *array, int lo, int hi, int maxnumber){

	//printf("Just entered kdtree\n");

	treenode *node = (treenode*)malloc(sizeof(treenode));

	//printf("Made the node\n");

	node->start = lo;
	node->end = hi;

	int medianIndex;
	int si;

	if(node->end - node->start > maxnumber) {

		medianIndex = (int) ceil(0.5*(node->end-node->start)) + node->start;
		//printf("The node end is %d, and node start is %d, the median index is %d\n",node->end, node->start, medianIndex);

		node->leftchild = kdtree(array,node->start,medianIndex,maxnumber);
	 	node->rightchild = kdtree(array,medianIndex,node->end,maxnumber);
	}

	else{

		//printf("Am leaf and here are my array elements\n");

		//for(int ii=node->start;ii<node->end;ii++) printf("Element %d, entry = %lf\n", ii, array[ii]);

		node->leftchild = NULL;
		node->rightchild = NULL;

	}

	return node;

}

void BuildBoundaries(float *array, treenode *node){

	int ii;
	double maxx=-1e6,minx=1e6;
	
	treenode *left, *right;

	if((node->leftchild==NULL && node->rightchild!=NULL) || (node->leftchild!=NULL && node->rightchild==NULL)) printf("something is wrong\n");

	else if(node->leftchild!=NULL && node->rightchild!=NULL){

		BuildBoundaries(array,node->leftchild);
		BuildBoundaries(array,node->rightchild);

		left = node->leftchild;
		right = node->rightchild;

		if(left->xmax > right->xmax) node->xmax = left->xmax;
		else node->xmax = right->xmax;

		if(left->xmin < right->xmin) node->xmin = left->xmin;
		else node->xmin = right->xmin;

		node->centre = (left->centre*(left->end - left->start) + right->centre*(right->end - right->start)) / (left->end - left->start + right->end - right->start);

	}
	else{

		node->centre = 0;

		for(ii=node->start;ii<node->end;ii++){

			if(array[ii] > maxx) maxx = array[ii];
			if(array[ii] < minx) minx = array[ii];

			node->centre = node->centre + array[ii];

		}

		node->xmax = maxx;
		node->xmin = minx;

		node->centre = node->centre / (node->end - node->start);


	}


}

treenode* BuildTree(float *array, int hi, int maxnumber){

	treenode *rootnode = (treenode*)malloc(sizeof(treenode));

	QuickSort(array,0,hi-1);
	//for(int ii=0;ii<hi;ii++) printf("%lf\n",array[ii]);
	rootnode = kdtree(array,0,hi,maxnumber);
	BuildBoundaries(array, rootnode);

	return rootnode;

}


void WalkTree(treenode *node, float *array, double position, double std, double error){

	int ii;

	double A,D;

	double opening_crit_dist;

	A = 1./sqrt(2*M_PI*std*std);

	D = 0.5/(std*std);

	opening_crit_dist = error*std*(node->end - node->start);

	//printf("The opening_crit_dist is %lf\n",opening_crit_dist);

	if(node->leftchild!=NULL && node->rightchild!=NULL){			

			if(position > node->xmin && position < node->xmax){

				//printf("Opening child: position was %lf and node centre was %lf and opening crit dist was %lf\n", position, node->centre, opening_crit_dist);
				WalkTree(node->leftchild,array,position,std,error);
				WalkTree(node->rightchild,array,position,std,error);
			}

			else if (fabs(position - node->centre) < opening_crit_dist){

				WalkTree(node->leftchild,array,position,std,error);
				WalkTree(node->rightchild,array,position,std,error);

			}

			else{

				//printf("The node is sufficient\n");

				answer = answer + (node->end - node->start)*A*exp(-D*(position-node->centre)*(position-node->centre));

			}

	}

	else{
		
		for(ii=node->start;ii<node->end;ii++){

			answer = answer + A*exp(-D*(position-array[ii])*(position-array[ii]));
			num_direct = num_direct + 1;

		}
	
	}


}




void Gaussian(double *x, int x_size, const float mean, double std, double *out){

	int ii=0;

	double A,D;

	A = 1./sqrt(2*M_PI*std*std);

	D = 0.5/(std*std);

	for(ii=0;ii<x_size;ii++){

		out[ii] = out[ii] + A*exp(-D*(x[ii]-mean)*(x[ii]-mean));

	}

}

double Bandwidth_calc(float *array, int a_size){

	int ii=0;

	double std,factor;
	double var=0;
	double mean=0;

	factor = pow((double)a_size,-0.2);

	for(ii=0;ii<a_size;ii++) mean = mean + array[ii];
	mean = mean/a_size;
	for(ii=0;ii<a_size;ii++) var = var + (array[ii]-mean)*(array[ii]-mean);
	var = var/a_size;

	std = factor*1.06*sqrt(var)/1.34;

	return std;

}

void KDE(double *x, int x_size, const float *array, int a_size, double bandwidth, double *output){

	int ii=0;

	for(ii=0;ii<a_size;ii++){

		Gaussian(x,x_size,array[ii],bandwidth,output);
	
	}

}

void DeleteTree(treenode *node){

	if(node->leftchild != NULL && node->rightchild != NULL){

		DeleteTree(node->leftchild);
		DeleteTree(node->rightchild);

		free(node->leftchild);
		free(node->rightchild);

	}

}

void approxKDE(double *x, int x_size, float *array, int a_size, double bandwidth, double error, double *output){

	treenode *rootnode = (treenode*)malloc(sizeof(treenode));

	//QuickSort(array,0,a_size);
	//printf("HERE and a_size is %d\n", a_size);


	struct timeval t1, t2;
	double elapsedTime;

	gettimeofday(&t1, NULL);

	rootnode = BuildTree(array, a_size, 8);

	gettimeofday(&t2, NULL);

	elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
	elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;

	printf("Building tree took %lf in ms\n", elapsedTime);



	gettimeofday(&t1, NULL);
	for(int ii=0;ii<x_size;ii++){

		answer = 0;
		num_direct = 0;

		WalkTree(rootnode,array,x[ii],bandwidth,error); 

		output[ii] = answer;

		//printf("At positon %lf, bandwidth of %lf, error of %lf, the output is %lf and number of direct eval is %d\n", x[ii], bandwidth,error,output[ii],num_direct);
	
	}


	gettimeofday(&t2, NULL);

	elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
	elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;

	printf("Walking tree took %lf in ms\n", elapsedTime);

	DeleteTree(rootnode);

}



void TwoPoint(const double *pos_x, const double *pos_y, int size, int nruns, const double *bounds, int x_size, double *DD, double *DR, double *RR){


	double xmin,xmax,ymin,ymax;

	xmin = bounds[0];
	xmax = bounds[1];
	
	ymin = bounds[2];
	ymax = bounds[3];

	int iDD=0;
	int iDR=0;
	int iRR=0;

	double rpos_x[size];
	double rpos_y[size];

	srand(1);	

	int DD_size=size*(size-1);
	int DR_size=nruns*size*size;
	int RR_size=nruns*size*(size-1);


	float *tot_DD = (float*)malloc(DD_size *sizeof(float));
	float *tot_RR = (float*)malloc(RR_size *sizeof(float));
	float *tot_DR = (float*)malloc(DR_size *sizeof(float));

	double x[x_size];

	double max_dist;

	max_dist = sqrt( (xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin) );

	for(int aa=0;aa<x_size;aa++){

		x[aa] = aa*max_dist/x_size;

	}

	for(int ii=0;ii<nruns;ii++){

		for(int aa=0; aa<size;aa++){
			rpos_x[aa] = (xmax-xmin)*(float)rand()/(float)(RAND_MAX) + xmin;
			rpos_y[aa] = (ymax-ymin)*(float)rand()/(float)(RAND_MAX) + ymin;

		}

		for(int jj=0;jj<size;jj++){

			for(int kk=0; kk<size; kk++){

				if(jj!=kk){

					if(ii==0){

						tot_DD[iDD] = sqrt( (pos_x[jj]-pos_x[kk])*(pos_x[jj]-pos_x[kk]) +  (pos_y[jj]-pos_y[kk])*(pos_y[jj]-pos_y[kk]));

						iDD = iDD + 1;

					}

					tot_RR[iRR] = sqrt( (rpos_x[jj]-rpos_x[kk])*(rpos_x[jj]-rpos_x[kk]) +  (rpos_y[jj]-rpos_y[kk])*(rpos_y[jj]-rpos_y[kk]));

					iRR = iRR + 1;

				}

				tot_DR[iDR] = sqrt( (pos_x[jj]-rpos_x[kk])*(pos_x[jj]-rpos_x[kk]) +  (pos_y[jj]-rpos_y[kk])*(pos_y[jj]-rpos_y[kk]));


				iDR = iDR + 1;

			}

		}

	}


	double bd;
	bd = Bandwidth_calc(tot_DR, DR_size);
	
	KDE(x,x_size,tot_DD,DD_size,bd,DD);


	struct timeval t1, t2;
	double elapsedTime;

	gettimeofday(&t1, NULL);

	KDE(x,x_size,tot_DR,DR_size,bd,DR);

	gettimeofday(&t2, NULL);

	elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
	elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;

	printf("The elapsed time in ms is %lf\n", elapsedTime);

	KDE(x,x_size,tot_RR,RR_size,bd,RR);

	for(int bb=0; bb<x_size;bb++){

		DD[bb] = DD[bb] / DD_size;
		DR[bb] = DR[bb] / DR_size;
		RR[bb] = RR[bb] / RR_size;

	}


}







void ApproxTwoPoint(const double *pos_x, const double *pos_y, int size, int nruns, const double *bounds, int x_size, double error, double *DD, double *DR, double *RR){


	double xmin,xmax,ymin,ymax;

	xmin = bounds[0];
	xmax = bounds[1];
	
	ymin = bounds[2];
	ymax = bounds[3];

	long int iDD=0;
	long int iDR=0;
	long int iRR=0;

	double rpos_x[size];
	double rpos_y[size];

	srand(1);	

	int DD_size=size*(size-1);
	long int DR_size;
	long int RR_size;

        DR_size = nruns*size;
        DR_size = DR_size*size;
        RR_size = nruns*size;
        RR_size = RR_size*(size-1);
	printf("before allocation\n");
	float *tot_DD = (float*)malloc(DD_size *sizeof(float));
	float *tot_RR = (float*)malloc(RR_size *sizeof(float));
	float *tot_DR = (float*)malloc(DR_size *sizeof(float));
	printf("after allocation\n");
	double max_dist;

	max_dist = sqrt( (xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin) );
        
        
        /*
        int i;
        double min_dist = 0;
        double xstart = 2;
        double xend = max_dist;
        double Lmin = 2;
        int bins;
        bins = (int)((log10((xend)/(xstart)))/(log10(1+(Lmin)/(xstart)))) + 1;
        Lmin = xstart*(pow((xend/xstart),(1.0/(bins-1))) - 1);
        double logxstart = log10(xstart);				
        double logxend = log10(xend);
        double logbinsize = (logxend - logxstart)/(bins-1);
        
        
        double *x;
        x = (double*)malloc(bins *sizeof(double));
        for (i = 0; i < bins; i++) x[i] = 0.0;
        for (i = 0; i < bins; i++)
        {
//                x[i] = pow(10, (logxstart + i*logbinsize));
                x[i] = logxstart + i*logbinsize;
        }
        */
        
        double xstart = 2;
        double xend = max_dist;
        double logxstart = log10(xstart);
        double logxend = log10(xend);
        double x[x_size];
	for(int aa=0;aa<x_size;aa++){

		x[aa] = logxstart + aa*(logxend-logxstart)/x_size;
        }
        
        printf("%d\n", size);
        printf("%d\n", nruns);
	printf("%ld\n", RR_size);
	int ii;
	for(ii=0;ii<nruns;ii++){
//		printf("%d\n", ii);
		for(int aa=0; aa<size;aa++){
			rpos_x[aa] = (xmax-xmin)*(float)rand()/(float)(RAND_MAX) + xmin;
			rpos_y[aa] = (ymax-ymin)*(float)rand()/(float)(RAND_MAX) + ymin;

		}

		for(int jj=0;jj<size;jj++){

			for(int kk=0; kk<size; kk++){

				if(jj!=kk){

					if(ii==0){

						tot_DD[iDD] = sqrt( (pos_x[jj]-pos_x[kk])*(pos_x[jj]-pos_x[kk]) +  (pos_y[jj]-pos_y[kk])*(pos_y[jj]-pos_y[kk]));
                                                
                                                tot_DD[iDD] = log10(tot_DD[iDD]);

						iDD = iDD + 1;

					}

					tot_RR[iRR] = sqrt( (rpos_x[jj]-rpos_x[kk])*(rpos_x[jj]-rpos_x[kk]) +  (rpos_y[jj]-rpos_y[kk])*(rpos_y[jj]-rpos_y[kk]));
                                        
                                        tot_RR[iRR] = log10(tot_RR[iRR]);

					iRR = iRR + 1;

				}
//				if (ii > 210)printf("%d\n", iDR);
				tot_DR[iDR] = sqrt( (pos_x[jj]-rpos_x[kk])*(pos_x[jj]-rpos_x[kk]) +  (pos_y[jj]-rpos_y[kk])*(pos_y[jj]-rpos_y[kk]));
                                
                                tot_DR[iDR] = log10(tot_DR[iDR]);


				iDR = iDR + 1;

			}

		}

	}


	double bd;
	bd = Bandwidth_calc(tot_DR, DR_size);
        printf("b = %lf\n", bd);


	KDE(x,x_size,tot_DD,DD_size,bd,DD);


	struct timeval t1, t2;
	double elapsedTime;

	gettimeofday(&t1, NULL);

	approxKDE(x,x_size,tot_DR,DR_size,bd,error,DR);


	gettimeofday(&t2, NULL);

	elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
	elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;

	printf("The elapsed time in ms is %lf\n", elapsedTime);


	approxKDE(x,x_size,tot_RR,RR_size,bd,error,RR);

	for(int bb=0; bb<x_size;bb++){

		DD[bb] = DD[bb] / DD_size;
		DR[bb] = DR[bb] / DR_size;
		RR[bb] = RR[bb] / RR_size;

	}

	//for(int ii=0;ii<DD_size;ii++) free(tot_DD[ii]);
	//for(int ii=0;ii<DR_size;ii++) free(tot_DR[ii]);	
	//for(int ii=0;ii<RR_size;ii++) free(tot_RR[ii]);

	free(tot_DD);
	free(tot_DR);
	free(tot_RR);

}



