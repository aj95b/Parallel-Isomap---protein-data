/***@Author:
             Arpita Joshi
    Parallelizing the loops of IsomapII code of TdSL using OpenMP                 ***/

#include<omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<complex.h>
#include<stdbool.h>


//#include"dijkstra.cpp"
double* data;
int columnCount, rowCount;
int k;
int dims = 5;
// #ifdef __cplusplus
//extern "C"
//#endif
//void mexFunction(double*, int, int);
//extern "C" {int compare(const void *, const void *); }
//extern "C" {int IndexOf(const double [], size_t , double);} 

int main(int argc, char* argv[])

{
  k= atoi(argv[2]);
  size_t count=100000000000;
  double* distance(double*, double*, int,int,int,double*,double*,double*,double*,double*,double*);
  void knn(double*);
  double* Transpose(double*,int*,int*);
  double* MergeSort(double*, int);
  void Merge(double*, int, int);
  void floyd_warshall(double*,int);
  bool connected(int, int, int*);
  int find(int,int*);
  int unite(int, int, int*, double*,int);
  void mds(double*,int*,int,int);
  //void embed(double*, double*, double*,int);
  void embed(double*, double*, int);
  extern void mexFunction(double*,int*, int,int);

  //char *line = malloc(1000000000000);
  //char *line = new char[100000000];
  char *line = (char *)malloc(1000000000000);
  FILE *file= fopen(argv[1], "r");
  double* data=(double*) malloc(100000000*sizeof(double));
//Read number of columns and number of rows
  getline(&line, &count, file);
  int read=-1,cur=0; columnCount=0;
  while(sscanf(line+cur, "%lf%n", &data[columnCount],&read)==1)
    {cur+=read;columnCount=columnCount+1;}
  rowCount=1;
  while(getline(&line, &count, file)!=-1) {rowCount=rowCount+1;}
  rewind(file);
//Reinitialize array using the number of rows and number of columns
  free(data);
  data=(double*) malloc(columnCount*rowCount*sizeof(double));
//Read file and store values
  int i=0;
  while(getline(&line, &count, file)!=-1){
    read=-1,cur=0;
    while(sscanf(line+cur, "%lf%n", &data[i],&read)==1){
            cur+=read;i=i+1;
     }
    }
 free(line);
 fclose(file);
 knn(data);
 //mexFunction(neighborhood_graph,rowCount,k);
 /*int rowCount=20; 
 FILE *fp= fopen("spt", "r");
 double* spt=(double*) malloc(rowCount*rowCount*sizeof(double));
 while(getline(&line, &count, fp)!=-1)
 {
    read=-1,cur=0;
    while(sscanf(line+cur, "%lf%n", &spt[i],&read)==1)
    {cur+=read;i=i+1;}
  }
 fclose(fp); 
printf("%lf\n",spt[i]);

i=0;
FILE *fp2= fopen("vec", "r");
 double* vec=(double*) malloc(rowCount*dims*sizeof(double));
 while(getline(&line, &count, fp2)!=-1)
 {
    read=-1,cur=0;
    while(sscanf(line+cur, "%lf%n", &vec[i],&read)==1)
    {cur+=read;i=i+1;}
  }
 fclose(fp2); 
printf("%lf\n",vec[8]);

i=0;
 FILE *fp1= fopen("val", "r");
 double* val=(double*) malloc(dims*sizeof(double));
 while(getline(&line, &count, fp1)!=-1)
 {
    read=-1,cur=0;
    while(sscanf(line+cur, "%lf%n", &val[i],&read)==1)
    {cur+=read;i=i+1;}
  }
 fclose(fp1); 
printf("%lf\n",val[i-1]);


 
  embed(spt,vec,val,rowCount);*/

}

int compare(const void *p, const void *q) 
{
    double x = *(const double *)p;
    double y = *(const double *)q;
   if (x < y)
        return -1;  // Return -1 if you want ascending, 1 if you want descending order. 
    else if (x > y)
        return 1;   // Return 1 if you want ascending, -1 if you want descending order. 

    return 0;
 }

int compare1(const void *p, const void *q) 
{
    double x = *(const double *)p;
    double y = *(const double *)q;
   if (x > y)
        return -1;  // Return -1 if you want ascending, 1 if you want descending order. 
    else if (x < y)
        return 1;   // Return 1 if you want ascending, -1 if you want descending order. 

    return 0;
 }


 double* Transpose(double* mat, int *rows, int *cols) 
 {
       double* temp;
       int i,j;
       temp = (double*)malloc((*rows)*(*cols)*sizeof(double));
       if(temp){
         for(i=0; i<*rows;i++)
           for(j=0; j<*cols; j++)
              *(temp + (*rows)*j +i) = *(mat+(*cols)*i+j);
        
                 }
        return temp;
    }

 double* distance(double* mat, double* X, int b,int rowCnt, int columnCnt,double*d,double*A,double*B,double*A_B,double*A_A,double*B_B)
  {
      int count=rowCnt; int i1;
      #pragma omp parallel for 
      for(i1=0; i1<columnCnt; i1++) {
         B[i1]=X[(i1*count)+b]*X[(i1*count)+b];
                                        }
        double sum1=0; int i2;
        #pragma omp parallel for reduction(+:sum1)
        for(i2=0; i2<columnCnt; i2++) {
         sum1 += B[i2];
                                          }
        int i3;
        #pragma omp parallel for
        for(i3=0; i3<rowCnt; i3++){
         B_B[i3]=sum1;
                                      }

        
        int i4; int count2=0;
        //#pragma omp parallel for
        for(i4=0; i4<rowCnt; i4++)
        { 
         double sum = 0; 
         int count1 = 0; int j;
         for(j=0; j<columnCnt; j++)
         {
           sum=sum+(mat[j+count2]*X[count1+b]);
           count1+=rowCnt;
          }
         count2+=columnCnt;
         A_B[i4]=sum;
          }
        int i;
       #pragma omp parallel for
       for( i=0; i<rowCnt; i++)
       {
        double e = A_A[i]+B_B[i]-2*A_B[i];
        
        if(e<0)
        {
        e=sqrtl(fabs(e));
        d[i]=0;
          }
        else d[i]=sqrtl(fabs(A_A[i]+B_B[i]-2*A_B[i]));
         }
        return d;
       
   }

  
  
  int IndexOf(const double a[], size_t size, double value) 
  {
  size_t index = 0;
  while ( index < size && a[index] != value ) ++index;
   return ( index == size ? -1 : index );                                                 
   }

 int find(int p,int *parent) 
  {
        while (p != parent[p])
          p = parent[p];
        return p;
    }

 bool connected(int p, int q, int *parent) 
 {
       if(find(p,parent) == find(q,parent))
         return true;
       else
         return false;
   }
    

    int unite(int p, int q, int* parent, double* size,int connected_components) 
   {
        int rootP = find(p,parent);
        int rootQ = find(q,parent);
        if (rootP == rootQ) 
         return connected_components;
        if (size[rootP] < size[rootQ]) {
          parent[rootP] = rootQ;
          size[rootQ] += size[rootP];
                                        }
        else {
         parent[rootQ] = rootP;
         size[rootP] += size[rootQ];
              }
        connected_components--;
        return connected_components;
     }

/**************************Writing residual variance and 3-d embedding files************************/

   //void embed(double* D, double* vec, double* val,int a)
    void embed(double* vec, double* val, int a)
    {
      
      double* R = (double*) malloc(dims*sizeof(double));
       

      for(int di=0; di<dims; di++)
     {
         int di_1 = di+1;
         double* d = (double*) malloc(a*sizeof(double));
         double* A_B = (double*) malloc(a*sizeof(double));
         double* A_A = (double*) malloc(a*sizeof(double));  
         double* B_B = (double*) malloc(a*sizeof(double));
         double* A = (double*) malloc(a*di_1*sizeof(double));
         double* B = (double*) malloc(di_1*sizeof(double));
         double* Y_coords = (double*)malloc(a*di_1*sizeof(double));
       
        int i; 
         #pragma omp parallel for
         for( i=0; i<a; i++)
         {
           double z; //int j;
           for(int j=0; j<di_1; j++)
            {
             z = val[j];
             
             if(z>0)
              *(Y_coords+di_1*i+j)=vec[dims*i+j]*sqrt(z);
             
             else
              *(Y_coords+di_1*i+j)=vec[dims*i+j]*sqrt(fabs(z))*(-1); 
            //printf("%lf\n",vec[a*i+j]); 
              }
           
           }
 
 if(di==2)
 {
   printf("Writing to an Embedding file...\n");
   FILE *fp5 = fopen("Embedding","w");
   for(int i=0; i<a; i++)
   {
    
     for(int j=0; j<(di_1); j++)
     {
       fprintf(fp5,"%lf\t",*(Y_coords+di_1*i+j));
      
      }
   fprintf(fp5,"\n");
   
   }
  
   fclose(fp5);

   }
     
      /*double* X = Transpose(Y_coords, &a, &di_1);
      int i1;
      #pragma omp parallel for
      for(i1=0; i1<a*(di+1); i1++)
         A[i1]=X[i1]*X[i1];

      int i2;
      //#pragma omp parallel for
      for(i2=0; i2<a; i2++) 
     {
       double sum=0;
       int count=0; 
      
       for(int j=0; j<(di+1); j++)
        { 
          sum+=A[count+i2];
          count+=a;
         }
       A_A[i2]= sum;      
      }

         double* dist = (double*) malloc(a*a*sizeof(double));

          if(dist == NULL)
          printf("Insufficient Memory");

         double* temp_dist = (double*) malloc(a*sizeof(double));
         int count1=0;
  
         for(int i=0; i<a; i++)
          {
          temp_dist=distance(Y_coords,X,i,a,(di_1),d,A,B,A_B,A_A,B_B);
           for(int j=0; j<a; j++)
             dist[count1+j]=temp_dist[j];
          count1+=a;
            }
     free(Y_coords); 
// free(temp_dist); free(d); free(A); free(B); free(A_B); free(A_A); free(B_B);

      
      FILE *fp2 = fopen("dist", "w");  
         for(int i=0; i<a; i++){ 
           for(int j=0; j<a; j++)
            fprintf(fp2,"%lf \t",dist[i]);   
          fprintf(fp2,"\n");  
                                  
         }
     fclose(fp2); 
     free(dist);
     system("Rscript corcoef.r");
     
     FILE *fp3 = fopen("corcoef", "r");
     double c; fscanf(fp3, "%lf", &c);
     fclose(fp3);
     printf("No error\n");
     double c; 
     double var1=0; double var2=0; double avg1=0; double avg2=0; double sum=0;
     //#pragma omp parallel for
     for(int i=0;i<a;i++)
      for(int j=0; j<a; j++)
     {
        avg1 += *(D+i*a+j);
        avg2 += *(dist+i*a+j);
      }
      avg1=avg1/a; avg2=avg2/a;

    //#pragma omp parallel for 
    for(int i=0; i<a; i++)
      for(int j=0; j<a; j++)
       {
        var1 += (*(D+i*a+j)-avg1)*(*(D+i*a+j)-avg1);
        var2 += (*(dist+i*a+j)-avg2)*(*(dist+i*a+j)-avg2);
        }

    var1=var1/a; var2=var2/a; double sd1 = sqrt(var1); double sd2 = sqrt(var2);
//printf("%lf\n",var1); printf("%lf\n",var2);
    for(int i=0; i<a; i++)
      for(int j=0; j<a; j++)
       sum += (((*(D+i*a+j)-avg1)*(*(dist+i*a+j)-avg2))/(sd1*sd2));

//printf("%lf\n",sum);
  //c=a-sum;  
c=sum/(a-1);
 c=1-(c*c);

     R[di]=c;*/
     
  }

  
  //free(D); 
  /*free(vec); free(val); 
  system("rm dist"); system("rm D");
  FILE *fp4 = fopen("Residual_Variance","w");
    for(int i=0; i<dims;i++)
      fprintf(fp4,"%lf\t",R[i]);
  fclose(fp4);*/
  
  }

/**********************Run classical multi-dimensional scaling on largest component*******************************/

    void mds(double * shortest_path_tree, int *Y_index, int a, int N) 
    // void mds(double * D, int *Y_index, int a, int N) 
  {
      double* D = (double*) calloc(a*a,sizeof(double)); 
      double* D_sqrd=(double*) calloc(a*a,sizeof(double));   
      
      int i1;
      #pragma omp parallel for collapse(2)
      for( i1=0; i1<a; i1++)
       for(int j=0; j<a; j++)
        {
         *(D+(a)*i1+j)=*(shortest_path_tree+((N)*(Y_index[i1]))+Y_index[j]);
          *(D_sqrd+(a)*i1+j)=(-0.5)*(*(D+(a)*i1+j))*(*(D+(a)*i1+j));
         }
FILE *fp10=fopen("D", "w");
      for (int i = 0; i < a; i++)
       {
        for (int j = 0; j < a; j++)
        {
            
                fprintf (fp10,"%lf\t", *(D+(a)*i+j));
        }
        fprintf(fp10,"\n");
       }
   fclose(fp10);   
      //free(shortest_path_tree);
      
      /*double* D_sqrd=(double*) calloc(a*a,sizeof(double));
      int i2;
      #pragma omp parallel for collapse(2)
      for(i2=0; i2<a; i2++)
        for(int j=0; j<a; j++)
         *(D_sqrd+(a)*i2+j)=(-0.5)*(*(D+(a)*i2+j))*(*(D+(a)*i2+j));

     double* ones_frac =(double*)malloc(a*sizeof(double)); int i3; double* ones=(double*)malloc(a*sizeof(double));
     #pragma omp parallel for
      for(i3=0; i3<a; i3++){
         ones_frac[i3]=(double)1/a;
         ones[i3]=1;
                           }*/
    

     double* D_sqrd_sum=(double*)malloc(a*sizeof(double));
     int i10;
     
     #pragma omp parallel for
      for(i10=0; i10<a; i10++)
      {
        double sum=0;
        for(int j=0; j<a; j++)
          sum+=*(D_sqrd+(a)*i10+j);
        D_sqrd_sum[i10]=sum;
       }
     
     /*double* D_sqrd_sum_frac=(double*)malloc(a*sizeof(double));
     int i11;
     #pragma omp parallel for
      for(i11=0; i11<a; i11++)
      {
        double sum=0;
        for(int j=0; j<a; j++)
          sum+=*(D_sqrd+(a)*i11+j);
        D_sqrd_sum_frac[i11]=sum/a;
       }
     
    

     double * mat1=(double*)malloc(a*a*sizeof(double));
     
     int i5;
     #pragma omp parallel for 
     for( i5=0; i5<a; i5++)
       for(int j=0; j<a; j++)
         *(mat1+(a)*i5+j)=(-0.5)*(*(D_sqrd_sum+i5))*((double)1/a);
      //free(ones_frac);

      

     double * mat2=(double*)malloc(a*a*sizeof(double));
     int i6;
     #pragma omp parallel for collapse(2)
     for( i6=0; i6<a; i6++)
       for(int j=0; j<a; j++)
         *(mat2+(a)*i6+j)=(-0.5)*D_sqrd_sum_frac[j];
      //free(ones);*/
     
     int i9;
     double A=0;
    #pragma omp parallel for reduction(+:A)
     for(i9=0; i9<a; i9++)
        A+=D_sqrd_sum[i9];
     A=(-0.5)*(double)(A/(a*a));
     //free(D_sqrd_sum_frac); free(D_sqrd_sum);

     printf("Performing Multi-dimensional Scaling\n"); 
     
     double * main_mat=(double*)malloc(a*a*sizeof(double));
     int i7;
    
    #pragma omp parallel for collapse(2)
       for(i7=0; i7<a; i7++)
          for(int j=0; j<a; j++)
             *(main_mat+(a)*i7+j) = *(D_sqrd+(a)*i7+j)+((-0.5)*D_sqrd_sum[i7]*((double)1/a))+((-0.5)*D_sqrd_sum[j]*((double)1/a))+A;
     //free(mat1); free(mat2); 
   free(D_sqrd);  free(D_sqrd_sum);
     
    FILE *fp = fopen("main_mat", "w"); 
    if(fp == NULL)
     printf("Can not open file\n");
        
         for(int i8=0; i8<a; i8++)
        {
           for(int j=0; j<a; j++) 
            fprintf(fp,"%f\t",*(main_mat+(a)*i8+j));    
          fprintf(fp,"\n");  
        }
     fclose(fp);

    system("Rscript eigs.r");   
    
    system("rm main_mat");
    
    free(main_mat);
     
    double* val = (double*) calloc(a,sizeof(double));
  
    double * vec = (double*) calloc(a*a,sizeof(double));
    
      
 /*float zmax,emax; int i,j;
   
 
 double* e = (double*) calloc(a,sizeof(double)); 
 for(int l=0; l<dims; l++)
  {
    printf("Calculating Eigen vectors and values...\n");
    double* x = (double*) calloc(a,sizeof(double));
    double* z = (double*) calloc(a,sizeof(double));

    x[l]=1;
    do
   
    {
        for(i=0; i<a; i++)
        {
            z[i]=0;
            for(j=0; j<a; j++)
            {
                z[i]=z[i]+(*(main_mat+(a)*i+j))*x[j];
            }
         }
        zmax=fabs(z[0]);
        
        for(i=0; i<a; i++)
        {
            if((fabs(z[i]))>zmax)
                zmax=fabs(z[i]);
        }
        //#pragma omp parallel for
        for(i=0; i<a; i++)
        {
            z[i]=z[i]/zmax;
        }
       //printf("Calculating Eigen vectors and values...\n");
        //#pragma omp parallel for
        for(i=0; i<a; i++)
        {
            e[i]=0;
            e[i]=fabs((fabs(z[i]))-(fabs(x[i])));
        }
        emax=e[0];
        
        for(i=1; i<a; i++)
        {
            if(e[i]>emax)
                emax=e[i];
        }
        //printf("%f\n",emax);
        for(i=0; i<a; i++)
        {
            x[i]=z[i];
        }
    }
    while(emax>0.001);
  
  val[l]=zmax;
  //#pragma omp parallel for
  for(int i=0; i<a; i++)
   *(vec+(a)*i+l)=z[i];

  free(x); free(z); 
  }
 free(e); free(main_mat);
 
  double * val1 = (double*) calloc(dims,sizeof(double));
  memmove(val1, val,dims*sizeof(double));
  qsort(val, dims, sizeof(val[0]), compare1);
   
 FILE *fp1 = fopen("val","w");  
 for(int i=0; i<dims; i++)
   fprintf(fp1,"%lf\t",*(val+i));


int * val_ind = (int*) calloc(dims,sizeof(int));
for(int i=0; i<dims; i++)
{
  val_ind[i] = IndexOf(val1,dims,val[i]);
  printf("%lf\n",val[i]);
  //printf("%d\n",val_ind[i]);
 }
double * vec1 = (double*) calloc(a*dims,sizeof(double));*/


/*for(int i=0; i<a; i++)
  for(int j=0; j<dims; j++)
   vec1[dims*i+j]=vec[a*i+j];

printf("No error\n");

FILE *fp2 = fopen("vec1", "w"); 
    if(fp2 == NULL)
     printf("Can not open file\n");
        
         for(int i=0; i=a; i++)
         {
           for(int j=0; j<dims; j++) 
            fprintf(fp2,"%lf\t",*(vec1+(dims)*i+j));    
          fprintf(fp2,"\n");  
          }
     fclose(fp2); 

 FILE *fp3 = fopen("vec", "w"); 
    if(fp3 == NULL)
     printf("Can not open file\n");
        
         for(int i=0; i<dims; i++)
         {
           for(int j=0; j<a; j++) 
           {
            int l = val_ind[i];
            fprintf(fp3,"%lf\n",vec[(a)*j+l]); 
            }   
          fprintf(fp3,"\t");  
          }
     fclose(fp3);  
FILE *fp2 = fopen("vec1", "w");    
    for(int i=0; i<a; i++)
     {
        
        for(int j=0; j<dims; j++) 
         {
           int l = val_ind[j];
           vec1[dims*i+j] = vec[a*i+l];
           fprintf(fp2,"%lf\t",vec1[(dims)*i+j]);  
          }  
        fprintf(fp2,"\n");
         
      }
fclose(fp2);
   free(vec);     

   //printf("No error\n");*/

 FILE *fp1=fopen("val", "r");
   size_t count1=100000000;
   char *line1 = malloc(1000000000);
   int i=0;
    while(getline(&line1, &count1, fp1)!=-1){
    int read=-1,cur=0;
    while(sscanf(line1+cur, "%lf%n", &val[i],&read)==1)
    {cur+=read;i=i+1;}
                                            }
    free(line1);
    fclose(fp1);
    system("rm val");
FILE *fp2=fopen("vec", "r");
   size_t count2=100000000;
   char *line2 = malloc(1000000000);
   //int i=0;
    i=0;
    while(getline(&line2, &count2, fp2)!=-1){
    int read1=-1,cur1=0;
    while(sscanf(line2+cur1, "%lf%n", &vec[i],&read1)==1)
    {cur1+=read1;i=i+1;}
                                            }
    free(line2);
    fclose(fp2);
    system("rm vec");
    free(D);
   //embed(D,vec,val,a);
   embed(vec,val,a);

  }

/**************************************Finding shortest pair-wise paths between vertices***************************************/ 
   
  void floyd_warshall(double * shortest_path_tree, int N)
  {
     
    /* double* shortest_path_tree=(double*)malloc(N*(N)*sizeof(double));
     #pragma omp parallel for collapse(2)
     for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
         *(shortest_path_tree+(N)*i+j) = *(neighborhood_graph+(N)*i+j);
    free(neighborhood_graph);
     #pragma omp parallel for collapse(2)
     for(int i=0; i<N; i++)
         for(int j=0; j<N; j++){
          if(*(shortest_path_tree+(N)*i+j)!=0)
            *(shortest_path_tree+(N)*j+i) = *(shortest_path_tree+(N)*i+j);
          if(i!=j && *(shortest_path_tree+(N)*i+j)==0)
           *(shortest_path_tree+(N)*i+j) = INFINITY;            
                                      }*/

     printf("Computing shortest paths...\n");
     
   /*FILE *fp=fopen("graph", "w");
      for (int i = 0; i < rowCount/1.5; i++)
       {
        for (int j = 0; j < rowCount/1.5; j++)
        {
            if (*(shortest_path_tree+(rowCount)*i+j) == 0)
              fprintf(fp,"%7s\t", "NA");
            else
                fprintf (fp,"%lf\t", *(shortest_path_tree+(rowCount)*i+j));
        }
        fprintf(fp,"\n");
       }
   fclose(fp);*/
  
   system("Rscript Floyd.r"); 
   system("rm ng");
   system("sed -i -- 's/NA/INF/g' *spt* ");
   FILE *fp1=fopen("spt", "r");
   size_t count1=100000000;
   char *line1 = malloc(1000000000);
   int i=0;
    while(getline(&line1, &count1, fp1)!=-1){
    int read=-1,cur=0;
    while(sscanf(line1+cur, "%lf%n", &shortest_path_tree[i],&read)==1)
    {cur+=read;i=i+1;}
                                            }
    free(line1);
    fclose(fp1);
    system("rm spt");

  
   /*  for(int a=0; a<N; a++)
     {
     if(a%100==0) printf("Iteration no. %d for Floyd\n", a);
       #pragma omp parallel for
       for(int b=0; b<N; b++)
         //#pragma omp parallel
         for(int c=0; c<N; c++)
             if (*(shortest_path_tree+(N)*b+a) + *(shortest_path_tree+(N)*a+c) < *(shortest_path_tree+(N)*b+c)) 
                *(shortest_path_tree+(N)*b+c) = *(shortest_path_tree+(N)*b+a) + *(shortest_path_tree+(N)*a+c);
        }    


  FILE *fp2=fopen("spt", "w");
      for (int i = 0; i < rowCount; i++)
       {
        for (int j = 0; j < rowCount; j++)
        {
            
                fprintf (fp2,"%lf\t", *(shortest_path_tree+(rowCount)*i+j));
        }
        fprintf(fp2,"\n");
       }
   fclose(fp2);*/
   printf("Finding number of connected components...\n");
    int* vertices=(int*)malloc((N)*N*sizeof(int));
    int l=0;
     //#pragma omp parallel
     for(int i=0; i<N; i++)
       //#pragma omp parallel
       for(int j=i+1; j<N; j++)
        if((*(shortest_path_tree+(N)*i+j)!=INFINITY))
         {
         vertices[l]=i;  
         l++;
         vertices[l]=j;
         l++;
          }
         
    
    int* edges=(int*)malloc(l*sizeof(int));
    #pragma omp parallel for
    for(int i=0; i<l; i++){
        edges[i]=vertices[i];
                             }
     free(vertices);
    int* parent=(int*)malloc(N*sizeof(int));
    double* size=(double*)malloc(N*sizeof(double));                                                         
    #pragma omp parallel for
    for(int i = 0; i <N; i++) {
         parent[i] = i;
         size[i] = 1;
                                      }

     int connected_components = N;
     for(int i=0,j=1; i<l-1 && j<l; i+=2,j+=2)
      {
       if(connected(edges[i],edges[j],parent)==true)
         continue;
        connected_components=unite(edges[i],edges[j],parent,size,connected_components);
       }
     free(edges);     

    double* size1=(double*)malloc(N*sizeof(double));
    memmove(size1, size,(N)*sizeof(double));
    qsort(size1, N, sizeof(size1[0]), compare);
    int size_lrgst_comp=size1[N-1];
    int root_lrgst_comp=parent[IndexOf((double*)size,N,(double)size_lrgst_comp)];     
    free(size1);
    free(size);

    int *Y_index=(int*)malloc(size_lrgst_comp*sizeof(int));
     
     for(int i=0,j=0; i<N; i++){
      if(j==size_lrgst_comp) break;
      if(parent[i]==root_lrgst_comp){
        Y_index[j]=i;
        j++;
                                     }
                                        }
    FILE *fp11=fopen("Indices", "w");
    for(int i=0; i<size_lrgst_comp;i++)
     fprintf (fp11,"%d\t", *(Y_index+i));
   fclose(fp11); 

     printf("Number of connected components:%d\n",connected_components);
     printf("Size of the largest component:%d\n",size_lrgst_comp);
     mds(shortest_path_tree,Y_index,size_lrgst_comp,N);
    } 



/*************************************Find k nearest neighbours and build the neighbourhood graph*************************************/
  
   void knn(double* data) 
   {
       
       int* Di=(int*)malloc(rowCount*(k)*sizeof(double));
       int* Dj=(int*)malloc(rowCount*(k)*sizeof(double));
       double* Ds=(double*)malloc(rowCount*(k)*sizeof(double));
       //double* dist=(double*)malloc(rowCount*sizeof(double));
       //double* c=(double*)malloc(rowCount*sizeof(double));
       //int* b=(int*)malloc(rowCount*sizeof(double));
     
      double* d = (double*) malloc(rowCount*sizeof(double));
      double* X = Transpose(data, &rowCount, &columnCount);
      double* A = (double*) malloc(rowCount*columnCount*sizeof(double));
      double* B = (double*) malloc(columnCount*sizeof(double));
      double* A_B = (double*) malloc(rowCount*sizeof(double));
      double* A_A = (double*) malloc(rowCount*sizeof(double));  
      double* B_B = (double*) malloc(rowCount*sizeof(double));
 
      int i1;
      #pragma omp parallel for
      for( i1=0; i1<rowCount*columnCount; i1++)
         A[i1]=X[i1]*X[i1];

      int i2;
      #pragma omp parallel for
      for(i2=0; i2<rowCount; i2++) 
     {
       double sum=0;
       int count=0;
       for(int j=0; j<columnCount; j++)
       { 
         sum+=A[count+i2];
         count+=rowCount;
        }
       A_A[i2]= sum;      
      }
      
       int i,j,counter; double* dist; double* c; int* b;
       #pragma omp parallel private(counter,i,j,dist,b,c)
    {
        dist=(double*)calloc(rowCount,sizeof(double));
        c=(double*)calloc(rowCount,sizeof(double));
        b = (int*) calloc(rowCount, sizeof(int)); counter=k;
     #pragma omp for
       for(i=0; i<rowCount; i++) 
       {
          
          if(i%100==0) printf("Iteration no. %d for knn\n", i);
          #pragma omp critical
        {
          dist=distance(data,X,i,rowCount,columnCount,d,A,B,A_B,A_A,B_B);
          memmove(c, dist,rowCount*sizeof(double));
          qsort(c, rowCount, sizeof(c[0]), compare);
        }
          for(j=0;j<rowCount;j++)
          {
              if(j!=0 && c[j]==c[j-1])
                b[j]=IndexOf(dist,rowCount,c[j])+1;
              else 
                b[j]=IndexOf(dist,rowCount,c[j]);
           for(int l=0; l<k; l++)
           {
             Di[l+(i*counter)]=i;
             if(l==0 && b[l]!=i) {Dj[l+(i*counter)]=b[l+1];  Ds[l+(i*counter)]=c[l+1];}
             Dj[l+(i*counter)]=b[l];  Ds[l+(i*counter)]=c[l];
            }
          }
       
        
        }
        
      }
           //free(b); free(c); free(dist);

     
     int N= rowCount;
     double * neighborhood_graph=(double*)calloc(N*N,sizeof(double));
    int cnt=0;                    
    //#pragma omp parallel                                  
    for(int i=0; i<N;i++)
    {
      for(int j=0; j<k; j++)
         *(neighborhood_graph+(N)*i+(Dj[j+cnt]))=Ds[j+cnt];
      cnt+=k;
     }

                                 
    free(Di); free(Dj); free(Ds);//free(d);free(A);free(B);free(A_B);free(A_A);free(B_B);
    free(data);
/*FILE *fp;
      fp = fopen("knn_out.txt", "w"); //Here, every time when the file is opened its contents are erased. 
       for(int i=0; i<rowCount*k; i++)
          fprintf(fp, "%d %d %f\n", Di[i], Dj[i], Ds[i]);
       fprintf(fp, "\n");
      fclose(fp);*/

  #pragma omp parallel for collapse(2)
     for(int i=0; i<N; i++)
         for(int j=0; j<N; j++)
         {
          if(*(neighborhood_graph+(N)*i+j)!=0)
            *(neighborhood_graph+(N)*j+i) = *(neighborhood_graph+(N)*i+j);
          if(i!=j && *(neighborhood_graph+(N)*i+j)==0)
           *(neighborhood_graph+(N)*i+j) = INFINITY;            
          }

 FILE *fp1=fopen("ng", "w");
      for (int i = 0; i < N; i++)
       {
        for (int j = 0; j < N; j++)
        {
            
                fprintf (fp1,"%lf\t", *(neighborhood_graph+(N)*i+j));
        }
        fprintf(fp1,"\n");
       }
   fclose(fp1);   

/*mexFunction(neighborhood_graph,N,k);

 int *Y_index=(int*)malloc(N*sizeof(int));
 for(int i=0; i<N; i++)
   Y_index[i] = i;
    //return neighborhood_graph;*/ 
     floyd_warshall(neighborhood_graph,N); 
  //mds(neighborhood_graph,Y_index,N,N); 
  
   }


        








