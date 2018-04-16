#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

double *box,*d;

double dist (int j, int n1, int i, int n2, double **x, double **y, double **z, double *d){
            
            double r;

            d[0]=x[j][n1]-x[i][n2];
            d[1]=y[j][n1]-y[i][n2];
            d[2]=z[j][n1]-z[i][n2];

            d[0] = d[0] - box[0]*( (int)(2.0*d[0]/box[0]) - (int)(d[0]/box[0]) );
            d[1] = d[1] - box[1]*( (int)(2.0*d[1]/box[1]) - (int)(d[1]/box[1]) );
            d[2] = d[2] - box[2]*( (int)(2.0*d[2]/box[2]) - (int)(d[2]/box[2]) );

            r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
            return r;
}



int main(int nvar,char **cvar){


  /*declaring variables and asking for inputs */
//---------------------------------------------------------------------------------

    FILE *rptr,*wptr1,*wptr2,*wptr3,*wptr4,*wptr5;
    int natom,nwater,ntotal,nmol,nframes,nskip,tstep,wn=0;
    int i,j,k,l,m,r,a,*bufftime,nregion;                        /* count1 is region1 and so on */
    int **count;
    double **x, **y, **z, **dp, **dplg, xx;
    double tdiff,cutoff1=2.6,cutoff2=5;          /* d1[0] is O-H for frame j  and d2[0] is for frame i */
    double *d1,*d2;
    char buf1[1000],buf2[1000], buf3[1000];
    int **region, *regionprev,*ifound,bufflimit,nmono=40;
    double *regiontot;

/*   give the name of the trajectory file after ./a.out    */
//----------------------------------------------------------------------------------

    if(nvar!=5){
      printf("  Usage: %s trajectoryfile %s abcboxsize \n",cvar[0],cvar[1]);
      exit(1);
    }
    rptr=fopen(cvar[1],"r");

    if(rptr==NULL){
        printf("  ERROR: could not open trajectory file %s\n",cvar[1]);
        exit(1);
    }

    box = (double *)malloc((3)*sizeof(double)); 

    sscanf(cvar[2], "%lf", &box[0]);
    sscanf(cvar[3], "%lf", &box[1]);
    sscanf(cvar[4], "%lf", &box[2]);
    printf(" Box size is %lf %lf %lf \n ",box[0],box[1],box[2]);

    printf("Enter the number of atoms in molecule, number of molecules, total number of water molecules present and number of regions \n");
    scanf("%d %d %d %d",&natom,&nmol,&nwater,&nregion);
   
    printf("Enter the total number of frames and skip value for reference frame and time betweem two frames \n");
    scanf("%d %d %lf",&nframes,&nskip,&tdiff);
    ntotal=natom*nmol+3*nwater;
    tstep=(int)(20.0/tdiff);
    bufflimit=(int)(0.5/tdiff);
    printf("Total number of atoms present in the system is %d bufflimit is %d \n",ntotal,bufflimit);

    wptr1=fopen("rotOside.dat","w");
    wptr2=fopen("rotC.dat","w");
    wptr3=fopen("rotB.dat","w");
    wptr4=fopen("regiontot.dat","w");
    wptr5=fopen("rotOmiddle.dat","w");



 /* allocate memory */
//----------------------------------------------------------------------------------
    x = (double **)malloc(nframes*sizeof(double *));
    for(i=0;i<nframes;i++) x[i] = (double *)malloc(ntotal*sizeof(double));
  
    y = (double **)malloc(nframes*sizeof(double *));
    for(i=0;i<nframes;i++) y[i] = (double *)malloc(ntotal*sizeof(double));
  
    z = (double **)malloc(nframes*sizeof(double *));
    for(i=0;i<nframes;i++) z[i] = (double *)malloc(ntotal*sizeof(double));


    d   = (double *)malloc((3)*sizeof(double)); 
    d1  = (double *)malloc((3)*sizeof(double)); 
    d2  = (double *)malloc((3)*sizeof(double)); 

    dplg = (double **)malloc(nregion*sizeof(double *));
    for(i=0;i<nregion;i++) dplg[i] = (double *)malloc(tstep*sizeof(double)); 

    dp   = (double **)malloc((nregion)*sizeof(double *));
    for(i=0;i<nregion;i++) dp[i] = (double *)malloc(tstep*sizeof(double)); 
    
    count = (int **)malloc(nregion*sizeof(int *));
    for(i=0;i<nregion;i++) count[i] = (int *)malloc(tstep*sizeof(int)); 

    region = (int **)malloc(nframes*sizeof(int *));
    for(i=0;i<nframes;i++) region[i] = (int *)malloc(ntotal*sizeof(int));

    regionprev = (int *)malloc((ntotal)*sizeof(int));
    for(i=0;i<ntotal;i++) regionprev[i]=-1;

    ifound = (int *)malloc((ntotal)*sizeof(int));
    for(i=0;i<ntotal;i++) ifound[i]=-1;

    regiontot = (double *)malloc((nregion)*sizeof(double));
    for(i=0;i<nregion;i++) regiontot[i]=0;

    bufftime = (int *)malloc((ntotal)*sizeof(int));
    for(i=0;i<ntotal;i++) bufftime[i]=0;


for(i=0;i<nregion;i++){ for(j=0;j<tstep;j++) count[i][j]=0;}
for(i=0;i<nregion;i++){ for(j=0;j<tstep;j++) dp[i][j]=0;}
for(i=0;i<nregion;i++){ for(j=0;j<tstep;j++) dplg[i][j]=0;}

printf("MEMORY ALLOCATION  FINISHED");


 /* start reading frames */
 //-----------------------------------------------------------------------------------
    for(j=0;j<nframes;j++){
       fscanf(rptr,"%d",&ntotal);

//first fgets is to read the end line which is left at the end of first line by scanf

       fgets(buf1,1000,rptr);
       fgets(buf1,1000,rptr);


       for(i=0;i<ntotal;i++){
          fscanf(rptr,"%s %lf %lf %lf",buf2,&x[j][i],&y[j][i],&z[j][i]);
       }

    }

printf("FRAME READING  FINISHED");


/* Deciding the region for each water molecule */
//---------------------------------------------------------------------------------------

 for(j=0;j<nframes;j++){
    for(i=0;i<ntotal;i++) ifound[i]=-1;
       for(k=natom*nmol;k<ntotal;k=k+3){
          for(l=k+1;l<k+3;l++){

             for(a=0;a<nmol*natom;a=a+natom){
                if(dist(j,l,j,a+1,x,y,z,d)<cutoff1){ region[j][l]=3; ifound[l]=1; break; }
             }
             if(ifound[l]!=1){
               for(a=0;a<nmol*natom;a=a+natom){
                  if(dist(j,l,j,a+10,x,y,z,d)<cutoff1 || dist(j,l,j,a+13,x,y,z,d)<cutoff1){   region[j][l]=0; ifound[l]=1; break; }  
               }
             } 
             if(ifound[l]!=1){
                for(a=0;a<nmol*natom;a=a+natom){
                   if(dist(j,k,j,a+0,x,y,z,d)<cutoff2 || dist(j,k,j,a+4,x,y,z,d)<cutoff2 || dist(j,k,j,a+5,x,y,z,d)<cutoff2 || dist(j,k,j,a+6,x,y,z,d)<cutoff2 || dist(j,k,j,a+15,x,y,z,d)<cutoff2 || dist(j,k,j,a+16,x,y,z,d)<cutoff2){   region[j][l]=1; ifound[l]=1; break; }  
                }
             } 
             if(ifound[l]!=1) region[j][l]=2;

/* Applying bufflimit to ignore small  fluctuations     */ 

             if(j>0){
               if(region[j][l]!= region[j-1][l]){
                 if(bufftime[l]<bufflimit && region[j][l]==regionprev[l]){
                   for(i=0;i<bufftime[l];i++){region[j-i][l]=regionprev[l];}
                 }
                 else{
                     bufftime[l]=0;
                     regionprev[l]=region[j-1][l];
                 }

               }
             }
             bufftime[l]++;
          }
       }
 }

printf("REGION CALCULATION FINISHED");


/* Calculating dipoles */
//---------------------------------------------------------------------------------------
 for(j=0;j<nframes;j=j+nskip){
    for(k=natom*nmol;k<ntotal;k=k+3){
       for(l=k+1;l<k+3;l++){ 

          xx = dist(j,k,j,l,x,y,z,d);
          for(m=0;m<3;m++) d1[m] = d[m];

          for(i=j;i<nframes && i<j+tstep;i++){

             xx = dist(i,k,i,l,x,y,z,d);
             for(m=0;m<3;m++) d2[m] = d[m];

             xx = (d1[0]*d2[0]+d1[1]*d2[1]+d1[2]*d2[2])/(sqrt(d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2])*sqrt(d2[0]*d2[0]+d2[1]*d2[1]+d2[2]*d2[2]));
    	     r=region[i][l];
   
             dp[r][i-j] += xx;
             dplg[r][i-j] += 0.5*(3.0 * xx*xx - 1.0);
             count[r][i-j]++;
          }
       }
    }
 }



printf("DIPOLE CALCULATION FINISHED \n");


/* printing the result */
//-------------------------------------------------------------------------------------------
for(i=0;i<tstep;i++)   fprintf(wptr1,"%lf \t %lf \t %lf \n",i*tdiff,(dp[0][i]/count[0][i]),(dplg[0][i]/count[0][i]));
for(i=0;i<tstep;i++)   fprintf(wptr2,"%lf \t %lf \t %lf \n",i*tdiff,(dp[1][i]/count[1][i]),(dplg[1][i]/count[1][i]));
for(i=0;i<tstep;i++)   fprintf(wptr3,"%lf \t %lf \t %lf \n",i*tdiff,(dp[2][i]/count[2][i]),(dplg[2][i]/count[2][i]));
for(i=0;i<tstep;i++)   fprintf(wptr5,"%lf \t %lf \t %lf \n",i*tdiff,(dp[3][i]/count[3][i]),(dplg[3][i]/count[3][i]));


   for(i=0;i<nframes;i++){
      for(j=natom*nmol;j<ntotal;j=j+3){
         regiontot[region[i][j]]++;
      }
   }

fprintf(wptr4,"regionOxygenmiddle is %lf regionOxygenside is %lf regionCarbon is %lf regionBulk is %lf \n ",regiontot[3]/(nframes),regiontot[0]/(nframes),regiontot[1]/(nframes),regiontot[2]/(nframes));




fclose(rptr);
fclose(wptr1);
fclose(wptr2);
fclose(wptr3);
fclose(wptr4);
fclose(wptr5);
return 0;

}




