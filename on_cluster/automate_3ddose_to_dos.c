#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

#define MAXPAR 10
#define FALSE 0
#define TRUE !FALSE

/* program to read .3ddose file from DOSXYZ and convert to .dos format written by VMC++ 

BUT unlike its namesake (minue 'automate'), this is gonna be automated st it'll run it on 200 files at a time (in much the same way as automateDosxyznrcSubmit.c)
*/

int fileExists(const char *path);

//char name_3ddose[100];
FILE *f3ddose,*dosFile;
int nx,ny,nz,i,j,k,addr;
float *xbound,*ybound,*zbound,*dose,max_dose,min_dose,*error,tmp,maxdose,mindose,thresh,meanerror;
int nreg;

main (int argc, char** argv)
{
  
  int numBeamlets = 323; //fixed for now (afaik this is the overall max?)
  int maxSubmit = 200;  //maximum number of jobs to submit to the cluster at a time


  //enter arguments when running like:
  // ./auto2ddt2dos 3ddose_file_base dos_file_base totnumbeams startingBeam startingBeamlet startingStf
  //startingStf is number of beamlet as given by the stf struct in matRad. Basically numbering w/out gaps
  //just either give it 1 if starting from scratch, or give it the value specified at the end of the last run

  if(argc != 7) {
      printf("You haven't given me the right number of  arguments. When running me, write: ./auto2ddt2dos 3ddose_file_base dos_file_base totnumbeams startingBeam startingBeamlet\n");
      return 0;
  }
  
  if(atoi(argv[4]) <= 0 || atoi(argv[4]) > atoi(argv[3])) {
      printf("the starting beam number you entered, n = %d, doesn't fall within the range of beams in use (total number of beams = %d)\n",atoi(argv[4]),atoi(argv[3]));
      return 0;
  }
  

  printf("total number of beams given = %d\n",atoi(argv[3]));
  printf("total number of beamlets = %d\n",numBeamlets);
  printf("starting beam = %d\n",atoi(argv[4]));
  printf("starting beamlet = %d\n",atoi(argv[5]));
  printf("starting beamlet numbered by stf struct in matRad = %d\n",atoi(argv[6]));
  
  
  //so for me files are named like filebaseBeamNBeamletM.3ddose
  //so take sprintf stuff from my other automated file
  //but overall I think I'll format this exactly like the autoDos file
  //including how it's to be run. Will make it so much simpler
  //I can do this from home so that's what I'll be doing now bye

  printf("Convert a 3ddose file to .dos\n");
  /*printf("Enter file name of 3ddose file (include extension):\n");
  scanf("%s",name_3ddose);*/
  
  int n, l;
  int m = atoi(argv[6]);
  int numSubmitted = 0;
  
  for(n = atoi(argv[4]); n <= atoi(argv[3]); n++) {
    for(l = atoi(argv[5]); l <= numBeamlets; l++) {
      char name_3ddose[256];
      sprintf(name_3ddose,"%sBeam%dBeamlet%d.3ddose",argv[1],n,l);
      printf("uhh %s\n",name_3ddose);
      
      //so that m, the beamlet as numbered in stf, resets when we go to the next beam
      if(l == 1) {
        m = 1;
      }
      
      //do the command if the file exists; otherwise go on to the next one
      if(fileExists(name_3ddose)==1) {
        numSubmitted = numSubmitted + 1;
        f3ddose=fopen(name_3ddose,"r");
        printf("File opened\n");
        fscanf(f3ddose,"%d %d %d\n",&nx,&ny,&nz);
        printf("The voxel dimensions in this file is: %d %d %d\n",nx,ny,nz);
      
        /*get all voxel boundaries in x-,y- and z-directions */
        xbound = (float *) malloc(sizeof(float)*(nx+1));
        ybound = (float *) malloc(sizeof(float)*(ny+1));
        zbound = (float *) malloc(sizeof(float)*(nz+1));
      
        for (i = 0; i <= nx; i++)
          {
            fscanf(f3ddose,"%f",&xbound[i]);
            //printf("x bound: %f\n",xbound[i]);
          }
      
        for (j = 0; j <= ny; j++)
          {
            fscanf(f3ddose,"%f",&ybound[i]);
            //printf("y bound: %f\n",ybound[i]);
          }
      
        for (k = 0; k <= nz; k++)
          {
            fscanf(f3ddose,"%f",&zbound[i]);
            //printf("z bound: %f\n",zbound[i]);
          }
      
      
        /* read dose */
        dose = (float *) malloc(nx*ny*nz*sizeof(float));  
        max_dose = -1.0E+15;
        min_dose = 1.0E+15;
      
        for (k=0;k<nz;k++)
          for (j=0;j<ny;j++)
            for (i=0;i<nx;i++)
      	{ 
      	  addr = i + j*nx + k*nx*ny;
      	  fscanf(f3ddose,"%f",&dose[addr]);
      	  //dose[addr]=(1.0e+10)*dose[addr];
      	  if (dose[addr]>=max_dose)
      	    max_dose=dose[addr];
      	  if (dose[addr]<=min_dose)
      	    min_dose=dose[addr];
      	}
        printf("Max dose %e Min dose %e\n", max_dose,min_dose);
      
        /* read error */
        error = (float *) malloc(nx*ny*nz*sizeof(float));
        for (k=0;k<nz;k++)
          for (j=0;j<ny;j++)
            for (i=0;i<nx;i++)
              {
                addr = i + j*nx + k*nx*ny;
                fscanf(f3ddose,"%f",&error[addr]);
      	  error[addr]=error[addr]*dose[addr];
      	}
      
        fclose(f3ddose);
      
        /* now write in .dos format */
        char name_dos[256];
        //filebaseMBeamNBeamletL.dos
        sprintf(name_dos,"%s%dBeam%dBeamlet%d.dos",argv[2],m,n,l);
        printf("writing converted dose to %s\n",name_dos);
        dosFile = fopen(name_dos,"wb");
        
        nreg = nx*ny*nz;
      
        fwrite(&nreg,sizeof(int),1,dosFile);
        printf("nreg %d\n",nreg);
      
        fwrite(dose,sizeof(float),nreg,dosFile);
        printf("dose read\n");
        fwrite(error,sizeof(float),nreg,dosFile);
        fclose(dosFile);
  
        m = m + 1;
      }
      
      if(numSubmitted == maxSubmit){
        printf("Last submitted for this run = beam %d beamlet %d.\n",n,l);
        if(l<numBeamlets) {  //so if we can keep in the same beam but go on to the enxt beamelt next time:
          printf("To continue this batch, tell me to start from beam %d, beamlet %d stf-beamlet %d.\n",n,l+1,m);
        } else {  //if that was the last beamlet in the beam:
          if(n == atoi(argv[3])) {  //if it's also the last beam -> we've reached the end:
            printf("All the jobs for all the beams and all the beamlets have been submitted!\nYay!\nDo your happy dance!\n");
          } else {  //if it was the last beamlet in that beam but more beam(s) await:
            printf("To continue this batch, tell me to start from beam %d, beamlet %d.\n",n+1,1);
          }
        }
        
        exit(0);
      } else if(n==atoi(argv[3]) && l == numBeamlets) {
        printf("Looks like we've reached the end of all the beams and beamlets!\n");
      }
    }
  }

};



//this function checks if the file exists:
int fileExists(const char *path)
{
    // Check for file existence
    if (access(path, F_OK) == -1)
        return 0;

    return 1;
}
