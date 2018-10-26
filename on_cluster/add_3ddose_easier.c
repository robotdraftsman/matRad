#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_NUMBER_OF_FIELDS 20

// same as similarly named file but reads in all phsp files at once

FILE *f3ddose;
int nx,ny,nz,i,j,k,nfield, ifield, addr;
float *temp[MAX_NUMBER_OF_FIELDS], *temps[MAX_NUMBER_OF_FIELDS], *weight,*sumfield, *sumerror, xbound[300],ybound[300],zbound[200];
char name_3ddose[100], num[10], name[100], field_filename[100];

float dx,dy,dz, output[MAX_NUMBER_OF_FIELDS],output_mc[MAX_NUMBER_OF_FIELDS],
  w_factor[MAX_NUMBER_OF_FIELDS],min_dose,max_dose,MU[MAX_NUMBER_OF_FIELDS],
  value,max_error;
int read_error,write_error,max_addr;

int main()

//make this work for field filenames being all identical except ending number, e.g. field1.3ddose, field2.3ddose...

{
  printf("Program to add n fields from DOSXYZ and produce one .3ddose file\n");

  printf("Enter the number of dose fields to be processed: ");
  scanf("%d", &nfield);

  weight = (float *) malloc(nfield*sizeof(float));

  //get base for file name (e.g. smolGroup)
  printf("Please enter the filename base (no extensions):");
  scanf("%s",field_filename);
  
  for (ifield=0;ifield<nfield;ifield++)
    {
      
      //printf("Enter field weight: ");
      //scanf("%f",&weight[ifield]); //instead just weight them all the same = 1/N
      weight[ifield] = (float)1/nfield;
      printf("The weight is %f\n",weight[ifield]);
      
      sprintf(name_3ddose, "%s%d.3ddose",field_filename,ifield+1);
      printf("File name of 3ddose file: %s\n",name_3ddose);
            while ((f3ddose=fopen(name_3ddose,"r"))==NULL)
	{
	  printf("File could not be opened\n");
	  printf("Reenter file name 3ddose file: ");
	  scanf("%s",field_filename);
	  sprintf(name_3ddose, "%s.d3d",field_filename);
	}
      
      /* read .3ddose file */
      fscanf(f3ddose,"%d %d %d\n",&nx,&ny,&nz);
      printf("nx = %d ny = %d nz = %d\n",nx,ny,nz);
      for (i=0;i<nx+1;i++) fscanf(f3ddose,"%f",&xbound[i]);
      for (i=0;i<ny+1;i++) fscanf(f3ddose,"%f",&ybound[i]);
      for (i=0;i<nz+1;i++) fscanf(f3ddose,"%f",&zbound[i]);

      temp[ifield] = (float *) malloc(nx*ny*nz*sizeof(float));
      temps[ifield] = (float *) malloc(nx*ny*nz*sizeof(float));
      for (i=0;i<nx*ny*nz;i++) 
	fscanf(f3ddose,"%e",(temp[ifield] +i));
      for (i=0;i<nx*ny*nz;i++) 
	fscanf(f3ddose,"%f",(temps[ifield] + i));
      fclose(f3ddose);
      printf("3ddose file for field nr %d read.\n",ifield);

      max_dose = -1.0E+15;
      min_dose = 1.0E+15;
      for (k=0;k<nz;k++)
	for (j=0;j<ny;j++)
	  for (i=0;i<nx;i++)
	    { 
	      addr = i + j*nx + k*nx*ny;
	      value = *(temp[ifield] + addr);
	      if (value>=max_dose)
		{
		  max_dose=value;
		  max_error=*(temps[ifield] + addr);
		}
	      if (value<=min_dose)
		min_dose=value;
	    }
      printf("Max dose %e +/- %f for field %d\n", 
	     max_dose,max_error,ifield+1);
      printf("***************************************\n");
    }
   
  sumfield = (float *) malloc(nx*ny*nz*sizeof(float));
  sumerror = (float *) malloc(nx*ny*nz*sizeof(float));
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      for (k=0;k<nz;k++)
	{ 
	  addr = i + j*nx + k*nx*ny;
	  *(sumfield + addr) = 0.0;
	  *(sumerror + addr) = 0.0;
	}

  for (ifield=0;ifield<nfield;ifield++)
    {
      for (k=0;k<nz;k++)
	for (j=0;j<ny;j++)
	  for (i=0;i<nx;i++)
	    { 
	      addr = i + j*nx + k*nx*ny;
	      *(sumfield + addr) += (*(temp[ifield] + addr))*(weight[ifield]);
	      *(sumerror + addr) += ((*(temps[ifield] + addr))*(*(temp[ifield] + addr)))*((*(temps[ifield] + addr))*(*(temp[ifield] + addr)))/(nfield*nfield);
	    }

      free(temp[ifield]);
    }

  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      for (k=0;k<nz;k++)
	{ 
	  addr = i + j*nx + k*nx*ny;

	  if((*(sumfield + addr))==0) (*(sumerror + addr))=0.999999881;
	  else {
	  *(sumerror + addr) = (sqrt(*(sumerror + addr)))/((*(sumfield+addr)));
	  }
	}  

  //find new maximum dose
  max_dose = -1.0E+15;
  min_dose = 1.0E+15;

  for (k=0;k<nz;k++)
    for (j=0;j<ny;j++)
      for (i=0;i<nx;i++)
	{ 
	  addr = i + j*nx + k*nx*ny;
	  value = *(sumfield + addr);
	  if (value>=max_dose)
	    {
	      max_dose=value;
	      max_error=*(sumerror + addr);
	    }
	  if (value<=min_dose)
	    min_dose=value;
	}
  printf("summed dose:\n");
  printf("Max dose %e +/- %f for field %d\n", 
	 max_dose,max_error,ifield+1);
  printf("***************************************\n");

  printf("Enter filename for sum field (no extension): ");
  scanf("%s",field_filename);
  sprintf(name_3ddose, "%s.3ddose",field_filename);
  f3ddose = fopen(name_3ddose,"w");
  fprintf(f3ddose,"%d %d %d\n",nx,ny,nz);
  for (i=0;i<nx+1;i++) fprintf(f3ddose,"%7.3f ",xbound[i]);
  fprintf(f3ddose,"\n");
  for (i=0;i<ny+1;i++) fprintf(f3ddose,"%7.3f ",ybound[i]);
  fprintf(f3ddose,"\n");
  for (i=0;i<nz+1;i++) fprintf(f3ddose,"%7.3f ",zbound[i]);
  fprintf(f3ddose,"\n");
      for (k=0;k<nz;k++)
	{
	for (j=0;j<ny;j++)
	  for (i=0;i<nx;i++)
	    { 
	      addr = i + j*nx + k*nx*ny;
	      if(*(sumfield + addr)==0) fprintf(f3ddose,"  0.");
	      else fprintf(f3ddose,"  %12.9e",*(sumfield+addr));
	    }
	fprintf(f3ddose,"\n");
	  }

      for (k=0;k<nz;k++)
	for (j=0;j<ny;j++)
	  for (i=0;i<nx;i++)
	    { 
	      addr = i + j*nx + k*nx*ny;
	      fprintf(f3ddose,"  %12.9f",*(sumerror+addr));
	    }      
      
  fclose(f3ddose);
  printf("Done.\n");
};



