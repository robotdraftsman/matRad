#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>

#define MAXPAR 10
#define FALSE 0
#define TRUE !FALSE


int fileExists(const char *path);
char getInfo(const char *prompt);

char command[600];

int main(int argc, char** argv) {

  /* note: this assumes that the filename is formatted like:
  filebaseBeamnBeamleti with n = beam number and i = beamlet number.
  i isn't beamlet numbered in beam, but the number that corresponds to the beamlet phsp file 
  */
  
  //so when calling fn, go like:
  // ./autoDos filebase totnumBeams startingBeam startingBeamlet
  // int numBeams = argv[2];
  int numBeamlets = 361; //this depends on bixel resolution and max/min beamlet locations. The maximum number of beamlets you can have for any given beam
  int maxSubmit = 100;  //maximum number of jobs to submit to the cluster at a time
  int numSubitted = 0; 

  if(argc != 5) {
    //note: filebase = the name of the file without beam/beamlet numbers, e.g. for 'inputsBeam3Beamlet5.egsinp' the base is 'inputs'
    //totNumOfBeams = all the beams used in the whole simulation
      printf("You haven't given me the right number of  arguments. When running me, write: ./autoDos filebase totNumOfBeams startingBeamNumber startingBeamletNumber\n");
    return 0;
  }

  if(atoi(argv[3]) <= 0 || atoi(argv[3]) > atoi(argv[2])) {
    printf("the starting beam number you entered, n = %d, doesn't fall within the range of beams in use (total number of beams = %d)\n",atoi(argv[3]),atoi(argv[2]));
    return 0;
  }

  printf("total number of beams given = %d\n",atoi(argv[2]));
  printf("maximum number of beamlets per beam = %d\n",numBeamlets);
  printf("starting beam = %d\n",atoi(argv[3]));
  printf("starting beamlet = %d\n",atoi(argv[4]));
  
  /* so after it submits whatever specified number of jobs to the cluster it'll stop and tell you what beam
  and beamlet number it just submitted, and give you the next beamlet number (if it's
  not on beamlet 361; if it is, it'll give you the next beam number, and beamlet number 1)
  so that when these jobs are done, you can go and run it again but starting from that number
  so that it'll submit the next whatever number. Keep doing that until you've done them all
  */
  
 int n, i;
  int numSubmitted = 0;
  
 int iStart = atoi(argv[4]);
  
  for(n = atoi(argv[3]); n <= atoi(argv[2]); n++) {
    for(i = iStart; i <= numBeamlets; i++) {
      char filename[256];
      sprintf(filename,"%sBeam%dBeamlet%d.egsinp",argv[1],n,i);
      // printf("%s\n",filename);
      
      //do the command if the file exists; otherwise go on to the next one
      if(fileExists(filename)==1) {
        //sprintf(command,"$HEN_HOUSE/scripts/run_user_code_batch dosxyznrc %s 700icru",file); //submit to cluster
        sprintf(command,"$HEN_HOUSE/scripts/run_user_code_batch dosxyznrc %s 700icru",filename);
        printf("commandline=%s\n",command);
        system(command);
        numSubmitted = numSubmitted+1;
      }
      
      if(numSubmitted == maxSubmit){
        printf("Last submitted for this run = beam %d beamlet %d.\n",n,i);
        if(i<numBeamlets) {  //so if we can keep in the same beam but go on to the enxt beamelt next time:
          printf("To continue this batch, tell me to start from beam %d, beamlet %d.\n",n,i+1);
        } else {  //if that was the last beamlet in the beam:
          if(n == atoi(argv[2])) {  //if it's also the last beam -> we've reached the end:
            printf("All the jobs for all the beams and all the beamlets have been submitted!\nYay!\nDo your happy dance!\n");
          } else {  //if it was the last beamlet in that beam but more beam(s) await:
            printf("To continue this batch, tell me to start from beam %d, beamlet %d.\n",n+1,1);
          }
        }
        
        exit(0);
      } else if(n==atoi(argv[2]) && i == numBeamlets) {
        printf("Looks like we've reached the end of all the beams and beamlets!\n");
      }
    }
    iStart = 1;
  }
  
  
  return 0;
}


//this function checks if the file exists:
int fileExists(const char *path)
{
    // Check for file existence
    if (access(path, F_OK) == -1)
        return 0;

    return 1;
}

char getInfo(const char *prompt) {
  printf(prompt);
  char *filebase = NULL;
  int read;
  unsigned int len;
  read = getline(&filebase, &len, stdin);
  if (-1 != read)
      puts(filebase);
  else
      printf("No line read...\n");

  //so this is our file name base:
  filebase[strlen(filebase)-1] = 0;
  printf("%s\n",filebase);
  
  return *filebase;
  //free(filebase);
}