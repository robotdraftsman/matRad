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
  //get the name of the files you want to run dosxyz on:
  /*printf("Gimme the name base of your files (no numbers or extensions): ");
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
  printf("%s\n",filebase); */

  /* note: this assumes that the filename is formatted like:
  filebaseBeamnBeamleti with n = beam number and i = beamlet number.
  i isn't beamlet numbered in beam, but the number that corresponds to the beamlet phsp file 
  */

  //loop over beams and beamlets:
  
  //so when calling fn, go like:
  // ./autoDos filebase totnumBeams startingBeam startingBeamlet
  // int numBeams = argv[2];
  int numBeamlets = 323; //fixed for now (afaik this is the overall max?)
  int maxSubmit = 200;  //maximum number of jobs to submit to the cluster at a time
  int numSubitted = 0; 
  /*int nat = 3;
  int * n = nat; //beam we start at
  int iat = 4;
  int * i = iat; //beamlet within the beam we start at*/

  //argv[0] is filebase, 1 is n, and 2 is i

  if(argc != 5) {
    printf("You haven't given me the right number of  arguments. When running me, write: ./autoDos filebase totNumOfBeams startingBeamNumber startingBeamletNumber\n");
    return 0;
  }

  if(atoi(argv[3]) <= 0 || atoi(argv[3]) > atoi(argv[2])) {
    printf("the starting beam number you entered, n = %d, doesn't fall within the range of beams in use (total number of beams = %d)\n",atoi(argv[3]),atoi(argv[2]));
    return 0;
  }

  printf("total number of beams given = %d\n",atoi(argv[2]));
  printf("total number of beamlets = %d\n",numBeamlets);
  printf("starting beam = %d\n",atoi(argv[3]));
  printf("starting beamlet = %d\n",atoi(argv[4]));
  
  
  //so now to make it such that you type in starting beam number and starting beamlet number before running it
  //if no numbers given, then ask for them in here
  //if wrong number of arguments, freak out and yell at the user
  //etc. etc.
  
  //cases: if too many arguments, if two few (but not 0), if 0.
  //If first two, complain/quit. If 0 prompt for starting point. If none of those, we have 2 arguments -> examine them.
  
  //no wait: if 3 args, first is filebase, and other two are starting point. If 2, just starting point.
  
  /* !!! wasting too much time here so I'll just run it with all 3 inputs going when you start it
  printf("start the if statements\n");
  if(argc == 4) {
    sprintf(filebase,argv[0]);
    n = atoi(argv[1]);
    i = atoi(argv[2]);
  } else if(argc == 3) {
    printf("doing NIII\n");
    n = atoi(argv[0]);
    i = atoi(argv[1]);
    
    filebase = getInfo("gimme the name of the file base: ");
    
  } else if(argc == 1) {
    filebase = getInfo("gimme the name of the file base: ");
    printf("ok %d\n",atoi(getInfo("Now give me the beam number to start from: ")));  //and after this works check that it's in the right range. Same for i
    i = atoi(getInfo("And may I please have the beamlet number from which to begin? "));
    printf("n: %d; i: %d\n",n,i);
  } else {
    printf("You have the wrong number of arguments. Options:\n 3: file base, starting beam, starting beamlet\n 2: starting beam, starting beamlet\n 0: provide all information during runtime");
  }
*/
  
  /*
  Here we gotta get the person to tell us what number beam (n) and beamlet (i) to start on
  And also make sure that 1 <= n <= 5 and 1 <= i <= 323
  */  
  
  /*char * argv[];
  printf("Which beam do you want to start on (0 <= n <= 5): ");
  read = argv[];*/
  

  
  /* so after it submits 200 jobs to the cluster it'll stop and tell you what beam
  and beamlet number it just submitted, and give you the next beamlet number (if it's
  not on beamlet 323; if it is, it'll give you the next beam number, and beamlet number 1)
  so that when these jobs are done, you can go and run it again but starting from that number
  so that it'll submit the next 200. Keep doing that until you've done them all
  */
  
  int n, i;
  int numSubmitted = 0;
  
  
  
  for(n = atoi(argv[3]); n <= atoi(argv[2]); n++) {
    for(i = atoi(argv[4]); i <= numBeamlets; i++) {
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
