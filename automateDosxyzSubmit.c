#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>

#define MAXPAR 10
#define FALSE 0
#define TRUE !FALSE


int main() {

sprintf(command,"$HEN_HOUSE/scripts/run_user_code_batch dosxyznrc smolGroup1 700icru"); /*submit to cluster*/
printf("commandline=%s\n",command);
system(command);

}