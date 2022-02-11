#include <stdio.h>

int main(){
    int i=0,j=0;
    while(i<=5){
        j = i++;
        printf("%d\n",j);
    }

    printf("\n\n");
    i=j=0;

    while(i<=5){
        j = ++i;
        printf("%d\n",j);
    }
}

//  Thus, we can see: i++ first gets assigned, then incremented
// but, ++i first gets incremented then assigned

// We really doesn't see any differences until we have something like
// `j = i++` or `j = ++i` in the code, where there is an assignment operator