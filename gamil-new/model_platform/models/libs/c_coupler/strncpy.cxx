#include <stdio.h>
#include <string.h>
int main() {
    char name[] = {"Chinanet"},dest[20]={};
    strncpy(dest,name,3);
    dest[3]='\0';
    printf("%s\n",dest);
}
