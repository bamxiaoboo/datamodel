#include <stdio.h>
#include <string.h>
#include <iostream>
using namespace std;

bool get_next_attr(char *attr, char **line)
{
	
    if ((*line)[0] == '\0') {
        (*line) = NULL;
        return false;
    }

    while ((*line)[0] == ' ' || (*line)[0] == '\t')
        (*line) ++;
    
    while ((*line)[0] != '\0' && (*line)[0] != '\t') {
        *attr = (*line)[0];
        attr ++;
        (*line) ++;
    }
    if ((*line)[0] == '\t')
        (*line) ++;

    if (*attr == ' ' || *attr == '\t') {
        while (*attr == ' ' || *attr == '\t') {
            attr --;
        }
        attr ++;
    }
    *attr = '\0';
    
    return true;
}

int main() {
        char* line1;
        line1= strdup("gamil\tUSTAR\tgamil_2D_decomp_phys\tgamil_grid\t0\tUSTAR\treal4");
        char attr1[64];
        bool a = get_next_attr(attr1,&line1);
        cout << attr1 << "||" << endl;
        return 0;
}

