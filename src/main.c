#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>

#include "oofit.h"

struct hash {
    int key;
    char value[50];
};

void fitstruct (struct oorow * record, int size, int parameter) {
    //fitting
    // figure out the number of concentrations of the recording
    int i;
    // ---- FIT BY FILE ----
    for (i=0;i<size;i++) {
        int n = ooarrsize(&record[i]);
        // define the arrays based on the size
        double dose_arr[n];
        double res_arr[n];
        ootoarr(&record[i], res_arr, dose_arr, n);
        // debugging
        /*
        int j;
        for (j=0;j<n;j++) {
            printf("%.7g ", res_arr[j]);
            printf("%.7g ", dose_arr[j]);
        }
        */
        double ret_arr[6];
        if (parameter == 1) {
            fit1d(res_arr, dose_arr, n, ret_arr);
            //pfit(ret_arr);
        }
        if (parameter == 2) {
            fit2d(res_arr, dose_arr, n, ret_arr);
            //pfit(ret_arr);
        }
        if (parameter == 3) {
            fit3d(res_arr, dose_arr, n, ret_arr);
            //pfit(ret_arr);
        }
        if (parameter == 4) {
            fit4d(res_arr, dose_arr, n, ret_arr);
            //pfit(ret_arr);
        }
        struct oorow * dptr = &record[i];
        dptr->logec50 = ret_arr[0];
        dptr->hillslope = ret_arr[1];
        dptr->ymin = ret_arr[2];
        dptr->ymax = ret_arr[3];
        dptr->equation = ret_arr[4];
        dptr->conv = ret_arr[5];
    }
}

void quitfile() /* write error message and quit */
{
    fprintf(stderr, "The file name that you entered was invalid\n");
    exit(1);
}

int main() {
    printf("|=======================================================|\n");
    printf("| oofit v 1.0 written by James Allen for sftlab & CFERV |\n");
    printf("|=======================================================|\n");
    // get the filename
    printf("To get started, please choose a file.\n");
    printf("Below are the loaded .oo files:\n");
    printf("------------------------------\n");
    DIR * dir;
    struct dirent * ent;
    dir = opendir("../dir");
    int fcount = 0;
    while ((ent = readdir (dir)) != NULL) {
        if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
            continue;
        } else {
            //printf("%s \n", ent->d_name);
            fcount++;
        }
    }
    rewinddir(dir);
    struct hash dhash[fcount];
    int i = 0;
    char dirstr[8] = "../dir/";
    while ((ent = readdir (dir)) != NULL) {
        if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
            continue;
        } else {
            //printf("%s \n", ent->d_name);
            dhash[i].key = i;
            strcpy(dhash[i].value, dirstr);
            strcat(dhash[i].value, ent->d_name);
            printf("%-5d : %-20s\n", i, ent->d_name);
            i++;
        }
    }
    closedir(dir);
    int fint;
    int result = scanf("%d", &fint);
    if (result == EOF) {
    /* ... you're not going to get any input ... */
    }
    if (result == 0) {
        while (fgetc(stdin) != '\n');
    }
    char fname[50];
    for (i=0; i<fcount; i++) {
        if (fint == dhash[i].key) {
            strcpy(fname, dhash[i].value);
            printf("%s\n", fname);
        }
    }
    FILE * file;
    file = fopen(fname, "r");
    if (file){
        fclose(file);
        printf("Parsing your file...");
    }else{
        quitfile();
    }
    // find the size of the file in rows
    int size = oosize(fname);
    // define the structure based on the size
    struct oorow record[size];
    // fill the structure with data
    ooparse(fname, record);

    int param;
    printf("Please enter an INTEGER parameter setting. Options below:\n");
    printf("1 = BOTTOM constrained to 0, TOP constrained to 100\n");
    printf("2 = BOTTOM unconstrained,    TOP constrained to 100\n");
    printf("3 = BOTTOM constrained to 0, TOP unconstrained\n");
    printf("4 = BOTTOM unconstrained,    TOP unconstrained\n");
    printf("Equation: ");
    result = scanf("%d", &param);

    if (result == EOF) {
    /* ... you're not going to get any input ... */
    }
    if (result == 0) {
        while (fgetc(stdin) != '\n');
    }
    fitstruct(record, size, param);
    printf("==========================================================\n");
    printf("| %-20s | %-7s | %-5s | %-5s | %-5s |\n", "Construct", "logEC50", "nH", "Ymin", "Ymax");
    printf("|----------------------|---------|-------|-------|-------|\n");
    for (i=0; i<size; i++) {
        char construct[1024];
        char dash[2] = "/";
        strcpy(construct, record[i].glun1);
        strncat(construct, dash, 2);
        strncat(construct, record[i].glun2, 20);
        printf("| %-20s | %-7.3g | %-5.3g | %-5.3g | %-5.3g |\n", construct, record[i].logec50, record[i].hillslope, record[i].ymin, record[i].ymax);
    }
    printf("==========================================================\n");
    return 0;
}
