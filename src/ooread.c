#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "oofit.h"

// create a function that reads the .oo file

char *zStrtok(char *str, const char *delim) {
    // https://github.com/fnoyanisi/zString
    // thanks for this function
    static char *static_str=0;      /* var to store last address */
    int index=0, strlength=0;           /* integers for indexes */
    int found = 0;                  /* check if delim is found */

    /* delimiter cannot be NULL
    * if no more char left, return NULL as well
    */
    if (delim==0 || (str == 0 && static_str == 0))
        return 0;

    if (str == 0)
        str = static_str;

    /* get length of string */
    while(str[strlength])
        strlength++;

    /* find the first occurance of delim */
    for (index=0;index<strlength;index++)
        if (str[index]==delim[0]) {
            found=1;
            break;
        }

    /* if delim is not contained in str, return str */
    if (!found) {
        static_str = 0;
        return str;
    }

    /* check for consecutive delimiters
    *if first char is delim, return delim
    */
    if (str[0]==delim[0]) {
        static_str = (str + 1);
        return (char *)delim;
    }

    /* terminate the string
    * this assignmetn requires char[], so str has to
    * be char[] rather than *char
    */
    str[index] = '\0';

    /* save the rest of the string */
    if ((str + index + 1)!=0)
        static_str = (str + index + 1);
    else
        static_str = 0;

        return str;
}

const char* getfield(char* line, int num) {
    const char* tok;
    for (tok = strtok(line, ",");
            tok && *tok;
            tok = strtok(NULL, ",\n"))
    {
        if (!--num)
            return tok;
    }
    return NULL;
}



int oosize(char * filename) {
    FILE * stream = fopen(filename, "r");
    if (stream == NULL) {
        printf("Error: file not found");
    }
    int lcount = 0;
    char line[1024];
    char * oo = "OOCYTE";
    int strsz = 6;
    while (fgets(line, 1024, stream)) {
        char * tmp = strdup(line);
        if (strncmp(getfield(tmp,1), oo, strsz) == 0) {
            lcount ++;
        }
        free(tmp);
    }
    return lcount;
}
void ooparse(char * filename, struct oorow * dptr) {
    // given a file name and a structure will parse the contents of the file
    // into the structure. No return value, just fills the structure.
   FILE * stream = fopen(filename, "r");
    if (stream == NULL) {
        printf("Error: file not found");
    }
    int lcount = 0;
    char line[1024];
    char * oo = "OOCYTE";
    int strsz = 6;
    while (fgets(line, 1024, stream)) {
        char * tmp2 = strdup(line);
        char * tmp3 = strdup(line);
        if (strncmp(getfield(tmp2,1), oo, strsz) == 0) {
            char *token = NULL;
            int n_tokens = 0;
            // Split the string into tokens delimited by spaces and commas
            token = zStrtok (tmp3,",");
            while (token != NULL) {
                // Different than strtok, splits the string even if null
                token = zStrtok (NULL, ",");
                char *eptr;
                switch(n_tokens) {
                    case 0:
                        strcpy(dptr->file, token);
                    case 1:
                        strcpy(dptr->glun1, token);
                    case 2:
                        strcpy(dptr->glun2, token);
                    case 3:
                        dptr->i = strtod(token, &eptr);
                    case 4:
                        if (strcmp(token, ",") > 1) {
                            dptr->sconc = strtod(token, &eptr);
                        } else {
                            dptr->sconc = NAN;
                        }
                    case 5:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm10.res = strtod(token, &eptr);
                        } else {
                            dptr->logm10.res = NAN;
                        }
                        dptr->logm10.dose = log10(0.0000000001);
                    case 6:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm9p5.res = strtod(token, &eptr);
                        } else {
                            dptr->logm9p5.res = NAN;
                        }
                        dptr->logm9p5.dose = log10(0.0000000003);
                    case 7:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm9.res = strtod(token, &eptr);
                        } else {
                            dptr->logm9.res = NAN;
                        }
                        dptr->logm9.dose = log10(0.000000001);
                    case 8:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm8p5.res = strtod(token, &eptr);
                        } else {
                            dptr->logm8p5.res = NAN;
                        }
                        dptr->logm8p5.dose = log10(0.000000003);
                    case 9:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm8.res = strtod(token, &eptr);
                        } else {
                            dptr->logm8.res = NAN;
                        }
                        dptr->logm8.dose = log10(0.00000001);
                    case 10:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm7p5.res = strtod(token, &eptr);
                        } else {
                            dptr->logm7p5.res = NAN;
                        }
                        dptr->logm7p5.dose = log10(0.00000003);
                    case 11:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm7.res = strtod(token, &eptr);
                        } else {
                            dptr->logm7.res = NAN;
                        }
                        dptr->logm7.dose = log10(0.0000001);
                    case 12:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm6p5.res = strtod(token, &eptr);
                        } else {
                            dptr->logm6p5.res = NAN;
                        }
                        dptr->logm6p5.dose = log10(0.0000003);
                    case 13:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm6.res = strtod(token, &eptr);
                        } else {
                            dptr->logm6.res = NAN;
                        }
                        dptr->logm6.dose = log10(0.000001);
                    case 14:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm5p5.res = strtod(token, &eptr);
                        } else {
                            dptr->logm5p5.res = NAN;
                        }
                        dptr->logm5p5.dose = log10(0.000003);
                        
                    case 15:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm5.res = strtod(token, &eptr);
                        } else {
                            dptr->logm5.res = NAN;
                        }
                        dptr->logm5.dose = log10(0.00001);
                    case 16:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm4p5.res = strtod(token, &eptr);
                        } else {
                            dptr->logm4p5.res = NAN;
                        }
                        dptr->logm4p5.dose = log10(0.00003);
                    case 17:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm4.res = strtod(token, &eptr);
                        } else {
                            dptr->logm4.res = NAN;
                        }
                        dptr->logm4.dose = log10(0.0001);
                    case 18:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm3p5.res = strtod(token, &eptr);
                        } else {
                            dptr->logm3p5.res = NAN;
                        }
                        dptr->logm3p5.dose = log10(0.0003);
                    case 19:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm3.res = strtod(token, &eptr);
                        } else {
                            dptr->logm3.res = NAN;
                        }
                        dptr->logm3.dose = log10(0.001);
                    case 20:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm2p5.res = strtod(token, &eptr);
                        } else {
                            dptr->logm2p5.res = NAN;
                        }
                        dptr->logm2p5.dose = log10(0.003);
                    case 21:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm2.res = strtod(token, &eptr);
                        } else {
                            dptr->logm2.res = NAN;
                        }
                        dptr->logm2.dose = log10(0.01);
                    case 22:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm1p5.res = strtod(token, &eptr);
                        } else {
                            dptr->logm1p5.res = NAN;
                        }
                        dptr->logm1p5.dose = log10(0.03);
                    case 23:
                        if (strcmp(token, ",") > 1) {
                            dptr->logm1.res = strtod(token, &eptr);
                        } else {
                            dptr->logm1.res = NAN;
                        }
                        dptr->logm1.dose = log10(0.1);
                    case 24:
                        strcpy(dptr->notes, token);
                    case 25:
                        strcpy(dptr->dbinfo, token);
                    case 26:
                        if (strcmp(token, ",") > 1) {
                            dptr->logec50 = strtod(token, &eptr);
                        } else {
                            dptr->logec50 = NAN;
                        }
                    case 27:
                        if (strcmp(token, ",") > 1) {
                            dptr->hillslope = strtod(token, &eptr);
                        } else {
                            dptr->hillslope = NAN;
                        }
                    case 28:
                        if (strcmp(token, ",") > 1) {
                            dptr->ymin = strtod(token, &eptr);
                        } else {
                            dptr->ymin = NAN;
                        }
                    case 29:
                        if (strcmp(token, ",") > 1) {
                            dptr->ymax = strtod(token, &eptr);
                        } else {
                            dptr->ymax = NAN;
                        }
                    case 30:
                        if (strcmp(token, ",") > 1) {
                            dptr->equation = strtod(token, &eptr);
                        } else {
                            dptr->equation = NAN;
                        }
                    case 31:
                        if (strcmp(token, ",") > 1) {
                            dptr->conv = strtod(token, &eptr);
                        } else {
                            dptr->conv = NAN;
                        }
                    case 32:
                        strcpy(dptr->assay, token);
                    case 33:
                        strcpy(dptr->date_rec, token);
                    case 34:
                        strcpy(dptr->date_inj, token);
                    case 35:
                        if (strcmp(token, ",") > 1) {
                            dptr->vhold = strtod(token, &eptr);
                        } else {
                            dptr->vhold = NAN;
                        }
                    case 36:
                        strcpy(dptr->coagonist, token);
                    case 37:
                        if (strcmp(token, ",") > 1) {
                            dptr->phsol = strtod(token, &eptr);
                        } else {
                            dptr->phsol = NAN;
                        }
                    case 38:
                        strcpy(dptr->drug_id, token);
                    case 39:
                        strcpy(dptr->rig, token);
                    case 40:
                        strcpy(dptr->initials, token);
                    case 41:
                        strcpy(dptr->experiment_info, token);
                }//end case
                n_tokens++;
            }//end while
        // pointer to a structure, the next iteration will now point to the next
        // row in the array of structures.
        dptr++;
        }
    }
    fclose(stream);
}

