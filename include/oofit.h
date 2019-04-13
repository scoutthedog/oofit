struct drp {
    double dose;
    double res;
};
struct oorow {
    char file[30];
    char glun1[30];
    char glun2[30];
    double i;
    double sconc;
    struct drp logm10;
    struct drp logm9p5;
    struct drp logm9;
    struct drp logm8p5;
    struct drp logm8;
    struct drp logm7p5;
    struct drp logm7;
    struct drp logm6p5;
    struct drp logm6;
    struct drp logm5p5;
    struct drp logm5;
    struct drp logm4p5;
    struct drp logm4;
    struct drp logm3p5;
    struct drp logm3;
    struct drp logm2p5;
    struct drp logm2;
    struct drp logm1p5;
    struct drp logm1;
    char notes[1024];
    char dbinfo[4];
    double logec50;
    double hillslope;
    double ymin;
    double ymax;
    double equation;
    double conv;
    char assay[10];
    char date_inj[1024];
    char date_rec[1024];
    double vhold;
    char coagonist[1024];
    double phsol;
    char drug_id[1024];
    char rig[1024];
    char initials[1024];
    char experiment_info[1024];
};

double * fit1(double * res, double * dose, int n, double *lb, double * ub, double * ret_arr);
double * fit1d(double * res, double * dose, int n, double * ret_arr);
double * fit2(double * res, double * dose, int n, double * lb, double * ub, double * ret_arr);
double * fit2d(double * res, double * dose, int n, double * ret_arr);
double * fit3(double * res, double * dose, int n, double * lb, double * ub, double * ret_arr);
double * fit3d(double * res, double * dose, int n, double * ret_arr);
double * fit4(double * res, double * dose, int n, double * lb, double * ub, double * ret_arr);
double * fit4d(double * res, double * dose, int n, double * ret_arr);
void pfit(double * ret_arr);
int oosize(char * filename);
int ooarrsize(struct oorow * record);
void * ootoarr(struct oorow * record, double * res, double * dose, int n);
void ooparse(char * filename, struct oorow * dptr);



