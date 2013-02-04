//header file for probterms

//utility
void addIndivsToAFs(int **geno1, int **geno2, double ***AFs, double ***newAFs, int npops, int *nalleles, int **npeople, int addpop);

//probabilities
void pGenosGivenK(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, double *pmatches, double theta, double *ks);
void pGenoGivenK(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *pmatches, double theta, double *ks);

//derivatives of probabilities
void deltapGenoGivenK(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *deltaps, double theta, double *ks, int allelenum);

//derivatives of log likelihood ratios
void delta_ln_LR_given_l_i(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *dlnLRli, double theta, double *ksP, int allelenum);

//functions of derivatives of log likelihood ratios
void sum_square_dlnLRl(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *sumsqs, double theta, double *ksP);
void sum_prod_dlnLRl(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *sums, double theta, double *ksP);

//variance of log liklihood ratios
void var_ln_LR_given_l(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *vars, double theta, double *ksP);
void var_ln_LR(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, double *vars, double theta, double *ksP);

//likelihood ratios
void ln_LR(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, double *LRs, double theta, double *ksp);
void LR(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, double *LRs, double theta, double *ksp);

//confidence intervals
void CIs(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, double *lci, double *uci, double theta, double *ksP, double z);

//Y hap functions
int checkHapMatch(int *hap1, int *hap2, int nloci);
int countHaps(int *hap, int ***refHaps, int popi, int nindivs, int nloci);
void YCLhap(int *hap1, int *hap2, int ***refHaps, int npops, int *nindivs, int nloci, double *YCLs, double theta);
void YCL(int hap1, int hap2, double **hapFreqs, int npops, int *nindivs, double *YCLs);
void YLR(int hap1, int hap2, double **hapFreqs, int npops, int *nindivs, double *YLRs);
void calcHapFreqs(int ***refHaps, int npops, int *nindivs, int nloci, double **hapFreqs);
