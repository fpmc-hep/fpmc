// Computes different helicity amplitudes as defined in 
// Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

int limits(double sred, double tred);

void Mxxxx_fermion(double x, double y, double * re, double * im);

void Mpppp_fermion(double sred, double tred, double *re, double *im, int exclude_loops);

void Mpmmp_fermion(double sred, double tred, double *re, double *im, int exclude_loops);

void Mpmpm_fermion(double sred, double tred, double *re, double *im, int exclude_loops);

void Mpppm_fermion(double sred, double tred, double * re, double * im, int exclude_loops);

void Mppmm_fermion(double sred, double tred, double * re, double * im, int exclude_loops);


void Mxxxx_vector(double x, double y, double * re, double * im);

void Mpppp_vector(double sred, double tred, double *re, double *im, int exclude_loops);

void Mpmmp_vector(double sred, double tred, double *re, double *im, int exclude_loops);

void Mpmpm_vector(double sred, double tred, double *re, double *im, int exclude_loops);

void Mpppm_vector(double sred, double tred, double * re, double * im, int exclude_loops);

void Mppmm_vector(double sred, double tred, double * re, double * im, int exclude_loops);


void Mxxxx_spin0even(double x, double y, double m, double f0, double w_const, double a2, double * re, double * im);

void Mpppp_spin0even(double x, double y, double m, double f0, double w_const, double a2, double * re, double * im);

void Mpmmp_spin0even(double x, double y, double m, double f0, double w_const, double a2, double * re, double * im);

void Mpmpm_spin0even(double x, double y, double m, double f0, double w_const, double a2, double * re, double * im);

void Mppmm_spin0even(double x, double y, double m, double f0, double w_const, double a2, double * re, double * im);

void Mpppm_spin0even(double x, double y, double m, double f0, double w_const, double a2, double * re, double * im);


void Mpppp_eft(double zeta1, double zeta2, double s, double t, double *re, double *im);
 
void Mpmmp_eft(double zeta1, double zeta2, double s, double t, double *re, double *im);
 
void Mpmpm_eft(double zeta1, double zeta2, double s, double t, double *re, double *im);
 
void Mppmm_eft(double zeta1, double zeta2, double s, double t, double *re, double *im);
 
void Mpppm_eft(double zeta1, double zeta2, double s, double t, double *re, double *im);


double width_gen(double s, double m, double f0, double w_const, double a2);

void MZZxxxx_spin0even(double x, double y, double m, double f0, double f0Z, double w_const, double a2, double * re, double * im);

void MZZpppp_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im);

void MZZpmmp_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im);

void MZZpmpm_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im);

void MZZppmm_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im);

void MZZpppm_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im);

void MZZpp00_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im);


void MZgxxxx_spin0even(double x, double y, double m, double f0, double f0Zg, double w_const, double a2, double * re, double * im);

void MZgpppp_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im);

void MZgpmmp_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im);

void MZgpmpm_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im);

void MZgppmm_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im);

void MZgpppm_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im);

void MZgpp00_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im);


void MWWxxxx_spin0even(double x, double y, double m, double f0, double f0W, double w_const, double a2, double * re, double * im);

void MWWpppp_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im);

void MWWpmmp_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im);

void MWWpmpm_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im);

void MWWppmm_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im);

void MWWpppm_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im);

void MWWpp00_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im);


void Mxxxx_spin2(double x, double y, double m, double f0, double w_const, double a2, double * re, double * im);

void Mpppp_spin2(double s, double t, double m, double f0, double w_const, double a2, double * re, double * im);

void Mpmmp_spin2(double s, double t, double m, double f0, double w_const, double a2, double * re, double * im);

void Mpmpm_spin2(double s, double t, double m, double f0, double w_const, double a2, double * re, double * im);

void Mppmm_spin2(double s, double t, double m, double f0, double w_const, double a2, double * re, double * im);

void Mpppm_spin2(double s, double t, double m, double f0, double w_const, double a2, double * re, double * im);


void Mggxxxx_spin0even(double x, double y, double m, double f0, double f0g, double w_const, double a2, double * re, double * im);

void Mggpppp_spin0even(double s, double t, double m, double f0, double f0g, double w_const, double a2, double *re, double *im);

void Mggpmmp_spin0even(double s, double t, double m, double f0, double f0g, double w_const, double a2, double *re, double *im);

void Mggpmpm_spin0even(double s, double t, double m, double f0, double f0g, double w_const, double a2, double *re, double *im);

void Mggppmm_spin0even(double s, double t, double m, double f0, double f0g, double w_const, double a2, double *re, double *im);

void Mggpppm_spin0even(double s, double t, double m, double f0, double f0g, double w_const, double a2, double *re, double *im);


void MhhOHxxxx_spin0even(double x, double y, double m, double f0, double f0H, double w_const, double a2, double * re, double * im);

void MhhOHpppp_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);

void MhhOHpmmp_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);

void MhhOHpmpm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);

void MhhOHppmm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);

void MhhOHpppm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);

void MhhOHpp00_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);


void Mpppp_eft_AZ(double zeta1, double s, double t, double *re, double *im);
 
void Mpmmp_eft_AZ(double zeta1, double s, double t, double *re, double *im);
 
void Mpmpm_eft_AZ(double zeta1, double s, double t, double *re, double *im);
 
void Mppmm_eft_AZ(double zeta1, double s, double t, double *re, double *im);
 
void Mpppm_eft_AZ(double zeta1, double s, double t, double *re, double *im);

void Mppmp_eft_AZ(double zeta1, double s, double t, double *re, double *im);


void MZZOHxxxx_spin0even(double x, double y, double m, double f0, double f0H, double w_const, double a2, double * re, double * im);

void MZZOHpppp_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);

void MZZOHpmmp_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);

void MZZOHpmpm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);

void MZZOHppmm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);

void MZZOHpppm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);

void MZZOHpp00_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im);
