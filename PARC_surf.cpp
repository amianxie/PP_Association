/* c++ P-p association rate constant prediction main program: code/
 
//pre-processor directive
/* #include #define #ifdef #endif*/
#include <iostream> //C++ newly include
#include <cstdio>   //prconst intf[ scanf[ getchar[ gets[ getch[ getche[
#include <cstdlib>  //system[ rand[ srand[
#include <math.h>   //sin[ tan[ sqrt[ exp[ log[
#include <time.h>   //time[ clock[
#include <cstring>  //strcat[ strcpy[ strlen[
#include <fstream>  //
#include <sstream>
#include <iomanip>  //

//external declaration
using namespace std;

    FILE *ifi, *ofi, *sfi, *lfi, *efi, *cfi, *tfi, *ffi; 
          int chainA=864;//839;//2237;//275;
          int chainD=693;//699;//886;//200;//missing residues
          int chainX=108;
          int chainY=87;
    const double dist_cutoff=5.0;//4.0; //7.5
    const int nof_res=4;
    const int trial_max=1000000;
    const int trj_max=1; //100
    const double PI=3.1415926;
    const double cell_range_x = 100.0;   // Unit: A
    const double cell_range_y = 100.0;  // Unit: A
    const double cell_range_z = 100.0; // Unit: A
    const double dist_constrain = cell_range_x/2.; //50.0 //25.0;
    const double clash_cutoff=5.0;//1.0; //3.0
    const double bond_thetapd=180.0,bond_thetapd_cutoff=30.0;
    const double bond_thetaot=90.0,bond_thetaot_cutoff=30.0;
    const int simu_step=1; //1001 //test-part !!can't be 1 or won't run
    const double temp_ini=1.0,temp_fin=1.0;
          double simu_temp=1.;
          double wel=1.;
          double whp=1.;
          double xi=8.2;
          double time_step=1.; // unit: nanosecond
          double distance_step=10.0*pow(time_step,0.5); // unit: Angstrom
          double MHC_D=11.5*time_step; //  unit: Angstrom square per nanosecond MHC ==> barnase
          double TCR_D=12.2*time_step; //  unit: Angstrom square per nanosecond TCR ==> barstar
          double MHC_rot_D=4.3*time_step;  //unit: degree per ns
          double TCR_rot_D=5.15*time_step; //unit: degree per ns
    const double clash_ene=20.0;//20.0;
    const double go_ene=-5.0;//-5.0;
    const double well_width=2.0;
    const double k_theta_ot=0.1;
    const double k_theta_pd=0.1;
    const double k_dist_constrain=20.0;
    const double RMSD_cutoff=4.0;
    
    int i=0,j=0,k=0;
    int nof_TCR=0,nof_MHC=0;
    double dist=0.;
    
    int nof_pair=0;
    int native_cutoff=-2;
    int base=0;
    double native_rmsd = 10.;
    double cm1_a_x=0.;
    double cm1_a_y=0.;
    double cm1_a_z=0.;
    double cm1_b_x=0.;
    double cm1_b_y=0.;
    double cm1_b_z=0.;
    double dist_cm1_IJ=0.;
    double dist_temp=0.,dist_vect=0.;
    double vect_x=0.,vect_y=0.,vect_z=0.;
    double RB_TCR_x[nof_res];
    double RB_TCR_y[nof_res];
    double RB_TCR_z[nof_res];
    double RB_MHC_x[nof_res];
    double RB_MHC_y[nof_res];
    double RB_MHC_z[nof_res];
    double RB_TCR_x_0[nof_res];
    double RB_TCR_y_0[nof_res];
    double RB_TCR_z_0[nof_res];
    double RB_MHC_x_0[nof_res];
    double RB_MHC_y_0[nof_res];
    double RB_MHC_z_0[nof_res];
    double RB_TCR_x_rec[nof_res];
    double RB_TCR_y_rec[nof_res];
    double RB_TCR_z_rec[nof_res];
    double RB_MHC_x_rec[nof_res];
    double RB_MHC_y_rec[nof_res];
    double RB_MHC_z_rec[nof_res];
    double rand_x,rand_y,rand_z;
    //double theta=0.,phi=0.,psi=0.,phai=0.;
    double t[3][3];
    int n_trial=0;
    int nof_trj=0;
    double vector[2][3]; //from coordinates - vector - theta //point_x[2],point_y[2],point_z[2];
    double theta_pd=0.,theta_ot=0.;
    int mc_time_step=0;
    double Prob_Diff=0.;
    int pair_index=0;
    int atompr[1][1000];
    double xa[1000],ya[1000],za[1000],xb[1000],yb[1000],zb[1000];
    int npair=0,npair_a=0,npair_b=0;
    double RMSD=0.;
    double sum_diff=0.,sum_diff2=0., ego = 0., edp = 0., elj = 0., ehp=0., eele=0., eclash=0.;
    double current_simu_time=0.;
    int tp=0;
    double tn = 0.;
    int output_str_step=10, output_ene_step=1;//1000; //10 //test-part
    bool clash_flag=false, output_ene_flag=true, output_str_flag=false, comp_rmsd=true;
    int mem[17];
    double trj[8];
    //int index;

    void ini_TM(int nof_res, int nof_tcr, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3);
    void translation(int nof_res, int nof_tcr, double dist, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3); //translate single
    void rotation (int nof_res, int nof_tcr, double upper_limit, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3); //rotate single
    void rr(int nof_p, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3); //record-recover
    void BC(int nof_res, int nof_res1, double cell_range_x, double cell_range_y, double cell_range_z, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3);
    int random_i(int upper_limit);
    int interface_np(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *p4, double *p5, double *q1, double *q2, double *q3, double *q4, double *q5, bool *r, double *s); 
    double random_d(double upper_limit);
    double golike(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r, double *s); //GO-like potential    
    double ori(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r); //GO-like potential    
    double dppdlike(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r, double *s); //DPPD-like potential    
    double ljlike(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r, double *s); //Lennard-Jones-like potential    
    double ljgo(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r, double *s); //Lennard-Jones-like potential    
    double Eele(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *p4, double *q1, double *q2, double *q3, double *q4, bool *r); //Sichun Yang's potential    
    double HP(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *p4, double *q1, double *q2, double *q3, double *q4, bool *r); //account for hydrophobic residues  
    double ncl(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3); //nof clash
    double distance(double a, double b, double c, double x, double y, double z);
    double dp(double *p, double *q);
    double dot2theta(double *p, double *q);
    double rmsdcall(int x, int y, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, double *r1, double *r2, double *r3, double *s1, double *s2, double *s3);
    double rmsdcompute(int n, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3);
    extern "C"
    {
        void rotlsq_(double *xo, double *yo, double *zo, int &nof_atom1, double *x2, double *y2, double *z2, int &nof_atom2, int *p, int &nof_atom3);
        //void cgsas_(int &nof_TCR, double *p1, double *p2, double *p3, char *p4, int &nof_MHC, double *q1, double *q2, double *q3, char *q4, double &ewu);
    }

    //         read structure
    int main (int argc, char *argv[])
    {
        double rmsd=0.0;
        chainA   = atoi(argv[5]);//atom
        chainD   = atoi(argv[6]);
        chainX   = atoi(argv[7]);//res.
        chainY   = atoi(argv[8]);
        MHC_D    = atof(argv[9])*time_step;
        TCR_D    = atof(argv[10])*time_step;
        MHC_rot_D= atof(argv[11])*time_step*2.;
        TCR_rot_D= atof(argv[12])*time_step*2.;
        wel      = atof(argv[13]);
        whp      = atof(argv[14]);
        xi       = atof(argv[15]);
        simu_temp= atof(argv[16]);
        
        char MHC_na[chainA][3];
        char TCR_na[chainD][3];
        char MHC_at[chainA][3];
        char TCR_at[chainD][3];
         
        double MHC_x_00[chainA];
        double MHC_y_00[chainA];
        double MHC_z_00[chainA];
        int    MHC_n_00[chainA]; //rec res no.
        double MHC_x_0[chainX][2]; //ini coor
        double MHC_y_0[chainX][2];
        double MHC_z_0[chainX][2];
        double MHC_x[chainX][2]; //tmp coor
        double MHC_y[chainX][2];
        double MHC_z[chainX][2];
        double MHC_q[chainX];  //ele
        double MHC_p[chainX];  //hp
        bool   MHC_i[chainX]; //interface
        double MHC_x_rec[chainX][2];
        double MHC_y_rec[chainX][2];
        double MHC_z_rec[chainX][2];
        double TCR_x_00[chainD];
        double TCR_y_00[chainD];
        double TCR_z_00[chainD];
        int    TCR_n_00[chainD];
        double TCR_x_0[chainY][2];
        double TCR_y_0[chainY][2];
        double TCR_z_0[chainY][2];
        double TCR_x[chainY][2];
        double TCR_y[chainY][2];
        double TCR_z[chainY][2];
        double TCR_q[chainY];
        double TCR_p[chainY];
        bool   TCR_i[chainY];
        double TCR_x_rec[chainY][2];
        double TCR_y_rec[chainY][2];
        double TCR_z_rec[chainY][2];
        bool pair_list[chainX][chainY];
        double pair_dist[chainX][chainY];
        
        ifstream ifi(argv[1],ios::in); 
        ofstream sfi(argv[3],ios::out|ios::app);
        string buf;
        
        if(!ifi)
        {
            cout << "unable to open ifi\n";
        }
        tp=-1;
        for(j=0;j<chainA;j++)
        {
            getline(ifi,buf);
            MHC_x_00[j] = atof(buf.substr(30,8).c_str());
            MHC_y_00[j] = atof(buf.substr(38,8).c_str());
            MHC_z_00[j] = atof(buf.substr(46,8).c_str());
            strcpy(MHC_na[j],buf.substr(17,3).c_str());
            strcpy(MHC_at[j],buf.substr(13,3).c_str());
            MHC_n_00[j] = tp;
            if (MHC_at[j][0] == 'N' and MHC_at[j][1] == ' ')
            {
                tn=0.; //
                tp++;  //
                MHC_q[tp]=0.;
		MHC_p[tp]=0.;
                MHC_n_00[j] = tp;
            }
            else if (MHC_at[j][0] == 'C' and MHC_at[j][1] == 'A')
            {
                MHC_x_0[tp][0]=MHC_x_00[j];
                MHC_y_0[tp][0]=MHC_y_00[j];
                MHC_z_0[tp][0]=MHC_z_00[j];
                tn=0.;
                MHC_q[tp]=0.;
		MHC_p[tp]=0.;
                MHC_x_0[tp][1]=MHC_x_00[j];
                MHC_y_0[tp][1]=MHC_y_00[j];
                MHC_z_0[tp][1]=MHC_z_00[j];
            }
            else if(MHC_na[j][0] == 'A' and MHC_na[j][1] == 'S')// and MHC_na[j][2] == 'P')
            {
                if(MHC_at[j][1] == 'D')
                {
                    tn++;
                    MHC_x_0[tp][1]=(MHC_x_0[tp][1]*(tn-1.)+MHC_x_00[j])/tn;
                    MHC_y_0[tp][1]=(MHC_y_0[tp][1]*(tn-1.)+MHC_y_00[j])/tn;
                    MHC_z_0[tp][1]=(MHC_z_0[tp][1]*(tn-1.)+MHC_z_00[j])/tn;
                    if (MHC_na[j][2] == 'P') MHC_q[tp]=-1.;
		    //MHC_p[tp]=0.;
                }
            }
            else if(MHC_na[j][0] == 'G' and MHC_na[j][1] == 'L')// and MHC_na[j][2] == 'U')
            {
                if(MHC_at[j][1] == 'E')//(MHC_at[j][0] == 'O' and MHC_at[j][1] == 'E' )
                {
                    tn++;
                    MHC_x_0[tp][1]=(MHC_x_0[tp][1]*(tn-1.)+MHC_x_00[j])/tn;
                    MHC_y_0[tp][1]=(MHC_y_0[tp][1]*(tn-1.)+MHC_y_00[j])/tn;
                    MHC_z_0[tp][1]=(MHC_z_0[tp][1]*(tn-1.)+MHC_z_00[j])/tn;
                    if (MHC_na[j][2] == 'U') MHC_q[tp]=-1.;
		    //MHC_p[tp]=0.;
                }
            }
            else if(MHC_na[j][0] == 'L' and MHC_na[j][1] == 'Y' and MHC_na[j][2] == 'S')
            {
                if(MHC_at[j][1] == 'Z')//(MHC_at[j][0] == 'N' and MHC_at[j][1] == 'Z' )
                {
                    MHC_x_0[tp][1]=MHC_x_00[j];
                    MHC_y_0[tp][1]=MHC_y_00[j];
                    MHC_z_0[tp][1]=MHC_z_00[j];
                    MHC_q[tp]=1.;
		    //MHC_p[tp]=0.;
                }
            }
            else if(MHC_na[j][0] == 'A' and MHC_na[j][1] == 'R' and MHC_na[j][2] == 'G')
            {
                if(MHC_at[j][1] == 'H')//(MHC_at[j][0] == 'N' and MHC_at[j][1] == 'H' )
                {
                    tn++;
                    MHC_x_0[tp][1]=(MHC_x_0[tp][1]*(tn-1.)+MHC_x_00[j])/tn;
                    MHC_y_0[tp][1]=(MHC_y_0[tp][1]*(tn-1.)+MHC_y_00[j])/tn;
                    MHC_z_0[tp][1]=(MHC_z_0[tp][1]*(tn-1.)+MHC_z_00[j])/tn;
                    MHC_q[tp]=1.;
		    //MHC_p[tp]=0.; 
                }
            }
            else if(MHC_na[j][0] == 'H' and MHC_na[j][1] == 'I' and MHC_na[j][2] == 'S')
            {
                if(MHC_at[j][1] == 'G' or MHC_at[j][1] == 'D' or MHC_at[j][1] == 'E')
                {
                    tn++;
                    MHC_x_0[tp][1]=(MHC_x_0[tp][1]*(tn-1.)+MHC_x_00[j])/tn;
                    MHC_y_0[tp][1]=(MHC_y_0[tp][1]*(tn-1.)+MHC_y_00[j])/tn;
                    MHC_z_0[tp][1]=(MHC_z_0[tp][1]*(tn-1.)+MHC_z_00[j])/tn;
                    MHC_q[tp]=.5;
		    //MHC_p[tp]=0.;
                }
            }
            else if(MHC_na[j][0] == 'P' and MHC_na[j][1] == 'H' and MHC_na[j][2] == 'E')
            {
                if(MHC_at[j][1] == 'G' or MHC_at[j][1] == 'D' or MHC_at[j][1] == 'E' or MHC_at[j][1] == 'Z')
                {
                    tn++;
                    MHC_x_0[tp][1]=(MHC_x_0[tp][1]*(tn-1.)+MHC_x_00[j])/tn;
                    MHC_y_0[tp][1]=(MHC_y_0[tp][1]*(tn-1.)+MHC_y_00[j])/tn;
                    MHC_z_0[tp][1]=(MHC_z_0[tp][1]*(tn-1.)+MHC_z_00[j])/tn;
                    MHC_p[tp]=-2.8;
                }
            }
            else if(MHC_na[j][0] == 'T' and MHC_na[j][1] == 'R' and MHC_na[j][2] == 'P')
            {
                if(MHC_at[j][1] == 'G' or MHC_at[j][1] == 'D' or MHC_at[j][1] == 'E' or MHC_at[j][1] == 'Z' or MHC_at[j][1] == 'H')
                {
                    tn++;
                    MHC_x_0[tp][1]=(MHC_x_0[tp][1]*(tn-1.)+MHC_x_00[j])/tn;
                    MHC_y_0[tp][1]=(MHC_y_0[tp][1]*(tn-1.)+MHC_y_00[j])/tn;
                    MHC_z_0[tp][1]=(MHC_z_0[tp][1]*(tn-1.)+MHC_z_00[j])/tn;
		    //MHC_p[tp]=0.;
                }
            }
            else if(MHC_na[j][0] == 'T' and MHC_na[j][1] == 'Y' and MHC_na[j][2] == 'R')
            {
                if(MHC_at[j][1] == 'H')
                {
                    MHC_x_0[tp][1]=MHC_x_00[j];
                    MHC_y_0[tp][1]=MHC_y_00[j];
                    MHC_z_0[tp][1]=MHC_z_00[j];
		    //MHC_p[tp]=0.;
                }
            }
            else if(MHC_na[j][0] == 'S' and MHC_na[j][1] == 'E' and MHC_na[j][2] == 'R')
            {
                if(MHC_at[j][1] == 'G')
                {
                    MHC_x_0[tp][1]=MHC_x_00[j];
                    MHC_y_0[tp][1]=MHC_y_00[j];
                    MHC_z_0[tp][1]=MHC_z_00[j];
		    //MHC_p[tp]=0.;
                }
            }
            else if(MHC_na[j][0] == 'C' and MHC_na[j][1] == 'Y' and MHC_na[j][2] == 'S')
            {
                if(MHC_at[j][1] == 'G')
                {
                    MHC_x_0[tp][1]=MHC_x_00[j];
                    MHC_y_0[tp][1]=MHC_y_00[j];
                    MHC_z_0[tp][1]=MHC_z_00[j];
                    MHC_p[tp]=-2.5;
                }
            }
            else if(MHC_na[j][0] == 'T' and MHC_na[j][1] == 'H' and MHC_na[j][2] == 'R')
            {
                if(MHC_at[j][1] == 'G')
                {
                    tn++;
                    MHC_x_0[tp][1]=(MHC_x_0[tp][1]*(tn-1.)+MHC_x_00[j])/tn;
                    MHC_y_0[tp][1]=(MHC_y_0[tp][1]*(tn-1.)+MHC_y_00[j])/tn;
                    MHC_z_0[tp][1]=(MHC_z_0[tp][1]*(tn-1.)+MHC_z_00[j])/tn;
		    //MHC_p[tp]=0.;
                }
            }
            else if(MHC_na[j][0] == 'V' and MHC_na[j][1] == 'A' and MHC_na[j][2] == 'L')
            {
                if(MHC_at[j][1] == 'G')
                {
                    tn++;
                    MHC_x_0[tp][1]=(MHC_x_0[tp][1]*(tn-1.)+MHC_x_00[j])/tn;
                    MHC_y_0[tp][1]=(MHC_y_0[tp][1]*(tn-1.)+MHC_y_00[j])/tn;
                    MHC_z_0[tp][1]=(MHC_z_0[tp][1]*(tn-1.)+MHC_z_00[j])/tn;
                    MHC_p[tp]=-4.2;
                }
            }
            else if(MHC_na[j][0] == 'L' and MHC_na[j][1] == 'E' and MHC_na[j][2] == 'U')
            {
                if(MHC_at[j][1] == 'D')
                {
                    tn++;
                    MHC_x_0[tp][1]=(MHC_x_0[tp][1]*(tn-1.)+MHC_x_00[j])/tn;
                    MHC_y_0[tp][1]=(MHC_y_0[tp][1]*(tn-1.)+MHC_y_00[j])/tn;
                    MHC_z_0[tp][1]=(MHC_z_0[tp][1]*(tn-1.)+MHC_z_00[j])/tn;
                    MHC_p[tp]=-3.8;
                }
            }
            else if(MHC_na[j][0] == 'M' and MHC_na[j][1] == 'E' and MHC_na[j][2] == 'T')
            {
                if(MHC_at[j][1] == 'E')
                {
                    MHC_x_0[tp][1]=MHC_x_00[j];
                    MHC_y_0[tp][1]=MHC_y_00[j];
                    MHC_z_0[tp][1]=MHC_z_00[j];
                    MHC_p[tp]=-1.9;
                }
            }
            else if(MHC_na[j][0] == 'I' and MHC_na[j][1] == 'L' and MHC_na[j][2] == 'E')
            {
                if(MHC_at[j][1] == 'D')
                {
                    MHC_x_0[tp][1]=MHC_x_00[j];
                    MHC_y_0[tp][1]=MHC_y_00[j];
                    MHC_z_0[tp][1]=MHC_z_00[j];
                    MHC_p[tp]=-4.5;
                }
            }
            else if(MHC_na[j][0] == 'A' and MHC_na[j][1] == 'L' and MHC_na[j][2] == 'A')
            {
                if(MHC_at[j][1] == 'B')
                {
                    MHC_x_0[tp][1]=MHC_x_00[j];
                    MHC_y_0[tp][1]=MHC_y_00[j];
                    MHC_z_0[tp][1]=MHC_z_00[j];
                    MHC_p[tp]=-1.8;
                }
            }
            else if(MHC_na[j][0] == 'P' and MHC_na[j][1] == 'R' and MHC_na[j][2] == 'O')
            {
                if(MHC_at[j][1] == 'G' or MHC_at[j][1] == 'D' or MHC_at[j][1] == 'B')
                {
                    tn++;
                    MHC_x_0[tp][1]=(MHC_x_0[tp][1]*(tn-1.)+MHC_x_00[j])/tn;
                    MHC_y_0[tp][1]=(MHC_y_0[tp][1]*(tn-1.)+MHC_y_00[j])/tn;
                    MHC_z_0[tp][1]=(MHC_z_0[tp][1]*(tn-1.)+MHC_z_00[j])/tn;
                    //MHC_p[tp]=0.;
                }
            }
        }
        
        tp=-1;
        for(j=0;j<chainD;j++)
        {
            getline (ifi,buf);
            TCR_x_00[j] = atof(buf.substr(30,8).c_str());
            TCR_y_00[j] = atof(buf.substr(38,8).c_str());
            TCR_z_00[j] = atof(buf.substr(46,8).c_str());
            strcpy(TCR_na[j],buf.substr(17,3).c_str());
            strcpy(TCR_at[j],buf.substr(13,3).c_str());
            TCR_n_00[j] = tp;
            if (TCR_at[j][0] == 'N' and TCR_at[j][1] == ' ')
            {
                tn = 0.;
                tp++;
                TCR_q[tp]=0.;
		TCR_p[tp]=0.;
                TCR_n_00[j] = tp;
            }
            else if (TCR_at[j][0] == 'C' and TCR_at[j][1] == 'A')
            {
                TCR_x_0[tp][0]=TCR_x_00[j];
                TCR_y_0[tp][0]=TCR_y_00[j];
                TCR_z_0[tp][0]=TCR_z_00[j];
                tn=0.;
                TCR_q[tp]=0.;
		TCR_p[tp]=0.;
                TCR_x_0[tp][1]=TCR_x_00[j];
                TCR_y_0[tp][1]=TCR_y_00[j];
                TCR_z_0[tp][1]=TCR_z_00[j];
            }
            else if(TCR_na[j][0] == 'A' and TCR_na[j][1] == 'S')// and TCR_na[j][2] == 'P')
            {
                if(TCR_at[j][1] == 'D')
                {
                    tn++;
                    TCR_x_0[tp][1]=(TCR_x_0[tp][1]*(tn-1.)+TCR_x_00[j])/tn;
                    TCR_y_0[tp][1]=(TCR_y_0[tp][1]*(tn-1.)+TCR_y_00[j])/tn;
                    TCR_z_0[tp][1]=(TCR_z_0[tp][1]*(tn-1.)+TCR_z_00[j])/tn;
                    if (TCR_na[j][2] == 'P') TCR_q[tp]=-1.;
		    //TCR_p[tp]=1.;
                }
            }
            else if(TCR_na[j][0] == 'G' and TCR_na[j][1] == 'L')// and TCR_na[j][2] == 'U')
            {
                if(TCR_at[j][1] == 'E')//(TCR_at[j][0] == 'O' and TCR_at[j][1] == 'E' )
                {
                    tn++;
                    TCR_x_0[tp][1]=(TCR_x_0[tp][1]*(tn-1.)+TCR_x_00[j])/tn;
                    TCR_y_0[tp][1]=(TCR_y_0[tp][1]*(tn-1.)+TCR_y_00[j])/tn;
                    TCR_z_0[tp][1]=(TCR_z_0[tp][1]*(tn-1.)+TCR_z_00[j])/tn;
                    if (TCR_na[j][2] == 'U') TCR_q[tp]=-1.;
		    //TCR_p[tp]=1.;
                }
            }
            else if(TCR_na[j][0] == 'L' and TCR_na[j][1] == 'Y' and TCR_na[j][2] == 'S')
            {
                if(TCR_at[j][1] == 'Z')//(TCR_at[j][0] == 'N' and TCR_at[j][1] == 'Z' )
                {
                    TCR_x_0[tp][1]=TCR_x_00[j];
                    TCR_y_0[tp][1]=TCR_y_00[j];
                    TCR_z_0[tp][1]=TCR_z_00[j];
                    TCR_q[tp]=1.;
		    //TCR_p[tp]=1.;
                }
            }
            else if(TCR_na[j][0] == 'A' and TCR_na[j][1] == 'R' and TCR_na[j][2] == 'G')
            {
                if(TCR_at[j][1] == 'H')//(TCR_at[j][0] == 'N' and TCR_at[j][1] == 'H' )
                {
                    tn++;
                    TCR_x_0[tp][1]=(TCR_x_0[tp][1]*(tn-1.)+TCR_x_00[j])/tn;
                    TCR_y_0[tp][1]=(TCR_y_0[tp][1]*(tn-1.)+TCR_y_00[j])/tn;
                    TCR_z_0[tp][1]=(TCR_z_0[tp][1]*(tn-1.)+TCR_z_00[j])/tn;
                    TCR_q[tp]=1.;
		    //TCR_p[tp]=1.; 
                }
            }
            else if(TCR_na[j][0] == 'H' and TCR_na[j][1] == 'I' and TCR_na[j][2] == 'S')
            {
                if(TCR_at[j][1] == 'G' or TCR_at[j][1] == 'D' or TCR_at[j][1] == 'E')
                {
                    tn++;
                    TCR_x_0[tp][1]=(TCR_x_0[tp][1]*(tn-1.)+TCR_x_00[j])/tn;
                    TCR_y_0[tp][1]=(TCR_y_0[tp][1]*(tn-1.)+TCR_y_00[j])/tn;
                    TCR_z_0[tp][1]=(TCR_z_0[tp][1]*(tn-1.)+TCR_z_00[j])/tn;
                    TCR_q[tp]=.5;
		    //TCR_p[tp]=1.;
                }
            }
            else if(TCR_na[j][0] == 'P' and TCR_na[j][1] == 'H' and TCR_na[j][2] == 'E')
            {
                if(TCR_at[j][1] == 'G' or TCR_at[j][1] == 'D' or TCR_at[j][1] == 'E' or TCR_at[j][1] == 'Z')
                {
                    tn++;
                    TCR_x_0[tp][1]=(TCR_x_0[tp][1]*(tn-1.)+TCR_x_00[j])/tn;
                    TCR_y_0[tp][1]=(TCR_y_0[tp][1]*(tn-1.)+TCR_y_00[j])/tn;
                    TCR_z_0[tp][1]=(TCR_z_0[tp][1]*(tn-1.)+TCR_z_00[j])/tn;
                    TCR_p[tp]=-2.8;
                }
            }
            else if(TCR_na[j][0] == 'T' and TCR_na[j][1] == 'R' and TCR_na[j][2] == 'P')
            {
                if(TCR_at[j][1] == 'G' or TCR_at[j][1] == 'D' or TCR_at[j][1] == 'E' or TCR_at[j][1] == 'Z' or TCR_at[j][1] == 'H')
                {
                    tn++;
                    TCR_x_0[tp][1]=(TCR_x_0[tp][1]*(tn-1.)+TCR_x_00[j])/tn;
                    TCR_y_0[tp][1]=(TCR_y_0[tp][1]*(tn-1.)+TCR_y_00[j])/tn;
                    TCR_z_0[tp][1]=(TCR_z_0[tp][1]*(tn-1.)+TCR_z_00[j])/tn;
		    //TCR_p[tp]=1.;
                }
            }
            else if(TCR_na[j][0] == 'T' and TCR_na[j][1] == 'Y' and TCR_na[j][2] == 'R')
            {
                if(TCR_at[j][1] == 'H')
                {
                    TCR_x_0[tp][1]=TCR_x_00[j];
                    TCR_y_0[tp][1]=TCR_y_00[j];
                    TCR_z_0[tp][1]=TCR_z_00[j];
		    //TCR_p[tp]=1.;
                }
            }
            else if(TCR_na[j][0] == 'S' and TCR_na[j][1] == 'E' and TCR_na[j][2] == 'R')
            {
                if(TCR_at[j][1] == 'G')
                {
                    TCR_x_0[tp][1]=TCR_x_00[j];
                    TCR_y_0[tp][1]=TCR_y_00[j];
                    TCR_z_0[tp][1]=TCR_z_00[j];
		    //TCR_p[tp]=1.;
                }
            }
            else if(TCR_na[j][0] == 'C' and TCR_na[j][1] == 'Y' and TCR_na[j][2] == 'S')
            {
                if(TCR_at[j][1] == 'G')
                {
                    TCR_x_0[tp][1]=TCR_x_00[j];
                    TCR_y_0[tp][1]=TCR_y_00[j];
                    TCR_z_0[tp][1]=TCR_z_00[j];
                    TCR_p[tp]=-2.5;
                }
            }
            else if(TCR_na[j][0] == 'T' and TCR_na[j][1] == 'H' and TCR_na[j][2] == 'R')
            {
                if(TCR_at[j][1] == 'G')
                {
                    tn++;
                    TCR_x_0[tp][1]=(TCR_x_0[tp][1]*(tn-1.)+TCR_x_00[j])/tn;
                    TCR_y_0[tp][1]=(TCR_y_0[tp][1]*(tn-1.)+TCR_y_00[j])/tn;
                    TCR_z_0[tp][1]=(TCR_z_0[tp][1]*(tn-1.)+TCR_z_00[j])/tn;
		    //TCR_p[tp]=1.;
                }
            }
            else if(TCR_na[j][0] == 'V' and TCR_na[j][1] == 'A' and TCR_na[j][2] == 'L')
            {
                if(TCR_at[j][1] == 'G')
                {
                    tn++;
                    TCR_x_0[tp][1]=(TCR_x_0[tp][1]*(tn-1.)+TCR_x_00[j])/tn;
                    TCR_y_0[tp][1]=(TCR_y_0[tp][1]*(tn-1.)+TCR_y_00[j])/tn;
                    TCR_z_0[tp][1]=(TCR_z_0[tp][1]*(tn-1.)+TCR_z_00[j])/tn;
                    TCR_p[tp]=-4.2;
                }
            }
            else if(TCR_na[j][0] == 'L' and TCR_na[j][1] == 'E' and TCR_na[j][2] == 'U')
            {
                if(TCR_at[j][1] == 'D')
                {
                    tn++;
                    TCR_x_0[tp][1]=(TCR_x_0[tp][1]*(tn-1.)+TCR_x_00[j])/tn;
                    TCR_y_0[tp][1]=(TCR_y_0[tp][1]*(tn-1.)+TCR_y_00[j])/tn;
                    TCR_z_0[tp][1]=(TCR_z_0[tp][1]*(tn-1.)+TCR_z_00[j])/tn;
                    TCR_p[tp]=-3.8;
                }
            }
            else if(TCR_na[j][0] == 'M' and TCR_na[j][1] == 'E' and TCR_na[j][2] == 'T')
            {
                if(TCR_at[j][1] == 'E')
                {
                    TCR_x_0[tp][1]=TCR_x_00[j];
                    TCR_y_0[tp][1]=TCR_y_00[j];
                    TCR_z_0[tp][1]=TCR_z_00[j];
                    TCR_p[tp]=-1.9;
                }
            }
            else if(TCR_na[j][0] == 'I' and TCR_na[j][1] == 'L' and TCR_na[j][2] == 'E')
            {
                if(TCR_at[j][1] == 'D')
                {
                    TCR_x_0[tp][1]=TCR_x_00[j];
                    TCR_y_0[tp][1]=TCR_y_00[j];
                    TCR_z_0[tp][1]=TCR_z_00[j];
                    TCR_p[tp]=-4.5;
                }
            }
            else if(TCR_na[j][0] == 'A' and TCR_na[j][1] == 'L' and TCR_na[j][2] == 'A')
            {
                if(TCR_at[j][1] == 'B')
                {
                    TCR_x_0[tp][1]=TCR_x_00[j];
                    TCR_y_0[tp][1]=TCR_y_00[j];
                    TCR_z_0[tp][1]=TCR_z_00[j];
                    TCR_p[tp]=-1.8;
                }
            }
            else if(TCR_na[j][0] == 'P' and TCR_na[j][1] == 'R' and TCR_na[j][2] == 'O')
            {
                if(TCR_at[j][1] == 'G' or TCR_at[j][1] == 'D' or TCR_at[j][1] == 'B')
                {
                    tn++;
                    TCR_x_0[tp][1]=(TCR_x_0[tp][1]*(tn-1.)+TCR_x_00[j])/tn;
                    TCR_y_0[tp][1]=(TCR_y_0[tp][1]*(tn-1.)+TCR_y_00[j])/tn;
                    TCR_z_0[tp][1]=(TCR_z_0[tp][1]*(tn-1.)+TCR_z_00[j])/tn;
                    //TCR_p[tp]=0.;
                }
            }
        }
        
        ifi.close();
        
        
        srand(time(NULL)); //srand
        //          check interacting residues
        nof_MHC=chainA;
        nof_TCR=chainD;
        
        nof_pair=0;
        nof_pair = interface_np(chainY, chainX,&TCR_x_0[0][0],&TCR_y_0[0][0],&TCR_z_0[0][0],&TCR_q[0],&TCR_p[0],&MHC_x_0[0][0],&MHC_y_0[0][0],&MHC_z_0[0][0],&MHC_q[0],&MHC_p[0], &pair_list[0][0],&pair_dist[0][0]);
        cout << "Npair " << nof_pair << " " << native_cutoff << endl;
        //          center of mass of TCR and MHC
        cm1_a_x=0;
        cm1_a_y=0;
        cm1_a_z=0;
        cm1_b_x=0;
        cm1_b_y=0;
        cm1_b_z=0;
        for (i=0;i<chainY;i++)
        {
           cm1_a_x+=TCR_x_0[i][0];
           cm1_a_y+=TCR_y_0[i][0];
           cm1_a_z+=TCR_z_0[i][0];
        }
        for (i=0;i<chainX;i++)
        {
           cm1_b_x+=MHC_x_0[i][0];
           cm1_b_y+=MHC_y_0[i][0];
           cm1_b_z+=MHC_z_0[i][0];
        }
        cm1_a_x/=((double)(chainY));
        cm1_a_y/=((double)(chainY));
        cm1_a_z/=((double)(chainY));
        cm1_b_x/=((double)(chainX));
        cm1_b_y/=((double)(chainX));
        cm1_b_z/=((double)(chainX));
        dist_cm1_IJ=distance(cm1_a_x,cm1_a_y,cm1_a_z,cm1_b_x,cm1_b_y,cm1_b_z); //distance of 2 geometry centers
        //          calculate the X, Y and Z axes of TCR and MHC
      
        RB_MHC_x_0[0]=cm1_b_x;
        RB_MHC_y_0[0]=cm1_b_y;
        RB_MHC_z_0[0]=cm1_b_z;
        
        RB_TCR_x_0[0]=cm1_a_x;
        RB_TCR_y_0[0]=cm1_a_y;
        RB_TCR_z_0[0]=cm1_a_z;
        
        RB_MHC_x_0[3]=(cm1_b_x+cm1_a_x)/2.;
        RB_MHC_y_0[3]=(cm1_b_y+cm1_a_y)/2.;
        RB_MHC_z_0[3]=(cm1_b_z+cm1_a_z)/2.;
        
        RB_TCR_x_0[3]=(cm1_b_x+cm1_a_x)/2.;
        RB_TCR_y_0[3]=(cm1_b_y+cm1_a_y)/2.;
        RB_TCR_z_0[3]=(cm1_b_z+cm1_a_z)/2.;
        dist_temp=distance(RB_MHC_x_0[3],RB_MHC_y_0[3],RB_MHC_z_0[3],RB_MHC_x_0[0],RB_MHC_y_0[0],RB_MHC_z_0[0]);
        
        RB_MHC_x_0[1]=cm1_b_x+(RB_MHC_y_0[3]-RB_MHC_y_0[0])*dist_cm1_IJ/(2.*dist_temp);
        RB_MHC_y_0[1]=cm1_b_y+(RB_MHC_x_0[0]-RB_MHC_x_0[3])*dist_cm1_IJ/(2.*dist_temp);
        RB_MHC_z_0[1]=cm1_b_z;
        
        RB_TCR_x_0[2]=cm1_a_x-(RB_MHC_y_0[3]-RB_MHC_y_0[0])*dist_cm1_IJ/(2.*dist_temp);
        RB_TCR_y_0[2]=cm1_a_y-(RB_MHC_x_0[0]-RB_MHC_x_0[3])*dist_cm1_IJ/(2.*dist_temp);
        RB_TCR_z_0[2]=cm1_a_z;
        
        dist_temp=distance(RB_MHC_x_0[1],RB_MHC_y_0[1],RB_MHC_z_0[1],RB_MHC_x_0[0],RB_MHC_y_0[0],RB_MHC_z_0[0]);
        
        vect_x=(RB_MHC_y_0[1]-cm1_b_y)*(RB_MHC_z_0[3]-cm1_b_z)-(RB_MHC_z_0[1]-cm1_b_z)*(RB_MHC_y_0[3]-cm1_b_y);
        vect_y=(RB_MHC_z_0[1]-cm1_b_z)*(RB_MHC_x_0[3]-cm1_b_x)-(RB_MHC_x_0[1]-cm1_b_x)*(RB_MHC_z_0[3]-cm1_b_z);
        vect_z=(RB_MHC_x_0[1]-cm1_b_x)*(RB_MHC_y_0[3]-cm1_b_y)-(RB_MHC_y_0[1]-cm1_b_y)*(RB_MHC_x_0[3]-cm1_b_x);
        dist_vect=sqrt(pow(vect_x,2)+pow(vect_y,2)+pow(vect_z,2));
        
        RB_MHC_x_0[2]=cm1_b_x+vect_x*dist_cm1_IJ/(2.*dist_vect);    
        RB_MHC_y_0[2]=cm1_b_y+vect_y*dist_cm1_IJ/(2.*dist_vect);     
        RB_MHC_z_0[2]=cm1_b_z+vect_z*dist_cm1_IJ/(2.*dist_vect);
        
        RB_TCR_x_0[1]=cm1_a_x+vect_x*dist_cm1_IJ/(2.*dist_vect);      
        RB_TCR_y_0[1]=cm1_a_y+vect_y*dist_cm1_IJ/(2.*dist_vect);      
        RB_TCR_z_0[1]=cm1_a_z+vect_z*dist_cm1_IJ/(2.*dist_vect);    
        
        dist_temp=distance(RB_MHC_x_0[2],RB_MHC_y_0[2],RB_MHC_z_0[2],RB_MHC_x_0[0],RB_MHC_y_0[0],RB_MHC_z_0[0]);
        dist_temp=distance(RB_MHC_x_0[3],RB_MHC_y_0[3],RB_MHC_z_0[3],RB_TCR_x_0[3],RB_TCR_y_0[3],RB_TCR_z_0[3]);
        for(i=0;i<nof_res;i++) //centerize
        {
            RB_MHC_x_0[i]-=(cm1_b_x+cm1_a_x)/2.;
            RB_MHC_y_0[i]-=(cm1_b_y+cm1_a_y)/2.;
            RB_MHC_z_0[i]-=(cm1_b_z+cm1_a_z)/2.;
            RB_TCR_x_0[i]-=(cm1_b_x+cm1_a_x)/2.;
            RB_TCR_y_0[i]-=(cm1_b_y+cm1_a_y)/2.;
            RB_TCR_z_0[i]-=(cm1_b_z+cm1_a_z)/2.;;
        }
        for(i=0;i<chainX;i++)
        {
            MHC_x_0[i][0]-=(cm1_b_x+cm1_a_x)/2.;
            MHC_y_0[i][0]-=(cm1_b_y+cm1_a_y)/2.;
            MHC_z_0[i][0]-=(cm1_b_z+cm1_a_z)/2.;
            MHC_x_0[i][1]-=(cm1_b_x+cm1_a_x)/2.;
            MHC_y_0[i][1]-=(cm1_b_y+cm1_a_y)/2.;
            MHC_z_0[i][1]-=(cm1_b_z+cm1_a_z)/2.;
        }
        for(i=0;i<chainY;i++)
        {
            TCR_x_0[i][0]-=(cm1_b_x+cm1_a_x)/2.;
            TCR_y_0[i][0]-=(cm1_b_y+cm1_a_y)/2.;
            TCR_z_0[i][0]-=(cm1_b_z+cm1_a_z)/2.;
            TCR_x_0[i][1]-=(cm1_b_x+cm1_a_x)/2.;
            TCR_y_0[i][1]-=(cm1_b_y+cm1_a_y)/2.;
            TCR_z_0[i][1]-=(cm1_b_z+cm1_a_z)/2.;;
        }
        
        //      a number of different trials, for each trial, generate initial conformation
        nof_trj=0;
        for (n_trial=0;n_trial<trial_max;n_trial++)//trial_max = 1M, 
        {
            if (nof_trj >= trj_max) break; //trj_max = 100
            
            /*ofstream//sfi(argv[2],ios::out|ios::app);    //io//
            ofstream lfi(argv[3],ios::out|ios::app);   //io//
            ofstream efi(argv[4],ios::out|ios::app);  //io//
            ofstream cfi(argv[5],ios::out|ios::app); //io//
            ofstream tfi(argv[6],ios::out|ios::app);//io//*/
            mem[15]=0; //re-association flag
            mem[16]=0; //association assignment flag
            rr(nof_res, &RB_MHC_x_0[0],&RB_MHC_y_0[0],&RB_MHC_z_0[0],&RB_MHC_x[0],&RB_MHC_y[0],&RB_MHC_z[0]);//recover-native
            rr(chainX*2, &MHC_x_0[0][0],&MHC_y_0[0][0],&MHC_z_0[0][0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0]);//recover-native
            rr(nof_res, &RB_TCR_x_0[0],&RB_TCR_y_0[0],&RB_TCR_z_0[0],&RB_TCR_x[0],&RB_TCR_y[0],&RB_TCR_z[0]);//recover-native
            rr(chainY*2, &TCR_x_0[0][0],&TCR_y_0[0][0],&TCR_z_0[0][0],&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0]);//recover-native
            
            /*ini_TM(nof_res, chainY*2, &RB_TCR_x[0], &RB_TCR_y[0], &RB_TCR_z[0],&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0]);//test-part
	    ini_TM(nof_res, chainX*2, &RB_MHC_x[0], &RB_MHC_y[0], &RB_MHC_z[0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0]);//test-part
            BC(nof_res, chainY*2, cell_range_x, cell_range_y, cell_range_z, &RB_TCR_x[0], &RB_TCR_y[0], &RB_TCR_z[0],&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0]);
            BC(nof_res, chainX*2, cell_range_x, cell_range_y, cell_range_z, &RB_MHC_x[0], &RB_MHC_y[0], &RB_MHC_z[0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0]);    //*/
                    
	    //    the initial deviated TCR have to satisfy constrain
            tp=0;
            /*clash_flag=false;
            for(i=0;i<chainY;i++)
            {
                for (j=0;j<chainX;j++)
                {
                    dist=distance(TCR_x[i][1],TCR_y[i][1],TCR_z[i][1],MHC_x[j][1],MHC_y[j][1],MHC_z[j][1]); 
                    if(dist<=clash_cutoff) //2.0 or 5.0?
                    {
                        clash_flag=true;
                        tp++;
                        break;
                    }
                    dist=distance(TCR_x[i][1],TCR_y[i][1],TCR_z[i][1],MHC_x[j][0],MHC_y[j][0],MHC_z[j][0]); 
                    if(dist<=clash_cutoff) //2.0 or 5.0?
                    {
                        clash_flag=true;
                        tp++;
                        break;
                    }
                    dist=distance(TCR_x[i][0],TCR_y[i][0],TCR_z[i][0],MHC_x[j][1],MHC_y[j][1],MHC_z[j][1]); 
                    if(dist<=clash_cutoff) //2.0 or 5.0?
                    {
                        clash_flag=true;
                        tp++;
                        break;
                    }
                    dist=distance(TCR_x[i][0],TCR_y[i][0],TCR_z[i][0],MHC_x[j][0],MHC_y[j][0],MHC_z[j][0]); 
                    if(dist<=clash_cutoff) //2.0 or 5.0?
                    {
                        clash_flag=true;
                        tp++;
                        break;
                    }
                }
                if (clash_flag) break;
            }
            if (clash_flag) continue; //test-part 
            // with certain distance and angular constrain
            // calculate distance of two reaction sites between TCR and MHC*/
            dist=distance(RB_MHC_x[3],RB_MHC_y[3],RB_MHC_z[3],RB_TCR_x[3],RB_TCR_y[3],RB_TCR_z[3]); 
            if (output_str_flag)
            {
                for(i=0;i<chainX;i++)
                {
                    sfi <<"ATOM"<<setw(7)<<i<<"  "<< "CA  "<<MHC_na[i][0]<<MHC_na[i][1]<<MHC_na[i][2]<<" A"<<setw(4)<<i<<"    "<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<MHC_x[i][0]<<setw(8)<<MHC_y[i][0]<<setw(8)<<MHC_z[i][0]<<"  1.00" << setprecision(2)<<setw(6)<<MHC_i[i]*100.<<endl;
                }
                sfi<<"TER"<<endl;
                
                for(i=0;i<chainY;i++)
                {
                    sfi <<"ATOM"<<setw(7)<<i<<"  "<< "CA  "<<TCR_na[i][0]<<TCR_na[i][1]<<TCR_na[i][2]<<" D"<<setw(4)<<i<<"    "<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<TCR_x_0[i][0]<<setw(8)<<TCR_y_0[i][0]<<setw(8)<<TCR_z_0[i][0]<<"  1.00" << setprecision(2)<<setw(6)<<TCR_i[i]*100.<<endl;
                }
                sfi<<"TER"<<endl;
                sfi<<"END"<<endl;
            }
            /*rr(nof_res, &RB_TCR_x[0],&RB_TCR_y[0],&RB_TCR_z[0],&RB_TCR_x_rec[0],&RB_TCR_y_rec[0],&RB_TCR_z_rec[0]);//record
            rr(chainY*2, &TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0],&TCR_x_rec[0][0],&TCR_y_rec[0][0],&TCR_z_rec[0][0]);//record
            rr(nof_res, &RB_MHC_x[0],&RB_MHC_y[0],&RB_MHC_z[0],&RB_MHC_x_rec[0],&RB_MHC_y_rec[0],&RB_MHC_z_rec[0]);//record
            rr(chainX*2, &MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0],&MHC_x_rec[0][0],&MHC_y_rec[0][0],&MHC_z_rec[0][0]);//record//*/
            rmsd = rmsdcall(chainY*2, chainX*2, &TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0], &TCR_x_0[0][0],&TCR_y_0[0][0],&TCR_z_0[0][0], &MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0], &MHC_x_0[0][0],&MHC_y_0[0][0],&MHC_z_0[0][0]);
            {
                nof_trj++;
                
                current_simu_time=0.0;
                sum_diff=0.; //
                for (mc_time_step=0;mc_time_step < simu_step;mc_time_step++) //max simu_step 10000
                {
                    //   translational movement 
                    if (mem[16])
                    {    
                    }
                    else
                    {
                        /*Prob_Diff=(6*TCR_D*time_step)/pow(distance_step,2); // TCR_D=0.01, time_step=250.0, distance_step=5.0
                        if(random_d(1.0)<Prob_Diff)
                        {
                            translation(nof_res, chainY*2, distance_step, &RB_TCR_x[0], &RB_TCR_y[0], &RB_TCR_z[0],&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0]); //distance_step=5.0//test-part
                        } 
                        Prob_Diff=(6*MHC_D*time_step)/pow(distance_step,2); // TCR_D=0.01, time_step=250.0, distance_step=5.0
                        if(random_d(1.0)<Prob_Diff)
                        {//cc>>     if smaller than the diffusion probability, move alone
                            translation(nof_res, chainX*2, distance_step, &RB_MHC_x[0], &RB_MHC_y[0], &RB_MHC_z[0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0]); //distance_step=5.0//test-part
                        } 
                        //   rotational movement
                        rotation(nof_res, chainY*2, TCR_rot_D, &RB_TCR_x[0], &RB_TCR_y[0], &RB_TCR_z[0], &TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0]); //TCR_rot_D=10.0//test-part
                        rotation(nof_res, chainX*2, MHC_rot_D, &RB_MHC_x[0], &RB_MHC_y[0], &RB_MHC_z[0], &MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0]); //TCR_rot_D=10.0//test-part
                        BC(nof_res, chainY*2, cell_range_x, cell_range_y, cell_range_z, &RB_TCR_x[0], &RB_TCR_y[0], &RB_TCR_z[0],&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0]);
                        BC(nof_res, chainX*2, cell_range_x, cell_range_y, cell_range_z, &RB_MHC_x[0], &RB_MHC_y[0], &RB_MHC_z[0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0]);    //*/
                        sum_diff2=0.;    
                        ego = golike(chainY, chainX,&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0],&pair_list[0][0],&pair_dist[0][0]);// ori(chainY*2, chainX*2,&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0],&pair_list[0][0]); //
                        eele = Eele(chainY, chainX,&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0],&TCR_q[0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0],&MHC_q[0],&pair_list[0][0]); //k0
                        ehp  =   HP(chainY, chainX,&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0],&TCR_p[0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0],&MHC_p[0],&pair_list[0][0]);
                        eclash = ncl(chainY*2, chainX*2,&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0]);
                        sum_diff2 = wel * eele + whp * ehp + eclash; //k0 eele + ehp + eclash;
                        dist=distance(RB_MHC_x[3],RB_MHC_y[3],RB_MHC_z[3],RB_TCR_x[3],RB_TCR_y[3],RB_TCR_z[3]); 
                                    
                        // calculate two angle
                        // perpendicular angle
                        theta_pd=0.;
                        vector[0][0]=RB_TCR_x[3]-RB_TCR_x[0]; 
                        vector[0][1]=RB_TCR_y[3]-RB_TCR_y[0];
                        vector[0][2]=RB_TCR_z[3]-RB_TCR_z[0];
                        vector[1][0]=RB_MHC_x[0]-RB_MHC_x[3]; 
                        vector[1][1]=RB_MHC_y[0]-RB_MHC_y[3]; 
                        vector[1][2]=RB_MHC_z[0]-RB_MHC_z[3]; 
                        theta_pd=dot2theta(&vector[0][0], &vector[1][0]);
                        // orientational angle
                        theta_ot=0.;
                        vector[0][0]=RB_TCR_x[1]-RB_TCR_x[0]; 
                        vector[0][1]=RB_TCR_y[1]-RB_TCR_y[0];
                        vector[0][2]=RB_TCR_z[1]-RB_TCR_z[0];
                        vector[1][0]=RB_MHC_x[0]-RB_MHC_x[1];
                        vector[1][1]=RB_MHC_y[0]-RB_MHC_y[1];
                        vector[1][2]=RB_MHC_z[0]-RB_MHC_z[1];
                        theta_ot=dot2theta(&vector[0][0], &vector[1][0]);
                        
                        if (1)//((sum_diff2 <= sum_diff) or ((eclash < 1.) and (exp(-(sum_diff2-sum_diff)/simu_temp) > random_d(1.0))))//((sum_diff2+0.1 <= sum_diff) or ( sum_diff2 <= sum_diff+5. and (0.02 > random_d(1.0)))) //((sum_diff2 <= sum_diff) or ((sum_diff2 > sum_diff) and (exp(-(sum_diff2-sum_diff)/simu_temp) > random_d(1.0))))// (sum_diff2-sum_diff)/10. > random_d(1.0))) //
                        {
                            sum_diff=sum_diff2;
                            current_simu_time+=time_step/1000; //time_step=250.0 ns
                            rr(nof_res, &RB_TCR_x[0],&RB_TCR_y[0],&RB_TCR_z[0],&RB_TCR_x_rec[0],&RB_TCR_y_rec[0],&RB_TCR_z_rec[0]);//record
                            rr(chainY*2, &TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0],&TCR_x_rec[0][0],&TCR_y_rec[0][0],&TCR_z_rec[0][0]);//record
                            rr(nof_res, &RB_MHC_x[0],&RB_MHC_y[0],&RB_MHC_z[0],&RB_MHC_x_rec[0],&RB_MHC_y_rec[0],&RB_MHC_z_rec[0]);//record
                            rr(chainX*2, &MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0],&MHC_x_rec[0][0],&MHC_y_rec[0][0],&MHC_z_rec[0][0]);//record
                            rmsd = rmsdcall(chainY*2, chainX*2, &TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0], &TCR_x_0[0][0],&TCR_y_0[0][0],&TCR_z_0[0][0], &MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0], &MHC_x_0[0][0],&MHC_y_0[0][0],&MHC_z_0[0][0]);
                        }
                        else
                        {
                            rr(nof_res, &RB_TCR_x_rec[0],&RB_TCR_y_rec[0],&RB_TCR_z_rec[0],&RB_TCR_x[0],&RB_TCR_y[0],&RB_TCR_z[0]);//recover
                            rr(chainY*2, &TCR_x_rec[0][0],&TCR_y_rec[0][0],&TCR_z_rec[0][0],&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0]);//recover
                            rr(nof_res, &RB_MHC_x_rec[0],&RB_MHC_y_rec[0],&RB_MHC_z_rec[0],&RB_MHC_x[0],&RB_MHC_y[0],&RB_MHC_z[0]);//recover
                            rr(chainX*2, &MHC_x_rec[0][0],&MHC_y_rec[0][0],&MHC_z_rec[0][0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0]);//recover
                            comp_rmsd = false;
                        }
                    }
                    if(output_str_flag)
		    {
                        if((mc_time_step % output_str_step) == 0) //10 steps /output conformation
                        {
                            for(i=0;i<chainX;i++)
                            {
                                sfi <<"ATOM"<<setw(7)<<i<<"  "<< "CA  "<<MHC_na[i][0]<<MHC_na[i][1]<<MHC_na[i][2]<<" A"<<setw(4)<<i<<"    "<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<MHC_x[i][0]<<setw(8)<<MHC_y[i][0]<<setw(8)<<MHC_z[i][0]<<"  1.00" << setprecision(2)<<setw(6)<<MHC_i[i]*100.<<endl;
                            }
                            sfi<<"TER"<<endl;
                            for(i=0;i<chainY;i++)
                            {
                                sfi <<"ATOM"<<setw(7)<<i<<"  "<<"CA  "<<TCR_na[i][0]<<TCR_na[i][1]<<TCR_na[i][2]<<" D"<<setw(4)<<i<<"    "<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<TCR_x[i][0]<<setw(8)<<TCR_y[i][0]<<setw(8)<<TCR_z[i][0]<<"  1.00" << setprecision(2)<<setw(6)<<TCR_i[i]*100.<<endl;
                            }
                            sfi<<"TER"<<endl;
                            sfi<<"END"<<endl;
                        }
                    }
                    
		    {
                        {
                            mem[11]++;
                            if (mem[15] < 1) mem[13]++;
                            if (mem[16])//association marker
                            {
                            }    
                            else if (sum_diff == sum_diff2)
                            {
                                trj[0]=rmsd;
                                trj[1]=ego;
                                trj[2]=ehp;
                                trj[3]=eele;
                                trj[4]=eclash;
                                trj[5]=dist;
                                trj[6]=theta_pd;
                                trj[7]=theta_ot;
                                
                                mem[12]++;
                                if (mem[15] < 1)
                                {
                                    mem[14]++;
                                    if (ego < native_cutoff)//-9)
                                    {
                                        mem[10]++; //
                                        mem[15]++; //
                                       
                                    }
                                }
                                if (mem[16])//association marker
                                {
                                }
                                else
                                {
                                    if(ego < native_cutoff and trj[0] < native_rmsd) //-9 parameters
                                    {
                                        mem[16]=1;
                                        rr(nof_res, &RB_MHC_x_0[0],&RB_MHC_y_0[0],&RB_MHC_z_0[0],&RB_MHC_x[0],&RB_MHC_y[0],&RB_MHC_z[0]);//recover-native
                                        rr(chainX*2, &MHC_x_0[0][0],&MHC_y_0[0][0],&MHC_z_0[0][0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0]);//recover-native
                                        rr(nof_res, &RB_TCR_x_0[0],&RB_TCR_y_0[0],&RB_TCR_z_0[0],&RB_TCR_x[0],&RB_TCR_y[0],&RB_TCR_z[0]);//recover-native
                                        rr(chainY*2, &TCR_x_0[0][0],&TCR_y_0[0][0],&TCR_z_0[0][0],&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0]);//recover-native
                                        trj[1] = golike(chainY, chainX,&TCR_x[0][0],&TCR_y[0][0],&TCR_z[0][0],&MHC_x[0][0],&MHC_y[0][0],&MHC_z[0][0],&pair_list[0][0],&pair_dist[0][0]);
                                    }
                                }
                            }
                            if(output_ene_flag and ((mc_time_step % output_ene_step) == 0))
                            {
                                ofstream efi(argv[4],ios::out|ios::app);
                                efi << " " << mc_time_step << " " << sum_diff << " " << sum_diff2 << " " << trj[5] << " " << trj[6] << " " << trj[7] << " " << trj[0] << " " << trj[1] << " " << trj[2] << " " << trj[3] << " " << trj[4] << " " << mem[10] << endl;
                                efi.close();
                            }
                            
                            if(mc_time_step == (simu_step-1))
                            {
                                if(ego < 0)
                                {
                                    mem[0]++;
                                    if(ego < -1)
                                    {
                                        mem[1]++;
                                        if(ego < -2)
                                        {
                                            mem[2]++;
                                            if(ego < -3)
                                            {
                                                mem[3]++;
                                                if(ego < -4)
                                                {
                                                    mem[4]++;
                                                }
                                            }
                                        }
                                    }
                                }
                                ofstream ofi(argv[2],ios::out|ios::app);
                                ofi<<nof_trj<<" "<<mem[10]<<" "<<mem[0]<<" "<<mem[1]<<" "<<mem[2]<<" "<<mem[3]<<" "<<mem[4]<<" "<<mem[11]<<" "<<mem[12]<<" "<<mem[13]<<" "<<mem[14]<<" "<<trj[0]<<" "<<trj[1]<< " " <<trj[2]<< " " <<trj[3]<< " " <<trj[4]<< endl;
                                ofi.close();
                            }
                        }
                    }
                }//end for mc_time
            }//end if
        }//end for nof_trial
        sfi.close();
    }
    
    
    void ini_TM(int nof_res, int nof_tcr, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3)
    {//c>>>>>   random position
        int k=0;
	
        rand_x=random_d(2.)*(dist_constrain)-dist_constrain; //25.
        rand_y=random_d(2.)*(dist_constrain)-dist_constrain;
        rand_z=random_d(2.)*(dist_constrain)-dist_constrain;
        
        for (k=0;k<nof_res;k++)
        {
            *(p1+k)+=rand_x;
            *(p2+k)+=rand_y;
            *(p3+k)+=rand_z;
        }
        for (k=0;k<nof_tcr;k++)
        {
            *(q1+k)+=rand_x;
            *(q2+k)+=rand_y;
            *(q3+k)+=rand_z;
        }    
        //c>>>>>   random orientation
        rotation(nof_res,nof_tcr,180.,p1,p2,p3,q1,q2,q3);   
    }
    void translation(int nof_res, int nof_tcr, double maxdist, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3) //translate single
    {
        int k=0;
        double theta=0.0, phi=0.0;
        double dist = random_d(1.0)*maxdist*2.; //2014/3/6
        theta=random_d(1.0)*PI;
        phi=random_d(2.0)*PI;
        
        for(k=0;k<nof_tcr;k++)
        {
            *(q1+k)+=dist*sin(theta)*cos(phi);
            *(q2+k)+=dist*sin(theta)*sin(phi);
            *(q3+k)+=dist*cos(theta);
        }
        for (k=0;k<nof_res;k++)
        {
            *(p1+k)+=dist*sin(theta)*cos(phi);
            *(p2+k)+=dist*sin(theta)*sin(phi);
            *(p3+k)+=dist*cos(theta);
        }
        //cin >> collision;
    }//end func translation	
    void rotation (int nof_res, int nof_tcr, double upper_limit, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3) //rotate single
    {
        int k=0;
        double theta=0.0,phi=0.0,psi=0.0;
        double ang[3][3];
        double tmp[3][nof_tcr];
        
        theta=(random_d(2.0)-1.)*PI*upper_limit/180.;
        phi=  (random_d(2.0)-1.)*PI*upper_limit/180.;
        psi=  (random_d(2.0)-1.)*PI*upper_limit/180.;
        
        ang[0][0]=cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi);
        ang[0][1]=-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi);
        ang[0][2]=sin(theta)*sin(phi);
        
        ang[1][0]=cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi);
        ang[1][1]=-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi);
        ang[1][2]=-sin(theta)*cos(phi);
        
        ang[2][0]=sin(psi)*sin(theta);
        ang[2][1]=cos(psi)*sin(theta);
        ang[2][2]=cos(theta);
        
        for (k=0;k<nof_res;k++)
        {
            tmp[0][k]=*(p1+k);
            tmp[1][k]=*(p2+k);
            tmp[2][k]=*(p3+k);
        }
        for (k=1;k<nof_res;k++)
        {
            *(p1+k)=ang[0][0]*(tmp[0][k]-*p1)+ang[0][1]*(tmp[1][k]-*p2)+ang[0][2]*(tmp[2][k]-*p3)+*p1;
            *(p2+k)=ang[1][0]*(tmp[0][k]-*p1)+ang[1][1]*(tmp[1][k]-*p2)+ang[1][2]*(tmp[2][k]-*p3)+*p2;
            *(p3+k)=ang[2][0]*(tmp[0][k]-*p1)+ang[2][1]*(tmp[1][k]-*p2)+ang[2][2]*(tmp[2][k]-*p3)+*p3;
        }
        for (k=0;k<nof_tcr;k++)
        {
            tmp[0][k]=*(q1+k);
            tmp[1][k]=*(q2+k);
            tmp[2][k]=*(q3+k);
        }
        for (k=0;k<nof_tcr;k++)
        {
            *(q1+k)=ang[0][0]*(tmp[0][k]-*p1)+ang[0][1]*(tmp[1][k]-*p2)+ang[0][2]*(tmp[2][k]-*p3)+*p1;
            *(q2+k)=ang[1][0]*(tmp[0][k]-*p1)+ang[1][1]*(tmp[1][k]-*p2)+ang[1][2]*(tmp[2][k]-*p3)+*p2;
            *(q3+k)=ang[2][0]*(tmp[0][k]-*p1)+ang[2][1]*(tmp[1][k]-*p2)+ang[2][2]*(tmp[2][k]-*p3)+*p3;
        }                                                
        //cin >> collision;
    }// end func rotation
    void rr (int nof_p, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3) //record-recover
    {
        int k=0;
        for (k=0;k<nof_p;k++)
        {
            *(q1+k)=*(p1+k);
            *(q2+k)=*(p2+k);
            *(q3+k)=*(p3+k);
        }
    }
    int np(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r, double *s) //nof pair  
    {   
        int i=0, j=0, k=0, np=0;
        double dist = 0.;
        for (i=0; i<nof_mhc;i++)
        {
            for(j=0;j<nof_tcr;j++)
            {
                *(r+i*nof_tcr+j)=false;
                *(s+i*nof_tcr+j)=0.;
            }
        }
        for(i=0;i<nof_mhc;i++)
        {
            for(j=0;j<nof_tcr;j++)
            {
                dist=distance(*(p1+j),*(p2+j),*(p3+j),*(q1+i),*(q2+i),*(q3+i));
                *(s+i*nof_tcr+j)=dist;
                if(dist <= dist_cutoff) //5.0 //7.5
                {
                    np++;
                    *(r+i*nof_tcr+j)=true;
                }
            }
        }
        return np;
    }//end func np
    int interface_np(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *p4, double *p5, double *q1, double *q2, double *q3, double *q4, double *q5, bool *r, double *s) //int noa_tcr, int noa_mhc, int nor_tcr, int nor_mhc, double *p1, double *p2, double *p3, int *p4, double *p5, double *p6, double *p7, double *q1, double *q2, double *q3, int *q4, double *q5, double *q6, double *q7, bool *r, double *s) //p1,2,3:x,y,z p4 res no. //p4,5,6,7:p(polar),q(ele),n(serial),i(inter)
    {
        int i=0, j=0, x=0, y=0, nm=0, nt=0, np=0, tmp=0;
        double dist = 0., scy=0., qi = 0., qj=0., vep = 0.00014163, Ds = 10.;
        for (i=0;i<nof_mhc;i++)
        {
            if (*(q4+i) == 0) continue;
            nm = i;
            for (j=0;j<nof_tcr;j++)
            {
                if (*(p4+j) == 0) continue;
                nt=j;
                dist=distance(*(p1+j*2+1),*(p2+j*2+1),*(p3+j*2+1),*(q1+i*2+1),*(q2+i*2+1),*(q3+i*2+1));  //sqrt((TCR_x_new(i)-MHC_x(j))**2+(TCR_y_new(i)-MHC_y(j))**2+(TCR_z_new(i)-MHC_z(j))**2)
                scy=(*(p4+j) * *(q4+i))/(4*PI*vep*(Ds*exp(dist/xi))*dist);
                if (scy < -1. and dist < 10.)//nat-1.0)//(scy < -2.0)
                {
                    np++;
                    *(r+nm*nof_tcr+nt)=true;
                    *(s+nm*nof_tcr+nt)=dist;
                    cout <<nt<<" "<<nm<<" "<< "dist--> " << dist << " scy: " << scy << endl;
                }
            }
        }
        for (i=0;i<nof_mhc;i++)
        {
            if (*(q5+i) == 0) continue;
            nm = i;
            for (j=0;j<nof_tcr;j++)
            {
                if (*(p5+j) == 0) continue;
                nt=j;
                dist=distance(*(p1+j*2+1),*(p2+j*2+1),*(p3+j*2+1),*(q1+i*2+1),*(q2+i*2+1),*(q3+i*2+1));  //sqrt((TCR_x_new(i)-MHC_x(j))**2+(TCR_y_new(i)-MHC_y(j))**2+(TCR_z_new(i)-MHC_z(j))**2)
                scy = *(q5+i) + *(p5+j);
                if (dist < 6) 
                {
                    np++;
                    *(r+nm*nof_tcr+nt)=true;
                    *(s+nm*nof_tcr+nt)=dist;
                    cout <<nt<<" "<<nm<<" "<< "dist--> " << dist << " scy: " << scy << endl;
                }
            }
        }//no hydrophobic count */
        cout << "nm=" << nm << ";nt=" << nt <<";" <<endl;
        return np;
    }
    
    double Eele(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *p4, double *q1, double *q2, double *q3, double *q4, bool *r) //Electrostatic-like potential    
    {
        int i=0, j=0, k=0;
        double dist = 0., tscy = 0., scy=0., qi = 0., qj=0., vep = 0.00014163, Ds = 10.; //, xi = 8.2; //vep = 8.854e-12 F/m 
        for (i=0;i<nof_mhc;i++)
        {
            if (*(q4+i) == 0) continue;
            for (j=0;j<nof_tcr;j++)
            {
                if (!*(r+i*nof_tcr+j)) continue; //All charge out of surface
                if (*(p4+j) == 0) continue;
                dist=distance(*(p1+j*2+1),*(p2+j*2+1),*(p3+j*2+1),*(q1+i*2+1),*(q2+i*2+1),*(q3+i*2+1));  //sqrt((TCR_x_new(i)-MHC_x(j))**2+(TCR_y_new(i)-MHC_y(j))**2+(TCR_z_new(i)-MHC_z(j))**2);
                scy+=(*(p4+j) * *(q4+i))/(4*PI*vep*(Ds*exp(dist/xi))*dist);
            }
        }
        return scy;
    }//end func Eele
    double HP(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *p4, double *q1, double *q2, double *q3, double *q4, bool *r) //account for hydrophobic residues  
    {
        int i=0, j=0, k=0;
        double dist = 0., hp=0.; 
        // hydrophobic interaction
        for (i=0;i<nof_mhc;i++)
        {
            for (j=0;j<nof_tcr;j++)
            {
                if (*(q4+i) == 0 & *(p4+j) == 0) continue;
                dist=distance(*(p1+j*2+1),*(p2+j*2+1),*(p3+j*2+1),*(q1+i*2+1),*(q2+i*2+1),*(q3+i*2+1));  //sqrt((TCR_x_new(i)-MHC_x(j))**2+(TCR_y_new(i)-MHC_y(j))**2+(TCR_z_new(i)-MHC_z(j))**2)
                if (dist < 6) hp+=(*(p4+j) + *(q4+i)); //7.5
            }
        }
        return hp;
    }//end func HP
    double golike(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r, double *s) //GO-like potential    
    {
        int i=0, j=0, k=0;
        double dist = 0., go_like=0.;
        for (i=0;i<nof_mhc;i++)
        {
            for (j=0;j<nof_tcr;j++)
            {
                if(*(r+i*nof_tcr+j)) //ori connection
                {
                    dist=distance(*(p1+j*2+1),*(p2+j*2+1),*(p3+j*2+1),*(q1+i*2+1),*(q2+i*2+1),*(q3+i*2+1));  //sqrt((TCR_x_new(i)-MHC_x(j))**2+(TCR_y_new(i)-MHC_y(j))**2+(TCR_z_new(i)-MHC_z(j))**2)
                    if(dist < (*(s+i*nof_tcr+j)+well_width)) //((dist > (*(s+i*nof_tcr+j)-well_width)) and (dist < (*(s+i*nof_tcr+j)+well_width)))//well_width = 1.
                    {
                        go_like+=-1.; //go_ene; //-5.
                    }//endif
                }
            }
        }
        return go_like;
    }//end function golike
    double ncl(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3) //nof clash  
    {
        int i=0, j=0, k=0;
        double dist = 0., p_clash=0., rmin=0., p_temp=0.;
        for (i=0;i<nof_mhc;i++)
        {
            for (j=0;j<nof_tcr;j++)
            {
                dist=distance(*(p1+j),*(p2+j),*(p3+j),*(q1+i),*(q2+i),*(q3+i));  //sqrt((TCR_x_new(i)-MHC_x(j))**2+(TCR_y_new(i)-MHC_y(j))**2+(TCR_z_new(i)-MHC_z(j))**2)
                if(i%2==0 and j%2==0)
                {
                    rmin=4.3; //5.//5.3;
                }
                else if(i%2==1 and j%2==1)
                {
                    rmin=2.5;
                }
                else
                {
                    rmin=3.1; //3.8 //4.;
                }
                p_temp=-1.*go_ene*(pow((rmin/dist),12)-2*pow((rmin/dist),6));
                if (p_temp>0) p_clash+=p_temp;
            }
        }
        return p_clash;
    }//end func ncl
    double ljlike(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r, double *s) //GO-like potential    
    {
        int i=0, j=0, k=0;
        double dist = 0., lj_like=0., rmin=2.9;
        for (i=0;i<nof_mhc;i++)
        {
            for (j=0;j<nof_tcr;j++)
            {
                dist=distance(*(p1+j),*(p2+j),*(p3+j),*(q1+i),*(q2+i),*(q3+i));  //sqrt((TCR_x_new(i)-MHC_x(j))**2+(TCR_y_new(i)-MHC_y(j))**2+(TCR_z_new(i)-MHC_z(j))**2)
                lj_like-=go_ene*(pow((rmin/dist),12)-2*pow((rmin/dist),6));
            }
        }
        return lj_like;
    }
    double ljgo(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r, double *s) //GO-like potential    
    {
        int i=0, j=0, k=0;
        double dist = 0., lj_like=0., rmin=2.9;
        for (i=0;i<nof_mhc;i++)
        {
            for (j=0;j<nof_tcr;j++)
            {
                dist=distance(*(p1+j),*(p2+j),*(p3+j),*(q1+i),*(q2+i),*(q3+i));  //sqrt((TCR_x_new(i)-MHC_x(j))**2+(TCR_y_new(i)-MHC_y(j))**2+(TCR_z_new(i)-MHC_z(j))**2)
                if(dist < clash_cutoff) //2.0 or 5.0?
                {
                    lj_like-=go_ene*(pow((rmin/dist),12)-2*pow((rmin/dist),6)); //20
                }//endif
                else if(*(r+i*nof_tcr+j)) //ori connection
                {
                    lj_like-=go_ene*(pow((rmin/dist),12)-2*pow((rmin/dist),6));
                }//*/
            }
        }
        return lj_like;
    }
    double dppdlike(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r, double *s) //GO-like potential    
    {
        int i=0, j=0, k=0;
        double dist = 0., d_like=0., tmp_ene = 0.;
        for (i=0;i<nof_mhc;i++)
        {
            for (j=0;j<nof_tcr;j++)
            {
                dist=distance(*(p1+j),*(p2+j),*(p3+j),*(q1+i),*(q2+i),*(q3+i));  //sqrt((TCR_x_new(i)-MHC_x(j))**2+(TCR_y_new(i)-MHC_y(j))**2+(TCR_z_new(i)-MHC_z(j))**2)
                if(dist < clash_cutoff) //2.0 or 5.0?
                {
                    d_like+=pow((dist-clash_cutoff),2)*100.; //20 
                }//endif
                else if(*(r+i*nof_tcr+j)) //ori connection
                {
                    tmp_ene=go_ene+pow((dist-*(s+i*nof_tcr+j)),2); //-5.
                    if (tmp_ene<0) d_like+=tmp_ene;  
                }
            }
        }
        return d_like;
    }
    double ori(int nof_tcr, int nof_mhc, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, bool *r) //GO-like potential    
    {
        int i=0, j=0, k=0;
        double dist = 0., go_like=0.;
        for (i=0;i<nof_mhc;i++)
        {
            for (j=0;j<nof_tcr;j++)
            {
                dist=distance(*(p1+j),*(p2+j),*(p3+j),*(q1+i),*(q2+i),*(q3+i));  //sqrt((TCR_x_new(i)-MHC_x(j))**2+(TCR_y_new(i)-MHC_y(j))**2+(TCR_z_new(i)-MHC_z(j))**2)
                if(dist < clash_cutoff) //2.0 or 5.0? //
                {
                    go_like+=0.01; //clash_ene; //20
                }//endif
                if(*(r+i*nof_tcr+j)) //ori connection
                {
                    if((dist > clash_cutoff) and (dist < dist_cutoff))//5.0 //well_width = 1.
                    {
                        go_like+=-1.; //go_ene; //-5.
                    }//endif
                }
            }
        }
        return go_like;
    }//end function ori
    double distance(double a, double b, double c, double x, double y, double z)
    {
        double dist=0.0;
        dist = sqrt(pow((a-x),2)+pow((b-y),2)+pow((c-z),2));
        return dist;
    }
    int random_i(int upper_limit)
    {
        int r = int((double) rand()/RAND_MAX*upper_limit);
        return r;
    }
    double random_d(double upper_limit)
    {
        double r = (double) rand()/RAND_MAX*upper_limit;
        return r;
    }
    double dp(double *p, double *q)
    {
        double dot_product=0.0;
        dot_product = (*p * *q + *(p+1) * *(q+1) + *(p+2) * *(q+2));
        return dot_product;
    }
    double dot2theta(double *p, double *q)
    {
        double theta=0.0;
        theta = acos( -1.0*dp(p,q)/(distance(*p, *(p+1), *(p+2), 0., 0., 0.)*distance(*q, *(q+1), *(q+2), 0., 0., 0.)))*180./PI;
        return theta;
    }
    double rmsdcall(int x, int y, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3, double *r1, double *r2, double *r3, double *s1, double *s2, double *s3)
    {
        int c=x+y;
        int atompr[c][2];
        double rmsd=0.0;
        double a1[c],a2[c],a3[c],b1[c],b2[c],b3[c];
        c=0;
        for (i=0;i<x;i++)
        {
            a1[c]=*(p1+i);
            a2[c]=*(p2+i);
            a3[c]=*(p3+i);
            b1[c]=*(q1+i);
            b2[c]=*(q2+i);
            b3[c]=*(q3+i);
            c++;
        }
        for (i=0;i<y;i++)
        {
            a1[c]=*(r1+i);
            a2[c]=*(r2+i);
            a3[c]=*(r3+i);
            b1[c]=*(s1+i);
            b2[c]=*(s2+i);
            b3[c]=*(s3+i);
            c++;
        }
        //for calling superimpose using least square
        for (i=0;i<c;i++)
        {
            atompr[i][0]=i+1;
            atompr[i][1]=i+1;
        }
        rotlsq_(&a1[0],&a2[0],&a3[0], c, &b1[0],&b2[0],&b3[0], c, &atompr[0][0], c);
        
        rmsd=rmsdcompute(c, &a1[0],&a2[0],&a3[0], &b1[0],&b2[0],&b3[0]);
        return rmsd;
    }
    double rmsdcompute(int n, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3)
    {
        int i;
        double rmsd=0.0;
        for (i=0;i<n;i++)
        {
            rmsd+=pow((*(p1+i)-*(q1+i)),2)+pow((*(p2+i)-*(q2+i)),2)+pow((*(p3+i)-*(q3+i)),2);
        }
        rmsd=sqrt(rmsd/(double)n);
        return rmsd;
    }
    
    void BC(int nof_res, int nof_res1, double cell_range_x, double cell_range_y, double cell_range_z, double *p1, double *p2, double *p3, double *q1, double *q2, double *q3)
    {
        int i=0;
        double PB_x=0.0, PB_y=0.0, PB_z=0.0;
        PB_x=cell_range_x*(double) int(*p1/cell_range_x*2.);//here may cause error in the future
        PB_y=cell_range_y*(double) int(*p2/cell_range_y*2.);
        PB_z=cell_range_z*(double) int(*p3/cell_range_z*2.);
        for(i=0;i<nof_res;i++)
        {
            *(p1+i)-=PB_x;
            *(p2+i)-=PB_y;
            *(p3+i)-=PB_z;
            //if (PB != 0) cin >> i;
        }
        for(i=0;i<nof_res1;i++)
        {
            *(q1+i)-=PB_x;
            *(q2+i)-=PB_y;
            *(q3+i)-=PB_z;
            //if (PB != 0) cin >> i;
        }
        
    }
