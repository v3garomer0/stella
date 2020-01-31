
#include "/home/gustavo/Documents/orsayData/libPerso.h"

const Int_t n_run = 10;
const Double_t eneBeam = 9.6; // MeV; nominal in 10.04


const Int_t n_ch_s3 = 24; // number of rings
Double_t mu_s3[n_ch_s3]; // mean strip position
Double_t mu_s3b[n_ch_s3], mu_low_s3b[n_ch_s3], mu_high_s3b[n_ch_s3]; // limits of strip means
Double_t bin_s3b[2*n_ch_s3]; // rad
const Double_t ene_low = 0.; // qdc chn
const Double_t ene_high = 15.;
const Int_t ene_bin = (ene_high - ene_low)*50;
const Double_t time_low = -420.; // ns
const Double_t time_high = 420.;
const Int_t time_bin = (time_high - time_low)/10;
const Int_t chn_s3b[25] = {76, 77, 78, 79, 80, 81, 82, 83, 85, 84, // converts plugs into channels (before January 2017) Jean Nippert
			   87, 86, 88, 89, 90, 91, 92, 93, 94, 95,
			   97, 96, 99, 98, 100};

const Double_t z_s3b = -5.6; // cm; z position
const Double_t r_i_s3 = 2.2/2.;
const Double_t r_a_s3 = 7.0/2.;
const Double_t d_s3 = 0.0886; // cm, strip pitch
const Double_t offs_s3 = 0.1368; // cm, center sensitive area from strip pitch in active area from the manual
const Double_t s_s3 = 0.001; // close gap between strips in angular distribution plot
const Double_t thetaMin_s3b = 180. + TMath::ATan (r_i_s3/z_s3b)*180./TMath::Pi (); // rad; range of kinematics lines; could also calc from bin_s3b[2*n_ch_s3] ..
const Double_t thetaMax_s3b = 180. + TMath::ATan (r_a_s3/z_s3b)*180./TMath::Pi ();

TH2F *h_polVsEne, *h_detVsDet;
TH1F *h_ene[25], *h_time[25];


const Int_t nPoints = 100; // .. in TGraph
TF1 *f_eneVsPol_a[8], *f_eneVsPol_p[11];
TGraph *g_polVsEne_a[8], *g_polVsEne_p[11];
const Double_t A1 = 12., A2 = 12., A3_p = 1.007825, A4_p = 22.989770, A3_a = 4.002603, A4_a = 19.992440; // AMU
const Double_t m = 931.5; // AMU -> MeV
const Double_t eneEx_a[8] = {0.0, 1.633674, 4.2477, 4.96651, 5.6214, 5.7877, 6.725, 7.004}; // MeV
const Double_t eneEx_p[11] = {0.0, 0.43999, 2.076011, 2.390732, 2.63985, 2.703500, 2.982060, 3.67760, 3.848069, 3.914239, 4.42964};

const Int_t nGran = 10e3; // granilarty scanning dead layers
const Double_t alu_density = 2.7; // g/cm3
const Double_t thick_bwd = 0.8 + 0.1; // µm   CHANGE THICKNESS OF FOIL HERE
const Double_t alu_bwd = thick_bwd*TMath::Power (10, -4)*alu_density; // surface density; g/cm2
const Double_t si_density = 2.329; // g/cm3
const Double_t thick_si = 0.5; // µm   same fwd and bwd according to Tuff thesis  
const Double_t si = thick_si*TMath::Power (10, -4)*si_density; // surface density; g/cm2

// energy loss parameters in aluminum foil from Guillaume
const Double_t ene_4alu_a[9] = {12.8, 9.8, 7.5, 5.9, 4.9, 3.7, 2.8, 2.0, 1.1}; // MeV
const Double_t parA_4alu_a[9] = {2011.9297231, 1961.7222252, 1890.5511761, 1824.2919473, 1766.45444027, 1672.35226723, 1474.69408107, 1349.07285, 27.8030303};
const Double_t parB_4alu_a[9] = {-0.72702125, -0.71707839, -0.7009390, -0.6833843, -0.66558456, -0.63033348, -0.53273324, -0.448489, -326.6075757576};
const Double_t parC_4alu_a = 1527.3412121212; // > 1.1 MeV
const Double_t ene_4si_a[8] = {12.5, 10.0, 7.0, 5.0, 3.5, 2.5, 1.5, 1.1}; // MeV
const Double_t parA_4si_a[8] = {2048.867961, 1990.501841, 1911.1161, 1817.25965, 1706.723733, 1497.2076761, 1361.289351, 1305.263};
const Double_t parB_4si_a[8] = {-0.725248, -0.7137736, -0.69588, -0.669928, -0.630125, -0.62242, -0.52441, -0.3134};
const Double_t ene_4alu_p[6] = {7.5, 5.5, 4.0, 2.35, 1.35, 0.8}; // MeV
const Double_t parA_4alu_p[6] = {195.76893592, 191.22, 184.577826565, 180.826916, 174.608207, 171.7649857};
const Double_t parB_4alu_p[6] = {-0.763213525, -0.751669, -0.730959, -0.714861, -0.6731153, -0.6180668};
const Double_t ene_4si_p[7] = {9.5, 7.5, 5.0, 3.0, 1.3, 0.8, 0.5}; // MeV
const Double_t parA_4si_p[7] = {202.220141, 199.028157, 193.736956, 186.793941, 178.973525, 175.269, 180.0021557};
const Double_t parB_4si_p[7] = {-0.766856, -0.7597, -0.74609, -0.7231, -0.682159, -0.62242, -0.517836};

TGraphErrors *gr_res_a, *gr_res_p;
TF1 *f_res_a, *f_res_p;
