

#include "event_p_s3f.h"


Double_t my_resS3 (Double_t *var, Double_t *par){
  // response function for effective relative energy S3 (one sigma)
  // reproduce the trend for the lowest energies
  // [1]: including effective interactions, fano factor etc.
  
  //  return par[0]/TMath::Sqrt (var[0]); // textbook solution with intrinsic efficiency
  return par[0] + par[1]/TMath::Power (var[0], 3); // reproduce alpha spec
}


Double_t ene_alpha (Double_t *var, Double_t *par){
  
  // par[0] = eneCms
  // par[1] = gam3
  // par[2] = alu thickness
  // par[3] = si thickness
  // var[0] = theta lab
  
  Double_t K = A1*A3_a*par[0]/(A2*(A1 + A2)*TMath::Power (par[1], 2));
  Double_t eneKin_a = 0.;
  if (2*par[1]*TMath::Cos (var[0])*TMath::Power (1 - TMath::Power (par[1]*TMath::Sin (var[0]), 2), 1./2) > TMath::Power (par[1]*TMath::Cos (var[0]), 2) + (1 - TMath::Power (par[1]*TMath::Sin (var[0]),2))) eneKin_a =  K*TMath::Power ((par[1]*TMath::Cos (var[0]) - TMath::Power (1 - TMath::Power (par[1]*TMath::Sin (var[0]), 2),1./2)), 2);
  else eneKin_a = K*TMath::Power ((par[1]*TMath::Cos (var[0]) + TMath::Power (1 - TMath::Power (par[1]*TMath::Sin (var[0]), 2), 1./2)), 2);

  Int_t reg_a;
  // absorption in aluminum protection foil
  for (Int_t i = 0; i < nGran; i++) // nGran steps in aluminum
    {
      reg_a = 0;
      for (Int_t j = 0; j < 9; j++) // find energy loss regime
  	{
  	  if (eneKin_a >= ene_4alu_a[j])
	    {
	      reg_a = j;
	      break;
	    }
  	}
      if (reg_a < 8) eneKin_a -=  1./nGran*(par[2]/TMath::Abs (TMath::Cos (var[0])))*parA_4alu_a[reg_a]*TMath::Power (eneKin_a, parB_4alu_a[reg_a]); // power law
      else if (reg_a == 8) eneKin_a -= 1./nGran*(par[2]/TMath::Abs (TMath::Cos (var[0])))*(parA_4alu_a[reg_a]*TMath::Power (eneKin_a, 2) + parB_4alu_a[reg_a]*eneKin_a + parC_4alu_a); // 2nd order polynomial law
    }

  // absorption in silicon dead layer DSSSD
  for (Int_t i = 0; i < nGran; i++) // nGran steps in aluminum
    {
      reg_a = 0;
      for (Int_t j = 0; j < 8; j++) // find energy loss regime
  	{
  	  if (eneKin_a >= ene_4si_a[j])
	    {
	      reg_a = j;
	      break;
	    }
  	}
      if (reg_a < 8) eneKin_a -=  1./nGran*(par[3]/TMath::Abs (TMath::Cos (var[0])))*parA_4si_a[reg_a]*TMath::Power (eneKin_a, parB_4si_a[reg_a]); // power law
    }
 
  return eneKin_a;
}


Double_t ene_proton (Double_t *var, Double_t *par){
  
  // par[0] = eneCms
  // par[1] = gam3
  // par[2] = alu thickness
  // par[3] = si thickness
  // var[0] = theta lab
  
  Double_t K = A1*A3_p*par[0]/(A2*(A1 + A2)*TMath::Power (par[1], 2));
  Double_t eneKin_p = 0.;
  if (2*par[1]*TMath::Cos (var[0])*TMath::Power (1 - TMath::Power (par[1]*TMath::Sin (var[0]), 2), 1./2) > TMath::Power (par[1]*TMath::Cos (var[0]), 2) + (1 - TMath::Power (par[1]*TMath::Sin (var[0]),2))) eneKin_p =  K*TMath::Power ((par[1]*TMath::Cos (var[0]) - TMath::Power (1 - TMath::Power (par[1]*TMath::Sin (var[0]), 2),1./2)), 2);
  else eneKin_p = K*TMath::Power ((par[1]*TMath::Cos (var[0]) + TMath::Power (1 - TMath::Power (par[1]*TMath::Sin (var[0]), 2), 1./2)), 2);

  Int_t reg_p;
  // absorption in aluminum protection foil
  for (Int_t i = 0; i < nGran; i++) // nGran steps in aluminum
    {
      reg_p = 0;
      for (Int_t j = 0; j < 6; j++) // find energy loss regime
  	{
  	  if (eneKin_p >= ene_4alu_p[j])
	    {
	      reg_p = j;
	      break;
	    }
  	}
      if (reg_p < 6) eneKin_p -=  1./nGran*(par[2]/TMath::Abs (TMath::Cos (var[0])))*parA_4alu_p[reg_p]*TMath::Power (eneKin_p, parB_4alu_p[reg_p]); // power law
    }

  // absorption in silicon dead layer DSSSD
  for (Int_t i = 0; i < nGran; i++) // nGran steps in aluminum
    {
      reg_p = 0;
      for (Int_t j = 0; j < 7; j++) // find energy loss regime
  	{
  	  if (eneKin_p >= ene_4si_p[j])
	    {
	      reg_p = j;
	      break;
	    }
  	}
      if (reg_p < 7) eneKin_p -=  1./nGran*(par[3]/TMath::Abs (TMath::Cos (var[0])))*parA_4si_p[reg_p]*TMath::Power (eneKin_p, parB_4si_p[reg_p]); // power law
    }
  
  return eneKin_p;
}
  
  
  
Int_t setBinWidth (){

  // generate histograms with variable bin width of  detectors
  // see downwritings 'DSSSD geometry'

  Int_t count_s3 = 0; 
  for (Int_t i = 0; i < n_ch_s3; i++)
    {
      mu_s3[i] = r_i_s3 + ((r_a_s3 - r_i_s3) - ((n_ch_s3 - 1)*d_s3))/2. + i*d_s3; // cm; mean radial position of the strip
      mu_s3f[i] = TMath::ATan (mu_s3[i]/z_s3f); // rad;
      mu_low_s3f[i]  = TMath::ATan ((mu_s3[i] - (d_s3 - s_s3)/2.)/z_s3f); // lower radial limit
      mu_high_s3f[i] = TMath::ATan ((mu_s3[i] + (d_s3 - s_s3)/2.)/z_s3f); // upper
      
      bin_s3f[count_s3] = mu_low_s3f[i]; // sort the limits for the binning of the histogram
      count_s3++;
      bin_s3f[count_s3] = mu_high_s3f[i];
      count_s3++;
    }

  // limiters of the strips/bins
  // for (Int_t i = 0; i < 2*n_ch_s3; i++) std::cout << " i (s3f) = " << i << ", limit = " << bin_s3f[i] << " rad" << std::endl;

  return 0;
}


Int_t getKinem (){

  // get kinematics ala bass1980.pdf
  const Double_t v_inf = TMath::Sqrt (2*eneBeam/(A1*931.5));
  const Double_t mu = (A1*A2/(A1 + A2))*m;
  const Double_t eneCms = TMath::Power (v_inf, 2)*mu/2;
  Double_t qVal_a[8], gam3_a[8], qVal_p[11], gam3_p[11];
  for (Int_t i = 0; i < 11; i++)
    {
      // protons
      qVal_p[i] = m*(A1 + A2 - A3_p - A4_p) - eneEx_p[i];
      gam3_p[i] = TMath::Sqrt (A1*A3_p/(A2*A4_p)*(eneCms/(eneCms + qVal_p[i])));

      f_eneVsPol_p[i] = new TF1 (Form ("f_eneVsPol_p[%i]", i), ene_proton, thetaMin_s3f/180*TMath::Pi (), thetaMin_s3f/180*TMath::Pi (), 4);
      f_eneVsPol_p[i]->SetParameters (eneCms, gam3_p[i], alu_fwd, si);
      
      if (i >= 8) continue;

      // alphas
      qVal_a[i] = m*(A1 + A2 - A3_a - A4_a) - eneEx_a[i];
      gam3_a[i] = TMath::Sqrt (A1*A3_a/(A2*A4_a)*(eneCms/(eneCms + qVal_a[i])));

      f_eneVsPol_a[i] = new TF1 (Form ("f_eneVsPol_a[%i]", i), ene_alpha, thetaMin_s3f/180*TMath::Pi (), thetaMin_s3f/180*TMath::Pi (), 4);
      f_eneVsPol_a[i]->SetParameters (eneCms, gam3_a[i], alu_fwd, si);
    }

  // rotate by 90 deg
  Double_t pol_s3f[nPoints], ene_a[8][nPoints], ene_p[11][nPoints];
  for (Int_t i = 0; i < 11; i++)
    {
      // protons
      for (Int_t j = 0; j < nPoints; j++)
  	{
	  if (i == 0) pol_s3f[j] = thetaMin_s3f*TMath::Pi ()/180 + (thetaMax_s3f - thetaMin_s3f)*TMath::Pi ()/(180.*nPoints)*j; // define angles only once
  	  ene_p[i][j] = f_eneVsPol_p[i]->Eval (pol_s3f[j]);
	}
      g_polVsEne_p[i] = new TGraph (nPoints, ene_p[i], pol_s3f);
      g_polVsEne_p[i]->SetLineColor (2);
      g_polVsEne_p[i]->SetLineWidth (3);
      
      if (i >= 8) continue;

      // alphas
      for (Int_t j = 0; j < nPoints; j++) ene_a[i][j] = f_eneVsPol_a[i]->Eval (pol_s3f[j]);
      g_polVsEne_a[i] = new TGraph (nPoints, ene_a[i], pol_s3f);
      g_polVsEne_a[i]->SetLineColor (1);
      g_polVsEne_a[i]->SetLineWidth (3);
    }
   
  return 0;
}


Int_t getRes (){
  // extract an effective energy resolution for the fit region/sigma
  // separately for alphas and protons
  // I did some fits by hand in different spectra to cover the energy region

  const Int_t nPoi_a = 4;
  const Int_t nPoi_p = 6;
  const Double_t ene_a[nPoi_a] = {5.7, 7.8, 11.1, 13.0}; // MeV
  const Double_t ene_p[nPoi_p] = {5.8, 6.4, 6.5, 6.7, 8.6, 9.1};
  const Double_t sig_a[nPoi_a] = {162., 138., 114., 121.}; // sigma in QDC chn.
  const Double_t sig_p[nPoi_p] = {82., 127., 56., 51., 82., 76.};
  const Double_t qdc_a[nPoi_a] = {6005., 8304., 11890., 13961.}; // QDC chn.
  const Double_t qdc_p[nPoi_p] = {6094., 6637., 6871., 7089., 9013., 9484.};
  
  const Double_t dEne_a[nPoi_a] = {0.02, 0.02, 0.02, 0.02}; // MeV
  const Double_t dEne_p[nPoi_p] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  const Double_t dSig_a[nPoi_a] = {16., 17., 7., 2.}; // sigma uncertainty in QDC chn.
  const Double_t dSig_p[nPoi_p] = {4., 5., 6., 4., 2., 3.};
  const Double_t dQDC_a[nPoi_a] = {5., 6., 9., 3.};  // uncertainty QDC chn.
  const Double_t dQDC_p[nPoi_p] = {3., 3., 6., 3., 2., 3.}; 
  Double_t res_a[nPoi_a], res_p[nPoi_p], dRes_a[nPoi_a], dRes_p[nPoi_p];

  for (Int_t i = 0; i < 6; i++)
    {
      res_p[i] = sig_p[i]/qdc_p[i];
      dRes_p[i] = TMath::Sqrt (TMath::Power (dSig_p[i]/sig_p[i], 2) + TMath::Power (dQDC_p[i]/qdc_p[i], 2)); // relative error
      dRes_p[i] *= res_p[i]; // absolute
      
      if (i >= 4) continue; // alphas
      res_a[i] = sig_a[i]/qdc_a[i];
      dRes_a[i] = TMath::Sqrt (TMath::Power (dSig_a[i]/sig_a[i], 2) + TMath::Power (dQDC_a[i]/qdc_a[i], 2)); // relative error
      dRes_a[i] *= res_a[i]; // absolute
    }
  gr_res_a = new TGraphErrors (nPoi_a, ene_a, res_a, dEne_a, dRes_a);
  gr_res_a->SetTitle ("Alphas: relative energy resolution");
  gr_res_p = new TGraphErrors (nPoi_p, ene_p, res_p, dEne_p, dRes_p);
  gr_res_p->SetTitle ("Protons: relative energy resolution");

  // fit function
  f_res_a = new TF1 ("f_res_a", my_resS3, 2., 14., 2);
  f_res_a->SetParLimits (0, 0., 1.);
  f_res_p = new TF1 ("f_res_p", my_resS3, 2., 14., 2);
  f_res_p->SetParLimits (0, 0., 1.);

  
  TCanvas *c_res = new TCanvas ("c_res", "c_res", 700, 500);
  c_res->Draw ();
  c_res->Divide (1, 2);
  
  c_res->cd (1)->SetTickx (); // alphas
  c_res->cd (1)->SetTicky ();
  gr_res_a->Draw ("AP");
  gr_res_a->Fit (f_res_a, "RQI");

  c_res->cd (2)->SetTickx (); // protons
  c_res->cd (2)->SetTicky ();
  gr_res_p->Draw ("AP");
  gr_res_p->Fit (f_res_p, "RQI");
  
   c_res->Print (Form ("r%03i_dsssdRes_event_p_s3f.pdf", n_run));  
  return 0;
}


Int_t fillExp (){

  // pull data
  Long64_t event_g, event_p, event_gp;
  Int_t mult_g, mult_p, mult_gp;
  std::vector <Int_t> *det_g = 0;
  std::vector <Int_t> *det_p = 0;
  std::vector <Long64_t> *time_g = 0;
  std::vector <Long64_t> *time_p = 0;
  std::vector <Double_t> *ene_g = 0;
  std::vector <Double_t> *ene_p = 0;
  TBranch *b_event_g, *b_mult_g, *b_det_g, *b_time_g, *b_ene_g, *b_event_p, *b_mult_p, *b_det_p, *b_time_p, *b_ene_p, *b_event_gp, *b_mult_gp;

  TChain *f_in = new TChain ("t_event");
  //  f_in->Add ("/home/mheine/Delete/goData.root");
  f_in->Add ("/home/gustavo/Documents/orsayData/r010_ene_p1_1.root");
  f_in->SetBranchAddress ("event_g", &event_g, &b_event_g);
  f_in->SetBranchAddress ("mult_g", &mult_g, &b_mult_g);
  f_in->SetBranchAddress ("det_g", &det_g, &b_det_g);
  f_in->SetBranchAddress ("time_g", &time_g, &b_time_g);
  f_in->SetBranchAddress ("ene_g", &ene_g, &b_ene_g);
  f_in->SetBranchAddress ("event_p", &event_p, &b_event_p);
  f_in->SetBranchAddress ("mult_p", &mult_p, &b_mult_p);
  f_in->SetBranchAddress ("det_p", &det_p, &b_det_p);
  f_in->SetBranchAddress ("time_p", &time_p, &b_time_p);
  f_in->SetBranchAddress ("ene_p", &ene_p, &b_ene_p);
  f_in->SetBranchAddress ("event_gp", &event_gp, &b_event_gp);
  f_in->SetBranchAddress ("mult_gp", &mult_gp, &b_mult_gp);

  
  setBinWidth ();  
  h_polVsEne = new TH2F ("h_polVsEne", "", ene_bin, ene_low, ene_high, 2*n_ch_s3 - 1, bin_s3f);
  h_detVsDet = new TH2F ("h_detVsDet", "", 26, -.5, 25.5, 26, -.5, 25.5);
  for (Int_t i = 0; i < 25; i++)
    {
      h_ene[i] = new TH1F (Form ("h_ene[%i]", i), Form ("det. %i", i + 1), ene_bin, ene_low, ene_high);
      h_time[i] = new TH1F (Form ("h_time[%i]", i), Form ("h_time[%i]", i), time_bin, time_low, time_high);
    }

  
  // data loop
  for (Long64_t i = 0; i < f_in->GetEntries (); i++)
    {
      f_in->GetEntry (i);
      if (i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << " " << std::setw (4) << (Long64_t) (i*100.0001)/f_in->GetEntries () << "% done" << std::endl;

      if (mult_p == 0) continue;
     
      
      for (Int_t j = 0; j < mult_p; j++)
	{
	  h_polVsEne->Fill (ene_p->at (j), mu_s3f[chn_s3f[det_p->at (j) - 1] - 1]);
	  h_ene[chn_s3f[det_p->at (j) - 1] - 1]->Fill (ene_p->at (j));
	}
      
      if (mult_p == 2)
	{
	  //	  h_detVsDet->Fill (det_p->at (0), det_p->at (1));
	  h_detVsDet->Fill (chn_s3f[det_p->at (0) - 1], chn_s3f[det_p->at (1) - 1]);
	  if (det_p->at (0) == det_p->at (1)) h_time[chn_s3f[det_p->at (0) - 1] - 1]->Fill (time_p->at (0) - time_p->at (1));
	}
    }
  
  return 0;
}


Int_t plotExp (){

  TCanvas *c_ene[6];
  TLine *li_a, *li_p;
  li_a = new TLine ();
  li_a->SetLineColor (1);
  li_a->SetLineStyle (1);
  li_a->SetLineWidth (4);
  li_p = new TLine ();
  li_p->SetLineColor (2);
  li_p->SetLineStyle (1);
  li_p->SetLineWidth (4);
  Double_t posX_a[25][8], posY_a[25][8], posX_p[25][11], posY_p[25][11], dPosX_a[25][8], dPosX_p[25][11]; // MeV; lines gating on peak regions
  for (Int_t i = 0; i < 6; i++) // chn. 25 is sync
    {
      c_ene[i] = new TCanvas (Form ("c_ene[%i]", i), Form ("c_ene[%i]", i), 700, 500);
      c_ene[i]->Divide (2, 2);

      for (Int_t j = 0; j < 4; j++)
	{
	  c_ene[i]->cd (j + 1)->SetTickx ();
	  c_ene[i]->cd (j + 1)->SetTicky ();

	  h_ene[4*i + j]->Draw ("HIST");

          for (Int_t k = 0; k < 11; k++) // mark peaks, set up fit parameters; protons
	    {
	      posX_p[4*i + j][k] = f_eneVsPol_p[k]->Eval (TMath::ATan ((mu_s3[4*i + j])/(z_s3f)));
	      posY_p[4*i + j][k] = h_ene[4*i + j]->GetBinContent (TMath::Abs ((ene_low - posX_p[4*i + j][k])*ene_bin/(ene_high - ene_low)));
	      dPosX_p[4*i + j][k] = TMath::Abs (posX_p[4*i + j][k]*f_res_p->Eval (f_eneVsPol_p[k]->Eval (TMath::ATan ((mu_s3[4*i + j])/(z_s3f)))));
	      li_p->DrawLine (posX_p[4*i + j][k] - dPosX_p[4*i + j][k], posY_p[4*i + j][k], posX_p[4*i + j][k] + dPosX_p[4*i + j][k], posY_p[4*i + j][k]);


	      if (k >= 8) continue; // alphas

	      posX_a[4*i + j][k] = f_eneVsPol_a[k]->Eval (TMath::ATan ((mu_s3[4*i + j])/(z_s3f)));
	      posY_a[4*i + j][k] = h_ene[4*i + j]->GetBinContent (TMath::Abs ((ene_low - posX_a[4*i + j][k])*ene_bin/(ene_high - ene_low)));
	      dPosX_a[4*i + j][k] = TMath::Abs (posX_a[4*i + j][k]*f_res_p->Eval (f_eneVsPol_a[k]->Eval (TMath::ATan ((mu_s3[4*i + j])/(z_s3f)))));
	      li_a->DrawLine (posX_a[4*i + j][k] - dPosX_a[4*i + j][k], posY_a[4*i + j][k], posX_a[4*i + j][k] + dPosX_a[4*i + j][k], posY_a[4*i + j][k]);
	    }
 
	  h_ene[4*i + j]->GetXaxis ()->SetRangeUser (2., 15.); // MeV
	  h_ene[4*i + j]->GetYaxis ()->SetRangeUser (0., 2.*h_ene[4*i + j]->GetBinContent (TMath::Abs ((ene_low - f_eneVsPol_a[1]->Eval (TMath::ATan ((mu_s3[4*i + j])/(z_s3f))))*ene_bin/(ene_high - ene_low))));
	}
     }
  c_ene[0]->Print (Form ("r%03i_eneDet_event_p_s3f.pdf(", n_run));  
  c_ene[1]->Print (Form ("r%03i_eneDet_event_p_s3f.pdf", n_run));  
  c_ene[2]->Print (Form ("r%03i_eneDet_event_p_s3f.pdf", n_run));  
  c_ene[3]->Print (Form ("r%03i_eneDet_event_p_s3f.pdf", n_run));  
  c_ene[4]->Print (Form ("r%03i_eneDet_event_p_s3f.pdf", n_run));  
  c_ene[5]->Print (Form ("r%03i_eneDet_event_p_s3f.pdf)", n_run));
  
  
  TCanvas *c_polVsEne = new TCanvas ("c_polVsEne", "c_polVsEne", 700, 500);
  c_polVsEne->Draw ();
  c_polVsEne->cd ();
  c_polVsEne->SetTickx ();
  c_polVsEne->SetTicky ();
  c_polVsEne->SetLogz ();
  
  h_polVsEne->Draw ("COLZ");

  for (Int_t i = 0; i < 11; i++)
    {
      g_polVsEne_p[i]->Draw ("SAME");

      if (i >= 8) continue;
      g_polVsEne_a[i]->Draw ("SAME");
    }
  c_polVsEne->Print (Form ("r%03i_polVsEne_event_p_s3f.pdf", n_run));  
  
  
  TCanvas *c_detVsDet = new TCanvas ("c_detVsDet", "c_detVsDet", 700, 500);
  c_detVsDet->Draw ();
  c_detVsDet->cd ();
  c_detVsDet->SetTickx ();
  c_detVsDet->SetTicky ();
  c_detVsDet->SetLogz ();

  h_detVsDet->GetXaxis ()->SetTitle ("det_{A}");
  h_detVsDet->GetYaxis ()->SetTitle ("det_{B}");
  h_detVsDet->Draw ("COLZ");
  c_detVsDet->Print (Form ("r%03i_detCorr_event_p_s3f.pdf", n_run));  

  
  TCanvas *c_time = new TCanvas ("c_time", "c_time", 700, 500);
  c_time->Divide (6, 4);
  for (Int_t i = 0; i < 24; i++) // chn. 25 is sync
    {
      c_time->cd (i + 1)->SetTickx ();
      c_time->cd (i + 1)->SetTicky ();

      h_time[i]->Draw ("HIST");
    }
  c_time->Print (Form ("r%03i_dT_sameDet_event_p_s3f.pdf", n_run));
  
  return 0;
}


void event_p_s3f (){

  // display calibrated data and compare to kinematics calculation

  gStyle->SetOptStat (0);

  getKinem ();
  fillExp ();

  getRes ();
  plotExp ();
}
