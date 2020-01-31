#include "/home/gustavo/Documents/orsayData/libPerso.h"


void event_gp (){


  const Int_t n_run = 10;
  gStyle->SetOptStat (1111110);

  
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
  //f_in->Add ("/home/mheine/Delete/goData.root");
  f_in->SetBranchAddress ("event_g" , &event_g , &b_event_g);
  f_in->SetBranchAddress ("mult_g"  , &mult_g  , &b_mult_g);
  f_in->SetBranchAddress ("det_g"   , &det_g   , &b_det_g);    // 1..40
  f_in->SetBranchAddress ("time_g"  , &time_g  , &b_time_g);
  f_in->SetBranchAddress ("ene_g"   , &ene_g   , &b_ene_g);
  f_in->SetBranchAddress ("event_p" , &event_p , &b_event_p);
  f_in->SetBranchAddress ("mult_p"  , &mult_p  , &b_mult_p);
  f_in->SetBranchAddress ("det_p"   , &det_p   , &b_det_p);    // 1..100
  f_in->SetBranchAddress ("time_p"  , &time_p  , &b_time_p);
  f_in->SetBranchAddress ("ene_p"   , &ene_p   , &b_ene_p);
  f_in->SetBranchAddress ("event_gp", &event_gp, &b_event_gp);
  f_in->SetBranchAddress ("mult_gp" , &mult_gp , &b_mult_gp);


  // declare objects to fill data
  const Double_t eneG_low = .05; // MeV
  const Double_t eneG_high = 3.7;
  const Int_t eneG_bin = (eneG_high - eneG_low)*200;
  const Double_t eneP1_low = 0.;
  const Double_t eneP1_high = 14.;
  const Int_t eneP1_bin = (eneP1_high - eneP1_low)*50;
  const Double_t eneP4_low = 0.;
  const Double_t eneP4_high = 6.;
  const Int_t eneP4_bin = (eneP4_high - eneP4_low)*50;
  const Double_t diff_low = -450.; // ns
  const Double_t diff_high = 450.;
  const Int_t diff_bin = (diff_high - diff_low)/4;

  TH1D *h_multG = new TH1D ("h_multG", "Gamma Multiplicities", 46, -0.5, 45.5); // deal with huge numbers
  TH1D *h_multP = new TH1D ("h_multP", "Particle Multiplicities", 26, -0.5, 25.5);
  TH1F *h_eneG[40], *h_eneG_gp[40], *h_eneP[100], *h_eneP_gp[100]; // energies
  TH1F *h_diffG[5][5], *h_diffP[4][4], *h_diffGP[5][4]; // timing
  TH2F *h_diffGP_eneG[5][4], *h_diffGP_eneP[5][4];
  
  for (Int_t i = 0; i < 100; i++) // energies
    {
      if (i < 25) // 0..24: S3F
	{
	  h_eneP[i] = new TH1F (Form ("h_eneP[%i]", i), Form ("input %i", i + 1), eneP1_bin, eneP1_low, eneP1_high);
	  h_eneP_gp[i] = new TH1F (Form ("h_eneP_gp[%i]", i), Form ("input %i", i + 1), eneP1_bin, eneP1_low, eneP1_high);
	  h_eneP_gp[i]->SetLineColor (2);
	}
      else if (i >=75) // 75..99: S3B
	{
	  h_eneP[i] = new TH1F (Form ("h_eneP[%i]", i), Form ("input %i", i + 1), eneP4_bin, eneP4_low, eneP4_high);
	  h_eneP_gp[i] = new TH1F (Form ("h_eneP_gp[%i]", i), Form ("input %i", i + 1), eneP4_bin, eneP4_low, eneP4_high);
	  h_eneP_gp[i]->SetLineColor (2);
	}

      if (i < 40)
	{
	  h_eneG[i] = new TH1F (Form ("h_eneG[%i]", i), Form ("input %i", i + 1), eneG_bin, eneG_low, eneG_high);
	  h_eneG_gp[i] = new TH1F (Form ("h_eneG_gp[%i]", i), Form ("input %i", i + 1), eneG_bin, eneG_low, eneG_high);
	  h_eneG_gp[i]->SetLineColor (2);
	}
    }
  for (Int_t i = 0; i < 5; i++) // timing
    {
      for (Int_t j = 0; j < 5; j++)
	{
	  h_diffG[i][j] = new TH1F (Form ("h_diffG[%i][%i]", i, j), Form ("mod %i vs. mod %i", i + 1, j + 1), diff_bin, diff_low, diff_high);
	  if (j < 4)
	    {
	      if (j == 0) // S3F
		{
		  h_diffGP[i][j] = new TH1F (Form ("h_diffGP[%i][%i]", i, j), Form ("mod %i vs. card %i", i + 1, j + 124), diff_bin, diff_low, diff_high);
		  h_diffGP_eneG[i][j] = new TH2F (Form ("h_diffGP_eneG[%i][%i]", i, j), Form ("mod %i vs. card %i", i + 1, j + 124), diff_bin, diff_low, diff_high, eneG_bin, eneG_low, eneG_high);
		  h_diffGP_eneP[i][j] = new TH2F (Form ("h_diffGP_eneP[%i][%i]", i, j), Form ("mod %i vs. card %i", i + 1, j + 124), diff_bin, diff_low, diff_high, eneP1_bin, eneP1_low, eneP1_high);
		  
		  if (i == 0 || i == 3) h_diffP[i][j] = new TH1F (Form ("h_diffP[%i][%i]", i, j), Form ("card %i vs. card %i", i + 124, j + 124), diff_bin, diff_low, diff_high);
		}
	      else if (j == 3) // S3B
		{
		  h_diffGP[i][j] = new TH1F (Form ("h_diffGP[%i][%i]", i, j), Form ("mod %i vs. card %i", i + 1, j + 124), diff_bin, diff_low, diff_high);
		  h_diffGP_eneG[i][j] = new TH2F (Form ("h_diffGP_eneG[%i][%i]", i, j), Form ("mod %i vs. card %i", i + 1, j + 124), diff_bin, diff_low, diff_high, eneG_bin, eneG_low, eneG_high);
		  h_diffGP_eneP[i][j] = new TH2F (Form ("h_diffGP_eneP[%i][%i]", i, j), Form ("mod %i vs. card %i", i + 1, j + 124), diff_bin, diff_low, diff_high, eneP4_bin, eneP4_low, eneP4_high);
		  
		  if (i == 0 || i == 3) h_diffP[i][j] = new TH1F (Form ("h_diffP[%i][%i]", i, j), Form ("card %i vs. card %i", i + 124, j + 124), diff_bin, diff_low, diff_high);
		}
	    }
	}
    }


  // data loop
  Int_t modA, modB, detP1_ID, detP4_ID;
  Double_t sumEneG, sumEneP1, sumEneP4;
  for (Long64_t i = 0; i < f_in->GetEntries (); i++)
    {
      f_in->GetEntry (i);
      if (i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << std::setw (4) << (Long64_t) (i*100.001)/f_in->GetEntries () << "% done" << std::endl;
      sumEneG = 0.;
      sumEneP1 = 0.;
      sumEneP4 = 0.;
      detP1_ID = 0;
      detP4_ID = 0;
      

      if (mult_g > 0) // gammas
      	{
      	  h_multG->Fill (mult_g); // gamma multiplicity
      	  if (mult_g == 2) // time offsets between two modules
      	    {
      	      modA = TMath::Max ((det_g->at (0) - 1)/8, (det_g->at (1) - 1)/8);
      	      modB = TMath::Min ((det_g->at (0) - 1)/8, (det_g->at (1) - 1)/8);
     	      h_diffG[modA][modB]->Fill (time_g->at (1) - time_g->at (0));
      	    }

      	  for (Int_t j = 0; j < mult_g; j++) // energies
      	    {
  	      sumEneG += ene_g->at (j);
  	      h_eneG[det_g->at (j) - 1]->Fill (ene_g->at (j));
      	      if (mult_gp > 0) h_eneG_gp[det_g->at (j) - 1]->Fill (ene_g->at (j));
      	    }
  	  h_eneG[0]->Fill (sumEneG); // sum energy
	  if (mult_gp > 0) h_eneG_gp[0]->Fill (sumEneG);
      	}

      
      if (mult_p > 0) // particles
      	{
      	  h_multP->Fill (mult_p); // particle multiplicity
      	  if (mult_p == 2) // time offsets between cards
      	    {
      	      modA = TMath::Max ((det_p->at (0) - 1)/25, (det_p->at (1) - 1)/25);
      	      modB = TMath::Min ((det_p->at (0) - 1)/25, (det_p->at (1) - 1)/25);

      	      h_diffP[modA][modB]->Fill (time_p->at (1) - time_p->at (0));
      	    }

      	  for (Int_t j = 0; j < mult_p; j++) // energies
      	    {
      	      h_eneP[det_p->at (j) - 1]->Fill (ene_p->at (j));
      	      if (det_p->at (j) <= 25) // S3F
		{
		  sumEneP1 += ene_p->at (j);
		  detP1_ID = j;
		}
      	      if (det_p->at (j) > 75) // S3B
		{
		  sumEneP4 += ene_p->at (j);
		  detP4_ID = j;
		}

      	      if (mult_gp > 0) h_eneP_gp[det_p->at (j) - 1]->Fill (ene_p->at (j));
      	    }

      	  if (sumEneP1 > 1.e-3) // sum energy: seed to shift with the kinematics
      	    {
      	      h_eneP[24]->Fill (sumEneP1);
      	      if (mult_gp > 0) h_eneP_gp[24]->Fill (sumEneP1);
      	    }
      	  if (sumEneP4 > 1.e-3)
      	    {
      	      h_eneP[99]->Fill (sumEneP4);
      	      if (mult_gp > 0) h_eneP_gp[99]->Fill (sumEneP4);
      	    }
      	}

      
      if (mult_gp > 0) // coincidence timing
      	{
  	  h_diffGP[(Int_t) (det_g->at (0) - 1)/8][(Int_t) (det_p->at (0) - 1)/25]->Fill (time_g->at (0) - time_p->at (0));
  	  h_diffGP_eneG[(Int_t) (det_g->at (0) - 1)/8][(Int_t) (det_p->at (0) - 1)/25]->Fill (time_g->at (0) - time_p->at (0), sumEneG);
	  if (sumEneP1 > 1.e-3) h_diffGP_eneP[(Int_t) (det_g->at (0) - 1)/8][(Int_t) (det_p->at (detP1_ID) - 1)/25]->Fill (time_g->at (0) - time_p->at (0), sumEneP1);
	  if (sumEneP4 > 1.e-3) h_diffGP_eneP[(Int_t) (det_g->at (0) - 1)/8][(Int_t) (det_p->at (detP4_ID) - 1)/25]->Fill (time_g->at (0) - time_p->at (0), sumEneP4);
      	}
    }


  // plot data
  TCanvas *c_mult = new TCanvas ("c_mult", "Detector Multiplicities", 700, 500); // multiplicity pattern
  c_mult->Draw ();
  c_mult->Divide (1, 2);
  
  c_mult->cd (1)->SetTickx (); // gammas
  c_mult->cd (1)->SetTicky ();
  c_mult->cd (1)->SetLogy ();
  h_multG->Draw ();
  c_mult->cd (2)->SetTickx (); // particles
  c_mult->cd (2)->SetTicky ();
  c_mult->cd (2)->SetLogy ();
  h_multP->Draw ();
  c_mult->Print (Form ("r%03i_event_gp_mult.pdf", n_run));

  
  TCanvas *c_eneG = new TCanvas ("c_eneG", "Gamma Energies", 700, 500); // gamma energies
  c_eneG->Draw ();
  c_eneG->Divide (8, 5);
  for (Int_t i = 0; i < 40; i++)
    {
      c_eneG->cd (i + 1)->SetTickx ();
      c_eneG->cd (i + 1)->SetTicky ();
      h_eneG[i]->Draw ();
      h_eneG_gp[i]->Draw ("SAME");
    }
  c_eneG->Print (Form ("r%03i_event_gp_eneG.pdf", n_run));

  
  TCanvas *c_eneP[4]; // particle energies
  for (Int_t i = 0; i < 4; i++)
    {
      if (i == 0 || i == 3)
  	{
  	  c_eneP[i]= new TCanvas (Form ("c_eneP[%i]", i), Form ("Particle Energies Card %i", i + 124), 700, 500);
  	  c_eneP[i]->Draw ();
  	  c_eneP[i]->Divide (5, 5);
  	  for (Int_t j = 0; j < 25; j++)
  	    {
  	      c_eneP[i]->cd (j + 1)->SetTickx ();
  	      c_eneP[i]->cd (j + 1)->SetTicky ();

  	      h_eneP[i*25 + j]->Draw ();
  	      h_eneP_gp[i*25 + j]->Draw ("SAME");
  	    }
  	}
    }
  c_eneP[0]->Print (Form ("r%03i_event_gp_eneP.pdf(", n_run));
  c_eneP[3]->Print (Form ("r%03i_event_gp_eneP.pdf)", n_run));

  
  TCanvas *c_diffG = new TCanvas ("c_diffG", "Gamma Module Time Difference", 700, 500); // gamma module time difference
  c_diffG->Draw ();
  c_diffG->Divide (5, 5);
  for (Int_t i = 0; i < 5; i++)
    {
      for (Int_t j = 0; j < 5; j++)
	{
	  c_diffG->cd (5*i + j + 1)->SetTickx ();
	  c_diffG->cd (5*i + j + 1)->SetTicky ();
	  c_diffG->cd (5*i + j + 1)->SetLogy ();

	  h_diffG[i][j]->Draw ();
	}
    }
  c_diffG->Print (Form ("r%03i_event_gp_diffG.pdf", n_run));

  
  TCanvas *c_diffP = new TCanvas ("c_diffP", "Particle Card Time Difference", 700, 500); // particle cards time difference
  c_diffP->Draw ();
  c_diffP->Divide (4, 4);
  for (Int_t i = 0; i < 4; i++)
    {
      for (Int_t j = 0; j < 4; j++)
	{
	  if (i == 1 || i == 2) continue;
	  if (j == 1 || j == 2) continue;
	  
	  c_diffP->cd (4*i + j + 1)->SetTickx ();
	  c_diffP->cd (4*i + j + 1)->SetTicky ();
	  c_diffP->cd (4*i + j + 1)->SetLogy ();

	  h_diffP[i][j]->Draw ();
	}
    }
  c_diffP->Print (Form ("r%03i_event_gp_diffG.pdf", n_run));

  
  TCanvas *c_diffGP = new TCanvas ("c_diffGP", "Gamma-Module Particle-Card Time Difference", 700, 500); // gamma particle time difference
  c_diffGP->Draw ();
  c_diffGP->Divide (5, 4);
  for (Int_t i = 0; i < 5; i++)
    {
      for (Int_t j = 0; j < 4; j++)
  	{
  	  if (j == 1 || j == 2) continue;
	  
  	  c_diffGP->cd (5*j + i + 1)->SetTickx ();
  	  c_diffGP->cd (5*j + i + 1)->SetTicky ();
  	  c_diffGP->cd (5*j + i + 1)->SetLogy ();

  	  h_diffGP[i][j]->Draw ();
  	}
    }
  c_diffGP->Print (Form ("r%03i_event_gp_diffGP.pdf", n_run));

  
  TCanvas *c_diffGP_eneP = new TCanvas ("c_diffGP_eneP", "Particle Energy versus Gamma-Module Particle-Card Time Difference", 700, 500); // particle energy versus gamma particle time difference
  c_diffGP_eneP->Draw ();
  c_diffGP_eneP->Divide (5, 4);
  for (Int_t i = 0; i < 5; i++)
    {
      for (Int_t j = 0; j < 4; j++)
	{
	  if (j == 1 || j == 2) continue;
	  
	  c_diffGP_eneP->cd (5*j + i + 1)->SetTickx ();
	  c_diffGP_eneP->cd (5*j + i + 1)->SetTicky ();

	  h_diffGP_eneP[i][j]->Draw ();
	}
    }
  c_diffGP_eneP->Print (Form ("r%03i_event_gp_diffGP_eneP.pdf", n_run));

  
  TCanvas *c_diffGP_eneG = new TCanvas ("c_diffGP_eneG", "Gamma Energy versus Gamma-Module Particle-Card Time Difference", 700, 500); // gamma energy versus gamma particle time difference
  c_diffGP_eneG->Draw ();
  c_diffGP_eneG->Divide (5, 4);
  for (Int_t i = 0; i < 5; i++)
    {
      for (Int_t j = 0; j < 4; j++)
	{
	  if (j == 1 || j == 2) continue;
	  
	  c_diffGP_eneG->cd (5*j + i + 1)->SetTickx ();
	  c_diffGP_eneG->cd (5*j + i + 1)->SetTicky ();

	  h_diffGP_eneG[i][j]->Draw ();
	}
    }
  c_diffGP_eneG->Print (Form ("r%03i_event_gp_diffGP_eneG.pdf", n_run));
  
}
