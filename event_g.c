

#include "/home/gustavo/Documents/orsayData/libPerso.h"


void event_g (){
  
  const Int_t n_run = 10;
  gStyle->SetOptStat (1111110);

  
  // pull data
  Long64_t event_g;
  Int_t mult_g;
  std::vector <Int_t> *det_g = 0;
  std::vector <Long64_t> *time_g = 0;
  std::vector <Double_t> *ene_g = 0;
  TBranch *b_event_g, *b_mult_g, *b_det_g, *b_time_g, *b_ene_g;

  TChain *f_in = new TChain ("t_event");
  //  f_in->Add ("/home/mheine/Delete/goData.root");
  f_in->Add ("/home/gustavo/Documents/orsayData/r010_ene_g_1.root");
  f_in->SetBranchAddress ("event_g" , &event_g , &b_event_g);
  f_in->SetBranchAddress ("mult_g"  , &mult_g  , &b_mult_g);
  f_in->SetBranchAddress ("det_g"   , &det_g   , &b_det_g);    // 1..40
  f_in->SetBranchAddress ("time_g"  , &time_g  , &b_time_g);
  f_in->SetBranchAddress ("ene_g"   , &ene_g   , &b_ene_g);


  // declare objects to fill data
  const Double_t eneG_low = .05; // MeV
  const Double_t eneG_high = 3.7;
  const Int_t eneG_bin = (eneG_high - eneG_low)*200;
  TH1F *h_eneG[40];

  TH1D *h_multG = new TH1D ("h_multG", "Gamma Multiplicities", 46, -0.5, 45.5); // deal with huge numbers
  for (Int_t i = 0; i < 40; i++) h_eneG[i] = new TH1F (Form ("h_eneG[%i]", i), Form ("input %i", i + 1), eneG_bin, eneG_low, eneG_high);


  // data loop
  Double_t sumEneG;
  for (Long64_t i = 0; i < f_in->GetEntries (); i++)
    {
      f_in->GetEntry (i);
      if (i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << std::setw (4) << (Long64_t) (i*100.001)/f_in->GetEntries () << "% done" << std::endl;

      if (mult_g < 1) continue;
      sumEneG = 0.;
      
      h_multG->Fill (mult_g);
      for (Int_t j = 0; j < mult_g; j++) // energies
	{
	  sumEneG += ene_g->at (j);
	  h_eneG[det_g->at (j) - 1]->Fill (ene_g->at (j));
	}
      h_eneG[0]->Fill (sumEneG);
    }


  // plot data
  TCanvas *c_mult = new TCanvas ("c_mult", "Detector Multiplicities", 700, 500); // multiplicity pattern
  c_mult->Draw ();
  c_mult->cd ()->SetTickx ();
  c_mult->cd ()->SetTicky ();
  c_mult->cd ()->SetLogy ();
  
  h_multG->Draw ();
  c_mult->Print (Form ("r%03i_event_g_mult.pdf", n_run));

  
  TCanvas *c_eneG = new TCanvas ("c_eneG", "Gamma Energies", 700, 500); // energies
  c_eneG->Draw ();
  c_eneG->Divide (8, 5);
  for (Int_t i = 0; i < 40; i++)
    {
      c_eneG->cd (i + 1)->SetTickx ();
      c_eneG->cd (i + 1)->SetTicky ();

      h_eneG[i]->GetXaxis ()->SetRangeUser (0.35, 2.0);
      h_eneG[i]->Draw ();
    }
  c_eneG->Print (Form ("r%03i_event_g_eneG.pdf", n_run));

}
