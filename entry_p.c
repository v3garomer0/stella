

#include "../Libs/libPerso.h"


void entry_p (){

  gStyle->SetOptStat (1111110);

  // pull data
  Long64_t entry_p, stamp_p, trig_p;
  Int_t plug_p, chan_p;
  TBranch *b_entry_p, *b_plug_p, *b_stamp_p, *b_trig_p, *b_chan_p;


  TChain *f_in = new TChain ("t_entry");
  //  f_in->Add ("../goData.root");                                                   // !!!!!!!!!!!
  //  f_in->Add ("../../goData_particle.root");
  f_in->Add ("/home/mheine/Delete/r010_raw_p4_1.root");
  f_in->SetBranchAddress ("entry_p", &entry_p, &b_entry_p);
  f_in->SetBranchAddress ("plug_p"  , &plug_p  , &b_plug_p);
  f_in->SetBranchAddress ("stamp_p", &stamp_p, &b_stamp_p);
  f_in->SetBranchAddress ("trig_p" , &trig_p , &b_trig_p);
  f_in->SetBranchAddress ("chan_p" , &chan_p , &b_chan_p);

  // fill data
  const Double_t chan_low = -60.e3;
  const Double_t chan_high = 0.;
  const Int_t chan_bin = (chan_high - chan_low)/20.;
  TH1F *h_input = new TH1F ("h_input", "", 129, -0.5, 128.5);
  TH2F *h_sta_ent = new TH2F ("h_sta_ent", "", 700, 0., 7000., 1e6, 0., 70.e9);
  TH1F *h_chan_p[4];
  for (Int_t i = 0; i < 4; i++) h_chan_p[i] = new TH1F (Form ("h_chan_p[%i]", i), Form ("h_chan_p[%i]", i), chan_bin, chan_low, chan_high);

  for (Long64_t i = 0; i < f_in->GetEntries (); i++)
    {
      f_in->GetEntry (i);

      //      std::cout << " " << i << " " << stamp_p << std::endl;
      h_input->Fill (plug_p);
      h_sta_ent->Fill (i, stamp_p);

      if (plug_p == 76) h_chan_p[0]->Fill (chan_p);
      else if (plug_p == 77) h_chan_p[1]->Fill (chan_p);
      else if (plug_p == 98) h_chan_p[2]->Fill (chan_p); // last to are swapped at beginning
      else h_chan_p[3]->Fill (chan_p);
    }


  // plot data
  TCanvas *c_ene = new TCanvas ("c_ene", "", 700, 500);
  c_ene->Draw ();
  c_ene->Divide (2, 2);
  for (Int_t pad = 0; pad < 4; pad++)
    {
      c_ene->cd (pad + 1)->SetTickx ();
      c_ene->cd (pad + 1)->SetTicky ();
      //      c_ene->cd (pad + 1)->SetLogy ();
      h_chan_p[pad]->Draw ();
    }

  // input channel distibution
  /* c_ene->SetTickx (); */
  /* c_ene->SetTicky (); */
  /* h_input->Draw (); */

  // stamp event correlation
  //  h_sta_ent->Draw ();
}
