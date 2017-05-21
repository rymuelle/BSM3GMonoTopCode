#define BSM3G_MonoTop_Analyzer_cxx
#include "BSM3G_MonoTop_Analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TMath.h"

void BSM3G_MonoTop_Analyzer::Loop(std::string outFileName)
{
//   In a ROOT session, you can do:
//      root> .L BSM3G_MonoTop_Analyzer.C
//      root> BSM3G_MonoTop_Analyzer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
	fChain->SetBranchStatus("*",0); 

	fChain->SetBranchStatus("Met_puppi_phi"             , 1 );  
   	fChain->SetBranchStatus("Met_puppi_pt"             , 1 );  
   	fChain->SetBranchStatus("Gen_Met"             , 1 );  
   	

	fChain->SetBranchStatus("Jet_AK8_puppi_pt"             , 1 );  
	fChain->SetBranchStatus("Jet_AK8_puppi_eta"            , 1 ); 
	fChain->SetBranchStatus("Jet_AK8_puppi_phi"            , 1 ); 
	fChain->SetBranchStatus("Jet_AK8_puppi_mass"            , 1 );  
	fChain->SetBranchStatus("Jet_AK8_bDiscriminator_pfCISVV2"            , 1 );  
	fChain->SetBranchStatus("Jet_AK8_puppi_tau1"            , 1 );  
	fChain->SetBranchStatus("Jet_AK8_puppi_tau2"            , 1 );  
	fChain->SetBranchStatus("Jet_AK8_puppi_tau3"            , 1 ); 
	fChain->SetBranchStatus("Jet_AK8_GEN_energy"            , 1 );


	fChain->SetBranchStatus("Jet_puppi_pt"             , 1 ); 
   	fChain->SetBranchStatus("Jet_puppi_eta"             , 1 );  
   	fChain->SetBranchStatus("Jet_puppi_phi"            , 1 ); 
   	fChain->SetBranchStatus("Jet_puppi_energy"            , 1 ); 
   	fChain->SetBranchStatus("Jet_puppi_bDiscriminator_pfCISVV2"            , 1 );  


	fChain->SetBranchStatus("Jet_AK8_subjet0_pt"            , 1 );
	fChain->SetBranchStatus("Jet_AK8_subjet0_eta"            , 1 );
	fChain->SetBranchStatus("Jet_AK8_subjet0_phi"            , 1 );
	fChain->SetBranchStatus("Jet_AK8_subjet0_CSVv2"            , 1 );
	fChain->SetBranchStatus("Jet_AK8_subjet0_mass"            , 1 );
 	fChain->SetBranchStatus("Jet_AK8_subjet0_energy"            , 1 ); 
 	fChain->SetBranchStatus("Jet_AK8_subjet1_pt"            , 1 );
	fChain->SetBranchStatus("Jet_AK8_subjet1_eta"            , 1 );
	fChain->SetBranchStatus("Jet_AK8_subjet1_phi"            , 1 );
	fChain->SetBranchStatus("Jet_AK8_subjet1_CSVv2"            , 1 );
	fChain->SetBranchStatus("Jet_AK8_subjet1_mass"            , 1 );
 	fChain->SetBranchStatus("Jet_AK8_subjet1_energy"            , 1 ); 



	fChain->SetBranchStatus("Muon_pt"            , 1 );
	fChain->SetBranchStatus("Muon_eta"            , 1 );
	fChain->SetBranchStatus("Muon_phi"            , 1 );
   fChain->SetBranchStatus("Muon_energy"            , 1 );
	fChain->SetBranchStatus("Muon_isoPU"            , 1 );

	fChain->SetBranchStatus("patElectron_pt"            , 1 );
	fChain->SetBranchStatus("patElectron_eta"            , 1 );
	fChain->SetBranchStatus("patElectron_phi"            , 1 );
   fChain->SetBranchStatus("patElectron_energy"            , 1 );
	fChain->SetBranchStatus("patElectron_isoPU"            , 1 );


   fChain->SetBranchStatus("Photon_pt"            , 1 );
   fChain->SetBranchStatus("Photon_eta"            , 1 );
   fChain->SetBranchStatus("Photon_phi"            , 1 );
   fChain->SetBranchStatus("Photon_energy"            , 1 );
   fChain->SetBranchStatus("Photon_PFChIso"            , 1 );
   fChain->SetBranchStatus("Photon_PFPhoIso"            , 1 );   
   fChain->SetBranchStatus("Photon_PFNeuIso"            , 1 );

	fChain->SetBranchStatus("Tau_pt"            , 1 );
	fChain->SetBranchStatus("Tau_phi"            , 1 );
	fChain->SetBranchStatus("Tau_eta"            , 1 );
   fChain->SetBranchStatus("Tau_energy"            , 1 );
	fChain->SetBranchStatus("Tau_chargedIsoPtSum"            , 1 );
	fChain->SetBranchStatus("Tau_neutralIsoPtSum"            , 1 );
	fChain->SetBranchStatus("Tau_puCorrPtSum"            , 1 );
	fChain->SetBranchStatus("Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits"            , 1 );



	TLorentzVector TL_leading_AK8_puppi;
	TLorentzVector TL_top_bjet;



   TFile * monoTopFile = new TFile(outFileName.c_str(),"RECREATE");
   TTree * monoTopTree = new TTree("monoTopTree","monoTopTree");


   //TLorentzVector leadingJet = TLorentzVector();
   //std::vector<float> px_v;
   //std::vector<float> py_v;
   float px, px_total;
   float py, py_total;
   int subjet_index;
   
   const float HT_cut_value                     = 110;
   const float MET_cut_value_minus_mu           = 110;

   const float leading_AK8Puppi_jet_pt_cut      = 250;
   const float leading_AK8Puppi_jet_eta_cut     = 2.5;

   const float btag_medium_wp                   = .68; //made up
   const float btag_loose_wp                   = .50; //made up
   const float mass_window[2]                   = {110, 210};
   const float tau3_tau2_wp                     = .6;
   const float QCD_delta_phi_cut_value          = 1.1;
   const float bjet_delta_r_cut_value           = 1.5;



   //cuts
   int   MET_selection_cut_minus_mu             =  -0;
   int   HT_selection_cut                       =  -0;
   int   leading_AK8Puppi_jet_cut               =  -0;
   int   top_jet_cut                            =  -0;
      int   top_jet_cut_bjet                       =  -0;
      int   top_jet_cut_mass                       =  -0;
      int   top_jet_cut_tau                       =  -0;
   int   QCD_delta_phi_cut                        =  -0;
   int   ttbar_bjet_cut                           =  -0;
   int   el_iso_cut                               =  -0;
   int   mu_iso_cut                               =  -0;
   int   tau_iso_cut                              =  -0;
   int   gamma_pt_cut                            =  -0;


   //values
   float MET_pt                        =  -100 ;
   float MET_pt_minus_mu               =  -100 ;
   float HT_pt                         =  -100 ;
   float HT_pt_minus_mu                =  -100 ;
   float leading_AK8Puppi_jet_pt       =  -100 ;
   float leading_AK8Puppi_jet_eta      =  -100 ;
   float leading_AK8Puppi_subjet0_CSVv2         =  -100 ;
   float leading_AK8Puppi_subjet1_CSVv2         =  -100 ;
   float leading_AK8Puppi_jet_mass              =  -100 ;
   float leading_AK8Puppi_jet_tau3_tau2         =  -100 ;
   float QCD_delta_phi_narrow_jets              =  -100 ;
   float bjet_delta_r_min                       =  -100 ;
   float el_iso_ratio                           =  -100 ;
   float mu_iso_ratio                           =  -100 ;
   float tau_iso_pt                          =  -100 ;
   float gamma_pt                        =  -100 ;


   //value branches 
   monoTopTree->Branch("HT_pt"     ,   &HT_pt      );
   monoTopTree->Branch("MET_pt_minus_mu"     ,   &MET_pt_minus_mu      );
   monoTopTree->Branch("MET_pt"     ,   &MET_pt      );
   monoTopTree->Branch("leading_AK8Puppi_jet_pt"     ,   &leading_AK8Puppi_jet_pt      );
   monoTopTree->Branch("leading_AK8Puppi_jet_eta"     ,   &leading_AK8Puppi_jet_eta      );
   monoTopTree->Branch("leading_AK8Puppi_subjet0_CSVv2"     ,   &leading_AK8Puppi_subjet0_CSVv2      );
   monoTopTree->Branch("leading_AK8Puppi_subjet1_CSVv2"     ,   &leading_AK8Puppi_subjet1_CSVv2      );
   monoTopTree->Branch("leading_AK8Puppi_jet_mass"     ,   &leading_AK8Puppi_jet_mass      );
   monoTopTree->Branch("leading_AK8Puppi_jet_tau3_tau2"     ,   &leading_AK8Puppi_jet_tau3_tau2      );
   monoTopTree->Branch("QCD_delta_phi_narrow_jets"     ,   &QCD_delta_phi_narrow_jets      );
   monoTopTree->Branch("bjet_delta_r_min"     ,   &bjet_delta_r_min      );
   monoTopTree->Branch("el_iso_ratio"     ,   &el_iso_ratio      );
   monoTopTree->Branch("mu_iso_ratio"     ,   &mu_iso_ratio      );
   monoTopTree->Branch("tau_iso_pt"     ,   &tau_iso_pt      );
   monoTopTree->Branch("gamma_pt"     ,   &gamma_pt      );
   


   //cut branches
   monoTopTree->Branch("MET_selection_cut_minus_mu"     ,   &MET_selection_cut_minus_mu      );
   monoTopTree->Branch("HT_selection_cut"     ,   &HT_selection_cut      );
   monoTopTree->Branch("leading_AK8Puppi_jet_cut"     ,   &leading_AK8Puppi_jet_cut      );
   monoTopTree->Branch("top_jet_cut"     ,   &top_jet_cut      );
   monoTopTree->Branch("top_jet_cut_bjet"     ,   &top_jet_cut_bjet      );
   monoTopTree->Branch("top_jet_cut_mass"     ,   &top_jet_cut_mass      );
   monoTopTree->Branch("top_jet_cut_tau"     ,   &top_jet_cut_tau      );
   monoTopTree->Branch("QCD_delta_phi_cut"     ,   &QCD_delta_phi_cut      );
   monoTopTree->Branch("ttbar_bjet_cut"     ,   &ttbar_bjet_cut      );
   monoTopTree->Branch("el_iso_cut"     ,   &el_iso_cut      );
   monoTopTree->Branch("mu_iso_cut"     ,   &mu_iso_cut      );
   monoTopTree->Branch("tau_iso_cut"     ,   &tau_iso_cut      );
   monoTopTree->Branch("gamma_pt_cut"     ,   &gamma_pt_cut      );






   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      MET_pt                        =  -100 ;
      MET_pt_minus_mu               =  -100 ;
      HT_pt                         =  -100 ;
      HT_pt_minus_mu                =  -100 ;
      leading_AK8Puppi_jet_pt       =  -100 ;
      leading_AK8Puppi_jet_eta      =  -100 ;
      leading_AK8Puppi_subjet0_CSVv2         =  -100 ;
      leading_AK8Puppi_subjet1_CSVv2         =  -100 ;
      leading_AK8Puppi_jet_mass              =  -100 ;
      leading_AK8Puppi_jet_tau3_tau2         =  -100 ;
      QCD_delta_phi_narrow_jets              =  -100 ;
      bjet_delta_r_min                       =  -100 ;
      el_iso_ratio                           =  -100 ;
      mu_iso_ratio                           =  -100 ;
      tau_iso_pt                          =  -100 ;
      gamma_pt                        =  -100 ; 

      MET_selection_cut_minus_mu             =  -0;
      HT_selection_cut                       =  -0;
      leading_AK8Puppi_jet_cut               =  -0;
      top_jet_cut                            =  -0;
      top_jet_cut_bjet                       =  -0;
      top_jet_cut_mass                       =  -0;
      top_jet_cut_tau                        =  -0;
      QCD_delta_phi_cut                      =  -0;
      ttbar_bjet_cut                         =  -0;
      el_iso_cut                               =  -0;
      mu_iso_cut                               =  -0;
      tau_iso_cut                              =  -0;
      gamma_pt_cut                            =  -0;



      ////////////////////////////////////////////////
      //HT miss 
      ////////////////////////////////////////////////

      //for (const pat::Jet &j : *puppijets) { 
       px_total  = 0;
      py_total = 0;

      for (int n=0; n<Jet_puppi_pt->size(); n++) {
         px_total = px_total + Jet_puppi_pt->at(n)*TMath::Cos(Jet_puppi_phi->at(n));
         py_total = py_total + Jet_puppi_pt->at(n)*TMath::Sin(Jet_puppi_phi->at(n));   
      }

      HT_pt = TMath::Sqrt(px_total*px_total + py_total*py_total);

      if(HT_pt > HT_cut_value) {
         HT_selection_cut = 1;
      }
      
      px_total  = 0;
      py_total = 0;

      ////////////////////////////////////////////////
      //MET miss minus mu
      ////////////////////////////////////////////////

      for (int n=0; n< Muon_pt->size(); n++) {
         px_total = px_total + Muon_pt->at(n)*TMath::Cos(Muon_phi->at(n));
         py_total = py_total + Muon_pt->at(n)*TMath::Sin(Muon_phi->at(n));   
     //    std::cout << n << " " << px_total << " " << py_total << " " << jetAK4CHS_pt->at(n)  <<std::endl;   
      }
     // std::cout << px_total << " " << py_total << " " << TMath::Sqrt(px_total*px_total + py_total*py_total)  <<std::endl;
      if(Met_puppi_pt){
         px = Met_puppi_pt*TMath::Cos(Met_puppi_phi);
         py = Met_puppi_pt*TMath::Sin(Met_puppi_phi);
         px_total =  px_total-px;
         py_total =  py_total-py;
         MET_pt_minus_mu = TMath::Sqrt(px_total*px_total + py_total*py_total);
         MET_pt = Met_puppi_pt; 
         if(MET_pt_minus_mu > MET_cut_value_minus_mu) {
            MET_selection_cut_minus_mu = 1;
         }
      }
      px_total  = 0;
      py_total = 0;

      ////////////////////////////////////////////////
      //leadingJet cut
      ////////////////////////////////////////////////

      if(Jet_AK8_puppi_pt->size() > 0){
         leading_AK8Puppi_jet_pt = Jet_AK8_puppi_pt->at(0);
         leading_AK8Puppi_jet_eta = Jet_AK8_puppi_eta->at(0);
         if(leading_AK8Puppi_jet_pt > leading_AK8Puppi_jet_pt_cut && TMath::Abs(leading_AK8Puppi_jet_eta) < leading_AK8Puppi_jet_eta_cut){
            leading_AK8Puppi_jet_cut = 1;
         }
      }

      ////////////////////////////////////////////////
      //topJet cut
      ////////////////////////////////////////////////


      if(Jet_AK8_puppi_pt->size() > 0){
      //btag cut
         if(Jet_AK8_subjet0_pt->size() > 0){
            leading_AK8Puppi_subjet0_CSVv2 = Jet_AK8_subjet0_CSVv2->at(0);
         }
 		if(Jet_AK8_subjet1_pt->size() > 0){
            leading_AK8Puppi_subjet1_CSVv2 = Jet_AK8_subjet1_CSVv2->at(0);
         }

         if(TMath::Max(leading_AK8Puppi_subjet0_CSVv2,leading_AK8Puppi_subjet1_CSVv2) > btag_medium_wp){
            top_jet_cut_bjet = 1;
         }

      //mass cut
         leading_AK8Puppi_jet_mass = Jet_AK8_puppi_mass->at(0);

         if(mass_window[0]< leading_AK8Puppi_jet_mass && leading_AK8Puppi_jet_mass< mass_window[1]){
            top_jet_cut_mass = 1;
         }


      //Tau cut
         leading_AK8Puppi_jet_tau3_tau2 = Jet_AK8_puppi_tau3->at(0)/Jet_AK8_puppi_tau2->at(0);
         if(leading_AK8Puppi_jet_tau3_tau2 > tau3_tau2_wp){
            top_jet_cut_tau = 1;
         }

         if( top_jet_cut_bjet*top_jet_cut_mass*top_jet_cut_tau > 0){
            top_jet_cut = 1;
         }
      }



      ////////////////////////////////////////////////
      //QCD delta met, jet cut
      ////////////////////////////////////////////////

      

      if(Met_puppi_pt > 0 and Jet_puppi_pt->size())
      {
         QCD_delta_phi_narrow_jets = 100;
         for(unsigned int n = 0; n < Jet_puppi_pt->size(); n++){
            float temp_delta_phi =  TMath::Abs(Met_puppi_phi - Jet_puppi_phi->at(n));
            if( Jet_puppi_phi->at(n) < 30){
               break;
            }
            if( Jet_puppi_pt->at(n) > 30 && TMath::Abs(Jet_puppi_eta->at(n)) < 4.5 ){
                QCD_delta_phi_narrow_jets = TMath::Min(QCD_delta_phi_narrow_jets, temp_delta_phi );
            }
         }
         if(QCD_delta_phi_narrow_jets > QCD_delta_phi_cut_value){
            QCD_delta_phi_cut = 1;
         }
      } 


      ////////////////////////////////////////////////
      //ttbar delta r bjet, fat jet
      ////////////////////////////////////////////////
     

      if(Jet_puppi_pt->size() > 0 && Jet_AK8_puppi_pt->size() > 0){

         for(unsigned int n = 0; n < Jet_puppi_pt->size(); n++){
            if( Jet_puppi_bDiscriminator_pfCISVV2->at(n) > btag_loose_wp){
               int not_subjet = 1;
               float delta_r = 100;

               //check subjets
         /*      if(subJet_puppi_size >= Jet_puppi_vSubjetIndex0->at(0) and Jet_puppi_vSubjetIndex0->at(0) >= 0){
                  subjet_index = Jet_puppi_vSubjetIndex0->at(0);
                  delta_r = TMath::Sqrt( TMath::Power( jetAK4CHS_phi->at(n) - subJet_puppi_Phi[subjet_index], 2  ) + TMath::Power( jetAK4CHS_Eta->at(n) - subJet_puppi_Eta[subjet_index], 2  ));
                  std::cout << "subjet0 " << delta_r <<std::endl;
                  if( delta_r < 0.4){
                     not_subjet = 0;
                  }
               }

               if(subJet_puppi_size >= Jet_puppi_vSubjetIndex1->at(0) and Jet_puppi_vSubjetIndex1->at(0) >= 0){
                  subjet_index = Jet_puppi_vSubjetIndex1->at(0);
                  delta_r = TMath::Sqrt( TMath::Power( jetAK4CHS_phi->at(n) - subJet_puppi_Phi[subjet_index], 2  ) + TMath::Power( jetAK4CHS_Eta->at(n) - subJet_puppi_Eta[subjet_index], 2  ));
                  std::cout << "subjet1 " << delta_r <<std::endl;
                  if( delta_r < 0.4){
                     not_subjet = 0;
                  }
               }*/

               if( not_subjet == 1){
                  delta_r = TMath::Sqrt( TMath::Power( Jet_puppi_phi->at(n) - Jet_AK8_puppi_phi->at(0), 2  ) + TMath::Power( Jet_puppi_eta->at(n) - Jet_AK8_puppi_eta->at(0), 2  ));
               }
              // std::cout<< " " << jetAK4CHS_phi->at(n) << " " << jetAK4CHS_Eta->at(n) << " " << delta_r << std::endl;
               float temp_delta_r = TMath::Abs(bjet_delta_r_min);
               if(delta_r > 0.8 ){
                  bjet_delta_r_min = TMath::Min(delta_r, temp_delta_r ) ;
               }
            }
         }

         if(bjet_delta_r_min > bjet_delta_r_cut_value || bjet_delta_r_min < 0){
            ttbar_bjet_cut = 1;
            //std::cout <<"what " << bjet_delta_r_min << std::endl;
         }

         //std::cout << ttbar_bjet_cut << " " << bjet_delta_r_min << std::endl;
      }
      ////////////////////////////////////////////////
      //iso cuts
      ////////////////////////////////////////////////       



      el_iso_cut = 1;
      for(unsigned int n = 0; n < patElectron_pt->size(); n++){
         if(patElectron_pt->at(n)>10){
            el_iso_ratio = patElectron_isoPU->at(n)/patElectron_energy->at(n);
            if(TMath::Abs(patElectron_eta->at(n) < 1.479 && el_iso_ratio > .126)){
               el_iso_cut = 0; 
               break;
            }
            if(TMath::Abs(patElectron_eta->at(n) > 1.479 && TMath::Abs(patElectron_eta->at(n)) > 2.5 && el_iso_ratio > .144)){
               el_iso_cut = 0;
               break;
            }
   
         }
      }
   
      mu_iso_cut = 1;
      for(unsigned int n = 0; n < Muon_pt->size(); n++){
         if(Muon_pt->at(n)>10){
            mu_iso_ratio = Muon_isoPU->at(n)/Muon_energy->at(n);
            if(mu_iso_ratio > .2){
               mu_iso_cut = 0;
               break;
            }
   
         }
      }
   
      tau_iso_cut = 1;
      for(unsigned int n = 0; n < Tau_pt->size(); n++){
          tau_iso_pt = (Tau_pt->at(n)- Tau_puCorrPtSum->at(n));
         if(tau_iso_pt > 5 &&  Tau_pt->at(n) >18){
            tau_iso_cut = 0;
            break;
         }
      }
      
      
      gamma_pt_cut = 1;
      if(Photon_pt->size() > 0){
         gamma_pt =  Photon_pt->at(0);
         if(gamma_pt > 15){
            gamma_pt_cut = 0;
         }
      }





         
      //std::cout << HT_pt << std::endl;


      monoTopTree->Fill();

      // if (Cut(ientry) < 0) continue;

     // leading_AK8_puppi.SetPtEtaPhiM(Jet_AK8_puppi_pt->at(0), Jet_AK8_puppi_eta->at(0), Jet_AK8_puppi_phi->at(0), Jet_AK8_puppi_mass->at(0)  )
      //std::cout << leading_AK8_puppi.Energy() << " " << Jet_AK8_GEN_energy->at(0) << std::endl;
      //std::cout << Jet_AK8_puppi_pt->at(0) << std::endl;
   }
   monoTopFile->cd();
   monoTopFile->Write();
   monoTopFile->Close();
}
