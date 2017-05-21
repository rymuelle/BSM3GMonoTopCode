//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 11 16:15:32 2017 by ROOT version 6.08/05
// from TTree BOOM/BOOM
// found on file: OutTree.root
//////////////////////////////////////////////////////////

#ifndef BSM3G_MonoTop_Analyzer_h
#define BSM3G_MonoTop_Analyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class BSM3G_MonoTop_Analyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<int>     *Trigger_decision;
   vector<string>  *Trigger_names;
   vector<int>     *Trigger_prescale;
   vector<double>  *Muon_pt;
   vector<double>  *Muon_eta;
   vector<double>  *Muon_phi;
   vector<double>  *Muon_energy;
   vector<double>  *Muon_charge;
   vector<bool>    *Muon_tight;
   vector<bool>    *Muon_loose;
   vector<bool>    *Muon_medium;
   vector<bool>    *Muon_soft;
   vector<bool>    *Muon_isHighPt;
   vector<bool>    *Muon_pf;
   vector<double>  *Muon_isoCharged;
   vector<double>  *Muon_isoSum;
   vector<double>  *Muon_isoCharParPt;
   vector<double>  *Muon_isTrackerMuon;
   vector<double>  *Muon_POGisGood;
   vector<double>  *Muon_chi2;
   vector<double>  *Muon_validHits;
   vector<double>  *Muon_validHitsInner;
   vector<double>  *Muon_matchedStat;
   vector<double>  *Muon_dxy_pv;
   vector<double>  *Muon_dxy_bs;
   vector<double>  *Muon_dz_pv;
   vector<double>  *Muon_dz_bs;
   vector<double>  *Muon_dxyError;
   vector<double>  *Muon_dzError;
   vector<double>  *Muon_ndof;
   vector<double>  *Muon_vtx;
   vector<double>  *Muon_vty;
   vector<double>  *Muon_vtz;
   vector<double>  *Muon_track_pt;
   vector<double>  *Muon_track_ptError;
   vector<bool>    *Muon_isGlobal;
   vector<double>  *Muon_TLayers;
   vector<double>  *Muon_isoNeutralHadron;
   vector<double>  *Muon_isoPhoton;
   vector<double>  *Muon_isoPU;
   vector<double>  *Muon_combinedIso;
   vector<double>  *Muon_trackRe_iso;
   vector<double>  *Muon_dB;
   vector<double>  *Muon_besttrack_pt;
   vector<double>  *Muon_besttrack_ptError;
   vector<double>  *Muon_tunePBestTrack_pt;
   vector<double>  *Muon_tunePBestTrackType;
   vector<double>  *Muon_chi2LocalPosition;
   vector<double>  *Muon_trkKink;
   vector<double>  *Muon_segmentCompatibility;
   vector<double>  *Muon_validFraction;
   vector<double>  *Muon_pixelLayersWithMeasurement;
   vector<double>  *Muon_qualityhighPurity;
   vector<double>  *Muon_track_PCAx_bs;
   vector<double>  *Muon_track_PCAy_bs;
   vector<double>  *Muon_track_PCAz_bs;
   vector<double>  *Muon_track_PCAx_pv;
   vector<double>  *Muon_track_PCAy_pv;
   vector<double>  *Muon_track_PCAz_pv;
   vector<double>  *Muon_trackFitErrorMatrix_00;
   vector<double>  *Muon_trackFitErrorMatrix_01;
   vector<double>  *Muon_trackFitErrorMatrix_02;
   vector<double>  *Muon_trackFitErrorMatrix_11;
   vector<double>  *Muon_trackFitErrorMatrix_12;
   vector<double>  *Muon_trackFitErrorMatrix_22;
   vector<double>  *patElectron_pt;
   vector<double>  *patElectron_eta;
   vector<double>  *patElectron_phi;
   vector<double>  *patElectron_energy;
   vector<double>  *patElectron_charge;
   vector<int>     *patElectron_isPassVeto;
   vector<int>     *patElectron_isPassLoose;
   vector<int>     *patElectron_isPassMedium;
   vector<int>     *patElectron_isPassTight;
   vector<int>     *patElectron_isPassHEEPId;
   vector<int>     *patElectron_passMV1wp1Id;
   vector<int>     *patElectron_passMV2wp1Id;
   vector<double>  *patElectron_isoChargedHadrons;
   vector<double>  *patElectron_isoNeutralHadrons;
   vector<double>  *patElectron_isoPhotons;
   vector<double>  *patElectron_isoPU;
   vector<int>     *patElectron_expectedMissingInnerHits;
   vector<int>     *patElectron_passConversionVeto;
   vector<double>  *patElectron_gsfTrack_dxy_pv;
   vector<double>  *patElectron_gsfTrack_dxy_bs;
   vector<double>  *patElectron_dxyError;
   vector<double>  *patElectron_gsfTrack_dz_pv;
   vector<double>  *patElectron_gsfTrack_dz_bs;
   vector<double>  *patElectron_gsfTrack_normChi2;
   vector<double>  *patElectron_gsfTrack_ndof;
   vector<double>  *patElectron_gsfTrack_vtx;
   vector<double>  *patElectron_gsfTrack_vty;
   vector<double>  *patElectron_gsfTrack_vtz;
   vector<double>  *patElectron_gsfTrack_PCAx_bs;
   vector<double>  *patElectron_gsfTrack_PCAy_bs;
   vector<double>  *patElectron_gsfTrack_PCAz_bs;
   vector<double>  *patElectron_gsfTrack_PCAx_pv;
   vector<double>  *patElectron_gsfTrack_PCAy_pv;
   vector<double>  *patElectron_gsfTrack_PCAz_pv;
   vector<double>  *patElectron_gsfTrackFitErrorMatrix_00;
   vector<double>  *patElectron_gsfTrackFitErrorMatrix_01;
   vector<double>  *patElectron_gsfTrackFitErrorMatrix_02;
   vector<double>  *patElectron_gsfTrackFitErrorMatrix_11;
   vector<double>  *patElectron_gsfTrackFitErrorMatrix_12;
   vector<double>  *patElectron_gsfTrackFitErrorMatrix_22;
   vector<double>  *Tau_eta;
   vector<double>  *Tau_phi;
   vector<double>  *Tau_pt;
   vector<double>  *Tau_energy;
   vector<double>  *Tau_charge;
   vector<int>     *Tau_decayMode;
   vector<double>  *Tau_nProngs;
   vector<double>  *Tau_chargedIsoPtSum;
   vector<double>  *Tau_neutralIsoPtSum;
   vector<double>  *Tau_puCorrPtSum;
   vector<double>  *Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;
   vector<int>     *Tau_decayModeFinding;
   vector<int>     *Tau_decayModeFindingNewDMs;
   vector<int>     *Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *Tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *Tau_byVLooseIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *Tau_byLooseIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *Tau_byMediumIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *Tau_byTightIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *Tau_byVTightIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *Tau_byVLooseIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *Tau_byLooseIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *Tau_byMediumIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *Tau_byTightIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *Tau_byVTightIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
   vector<int>     *Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   vector<int>     *Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;
   //vector<int>     *Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<int>     *Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<int>     *Tau_againstMuonLoose3;
   vector<int>     *Tau_againstMuonTight3;
   vector<int>     *Tau_againstElectronMVAVLooseMVA6;
   vector<int>     *Tau_againstElectronMVALooseMVA6;
   vector<int>     *Tau_againstElectronMVAMediumMVA6;
   vector<int>     *Tau_againstElectronMVATightMVA6;
   vector<double>  *Tau_leadChargedCandPt;
   vector<double>  *Tau_leadChargedCandEta;
   vector<double>  *Tau_leadChargedCandPhi;
   vector<double>  *Tau_leadChargedCandCharge;
   vector<double>  *Tau_leadChargedCandChi2;
   vector<double>  *Tau_leadChargedCandValidHits;
   vector<double>  *Tau_leadChargedCandDxy_pv;
   vector<double>  *Tau_leadChargedCandDxy_bs;
   vector<double>  *Tau_leadChargedCandDz_pv;
   vector<double>  *Tau_leadChargedCandDz_bs;
   vector<double>  *Tau_leadChargedCandDzError;
   vector<double>  *Tau_leadChargedCandDxyError;
   vector<double>  *Tau_leadChargedCandNdof;
   vector<double>  *Tau_leadChargedCandVtx;
   vector<double>  *Tau_leadChargedCandVty;
   vector<double>  *Tau_leadChargedCandVtz;
   vector<double>  *Tau_leadChargedCandTrack_pt;
   vector<double>  *Tau_leadChargedCandTrack_ptError;
   vector<double>  *Tau_leadChargedCandTrack_PCAx_bs;
   vector<double>  *Tau_leadChargedCandTrack_PCAy_bs;
   vector<double>  *Tau_leadChargedCandTrack_PCAz_bs;
   vector<double>  *Tau_leadChargedCandTrack_PCAx_pv;
   vector<double>  *Tau_leadChargedCandTrack_PCAy_pv;
   vector<double>  *Tau_leadChargedCandTrack_PCAz_pv;
   vector<double>  *Tau_leadChargedCandTrackFitErrorMatrix_00;
   vector<double>  *Tau_leadChargedCandTrackFitErrorMatrix_01;
   vector<double>  *Tau_leadChargedCandTrackFitErrorMatrix_02;
   vector<double>  *Tau_leadChargedCandTrackFitErrorMatrix_11;
   vector<double>  *Tau_leadChargedCandTrackFitErrorMatrix_12;
   vector<double>  *Tau_leadChargedCandTrackFitErrorMatrix_22;
   vector<double>  *Tau_defaultDxy;
   vector<double>  *Tau_defaultDxyError;
   vector<double>  *Tau_defaultDxySig;
   vector<double>  *Tau_defaultFlightLengthX;
   vector<double>  *Tau_defaultFlightLengthY;
   vector<double>  *Tau_defaultFlightLengthZ;
   vector<double>  *Tau_defaultFlightLengthSig;
   vector<double>  *Tau_default_PCAx_pv;
   vector<double>  *Tau_default_PCAy_pv;
   vector<double>  *Tau_default_PCAz_pv;
   vector<double>  *Jet_pt;
   vector<double>  *Jet_eta;
   vector<double>  *Jet_phi;
   vector<double>  *Jet_energy;
   vector<double>  *Jet_bDiscriminator;
   vector<double>  *Jet_bDiscriminator_CISVV2;
   vector<double>  *Jet_bDiscriminator_pfCISVV2;
   vector<int>     *Jet_partonFlavour;
   vector<double>  *Jet_neutralHadEnergyFraction;
   vector<double>  *Jet_neutralEmEmEnergyFraction;
   vector<double>  *Jet_chargedHadronEnergyFraction;
   vector<double>  *Jet_chargedEmEnergyFraction;
   vector<double>  *Jet_muonEnergyFraction;
   vector<int>     *Jet_numberOfConstituents;
   vector<int>     *Jet_chargedMultiplicity;
   vector<double>  *Jet_puppi_pt;
   vector<double>  *Jet_puppi_eta;
   vector<double>  *Jet_puppi_phi;
   vector<double>  *Jet_puppi_energy;
   vector<double>  *Jet_puppi_bDiscriminator;
   vector<double>  *Jet_puppi_bDiscriminator_CISVV2;
   vector<double>  *Jet_puppi_bDiscriminator_pfCISVV2;
   vector<int>     *Jet_puppi_partonFlavour;
   vector<double>  *Jet_puppi_neutralHadEnergyFraction;
   vector<double>  *Jet_puppi_neutralEmEmEnergyFraction;
   vector<double>  *Jet_puppi_chargedHadronEnergyFraction;
   vector<double>  *Jet_puppi_chargedEmEnergyFraction;
   vector<double>  *Jet_puppi_muonEnergyFraction;
   vector<int>     *Jet_puppi_numberOfConstituents;
   vector<int>     *Jet_puppi_chargedMultiplicity;
   vector<double>  *Jet_AK8_pt;
   vector<double>  *Jet_AK8_eta;
   vector<double>  *Jet_AK8_phi;
   vector<double>  *Jet_AK8_energy;
   vector<double>  *Jet_AK8_bDiscriminator;
   vector<double>  *Jet_AK8_bDiscriminator_CISVV2;
   vector<double>  *Jet_AK8_bDiscriminator_pfCISVV2;
   vector<int>     *Jet_AK8_partonFlavour;
   vector<double>  *Jet_AK8_neutralHadEnergyFraction;
   vector<double>  *Jet_AK8_neutralEmEmEnergyFraction;
   vector<double>  *Jet_AK8_chargedHadronEnergyFraction;
   vector<double>  *Jet_AK8_chargedEmEnergyFraction;
   vector<double>  *Jet_AK8_muonEnergyFraction;
   vector<int>     *Jet_AK8_numberOfConstituents;
   vector<int>     *Jet_AK8_chargedMultiplicity;
   vector<double>  *Jet_AK8_puppi_pt;
   vector<double>  *Jet_AK8_puppi_eta;
   vector<double>  *Jet_AK8_puppi_mass;
   vector<double>  *Jet_AK8_puppi_phi;
   vector<double>  *Jet_AK8_puppi_tau1;
   vector<double>  *Jet_AK8_puppi_tau2;
   vector<double>  *Jet_AK8_puppi_tau3;
   vector<double>  *Jet_AK8_GEN_pt;
   vector<double>  *Jet_AK8_GEN_phi;
   vector<double>  *Jet_AK8_GEN_eta;
   vector<double>  *Jet_AK8_GEN_mass;
   vector<double>  *Jet_AK8_GEN_energy;
   vector<double>  *Jet_AK8_GEN_parton;
   vector<double>  *Jet_AK8_GEN_mother;
   vector<double>  *Jet_AK8_subjet0_pt;
   vector<double>  *Jet_AK8_subjet0_phi;
   vector<double>  *Jet_AK8_subjet0_eta;
   vector<double>  *Jet_AK8_subjet0_mass;
   vector<double>  *Jet_AK8_subjet0_energy;
   vector<double>  *Jet_AK8_subjet0_CSVv2;
   vector<double>  *Jet_AK8_subjet1_pt;
   vector<double>  *Jet_AK8_subjet1_phi;
   vector<double>  *Jet_AK8_subjet1_eta;
   vector<double>  *Jet_AK8_subjet1_mass;
   vector<double>  *Jet_AK8_subjet1_energy;
   vector<double>  *Jet_AK8_subjet1_CSVv2;
   vector<double>  *Jet_mass;
   vector<double>  *Jet_electronEnergy;
   vector<double>  *Jet_photonEnergy;
   vector<double>  *UncorrJet_pt;
   vector<double>  *Jet_puppi_mass;
   vector<double>  *Jet_puppi_electronEnergy;
   vector<double>  *Jet_puppi_photonEnergy;
   vector<double>  *UncorrJet_puppi_pt;
   vector<double>  *Jet_AK8_mass;
   vector<double>  *Jet_AK8_electronEnergy;
   vector<double>  *Jet_AK8_photonEnergy;
   vector<double>  *UncorrJet_AK8_pt;
   vector<double>  *Gen_pt;
   vector<double>  *Gen_eta;
   vector<double>  *Gen_phi;
   vector<double>  *Gen_status;
   vector<double>  *Gen_pdg_id;
   vector<double>  *Gen_motherpdg_id;
   vector<double>  *Gen_energy;
   vector<double>  *Gen_vx;
   vector<double>  *Gen_vy;
   vector<double>  *Gen_vz;
   vector<double>  *Gen_charge;
   vector<double>  *Gen_numDaught;
   vector<double>  *Gen_numMother;
   vector<int>     *Gen_BmotherIndices;
   vector<int>     *Gen_BdaughtIndices;
   vector<int>     *Gen_BmotherIndex;
   Double_t        weightevt;
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           lumiBlock;
   Int_t           nObservedInTimePUVertices;
   Float_t         nTruePUInteractions;
   Int_t           nObservedOutOfTimePUVertices;
   Int_t           nObservedPlus1BXPUVertices;
   Int_t           nObservedMinus1BXPUVertices;
   Int_t           bestVertices;
   vector<float>   *pvertex_x;
   vector<float>   *pvertex_y;
   vector<float>   *pvertex_z;
   vector<float>   *pvertex_xError;
   vector<float>   *pvertex_yError;
   vector<float>   *pvertex_zError;
   vector<float>   *beamSpot_x0;
   vector<float>   *beamSpot_y0;
   vector<float>   *beamSpot_z0;
   vector<float>   *beamSpot_xWidth;
   vector<float>   *beamSpot_yWidth;
   Double_t        Met_type1PF_pt;
   Double_t        Met_type1PF_px;
   Double_t        Met_type1PF_py;
   Double_t        Met_type1PF_pz;
   Double_t        Met_type1PF_sumEt;
   Double_t        Met_type1PF_phi;
   Double_t        Met_type1PF_cov00;
   Double_t        Met_type1PF_cov01;
   Double_t        Met_type1PF_cov10;
   Double_t        Met_type1PF_cov11;
   Double_t        Met_puppi_pt;
   Double_t        Met_puppi_px;
   Double_t        Met_puppi_py;
   Double_t        Met_puppi_pz;
   Double_t        Met_puppi_sumEt;
   Double_t        Met_puppi_phi;
   Double_t        Met_NoHF_pt;
   Double_t        Met_NoHF_px;
   Double_t        Met_NoHF_py;
   Double_t        Met_NoHF_pz;
   Double_t        Met_NoHF_sumEt;
   Double_t        Met_NoHF_phi;
   Double_t        Gen_Met;
   Double_t        Met_type1PF_shiftedPtUp;
   Double_t        Met_type1PF_shiftedPtDown;
   vector<float>   *Photon_pt;
   vector<float>   *Photon_eta;
   vector<float>   *Photon_phi;
   vector<float>   *Photon_energy;
   vector<float>   *Photon_et;
   vector<float>   *Photon_HoverE;
   vector<float>   *Photon_phoR9;
   vector<float>   *Photon_SigmaIEtaIEta;
   vector<float>   *Photon_SigmaIPhiIPhi;
   vector<float>   *Photon_PFChIso;
   vector<float>   *Photon_PFPhoIso;
   vector<float>   *Photon_PFNeuIso;
   vector<int>     *Photon_EleVeto;
   vector<int>     *Photon_hasPixelSeed;

   // List of branches
   TBranch        *b_Trigger_decision;   //!
   TBranch        *b_Trigger_names;   //!
   TBranch        *b_Trigger_prescale;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_tight;   //!
   TBranch        *b_Muon_loose;   //!
   TBranch        *b_Muon_medium;   //!
   TBranch        *b_Muon_soft;   //!
   TBranch        *b_Muon_isHighPt;   //!
   TBranch        *b_Muon_pf;   //!
   TBranch        *b_Muon_isoCharged;   //!
   TBranch        *b_Muon_isoSum;   //!
   TBranch        *b_Muon_isoCharParPt;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_Muon_POGisGood;   //!
   TBranch        *b_Muon_chi2;   //!
   TBranch        *b_Muon_validHits;   //!
   TBranch        *b_Muon_validHitsInner;   //!
   TBranch        *b_Muon_matchedStat;   //!
   TBranch        *b_Muon_dxy_pv;   //!
   TBranch        *b_Muon_dxy_bs;   //!
   TBranch        *b_Muon_dz_pv;   //!
   TBranch        *b_Muon_dz_bs;   //!
   TBranch        *b_Muon_dxyError;   //!
   TBranch        *b_Muon_dzError;   //!
   TBranch        *b_Muon_ndof;   //!
   TBranch        *b_Muon_vtx;   //!
   TBranch        *b_Muon_vty;   //!
   TBranch        *b_Muon_vtz;   //!
   TBranch        *b_Muon_track_pt;   //!
   TBranch        *b_Muon_track_ptError;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_TLayers;   //!
   TBranch        *b_Muon_isoNeutralHadron;   //!
   TBranch        *b_Muon_isoPhoton;   //!
   TBranch        *b_Muon_isoPU;   //!
   TBranch        *b_Muon_combinedIso;   //!
   TBranch        *b_Muon_trackRe_iso;   //!
   TBranch        *b_Muon_dB;   //!
   TBranch        *b_Muon_besttrack_pt;   //!
   TBranch        *b_Muon_besttrack_ptError;   //!
   TBranch        *b_Muon_tunePBestTrack_pt;   //!
   TBranch        *b_Muon_tunePBestTrackType;   //!
   TBranch        *b_Muon_chi2LocalPosition;   //!
   TBranch        *b_Muon_trkKink;   //!
   TBranch        *b_Muon_segmentCompatibility;   //!
   TBranch        *b_Muon_validFraction;   //!
   TBranch        *b_Muon_pixelLayersWithMeasurement;   //!
   TBranch        *b_Muon_qualityhighPurity;   //!
   TBranch        *b_Muon_track_PCAx_bs;   //!
   TBranch        *b_Muon_track_PCAy_bs;   //!
   TBranch        *b_Muon_track_PCAz_bs;   //!
   TBranch        *b_Muon_track_PCAx_pv;   //!
   TBranch        *b_Muon_track_PCAy_pv;   //!
   TBranch        *b_Muon_track_PCAz_pv;   //!
   TBranch        *b_Muon_trackFitErrorMatrix_00;   //!
   TBranch        *b_Muon_trackFitErrorMatrix_01;   //!
   TBranch        *b_Muon_trackFitErrorMatrix_02;   //!
   TBranch        *b_Muon_trackFitErrorMatrix_11;   //!
   TBranch        *b_Muon_trackFitErrorMatrix_12;   //!
   TBranch        *b_Muon_trackFitErrorMatrix_22;   //!
   TBranch        *b_patElectron_pt;   //!
   TBranch        *b_patElectron_eta;   //!
   TBranch        *b_patElectron_phi;   //!
   TBranch        *b_patElectron_energy;   //!
   TBranch        *b_patElectron_charge;   //!
   TBranch        *b_patElectron_isPassVeto;   //!
   TBranch        *b_patElectron_isPassLoose;   //!
   TBranch        *b_patElectron_isPassMedium;   //!
   TBranch        *b_patElectron_isPassTight;   //!
   TBranch        *b_patElectron_isPassHEEPId;   //!
   TBranch        *b_patElectron_passMV1wp1Id;   //!
   TBranch        *b_patElectron_passMV2wp1Id;   //!
   TBranch        *b_patElectron_isoChargedHadrons;   //!
   TBranch        *b_patElectron_isoNeutralHadrons;   //!
   TBranch        *b_patElectron_isoPhotons;   //!
   TBranch        *b_patElectron_isoPU;   //!
   TBranch        *b_patElectron_expectedMissingInnerHits;   //!
   TBranch        *b_patElectron_passConversionVeto;   //!
   TBranch        *b_patElectron_gsfTrack_dxy_pv;   //!
   TBranch        *b_patElectron_gsfTrack_dxy_bs;   //!
   TBranch        *b_patElectron_dxyError;   //!
   TBranch        *b_patElectron_gsfTrack_dz_pv;   //!
   TBranch        *b_patElectron_gsfTrack_dz_bs;   //!
   TBranch        *b_patElectron_gsfTrack_normChi2;   //!
   TBranch        *b_patElectron_gsfTrack_ndof;   //!
   TBranch        *b_patElectron_gsfTrack_vtx;   //!
   TBranch        *b_patElectron_gsfTrack_vty;   //!
   TBranch        *b_patElectron_gsfTrack_vtz;   //!
   TBranch        *b_patElectron_gsfTrack_PCAx_bs;   //!
   TBranch        *b_patElectron_gsfTrack_PCAy_bs;   //!
   TBranch        *b_patElectron_gsfTrack_PCAz_bs;   //!
   TBranch        *b_patElectron_gsfTrack_PCAx_pv;   //!
   TBranch        *b_patElectron_gsfTrack_PCAy_pv;   //!
   TBranch        *b_patElectron_gsfTrack_PCAz_pv;   //!
   TBranch        *b_patElectron_gsfTrackFitErrorMatrix_00;   //!
   TBranch        *b_patElectron_gsfTrackFitErrorMatrix_01;   //!
   TBranch        *b_patElectron_gsfTrackFitErrorMatrix_02;   //!
   TBranch        *b_patElectron_gsfTrackFitErrorMatrix_11;   //!
   TBranch        *b_patElectron_gsfTrackFitErrorMatrix_12;   //!
   TBranch        *b_patElectron_gsfTrackFitErrorMatrix_22;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_energy;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_decayMode;   //!
   TBranch        *b_Tau_nProngs;   //!
   TBranch        *b_Tau_chargedIsoPtSum;   //!
   TBranch        *b_Tau_neutralIsoPtSum;   //!
   TBranch        *b_Tau_puCorrPtSum;   //!
   TBranch        *b_Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_Tau_decayModeFinding;   //!
   TBranch        *b_Tau_decayModeFindingNewDMs;   //!
   TBranch        *b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byVLooseIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_Tau_byLooseIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_Tau_byTightIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_Tau_byVTightIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_Tau_byVLooseIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_Tau_byLooseIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_Tau_byTightIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_Tau_byVTightIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   //TBranch        *b_Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_Tau_againstMuonLoose3;   //!
   TBranch        *b_Tau_againstMuonTight3;   //!
   TBranch        *b_Tau_againstElectronMVAVLooseMVA6;   //!
   TBranch        *b_Tau_againstElectronMVALooseMVA6;   //!
   TBranch        *b_Tau_againstElectronMVAMediumMVA6;   //!
   TBranch        *b_Tau_againstElectronMVATightMVA6;   //!
   TBranch        *b_Tau_leadChargedCandPt;   //!
   TBranch        *b_Tau_leadChargedCandEta;   //!
   TBranch        *b_Tau_leadChargedCandPhi;   //!
   TBranch        *b_Tau_leadChargedCandCharge;   //!
   TBranch        *b_Tau_leadChargedCandChi2;   //!
   TBranch        *b_Tau_leadChargedCandValidHits;   //!
   TBranch        *b_Tau_leadChargedCandDxy_pv;   //!
   TBranch        *b_Tau_leadChargedCandDxy_bs;   //!
   TBranch        *b_Tau_leadChargedCandDz_pv;   //!
   TBranch        *b_Tau_leadChargedCandDz_bs;   //!
   TBranch        *b_Tau_leadChargedCandDzError;   //!
   TBranch        *b_Tau_leadChargedCandDxyError;   //!
   TBranch        *b_Tau_leadChargedCandNdof;   //!
   TBranch        *b_Tau_leadChargedCandVtx;   //!
   TBranch        *b_Tau_leadChargedCandVty;   //!
   TBranch        *b_Tau_leadChargedCandVtz;   //!
   TBranch        *b_Tau_leadChargedCandTrack_pt;   //!
   TBranch        *b_Tau_leadChargedCandTrack_ptError;   //!
   TBranch        *b_Tau_leadChargedCandTrack_PCAx_bs;   //!
   TBranch        *b_Tau_leadChargedCandTrack_PCAy_bs;   //!
   TBranch        *b_Tau_leadChargedCandTrack_PCAz_bs;   //!
   TBranch        *b_Tau_leadChargedCandTrack_PCAx_pv;   //!
   TBranch        *b_Tau_leadChargedCandTrack_PCAy_pv;   //!
   TBranch        *b_Tau_leadChargedCandTrack_PCAz_pv;   //!
   TBranch        *b_Tau_leadChargedCandTrackFitErrorMatrix_00;   //!
   TBranch        *b_Tau_leadChargedCandTrackFitErrorMatrix_01;   //!
   TBranch        *b_Tau_leadChargedCandTrackFitErrorMatrix_02;   //!
   TBranch        *b_Tau_leadChargedCandTrackFitErrorMatrix_11;   //!
   TBranch        *b_Tau_leadChargedCandTrackFitErrorMatrix_12;   //!
   TBranch        *b_Tau_leadChargedCandTrackFitErrorMatrix_22;   //!
   TBranch        *b_Tau_defaultDxy;   //!
   TBranch        *b_Tau_defaultDxyError;   //!
   TBranch        *b_Tau_defaultDxySig;   //!
   TBranch        *b_Tau_defaultFlightLengthX;   //!
   TBranch        *b_Tau_defaultFlightLengthY;   //!
   TBranch        *b_Tau_defaultFlightLengthZ;   //!
   TBranch        *b_Tau_defaultFlightLengthSig;   //!
   TBranch        *b_Tau_default_PCAx_pv;   //!
   TBranch        *b_Tau_default_PCAy_pv;   //!
   TBranch        *b_Tau_default_PCAz_pv;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_bDiscriminator;   //!
   TBranch        *b_Jet_bDiscriminator_CISVV2;   //!
   TBranch        *b_Jet_bDiscriminator_pfCISVV2;   //!
   TBranch        *b_Jet_partonFlavour;   //!
   TBranch        *b_Jet_neutralHadEnergyFraction;   //!
   TBranch        *b_Jet_neutralEmEmEnergyFraction;   //!
   TBranch        *b_Jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jet_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_muonEnergyFraction;   //!
   TBranch        *b_Jet_numberOfConstituents;   //!
   TBranch        *b_Jet_chargedMultiplicity;   //!
   TBranch        *b_Jet_puppi_pt;   //!
   TBranch        *b_Jet_puppi_eta;   //!
   TBranch        *b_Jet_puppi_phi;   //!
   TBranch        *b_Jet_puppi_energy;   //!
   TBranch        *b_Jet_puppi_bDiscriminator;   //!
   TBranch        *b_Jet_puppi_bDiscriminator_CISVV2;   //!
   TBranch        *b_Jet_puppi_bDiscriminator_pfCISVV2;   //!
   TBranch        *b_Jet_puppi_partonFlavour;   //!
   TBranch        *b_Jet_puppi_neutralHadEnergyFraction;   //!
   TBranch        *b_Jet_puppi_neutralEmEmEnergyFraction;   //!
   TBranch        *b_Jet_puppi_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jet_puppi_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_puppi_muonEnergyFraction;   //!
   TBranch        *b_Jet_puppi_numberOfConstituents;   //!
   TBranch        *b_Jet_puppi_chargedMultiplicity;   //!
   TBranch        *b_Jet_AK8_pt;   //!
   TBranch        *b_Jet_AK8_eta;   //!
   TBranch        *b_Jet_AK8_phi;   //!
   TBranch        *b_Jet_AK8_energy;   //!
   TBranch        *b_Jet_AK8_bDiscriminator;   //!
   TBranch        *b_Jet_AK8_bDiscriminator_CISVV2;   //!
   TBranch        *b_Jet_AK8_bDiscriminator_pfCISVV2;   //!
   TBranch        *b_Jet_AK8_partonFlavour;   //!
   TBranch        *b_Jet_AK8_neutralHadEnergyFraction;   //!
   TBranch        *b_Jet_AK8_neutralEmEmEnergyFraction;   //!
   TBranch        *b_Jet_AK8_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jet_AK8_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_AK8_muonEnergyFraction;   //!
   TBranch        *b_Jet_AK8_numberOfConstituents;   //!
   TBranch        *b_Jet_AK8_chargedMultiplicity;   //!
   TBranch        *b_Jet_AK8_puppi_pt;   //!
   TBranch        *b_Jet_AK8_puppi_eta;   //!
   TBranch        *b_Jet_AK8_puppi_mass;   //!
   TBranch        *b_Jet_AK8_puppi_phi;   //!
   TBranch        *b_Jet_AK8_puppi_tau1;   //!
   TBranch        *b_Jet_AK8_puppi_tau2;   //!
   TBranch        *b_Jet_AK8_puppi_tau3;   //!
   TBranch        *b_Jet_AK8_GEN_pt;   //!
   TBranch        *b_Jet_AK8_GEN_phi;   //!
   TBranch        *b_Jet_AK8_GEN_eta;   //!
   TBranch        *b_Jet_AK8_GEN_mass;   //!
   TBranch        *b_Jet_AK8_GEN_energy;   //!
   TBranch        *b_Jet_AK8_GEN_parton;   //!
   TBranch        *b_Jet_AK8_GEN_mother;   //!
   TBranch        *b_Jet_AK8_subjet0_pt;   //!
   TBranch        *b_Jet_AK8_subjet0_phi;   //!
   TBranch        *b_Jet_AK8_subjet0_eta;   //!
   TBranch        *b_Jet_AK8_subjet0_mass;   //!
   TBranch        *b_Jet_AK8_subjet0_energy;   //!
   TBranch        *b_Jet_AK8_subjet0_CSVv2;   //!
   TBranch        *b_Jet_AK8_subjet1_pt;   //!
   TBranch        *b_Jet_AK8_subjet1_phi;   //!
   TBranch        *b_Jet_AK8_subjet1_eta;   //!
   TBranch        *b_Jet_AK8_subjet1_mass;   //!
   TBranch        *b_Jet_AK8_subjet1_energy;   //!
   TBranch        *b_Jet_AK8_subjet1_CSVv2;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_electronEnergy;   //!
   TBranch        *b_Jet_photonEnergy;   //!
   TBranch        *b_UncorrJet_pt;   //!
   TBranch        *b_Jet_puppi_mass;   //!
   TBranch        *b_Jet_puppi_electronEnergy;   //!
   TBranch        *b_Jet_puppi_photonEnergy;   //!
   TBranch        *b_UncorrJet_puppi_pt;   //!
   TBranch        *b_Jet_AK8_mass;   //!
   TBranch        *b_Jet_AK8_electronEnergy;   //!
   TBranch        *b_Jet_AK8_photonEnergy;   //!
   TBranch        *b_UncorrJet_AK8_pt;   //!
   TBranch        *b_Gen_pt;   //!
   TBranch        *b_Gen_eta;   //!
   TBranch        *b_Gen_phi;   //!
   TBranch        *b_Gen_status;   //!
   TBranch        *b_Gen_pdg_id;   //!
   TBranch        *b_Gen_motherpdg_id;   //!
   TBranch        *b_Gen_energy;   //!
   TBranch        *b_Gen_vx;   //!
   TBranch        *b_Gen_vy;   //!
   TBranch        *b_Gen_vz;   //!
   TBranch        *b_Gen_charge;   //!
   TBranch        *b_Gen_numDaught;   //!
   TBranch        *b_Gen_numMother;   //!
   TBranch        *b_Gen_BmotherIndices;   //!
   TBranch        *b_Gen_BdaughtIndices;   //!
   TBranch        *b_Gen_BmotherIndex;   //!
   TBranch        *b_weightevt;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_nObservedInTimePUVertices;   //!
   TBranch        *b_nTruePUInteractions;   //!
   TBranch        *b_nObservedOutOfTimePUVertices;   //!
   TBranch        *b_nObservedPlus1BXPUVertices;   //!
   TBranch        *b_nObservedMinus1BXPUVertices;   //!
   TBranch        *b_bestVertices;   //!
   TBranch        *b_pvertex_x;   //!
   TBranch        *b_pvertex_y;   //!
   TBranch        *b_pvertex_z;   //!
   TBranch        *b_pvertex_xError;   //!
   TBranch        *b_pvertex_yError;   //!
   TBranch        *b_pvertex_zError;   //!
   TBranch        *b_beamSpot_x0;   //!
   TBranch        *b_beamSpot_y0;   //!
   TBranch        *b_beamSpot_z0;   //!
   TBranch        *b_beamSpot_xWidth;   //!
   TBranch        *b_beamSpot_yWidth;   //!
   TBranch        *b_Met_type1PF_pt;   //!
   TBranch        *b_Met_type1PF_px;   //!
   TBranch        *b_Met_type1PF_py;   //!
   TBranch        *b_Met_type1PF_pz;   //!
   TBranch        *b_Met_type1PF_sumEt;   //!
   TBranch        *b_Met_type1PF_phi;   //!
   TBranch        *b_Met_type1PF_cov00;   //!
   TBranch        *b_Met_type1PF_cov01;   //!
   TBranch        *b_Met_type1PF_cov10;   //!
   TBranch        *b_Met_type1PF_cov11;   //!
   TBranch        *b_Met_puppi_pt;   //!
   TBranch        *b_Met_puppi_px;   //!
   TBranch        *b_Met_puppi_py;   //!
   TBranch        *b_Met_puppi_pz;   //!
   TBranch        *b_Met_puppi_sumEt;   //!
   TBranch        *b_Met_puppi_phi;   //!
   TBranch        *b_Met_NoHF_pt;   //!
   TBranch        *b_Met_NoHF_px;   //!
   TBranch        *b_Met_NoHF_py;   //!
   TBranch        *b_Met_NoHF_pz;   //!
   TBranch        *b_Met_NoHF_sumEt;   //!
   TBranch        *b_Met_NoHF_phi;   //!
   TBranch        *b_Gen_Met;   //!
   TBranch        *b_Met_type1PF_shiftedPtUp;   //!
   TBranch        *b_Met_type1PF_shiftedPtDown;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_energy;   //!
   TBranch        *b_Photon_et;   //!
   TBranch        *b_Photon_HoverE;   //!
   TBranch        *b_Photon_phoR9;   //!
   TBranch        *b_Photon_SigmaIEtaIEta;   //!
   TBranch        *b_Photon_SigmaIPhiIPhi;   //!
   TBranch        *b_Photon_PFChIso;   //!
   TBranch        *b_Photon_PFPhoIso;   //!
   TBranch        *b_Photon_PFNeuIso;   //!
   TBranch        *b_Photon_EleVeto;   //!
   TBranch        *b_Photon_hasPixelSeed;   //!

   BSM3G_MonoTop_Analyzer(TTree *tree=0);
   virtual ~BSM3G_MonoTop_Analyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string outFileName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BSM3G_MonoTop_Analyzer_cxx
BSM3G_MonoTop_Analyzer::BSM3G_MonoTop_Analyzer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("OutTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("OutTree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("OutTree.root:/TNT");
      dir->GetObject("BOOM",tree);

   }
   Init(tree);
}

BSM3G_MonoTop_Analyzer::~BSM3G_MonoTop_Analyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BSM3G_MonoTop_Analyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BSM3G_MonoTop_Analyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BSM3G_MonoTop_Analyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Trigger_decision = 0;
   Trigger_names = 0;
   Trigger_prescale = 0;
   Muon_pt = 0;
   Muon_eta = 0;
   Muon_phi = 0;
   Muon_energy = 0;
   Muon_charge = 0;
   Muon_tight = 0;
   Muon_loose = 0;
   Muon_medium = 0;
   Muon_soft = 0;
   Muon_isHighPt = 0;
   Muon_pf = 0;
   Muon_isoCharged = 0;
   Muon_isoSum = 0;
   Muon_isoCharParPt = 0;
   Muon_isTrackerMuon = 0;
   Muon_POGisGood = 0;
   Muon_chi2 = 0;
   Muon_validHits = 0;
   Muon_validHitsInner = 0;
   Muon_matchedStat = 0;
   Muon_dxy_pv = 0;
   Muon_dxy_bs = 0;
   Muon_dz_pv = 0;
   Muon_dz_bs = 0;
   Muon_dxyError = 0;
   Muon_dzError = 0;
   Muon_ndof = 0;
   Muon_vtx = 0;
   Muon_vty = 0;
   Muon_vtz = 0;
   Muon_track_pt = 0;
   Muon_track_ptError = 0;
   Muon_isGlobal = 0;
   Muon_TLayers = 0;
   Muon_isoNeutralHadron = 0;
   Muon_isoPhoton = 0;
   Muon_isoPU = 0;
   Muon_combinedIso = 0;
   Muon_trackRe_iso = 0;
   Muon_dB = 0;
   Muon_besttrack_pt = 0;
   Muon_besttrack_ptError = 0;
   Muon_tunePBestTrack_pt = 0;
   Muon_tunePBestTrackType = 0;
   Muon_chi2LocalPosition = 0;
   Muon_trkKink = 0;
   Muon_segmentCompatibility = 0;
   Muon_validFraction = 0;
   Muon_pixelLayersWithMeasurement = 0;
   Muon_qualityhighPurity = 0;
   Muon_track_PCAx_bs = 0;
   Muon_track_PCAy_bs = 0;
   Muon_track_PCAz_bs = 0;
   Muon_track_PCAx_pv = 0;
   Muon_track_PCAy_pv = 0;
   Muon_track_PCAz_pv = 0;
   Muon_trackFitErrorMatrix_00 = 0;
   Muon_trackFitErrorMatrix_01 = 0;
   Muon_trackFitErrorMatrix_02 = 0;
   Muon_trackFitErrorMatrix_11 = 0;
   Muon_trackFitErrorMatrix_12 = 0;
   Muon_trackFitErrorMatrix_22 = 0;
   patElectron_pt = 0;
   patElectron_eta = 0;
   patElectron_phi = 0;
   patElectron_energy = 0;
   patElectron_charge = 0;
   patElectron_isPassVeto = 0;
   patElectron_isPassLoose = 0;
   patElectron_isPassMedium = 0;
   patElectron_isPassTight = 0;
   patElectron_isPassHEEPId = 0;
   patElectron_passMV1wp1Id = 0;
   patElectron_passMV2wp1Id = 0;
   patElectron_isoChargedHadrons = 0;
   patElectron_isoNeutralHadrons = 0;
   patElectron_isoPhotons = 0;
   patElectron_isoPU = 0;
   patElectron_expectedMissingInnerHits = 0;
   patElectron_passConversionVeto = 0;
   patElectron_gsfTrack_dxy_pv = 0;
   patElectron_gsfTrack_dxy_bs = 0;
   patElectron_dxyError = 0;
   patElectron_gsfTrack_dz_pv = 0;
   patElectron_gsfTrack_dz_bs = 0;
   patElectron_gsfTrack_normChi2 = 0;
   patElectron_gsfTrack_ndof = 0;
   patElectron_gsfTrack_vtx = 0;
   patElectron_gsfTrack_vty = 0;
   patElectron_gsfTrack_vtz = 0;
   patElectron_gsfTrack_PCAx_bs = 0;
   patElectron_gsfTrack_PCAy_bs = 0;
   patElectron_gsfTrack_PCAz_bs = 0;
   patElectron_gsfTrack_PCAx_pv = 0;
   patElectron_gsfTrack_PCAy_pv = 0;
   patElectron_gsfTrack_PCAz_pv = 0;
   patElectron_gsfTrackFitErrorMatrix_00 = 0;
   patElectron_gsfTrackFitErrorMatrix_01 = 0;
   patElectron_gsfTrackFitErrorMatrix_02 = 0;
   patElectron_gsfTrackFitErrorMatrix_11 = 0;
   patElectron_gsfTrackFitErrorMatrix_12 = 0;
   patElectron_gsfTrackFitErrorMatrix_22 = 0;
   Tau_eta = 0;
   Tau_phi = 0;
   Tau_pt = 0;
   Tau_energy = 0;
   Tau_charge = 0;
   Tau_decayMode = 0;
   Tau_nProngs = 0;
   Tau_chargedIsoPtSum = 0;
   Tau_neutralIsoPtSum = 0;
   Tau_puCorrPtSum = 0;
   Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits = 0;
   Tau_decayModeFinding = 0;
   Tau_decayModeFindingNewDMs = 0;
   Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byVLooseIsolationMVArun2v1DBoldDMwLT = 0;
   Tau_byLooseIsolationMVArun2v1DBoldDMwLT = 0;
   Tau_byMediumIsolationMVArun2v1DBoldDMwLT = 0;
   Tau_byTightIsolationMVArun2v1DBoldDMwLT = 0;
   Tau_byVTightIsolationMVArun2v1DBoldDMwLT = 0;
   Tau_byVLooseIsolationMVArun2v1DBnewDMwLT = 0;
   Tau_byLooseIsolationMVArun2v1DBnewDMwLT = 0;
   Tau_byMediumIsolationMVArun2v1DBnewDMwLT = 0;
   Tau_byTightIsolationMVArun2v1DBnewDMwLT = 0;
   Tau_byVTightIsolationMVArun2v1DBnewDMwLT = 0;
   Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = 0;
   Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT = 0;
   Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT = 0;
   Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT = 0;
   Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT = 0;
   Tau_againstMuonLoose3 = 0;
   Tau_againstMuonTight3 = 0;
   Tau_againstElectronMVAVLooseMVA6 = 0;
   Tau_againstElectronMVALooseMVA6 = 0;
   Tau_againstElectronMVAMediumMVA6 = 0;
   Tau_againstElectronMVATightMVA6 = 0;
   Tau_leadChargedCandPt = 0;
   Tau_leadChargedCandEta = 0;
   Tau_leadChargedCandPhi = 0;
   Tau_leadChargedCandCharge = 0;
   Tau_leadChargedCandChi2 = 0;
   Tau_leadChargedCandValidHits = 0;
   Tau_leadChargedCandDxy_pv = 0;
   Tau_leadChargedCandDxy_bs = 0;
   Tau_leadChargedCandDz_pv = 0;
   Tau_leadChargedCandDz_bs = 0;
   Tau_leadChargedCandDzError = 0;
   Tau_leadChargedCandDxyError = 0;
   Tau_leadChargedCandNdof = 0;
   Tau_leadChargedCandVtx = 0;
   Tau_leadChargedCandVty = 0;
   Tau_leadChargedCandVtz = 0;
   Tau_leadChargedCandTrack_pt = 0;
   Tau_leadChargedCandTrack_ptError = 0;
   Tau_leadChargedCandTrack_PCAx_bs = 0;
   Tau_leadChargedCandTrack_PCAy_bs = 0;
   Tau_leadChargedCandTrack_PCAz_bs = 0;
   Tau_leadChargedCandTrack_PCAx_pv = 0;
   Tau_leadChargedCandTrack_PCAy_pv = 0;
   Tau_leadChargedCandTrack_PCAz_pv = 0;
   Tau_leadChargedCandTrackFitErrorMatrix_00 = 0;
   Tau_leadChargedCandTrackFitErrorMatrix_01 = 0;
   Tau_leadChargedCandTrackFitErrorMatrix_02 = 0;
   Tau_leadChargedCandTrackFitErrorMatrix_11 = 0;
   Tau_leadChargedCandTrackFitErrorMatrix_12 = 0;
   Tau_leadChargedCandTrackFitErrorMatrix_22 = 0;
   Tau_defaultDxy = 0;
   Tau_defaultDxyError = 0;
   Tau_defaultDxySig = 0;
   Tau_defaultFlightLengthX = 0;
   Tau_defaultFlightLengthY = 0;
   Tau_defaultFlightLengthZ = 0;
   Tau_defaultFlightLengthSig = 0;
   Tau_default_PCAx_pv = 0;
   Tau_default_PCAy_pv = 0;
   Tau_default_PCAz_pv = 0;
   Jet_pt = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_energy = 0;
   Jet_bDiscriminator = 0;
   Jet_bDiscriminator_CISVV2 = 0;
   Jet_bDiscriminator_pfCISVV2 = 0;
   Jet_partonFlavour = 0;
   Jet_neutralHadEnergyFraction = 0;
   Jet_neutralEmEmEnergyFraction = 0;
   Jet_chargedHadronEnergyFraction = 0;
   Jet_chargedEmEnergyFraction = 0;
   Jet_muonEnergyFraction = 0;
   Jet_numberOfConstituents = 0;
   Jet_chargedMultiplicity = 0;
   Jet_puppi_pt = 0;
   Jet_puppi_eta = 0;
   Jet_puppi_phi = 0;
   Jet_puppi_energy = 0;
   Jet_puppi_bDiscriminator = 0;
   Jet_puppi_bDiscriminator_CISVV2 = 0;
   Jet_puppi_bDiscriminator_pfCISVV2 = 0;
   Jet_puppi_partonFlavour = 0;
   Jet_puppi_neutralHadEnergyFraction = 0;
   Jet_puppi_neutralEmEmEnergyFraction = 0;
   Jet_puppi_chargedHadronEnergyFraction = 0;
   Jet_puppi_chargedEmEnergyFraction = 0;
   Jet_puppi_muonEnergyFraction = 0;
   Jet_puppi_numberOfConstituents = 0;
   Jet_puppi_chargedMultiplicity = 0;
   Jet_AK8_pt = 0;
   Jet_AK8_eta = 0;
   Jet_AK8_phi = 0;
   Jet_AK8_energy = 0;
   Jet_AK8_bDiscriminator = 0;
   Jet_AK8_bDiscriminator_CISVV2 = 0;
   Jet_AK8_bDiscriminator_pfCISVV2 = 0;
   Jet_AK8_partonFlavour = 0;
   Jet_AK8_neutralHadEnergyFraction = 0;
   Jet_AK8_neutralEmEmEnergyFraction = 0;
   Jet_AK8_chargedHadronEnergyFraction = 0;
   Jet_AK8_chargedEmEnergyFraction = 0;
   Jet_AK8_muonEnergyFraction = 0;
   Jet_AK8_numberOfConstituents = 0;
   Jet_AK8_chargedMultiplicity = 0;
   Jet_AK8_puppi_pt = 0;
   Jet_AK8_puppi_eta = 0;
   Jet_AK8_puppi_mass = 0;
   Jet_AK8_puppi_phi = 0;
   Jet_AK8_puppi_tau1 = 0;
   Jet_AK8_puppi_tau2 = 0;
   Jet_AK8_puppi_tau3 = 0;
   Jet_AK8_GEN_pt = 0;
   Jet_AK8_GEN_phi = 0;
   Jet_AK8_GEN_eta = 0;
   Jet_AK8_GEN_mass = 0;
   Jet_AK8_GEN_energy = 0;
   Jet_AK8_GEN_parton = 0;
   Jet_AK8_GEN_mother = 0;
   Jet_AK8_subjet0_pt = 0;
   Jet_AK8_subjet0_phi = 0;
   Jet_AK8_subjet0_eta = 0;
   Jet_AK8_subjet0_mass = 0;
   Jet_AK8_subjet0_energy = 0;
   Jet_AK8_subjet0_CSVv2 = 0;
   Jet_AK8_subjet1_pt = 0;
   Jet_AK8_subjet1_phi = 0;
   Jet_AK8_subjet1_eta = 0;
   Jet_AK8_subjet1_mass = 0;
   Jet_AK8_subjet1_energy = 0;
   Jet_AK8_subjet1_CSVv2 = 0;
   Jet_mass = 0;
   Jet_electronEnergy = 0;
   Jet_photonEnergy = 0;
   UncorrJet_pt = 0;
   Jet_puppi_mass = 0;
   Jet_puppi_electronEnergy = 0;
   Jet_puppi_photonEnergy = 0;
   UncorrJet_puppi_pt = 0;
   Jet_AK8_mass = 0;
   Jet_AK8_electronEnergy = 0;
   Jet_AK8_photonEnergy = 0;
   UncorrJet_AK8_pt = 0;
   Gen_pt = 0;
   Gen_eta = 0;
   Gen_phi = 0;
   Gen_status = 0;
   Gen_pdg_id = 0;
   Gen_motherpdg_id = 0;
   Gen_energy = 0;
   Gen_vx = 0;
   Gen_vy = 0;
   Gen_vz = 0;
   Gen_charge = 0;
   Gen_numDaught = 0;
   Gen_numMother = 0;
   Gen_BmotherIndices = 0;
   Gen_BdaughtIndices = 0;
   Gen_BmotherIndex = 0;
   pvertex_x = 0;
   pvertex_y = 0;
   pvertex_z = 0;
   pvertex_xError = 0;
   pvertex_yError = 0;
   pvertex_zError = 0;
   beamSpot_x0 = 0;
   beamSpot_y0 = 0;
   beamSpot_z0 = 0;
   beamSpot_xWidth = 0;
   beamSpot_yWidth = 0;
   Photon_pt = 0;
   Photon_eta = 0;
   Photon_phi = 0;
   Photon_energy = 0;
   Photon_et = 0;
   Photon_HoverE = 0;
   Photon_phoR9 = 0;
   Photon_SigmaIEtaIEta = 0;
   Photon_SigmaIPhiIPhi = 0;
   Photon_PFChIso = 0;
   Photon_PFPhoIso = 0;
   Photon_PFNeuIso = 0;
   Photon_EleVeto = 0;
   Photon_hasPixelSeed = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Trigger_decision", &Trigger_decision, &b_Trigger_decision);
   fChain->SetBranchAddress("Trigger_names", &Trigger_names, &b_Trigger_names);
   fChain->SetBranchAddress("Trigger_prescale", &Trigger_prescale, &b_Trigger_prescale);
   fChain->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_energy", &Muon_energy, &b_Muon_energy);
   fChain->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_tight", &Muon_tight, &b_Muon_tight);
   fChain->SetBranchAddress("Muon_loose", &Muon_loose, &b_Muon_loose);
   fChain->SetBranchAddress("Muon_medium", &Muon_medium, &b_Muon_medium);
   fChain->SetBranchAddress("Muon_soft", &Muon_soft, &b_Muon_soft);
   fChain->SetBranchAddress("Muon_isHighPt", &Muon_isHighPt, &b_Muon_isHighPt);
   fChain->SetBranchAddress("Muon_pf", &Muon_pf, &b_Muon_pf);
   fChain->SetBranchAddress("Muon_isoCharged", &Muon_isoCharged, &b_Muon_isoCharged);
   fChain->SetBranchAddress("Muon_isoSum", &Muon_isoSum, &b_Muon_isoSum);
   fChain->SetBranchAddress("Muon_isoCharParPt", &Muon_isoCharParPt, &b_Muon_isoCharParPt);
   fChain->SetBranchAddress("Muon_isTrackerMuon", &Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
   fChain->SetBranchAddress("Muon_POGisGood", &Muon_POGisGood, &b_Muon_POGisGood);
   fChain->SetBranchAddress("Muon_chi2", &Muon_chi2, &b_Muon_chi2);
   fChain->SetBranchAddress("Muon_validHits", &Muon_validHits, &b_Muon_validHits);
   fChain->SetBranchAddress("Muon_validHitsInner", &Muon_validHitsInner, &b_Muon_validHitsInner);
   fChain->SetBranchAddress("Muon_matchedStat", &Muon_matchedStat, &b_Muon_matchedStat);
   fChain->SetBranchAddress("Muon_dxy_pv", &Muon_dxy_pv, &b_Muon_dxy_pv);
   fChain->SetBranchAddress("Muon_dxy_bs", &Muon_dxy_bs, &b_Muon_dxy_bs);
   fChain->SetBranchAddress("Muon_dz_pv", &Muon_dz_pv, &b_Muon_dz_pv);
   fChain->SetBranchAddress("Muon_dz_bs", &Muon_dz_bs, &b_Muon_dz_bs);
   fChain->SetBranchAddress("Muon_dxyError", &Muon_dxyError, &b_Muon_dxyError);
   fChain->SetBranchAddress("Muon_dzError", &Muon_dzError, &b_Muon_dzError);
   fChain->SetBranchAddress("Muon_ndof", &Muon_ndof, &b_Muon_ndof);
   fChain->SetBranchAddress("Muon_vtx", &Muon_vtx, &b_Muon_vtx);
   fChain->SetBranchAddress("Muon_vty", &Muon_vty, &b_Muon_vty);
   fChain->SetBranchAddress("Muon_vtz", &Muon_vtz, &b_Muon_vtz);
   fChain->SetBranchAddress("Muon_track_pt", &Muon_track_pt, &b_Muon_track_pt);
   fChain->SetBranchAddress("Muon_track_ptError", &Muon_track_ptError, &b_Muon_track_ptError);
   fChain->SetBranchAddress("Muon_isGlobal", &Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_TLayers", &Muon_TLayers, &b_Muon_TLayers);
   fChain->SetBranchAddress("Muon_isoNeutralHadron", &Muon_isoNeutralHadron, &b_Muon_isoNeutralHadron);
   fChain->SetBranchAddress("Muon_isoPhoton", &Muon_isoPhoton, &b_Muon_isoPhoton);
   fChain->SetBranchAddress("Muon_isoPU", &Muon_isoPU, &b_Muon_isoPU);
   fChain->SetBranchAddress("Muon_combinedIso", &Muon_combinedIso, &b_Muon_combinedIso);
   fChain->SetBranchAddress("Muon_trackRe_iso", &Muon_trackRe_iso, &b_Muon_trackRe_iso);
   fChain->SetBranchAddress("Muon_dB", &Muon_dB, &b_Muon_dB);
   fChain->SetBranchAddress("Muon_besttrack_pt", &Muon_besttrack_pt, &b_Muon_besttrack_pt);
   fChain->SetBranchAddress("Muon_besttrack_ptError", &Muon_besttrack_ptError, &b_Muon_besttrack_ptError);
   fChain->SetBranchAddress("Muon_tunePBestTrack_pt", &Muon_tunePBestTrack_pt, &b_Muon_tunePBestTrack_pt);
   fChain->SetBranchAddress("Muon_tunePBestTrackType", &Muon_tunePBestTrackType, &b_Muon_tunePBestTrackType);
   fChain->SetBranchAddress("Muon_chi2LocalPosition", &Muon_chi2LocalPosition, &b_Muon_chi2LocalPosition);
   fChain->SetBranchAddress("Muon_trkKink", &Muon_trkKink, &b_Muon_trkKink);
   fChain->SetBranchAddress("Muon_segmentCompatibility", &Muon_segmentCompatibility, &b_Muon_segmentCompatibility);
   fChain->SetBranchAddress("Muon_validFraction", &Muon_validFraction, &b_Muon_validFraction);
   fChain->SetBranchAddress("Muon_pixelLayersWithMeasurement", &Muon_pixelLayersWithMeasurement, &b_Muon_pixelLayersWithMeasurement);
   fChain->SetBranchAddress("Muon_qualityhighPurity", &Muon_qualityhighPurity, &b_Muon_qualityhighPurity);
   fChain->SetBranchAddress("Muon_track_PCAx_bs", &Muon_track_PCAx_bs, &b_Muon_track_PCAx_bs);
   fChain->SetBranchAddress("Muon_track_PCAy_bs", &Muon_track_PCAy_bs, &b_Muon_track_PCAy_bs);
   fChain->SetBranchAddress("Muon_track_PCAz_bs", &Muon_track_PCAz_bs, &b_Muon_track_PCAz_bs);
   fChain->SetBranchAddress("Muon_track_PCAx_pv", &Muon_track_PCAx_pv, &b_Muon_track_PCAx_pv);
   fChain->SetBranchAddress("Muon_track_PCAy_pv", &Muon_track_PCAy_pv, &b_Muon_track_PCAy_pv);
   fChain->SetBranchAddress("Muon_track_PCAz_pv", &Muon_track_PCAz_pv, &b_Muon_track_PCAz_pv);
   fChain->SetBranchAddress("Muon_trackFitErrorMatrix_00", &Muon_trackFitErrorMatrix_00, &b_Muon_trackFitErrorMatrix_00);
   fChain->SetBranchAddress("Muon_trackFitErrorMatrix_01", &Muon_trackFitErrorMatrix_01, &b_Muon_trackFitErrorMatrix_01);
   fChain->SetBranchAddress("Muon_trackFitErrorMatrix_02", &Muon_trackFitErrorMatrix_02, &b_Muon_trackFitErrorMatrix_02);
   fChain->SetBranchAddress("Muon_trackFitErrorMatrix_11", &Muon_trackFitErrorMatrix_11, &b_Muon_trackFitErrorMatrix_11);
   fChain->SetBranchAddress("Muon_trackFitErrorMatrix_12", &Muon_trackFitErrorMatrix_12, &b_Muon_trackFitErrorMatrix_12);
   fChain->SetBranchAddress("Muon_trackFitErrorMatrix_22", &Muon_trackFitErrorMatrix_22, &b_Muon_trackFitErrorMatrix_22);
   fChain->SetBranchAddress("patElectron_pt", &patElectron_pt, &b_patElectron_pt);
   fChain->SetBranchAddress("patElectron_eta", &patElectron_eta, &b_patElectron_eta);
   fChain->SetBranchAddress("patElectron_phi", &patElectron_phi, &b_patElectron_phi);
   fChain->SetBranchAddress("patElectron_energy", &patElectron_energy, &b_patElectron_energy);
   fChain->SetBranchAddress("patElectron_charge", &patElectron_charge, &b_patElectron_charge);
   fChain->SetBranchAddress("patElectron_isPassVeto", &patElectron_isPassVeto, &b_patElectron_isPassVeto);
   fChain->SetBranchAddress("patElectron_isPassLoose", &patElectron_isPassLoose, &b_patElectron_isPassLoose);
   fChain->SetBranchAddress("patElectron_isPassMedium", &patElectron_isPassMedium, &b_patElectron_isPassMedium);
   fChain->SetBranchAddress("patElectron_isPassTight", &patElectron_isPassTight, &b_patElectron_isPassTight);
   fChain->SetBranchAddress("patElectron_isPassHEEPId", &patElectron_isPassHEEPId, &b_patElectron_isPassHEEPId);
   fChain->SetBranchAddress("patElectron_passMV1wp1Id", &patElectron_passMV1wp1Id, &b_patElectron_passMV1wp1Id);
   fChain->SetBranchAddress("patElectron_passMV2wp1Id", &patElectron_passMV2wp1Id, &b_patElectron_passMV2wp1Id);
   fChain->SetBranchAddress("patElectron_isoChargedHadrons", &patElectron_isoChargedHadrons, &b_patElectron_isoChargedHadrons);
   fChain->SetBranchAddress("patElectron_isoNeutralHadrons", &patElectron_isoNeutralHadrons, &b_patElectron_isoNeutralHadrons);
   fChain->SetBranchAddress("patElectron_isoPhotons", &patElectron_isoPhotons, &b_patElectron_isoPhotons);
   fChain->SetBranchAddress("patElectron_isoPU", &patElectron_isoPU, &b_patElectron_isoPU);
   fChain->SetBranchAddress("patElectron_expectedMissingInnerHits", &patElectron_expectedMissingInnerHits, &b_patElectron_expectedMissingInnerHits);
   fChain->SetBranchAddress("patElectron_passConversionVeto", &patElectron_passConversionVeto, &b_patElectron_passConversionVeto);
   fChain->SetBranchAddress("patElectron_gsfTrack_dxy_pv", &patElectron_gsfTrack_dxy_pv, &b_patElectron_gsfTrack_dxy_pv);
   fChain->SetBranchAddress("patElectron_gsfTrack_dxy_bs", &patElectron_gsfTrack_dxy_bs, &b_patElectron_gsfTrack_dxy_bs);
   fChain->SetBranchAddress("patElectron_dxyError", &patElectron_dxyError, &b_patElectron_dxyError);
   fChain->SetBranchAddress("patElectron_gsfTrack_dz_pv", &patElectron_gsfTrack_dz_pv, &b_patElectron_gsfTrack_dz_pv);
   fChain->SetBranchAddress("patElectron_gsfTrack_dz_bs", &patElectron_gsfTrack_dz_bs, &b_patElectron_gsfTrack_dz_bs);
   fChain->SetBranchAddress("patElectron_gsfTrack_normChi2", &patElectron_gsfTrack_normChi2, &b_patElectron_gsfTrack_normChi2);
   fChain->SetBranchAddress("patElectron_gsfTrack_ndof", &patElectron_gsfTrack_ndof, &b_patElectron_gsfTrack_ndof);
   fChain->SetBranchAddress("patElectron_gsfTrack_vtx", &patElectron_gsfTrack_vtx, &b_patElectron_gsfTrack_vtx);
   fChain->SetBranchAddress("patElectron_gsfTrack_vty", &patElectron_gsfTrack_vty, &b_patElectron_gsfTrack_vty);
   fChain->SetBranchAddress("patElectron_gsfTrack_vtz", &patElectron_gsfTrack_vtz, &b_patElectron_gsfTrack_vtz);
   fChain->SetBranchAddress("patElectron_gsfTrack_PCAx_bs", &patElectron_gsfTrack_PCAx_bs, &b_patElectron_gsfTrack_PCAx_bs);
   fChain->SetBranchAddress("patElectron_gsfTrack_PCAy_bs", &patElectron_gsfTrack_PCAy_bs, &b_patElectron_gsfTrack_PCAy_bs);
   fChain->SetBranchAddress("patElectron_gsfTrack_PCAz_bs", &patElectron_gsfTrack_PCAz_bs, &b_patElectron_gsfTrack_PCAz_bs);
   fChain->SetBranchAddress("patElectron_gsfTrack_PCAx_pv", &patElectron_gsfTrack_PCAx_pv, &b_patElectron_gsfTrack_PCAx_pv);
   fChain->SetBranchAddress("patElectron_gsfTrack_PCAy_pv", &patElectron_gsfTrack_PCAy_pv, &b_patElectron_gsfTrack_PCAy_pv);
   fChain->SetBranchAddress("patElectron_gsfTrack_PCAz_pv", &patElectron_gsfTrack_PCAz_pv, &b_patElectron_gsfTrack_PCAz_pv);
   fChain->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_00", &patElectron_gsfTrackFitErrorMatrix_00, &b_patElectron_gsfTrackFitErrorMatrix_00);
   fChain->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_01", &patElectron_gsfTrackFitErrorMatrix_01, &b_patElectron_gsfTrackFitErrorMatrix_01);
   fChain->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_02", &patElectron_gsfTrackFitErrorMatrix_02, &b_patElectron_gsfTrackFitErrorMatrix_02);
   fChain->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_11", &patElectron_gsfTrackFitErrorMatrix_11, &b_patElectron_gsfTrackFitErrorMatrix_11);
   fChain->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_12", &patElectron_gsfTrackFitErrorMatrix_12, &b_patElectron_gsfTrackFitErrorMatrix_12);
   fChain->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_22", &patElectron_gsfTrackFitErrorMatrix_22, &b_patElectron_gsfTrackFitErrorMatrix_22);
   fChain->SetBranchAddress("Tau_eta", &Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_phi", &Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_pt", &Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_energy", &Tau_energy, &b_Tau_energy);
   fChain->SetBranchAddress("Tau_charge", &Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_decayMode", &Tau_decayMode, &b_Tau_decayMode);
   fChain->SetBranchAddress("Tau_nProngs", &Tau_nProngs, &b_Tau_nProngs);
   fChain->SetBranchAddress("Tau_chargedIsoPtSum", &Tau_chargedIsoPtSum, &b_Tau_chargedIsoPtSum);
   fChain->SetBranchAddress("Tau_neutralIsoPtSum", &Tau_neutralIsoPtSum, &b_Tau_neutralIsoPtSum);
   fChain->SetBranchAddress("Tau_puCorrPtSum", &Tau_puCorrPtSum, &b_Tau_puCorrPtSum);
   fChain->SetBranchAddress("Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", &Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("Tau_decayModeFinding", &Tau_decayModeFinding, &b_Tau_decayModeFinding);
   fChain->SetBranchAddress("Tau_decayModeFindingNewDMs", &Tau_decayModeFindingNewDMs, &b_Tau_decayModeFindingNewDMs);
   fChain->SetBranchAddress("Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("Tau_byTightCombinedIsolationDeltaBetaCorr3Hits", &Tau_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("Tau_byVLooseIsolationMVArun2v1DBoldDMwLT", &Tau_byVLooseIsolationMVArun2v1DBoldDMwLT, &b_Tau_byVLooseIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("Tau_byLooseIsolationMVArun2v1DBoldDMwLT", &Tau_byLooseIsolationMVArun2v1DBoldDMwLT, &b_Tau_byLooseIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("Tau_byMediumIsolationMVArun2v1DBoldDMwLT", &Tau_byMediumIsolationMVArun2v1DBoldDMwLT, &b_Tau_byMediumIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("Tau_byTightIsolationMVArun2v1DBoldDMwLT", &Tau_byTightIsolationMVArun2v1DBoldDMwLT, &b_Tau_byTightIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("Tau_byVTightIsolationMVArun2v1DBoldDMwLT", &Tau_byVTightIsolationMVArun2v1DBoldDMwLT, &b_Tau_byVTightIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("Tau_byVLooseIsolationMVArun2v1DBnewDMwLT", &Tau_byVLooseIsolationMVArun2v1DBnewDMwLT, &b_Tau_byVLooseIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("Tau_byLooseIsolationMVArun2v1DBnewDMwLT", &Tau_byLooseIsolationMVArun2v1DBnewDMwLT, &b_Tau_byLooseIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("Tau_byMediumIsolationMVArun2v1DBnewDMwLT", &Tau_byMediumIsolationMVArun2v1DBnewDMwLT, &b_Tau_byMediumIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("Tau_byTightIsolationMVArun2v1DBnewDMwLT", &Tau_byTightIsolationMVArun2v1DBnewDMwLT, &b_Tau_byTightIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("Tau_byVTightIsolationMVArun2v1DBnewDMwLT", &Tau_byVTightIsolationMVArun2v1DBnewDMwLT, &b_Tau_byVTightIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT", &Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT, &b_Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT", &Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT, &b_Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
  // fChain->SetBranchAddress("Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT", &Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, &b_Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT);
//    fChain->SetBranchAddress("Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT", &Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, &b_Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT", &Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT, &b_Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("Tau_againstMuonLoose3", &Tau_againstMuonLoose3, &b_Tau_againstMuonLoose3);
   fChain->SetBranchAddress("Tau_againstMuonTight3", &Tau_againstMuonTight3, &b_Tau_againstMuonTight3);
   fChain->SetBranchAddress("Tau_againstElectronMVAVLooseMVA6", &Tau_againstElectronMVAVLooseMVA6, &b_Tau_againstElectronMVAVLooseMVA6);
   fChain->SetBranchAddress("Tau_againstElectronMVALooseMVA6", &Tau_againstElectronMVALooseMVA6, &b_Tau_againstElectronMVALooseMVA6);
   fChain->SetBranchAddress("Tau_againstElectronMVAMediumMVA6", &Tau_againstElectronMVAMediumMVA6, &b_Tau_againstElectronMVAMediumMVA6);
   fChain->SetBranchAddress("Tau_againstElectronMVATightMVA6", &Tau_againstElectronMVATightMVA6, &b_Tau_againstElectronMVATightMVA6);
   fChain->SetBranchAddress("Tau_leadChargedCandPt", &Tau_leadChargedCandPt, &b_Tau_leadChargedCandPt);
   fChain->SetBranchAddress("Tau_leadChargedCandEta", &Tau_leadChargedCandEta, &b_Tau_leadChargedCandEta);
   fChain->SetBranchAddress("Tau_leadChargedCandPhi", &Tau_leadChargedCandPhi, &b_Tau_leadChargedCandPhi);
   fChain->SetBranchAddress("Tau_leadChargedCandCharge", &Tau_leadChargedCandCharge, &b_Tau_leadChargedCandCharge);
   fChain->SetBranchAddress("Tau_leadChargedCandChi2", &Tau_leadChargedCandChi2, &b_Tau_leadChargedCandChi2);
   fChain->SetBranchAddress("Tau_leadChargedCandValidHits", &Tau_leadChargedCandValidHits, &b_Tau_leadChargedCandValidHits);
   fChain->SetBranchAddress("Tau_leadChargedCandDxy_pv", &Tau_leadChargedCandDxy_pv, &b_Tau_leadChargedCandDxy_pv);
   fChain->SetBranchAddress("Tau_leadChargedCandDxy_bs", &Tau_leadChargedCandDxy_bs, &b_Tau_leadChargedCandDxy_bs);
   fChain->SetBranchAddress("Tau_leadChargedCandDz_pv", &Tau_leadChargedCandDz_pv, &b_Tau_leadChargedCandDz_pv);
   fChain->SetBranchAddress("Tau_leadChargedCandDz_bs", &Tau_leadChargedCandDz_bs, &b_Tau_leadChargedCandDz_bs);
   fChain->SetBranchAddress("Tau_leadChargedCandDzError", &Tau_leadChargedCandDzError, &b_Tau_leadChargedCandDzError);
   fChain->SetBranchAddress("Tau_leadChargedCandDxyError", &Tau_leadChargedCandDxyError, &b_Tau_leadChargedCandDxyError);
   fChain->SetBranchAddress("Tau_leadChargedCandNdof", &Tau_leadChargedCandNdof, &b_Tau_leadChargedCandNdof);
   fChain->SetBranchAddress("Tau_leadChargedCandVtx", &Tau_leadChargedCandVtx, &b_Tau_leadChargedCandVtx);
   fChain->SetBranchAddress("Tau_leadChargedCandVty", &Tau_leadChargedCandVty, &b_Tau_leadChargedCandVty);
   fChain->SetBranchAddress("Tau_leadChargedCandVtz", &Tau_leadChargedCandVtz, &b_Tau_leadChargedCandVtz);
   fChain->SetBranchAddress("Tau_leadChargedCandTrack_pt", &Tau_leadChargedCandTrack_pt, &b_Tau_leadChargedCandTrack_pt);
   fChain->SetBranchAddress("Tau_leadChargedCandTrack_ptError", &Tau_leadChargedCandTrack_ptError, &b_Tau_leadChargedCandTrack_ptError);
   fChain->SetBranchAddress("Tau_leadChargedCandTrack_PCAx_bs", &Tau_leadChargedCandTrack_PCAx_bs, &b_Tau_leadChargedCandTrack_PCAx_bs);
   fChain->SetBranchAddress("Tau_leadChargedCandTrack_PCAy_bs", &Tau_leadChargedCandTrack_PCAy_bs, &b_Tau_leadChargedCandTrack_PCAy_bs);
   fChain->SetBranchAddress("Tau_leadChargedCandTrack_PCAz_bs", &Tau_leadChargedCandTrack_PCAz_bs, &b_Tau_leadChargedCandTrack_PCAz_bs);
   fChain->SetBranchAddress("Tau_leadChargedCandTrack_PCAx_pv", &Tau_leadChargedCandTrack_PCAx_pv, &b_Tau_leadChargedCandTrack_PCAx_pv);
   fChain->SetBranchAddress("Tau_leadChargedCandTrack_PCAy_pv", &Tau_leadChargedCandTrack_PCAy_pv, &b_Tau_leadChargedCandTrack_PCAy_pv);
   fChain->SetBranchAddress("Tau_leadChargedCandTrack_PCAz_pv", &Tau_leadChargedCandTrack_PCAz_pv, &b_Tau_leadChargedCandTrack_PCAz_pv);
   fChain->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_00", &Tau_leadChargedCandTrackFitErrorMatrix_00, &b_Tau_leadChargedCandTrackFitErrorMatrix_00);
   fChain->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_01", &Tau_leadChargedCandTrackFitErrorMatrix_01, &b_Tau_leadChargedCandTrackFitErrorMatrix_01);
   fChain->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_02", &Tau_leadChargedCandTrackFitErrorMatrix_02, &b_Tau_leadChargedCandTrackFitErrorMatrix_02);
   fChain->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_11", &Tau_leadChargedCandTrackFitErrorMatrix_11, &b_Tau_leadChargedCandTrackFitErrorMatrix_11);
   fChain->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_12", &Tau_leadChargedCandTrackFitErrorMatrix_12, &b_Tau_leadChargedCandTrackFitErrorMatrix_12);
   fChain->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_22", &Tau_leadChargedCandTrackFitErrorMatrix_22, &b_Tau_leadChargedCandTrackFitErrorMatrix_22);
   fChain->SetBranchAddress("Tau_defaultDxy", &Tau_defaultDxy, &b_Tau_defaultDxy);
   fChain->SetBranchAddress("Tau_defaultDxyError", &Tau_defaultDxyError, &b_Tau_defaultDxyError);
   fChain->SetBranchAddress("Tau_defaultDxySig", &Tau_defaultDxySig, &b_Tau_defaultDxySig);
   fChain->SetBranchAddress("Tau_defaultFlightLengthX", &Tau_defaultFlightLengthX, &b_Tau_defaultFlightLengthX);
   fChain->SetBranchAddress("Tau_defaultFlightLengthY", &Tau_defaultFlightLengthY, &b_Tau_defaultFlightLengthY);
   fChain->SetBranchAddress("Tau_defaultFlightLengthZ", &Tau_defaultFlightLengthZ, &b_Tau_defaultFlightLengthZ);
   fChain->SetBranchAddress("Tau_defaultFlightLengthSig", &Tau_defaultFlightLengthSig, &b_Tau_defaultFlightLengthSig);
   fChain->SetBranchAddress("Tau_default_PCAx_pv", &Tau_default_PCAx_pv, &b_Tau_default_PCAx_pv);
   fChain->SetBranchAddress("Tau_default_PCAy_pv", &Tau_default_PCAy_pv, &b_Tau_default_PCAy_pv);
   fChain->SetBranchAddress("Tau_default_PCAz_pv", &Tau_default_PCAz_pv, &b_Tau_default_PCAz_pv);
   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
   fChain->SetBranchAddress("Jet_bDiscriminator", &Jet_bDiscriminator, &b_Jet_bDiscriminator);
   fChain->SetBranchAddress("Jet_bDiscriminator_CISVV2", &Jet_bDiscriminator_CISVV2, &b_Jet_bDiscriminator_CISVV2);
   fChain->SetBranchAddress("Jet_bDiscriminator_pfCISVV2", &Jet_bDiscriminator_pfCISVV2, &b_Jet_bDiscriminator_pfCISVV2);
   fChain->SetBranchAddress("Jet_partonFlavour", &Jet_partonFlavour, &b_Jet_partonFlavour);
   fChain->SetBranchAddress("Jet_neutralHadEnergyFraction", &Jet_neutralHadEnergyFraction, &b_Jet_neutralHadEnergyFraction);
   fChain->SetBranchAddress("Jet_neutralEmEmEnergyFraction", &Jet_neutralEmEmEnergyFraction, &b_Jet_neutralEmEmEnergyFraction);
   fChain->SetBranchAddress("Jet_chargedHadronEnergyFraction", &Jet_chargedHadronEnergyFraction, &b_Jet_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jet_chargedEmEnergyFraction", &Jet_chargedEmEnergyFraction, &b_Jet_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jet_muonEnergyFraction", &Jet_muonEnergyFraction, &b_Jet_muonEnergyFraction);
   fChain->SetBranchAddress("Jet_numberOfConstituents", &Jet_numberOfConstituents, &b_Jet_numberOfConstituents);
   fChain->SetBranchAddress("Jet_chargedMultiplicity", &Jet_chargedMultiplicity, &b_Jet_chargedMultiplicity);
   fChain->SetBranchAddress("Jet_puppi_pt", &Jet_puppi_pt, &b_Jet_puppi_pt);
   fChain->SetBranchAddress("Jet_puppi_eta", &Jet_puppi_eta, &b_Jet_puppi_eta);
   fChain->SetBranchAddress("Jet_puppi_phi", &Jet_puppi_phi, &b_Jet_puppi_phi);
   fChain->SetBranchAddress("Jet_puppi_energy", &Jet_puppi_energy, &b_Jet_puppi_energy);
   fChain->SetBranchAddress("Jet_puppi_bDiscriminator", &Jet_puppi_bDiscriminator, &b_Jet_puppi_bDiscriminator);
   fChain->SetBranchAddress("Jet_puppi_bDiscriminator_CISVV2", &Jet_puppi_bDiscriminator_CISVV2, &b_Jet_puppi_bDiscriminator_CISVV2);
   fChain->SetBranchAddress("Jet_puppi_bDiscriminator_pfCISVV2", &Jet_puppi_bDiscriminator_pfCISVV2, &b_Jet_puppi_bDiscriminator_pfCISVV2);
   fChain->SetBranchAddress("Jet_puppi_partonFlavour", &Jet_puppi_partonFlavour, &b_Jet_puppi_partonFlavour);
   fChain->SetBranchAddress("Jet_puppi_neutralHadEnergyFraction", &Jet_puppi_neutralHadEnergyFraction, &b_Jet_puppi_neutralHadEnergyFraction);
   fChain->SetBranchAddress("Jet_puppi_neutralEmEmEnergyFraction", &Jet_puppi_neutralEmEmEnergyFraction, &b_Jet_puppi_neutralEmEmEnergyFraction);
   fChain->SetBranchAddress("Jet_puppi_chargedHadronEnergyFraction", &Jet_puppi_chargedHadronEnergyFraction, &b_Jet_puppi_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jet_puppi_chargedEmEnergyFraction", &Jet_puppi_chargedEmEnergyFraction, &b_Jet_puppi_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jet_puppi_muonEnergyFraction", &Jet_puppi_muonEnergyFraction, &b_Jet_puppi_muonEnergyFraction);
   fChain->SetBranchAddress("Jet_puppi_numberOfConstituents", &Jet_puppi_numberOfConstituents, &b_Jet_puppi_numberOfConstituents);
   fChain->SetBranchAddress("Jet_puppi_chargedMultiplicity", &Jet_puppi_chargedMultiplicity, &b_Jet_puppi_chargedMultiplicity);
   fChain->SetBranchAddress("Jet_AK8_pt", &Jet_AK8_pt, &b_Jet_AK8_pt);
   fChain->SetBranchAddress("Jet_AK8_eta", &Jet_AK8_eta, &b_Jet_AK8_eta);
   fChain->SetBranchAddress("Jet_AK8_phi", &Jet_AK8_phi, &b_Jet_AK8_phi);
   fChain->SetBranchAddress("Jet_AK8_energy", &Jet_AK8_energy, &b_Jet_AK8_energy);
   fChain->SetBranchAddress("Jet_AK8_bDiscriminator", &Jet_AK8_bDiscriminator, &b_Jet_AK8_bDiscriminator);
   fChain->SetBranchAddress("Jet_AK8_bDiscriminator_CISVV2", &Jet_AK8_bDiscriminator_CISVV2, &b_Jet_AK8_bDiscriminator_CISVV2);
   fChain->SetBranchAddress("Jet_AK8_bDiscriminator_pfCISVV2", &Jet_AK8_bDiscriminator_pfCISVV2, &b_Jet_AK8_bDiscriminator_pfCISVV2);
   fChain->SetBranchAddress("Jet_AK8_partonFlavour", &Jet_AK8_partonFlavour, &b_Jet_AK8_partonFlavour);
   fChain->SetBranchAddress("Jet_AK8_neutralHadEnergyFraction", &Jet_AK8_neutralHadEnergyFraction, &b_Jet_AK8_neutralHadEnergyFraction);
   fChain->SetBranchAddress("Jet_AK8_neutralEmEmEnergyFraction", &Jet_AK8_neutralEmEmEnergyFraction, &b_Jet_AK8_neutralEmEmEnergyFraction);
   fChain->SetBranchAddress("Jet_AK8_chargedHadronEnergyFraction", &Jet_AK8_chargedHadronEnergyFraction, &b_Jet_AK8_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jet_AK8_chargedEmEnergyFraction", &Jet_AK8_chargedEmEnergyFraction, &b_Jet_AK8_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jet_AK8_muonEnergyFraction", &Jet_AK8_muonEnergyFraction, &b_Jet_AK8_muonEnergyFraction);
   fChain->SetBranchAddress("Jet_AK8_numberOfConstituents", &Jet_AK8_numberOfConstituents, &b_Jet_AK8_numberOfConstituents);
   fChain->SetBranchAddress("Jet_AK8_chargedMultiplicity", &Jet_AK8_chargedMultiplicity, &b_Jet_AK8_chargedMultiplicity);
   fChain->SetBranchAddress("Jet_AK8_puppi_pt", &Jet_AK8_puppi_pt, &b_Jet_AK8_puppi_pt);
   fChain->SetBranchAddress("Jet_AK8_puppi_eta", &Jet_AK8_puppi_eta, &b_Jet_AK8_puppi_eta);
   fChain->SetBranchAddress("Jet_AK8_puppi_mass", &Jet_AK8_puppi_mass, &b_Jet_AK8_puppi_mass);
   fChain->SetBranchAddress("Jet_AK8_puppi_phi", &Jet_AK8_puppi_phi, &b_Jet_AK8_puppi_phi);
   fChain->SetBranchAddress("Jet_AK8_puppi_tau1", &Jet_AK8_puppi_tau1, &b_Jet_AK8_puppi_tau1);
   fChain->SetBranchAddress("Jet_AK8_puppi_tau2", &Jet_AK8_puppi_tau2, &b_Jet_AK8_puppi_tau2);
   fChain->SetBranchAddress("Jet_AK8_puppi_tau3", &Jet_AK8_puppi_tau3, &b_Jet_AK8_puppi_tau3);
   fChain->SetBranchAddress("Jet_AK8_GEN_pt", &Jet_AK8_GEN_pt, &b_Jet_AK8_GEN_pt);
   fChain->SetBranchAddress("Jet_AK8_GEN_phi", &Jet_AK8_GEN_phi, &b_Jet_AK8_GEN_phi);
   fChain->SetBranchAddress("Jet_AK8_GEN_eta", &Jet_AK8_GEN_eta, &b_Jet_AK8_GEN_eta);
   fChain->SetBranchAddress("Jet_AK8_GEN_mass", &Jet_AK8_GEN_mass, &b_Jet_AK8_GEN_mass);
   fChain->SetBranchAddress("Jet_AK8_GEN_energy", &Jet_AK8_GEN_energy, &b_Jet_AK8_GEN_energy);
   fChain->SetBranchAddress("Jet_AK8_GEN_parton", &Jet_AK8_GEN_parton, &b_Jet_AK8_GEN_parton);
   fChain->SetBranchAddress("Jet_AK8_GEN_mother", &Jet_AK8_GEN_mother, &b_Jet_AK8_GEN_mother);
   fChain->SetBranchAddress("Jet_AK8_subjet0_pt", &Jet_AK8_subjet0_pt, &b_Jet_AK8_subjet0_pt);
   fChain->SetBranchAddress("Jet_AK8_subjet0_phi", &Jet_AK8_subjet0_phi, &b_Jet_AK8_subjet0_phi);
   fChain->SetBranchAddress("Jet_AK8_subjet0_eta", &Jet_AK8_subjet0_eta, &b_Jet_AK8_subjet0_eta);
   fChain->SetBranchAddress("Jet_AK8_subjet0_mass", &Jet_AK8_subjet0_mass, &b_Jet_AK8_subjet0_mass);
   fChain->SetBranchAddress("Jet_AK8_subjet0_energy", &Jet_AK8_subjet0_energy, &b_Jet_AK8_subjet0_energy);
   fChain->SetBranchAddress("Jet_AK8_subjet0_CSVv2", &Jet_AK8_subjet0_CSVv2, &b_Jet_AK8_subjet0_CSVv2);
   fChain->SetBranchAddress("Jet_AK8_subjet1_pt", &Jet_AK8_subjet1_pt, &b_Jet_AK8_subjet1_pt);
   fChain->SetBranchAddress("Jet_AK8_subjet1_phi", &Jet_AK8_subjet1_phi, &b_Jet_AK8_subjet1_phi);
   fChain->SetBranchAddress("Jet_AK8_subjet1_eta", &Jet_AK8_subjet1_eta, &b_Jet_AK8_subjet1_eta);
   fChain->SetBranchAddress("Jet_AK8_subjet1_mass", &Jet_AK8_subjet1_mass, &b_Jet_AK8_subjet1_mass);
   fChain->SetBranchAddress("Jet_AK8_subjet1_energy", &Jet_AK8_subjet1_energy, &b_Jet_AK8_subjet1_energy);
   fChain->SetBranchAddress("Jet_AK8_subjet1_CSVv2", &Jet_AK8_subjet1_CSVv2, &b_Jet_AK8_subjet1_CSVv2);
   fChain->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_electronEnergy", &Jet_electronEnergy, &b_Jet_electronEnergy);
   fChain->SetBranchAddress("Jet_photonEnergy", &Jet_photonEnergy, &b_Jet_photonEnergy);
   fChain->SetBranchAddress("UncorrJet_pt", &UncorrJet_pt, &b_UncorrJet_pt);
   fChain->SetBranchAddress("Jet_puppi_mass", &Jet_puppi_mass, &b_Jet_puppi_mass);
   fChain->SetBranchAddress("Jet_puppi_electronEnergy", &Jet_puppi_electronEnergy, &b_Jet_puppi_electronEnergy);
   fChain->SetBranchAddress("Jet_puppi_photonEnergy", &Jet_puppi_photonEnergy, &b_Jet_puppi_photonEnergy);
   fChain->SetBranchAddress("UncorrJet_puppi_pt", &UncorrJet_puppi_pt, &b_UncorrJet_puppi_pt);
   fChain->SetBranchAddress("Jet_AK8_mass", &Jet_AK8_mass, &b_Jet_AK8_mass);
   fChain->SetBranchAddress("Jet_AK8_electronEnergy", &Jet_AK8_electronEnergy, &b_Jet_AK8_electronEnergy);
   fChain->SetBranchAddress("Jet_AK8_photonEnergy", &Jet_AK8_photonEnergy, &b_Jet_AK8_photonEnergy);
   fChain->SetBranchAddress("UncorrJet_AK8_pt", &UncorrJet_AK8_pt, &b_UncorrJet_AK8_pt);
   fChain->SetBranchAddress("Gen_pt", &Gen_pt, &b_Gen_pt);
   fChain->SetBranchAddress("Gen_eta", &Gen_eta, &b_Gen_eta);
   fChain->SetBranchAddress("Gen_phi", &Gen_phi, &b_Gen_phi);
   fChain->SetBranchAddress("Gen_status", &Gen_status, &b_Gen_status);
   fChain->SetBranchAddress("Gen_pdg_id", &Gen_pdg_id, &b_Gen_pdg_id);
   fChain->SetBranchAddress("Gen_motherpdg_id", &Gen_motherpdg_id, &b_Gen_motherpdg_id);
   fChain->SetBranchAddress("Gen_energy", &Gen_energy, &b_Gen_energy);
   fChain->SetBranchAddress("Gen_vx", &Gen_vx, &b_Gen_vx);
   fChain->SetBranchAddress("Gen_vy", &Gen_vy, &b_Gen_vy);
   fChain->SetBranchAddress("Gen_vz", &Gen_vz, &b_Gen_vz);
   fChain->SetBranchAddress("Gen_charge", &Gen_charge, &b_Gen_charge);
   fChain->SetBranchAddress("Gen_numDaught", &Gen_numDaught, &b_Gen_numDaught);
   fChain->SetBranchAddress("Gen_numMother", &Gen_numMother, &b_Gen_numMother);
   fChain->SetBranchAddress("Gen_BmotherIndices", &Gen_BmotherIndices, &b_Gen_BmotherIndices);
   fChain->SetBranchAddress("Gen_BdaughtIndices", &Gen_BdaughtIndices, &b_Gen_BdaughtIndices);
   fChain->SetBranchAddress("Gen_BmotherIndex", &Gen_BmotherIndex, &b_Gen_BmotherIndex);
   fChain->SetBranchAddress("weightevt", &weightevt, &b_weightevt);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("nObservedInTimePUVertices", &nObservedInTimePUVertices, &b_nObservedInTimePUVertices);
   fChain->SetBranchAddress("nTruePUInteractions", &nTruePUInteractions, &b_nTruePUInteractions);
   fChain->SetBranchAddress("nObservedOutOfTimePUVertices", &nObservedOutOfTimePUVertices, &b_nObservedOutOfTimePUVertices);
   fChain->SetBranchAddress("nObservedPlus1BXPUVertices", &nObservedPlus1BXPUVertices, &b_nObservedPlus1BXPUVertices);
   fChain->SetBranchAddress("nObservedMinus1BXPUVertices", &nObservedMinus1BXPUVertices, &b_nObservedMinus1BXPUVertices);
   fChain->SetBranchAddress("bestVertices", &bestVertices, &b_bestVertices);
   fChain->SetBranchAddress("pvertex_x", &pvertex_x, &b_pvertex_x);
   fChain->SetBranchAddress("pvertex_y", &pvertex_y, &b_pvertex_y);
   fChain->SetBranchAddress("pvertex_z", &pvertex_z, &b_pvertex_z);
   fChain->SetBranchAddress("pvertex_xError", &pvertex_xError, &b_pvertex_xError);
   fChain->SetBranchAddress("pvertex_yError", &pvertex_yError, &b_pvertex_yError);
   fChain->SetBranchAddress("pvertex_zError", &pvertex_zError, &b_pvertex_zError);
   fChain->SetBranchAddress("beamSpot_x0", &beamSpot_x0, &b_beamSpot_x0);
   fChain->SetBranchAddress("beamSpot_y0", &beamSpot_y0, &b_beamSpot_y0);
   fChain->SetBranchAddress("beamSpot_z0", &beamSpot_z0, &b_beamSpot_z0);
   fChain->SetBranchAddress("beamSpot_xWidth", &beamSpot_xWidth, &b_beamSpot_xWidth);
   fChain->SetBranchAddress("beamSpot_yWidth", &beamSpot_yWidth, &b_beamSpot_yWidth);
   fChain->SetBranchAddress("Met_type1PF_pt", &Met_type1PF_pt, &b_Met_type1PF_pt);
   fChain->SetBranchAddress("Met_type1PF_px", &Met_type1PF_px, &b_Met_type1PF_px);
   fChain->SetBranchAddress("Met_type1PF_py", &Met_type1PF_py, &b_Met_type1PF_py);
   fChain->SetBranchAddress("Met_type1PF_pz", &Met_type1PF_pz, &b_Met_type1PF_pz);
   fChain->SetBranchAddress("Met_type1PF_sumEt", &Met_type1PF_sumEt, &b_Met_type1PF_sumEt);
   fChain->SetBranchAddress("Met_type1PF_phi", &Met_type1PF_phi, &b_Met_type1PF_phi);
   fChain->SetBranchAddress("Met_type1PF_cov00", &Met_type1PF_cov00, &b_Met_type1PF_cov00);
   fChain->SetBranchAddress("Met_type1PF_cov01", &Met_type1PF_cov01, &b_Met_type1PF_cov01);
   fChain->SetBranchAddress("Met_type1PF_cov10", &Met_type1PF_cov10, &b_Met_type1PF_cov10);
   fChain->SetBranchAddress("Met_type1PF_cov11", &Met_type1PF_cov11, &b_Met_type1PF_cov11);
   fChain->SetBranchAddress("Met_puppi_pt", &Met_puppi_pt, &b_Met_puppi_pt);
   fChain->SetBranchAddress("Met_puppi_px", &Met_puppi_px, &b_Met_puppi_px);
   fChain->SetBranchAddress("Met_puppi_py", &Met_puppi_py, &b_Met_puppi_py);
   fChain->SetBranchAddress("Met_puppi_pz", &Met_puppi_pz, &b_Met_puppi_pz);
   fChain->SetBranchAddress("Met_puppi_sumEt", &Met_puppi_sumEt, &b_Met_puppi_sumEt);
   fChain->SetBranchAddress("Met_puppi_phi", &Met_puppi_phi, &b_Met_puppi_phi);
   fChain->SetBranchAddress("Met_NoHF_pt", &Met_NoHF_pt, &b_Met_NoHF_pt);
   fChain->SetBranchAddress("Met_NoHF_px", &Met_NoHF_px, &b_Met_NoHF_px);
   fChain->SetBranchAddress("Met_NoHF_py", &Met_NoHF_py, &b_Met_NoHF_py);
   fChain->SetBranchAddress("Met_NoHF_pz", &Met_NoHF_pz, &b_Met_NoHF_pz);
   fChain->SetBranchAddress("Met_NoHF_sumEt", &Met_NoHF_sumEt, &b_Met_NoHF_sumEt);
   fChain->SetBranchAddress("Met_NoHF_phi", &Met_NoHF_phi, &b_Met_NoHF_phi);
   fChain->SetBranchAddress("Gen_Met", &Gen_Met, &b_Gen_Met);
   fChain->SetBranchAddress("Met_type1PF_shiftedPtUp", &Met_type1PF_shiftedPtUp, &b_Met_type1PF_shiftedPtUp);
   fChain->SetBranchAddress("Met_type1PF_shiftedPtDown", &Met_type1PF_shiftedPtDown, &b_Met_type1PF_shiftedPtDown);
   fChain->SetBranchAddress("Photon_pt", &Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_eta", &Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_phi", &Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_energy", &Photon_energy, &b_Photon_energy);
   fChain->SetBranchAddress("Photon_et", &Photon_et, &b_Photon_et);
   fChain->SetBranchAddress("Photon_HoverE", &Photon_HoverE, &b_Photon_HoverE);
   fChain->SetBranchAddress("Photon_phoR9", &Photon_phoR9, &b_Photon_phoR9);
   fChain->SetBranchAddress("Photon_SigmaIEtaIEta", &Photon_SigmaIEtaIEta, &b_Photon_SigmaIEtaIEta);
   fChain->SetBranchAddress("Photon_SigmaIPhiIPhi", &Photon_SigmaIPhiIPhi, &b_Photon_SigmaIPhiIPhi);
   fChain->SetBranchAddress("Photon_PFChIso", &Photon_PFChIso, &b_Photon_PFChIso);
   fChain->SetBranchAddress("Photon_PFPhoIso", &Photon_PFPhoIso, &b_Photon_PFPhoIso);
   fChain->SetBranchAddress("Photon_PFNeuIso", &Photon_PFNeuIso, &b_Photon_PFNeuIso);
   fChain->SetBranchAddress("Photon_EleVeto", &Photon_EleVeto, &b_Photon_EleVeto);
   fChain->SetBranchAddress("Photon_hasPixelSeed", &Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
   Notify();
}

Bool_t BSM3G_MonoTop_Analyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BSM3G_MonoTop_Analyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BSM3G_MonoTop_Analyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BSM3G_MonoTop_Analyzer_cxx
