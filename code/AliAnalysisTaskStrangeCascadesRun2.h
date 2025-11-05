/**
 * @file AliAnalysisTaskStrangeCascadesRun2.h
 * @author Ishaan Ahuja (ishaanahuja0@gmail.com)
 * @date 14-02-2025
 * @brief Analyse Xi and Omega cascade candidates
 *
 * Thesis: Multi-strange particle production in p–Pb collisions at √sNN = 8.16 TeV
 * DOI: https://doi.org/10.17181/cwcde-g1z94
 *
 * @version 1.0
 * 
 * @copyright Copyright (c) 2025
 */

#ifndef AliAnalysisTaskStrangeCascadesRun2_H
#define AliAnalysisTaskStrangeCascadesRun2_H

#include "TString.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"

/**
 * @brief Analysis task for strange cascade (Xi and Omega) reconstruction in Run 2
 * 
 * This class implements an analysis task for the reconstruction and study of strange
 * cascade particles (Xi and Omega baryons) in ALICE Run 2 data. It inherits from
 * AliAnalysisTaskSE and provides functionality for:
 * 
 * - Topological and kinematic selections for cascade candidates
 * - PID selection for daughter tracks
 * - Event selection and pile-up rejection
 * - Support for both real data and Monte Carlo analysis
 * - Configurable binning for centrality, mass and pt spectra
 * - Cut variation studies with multiple selection criteria
 * 
 * The task can analyze both Xi (Ξ⁻ and Ξ⁺) and Omega (Ω⁻ and Ω⁺) cascades,
 * with configurable selection criteria for each particle species.
 * 
 * @author Ishaan Ahuja
 * @date 14-02-2025
 * 
 * @note This task is designed for ALICE Run 2 data analysis
 * 
 * Main features:
 * - Flexible cut configuration system with default and variable cuts
 * - Support for MC association studies
 * - Parametric cuts for bachelor-baryon cosine pointing angle
 * - Geometric track length cuts
 * - Multiple analysis ranges (low, mid, high pt)
 * 
 * @see AliAnalysisTaskSE
 * @see AliEventCuts
 * @see AliPIDResponse
 */
class AliAnalysisTaskStrangeCascadesRun2 : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskStrangeCascadesRun2();
    AliAnalysisTaskStrangeCascadesRun2(const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskStrangeCascadesRun2();

    enum particles
    {
        kXi,
        kOm,
        knumpart
    };
    enum signedparticles
    {
        kxip,
        kxim,
        komp,
        komm,
        ksignednumpart
    };
    enum ptInterval
    {
        kLow,
        kMid,
        kHigh,
        kNumPtInterval
    };
    enum cascEvCuts
    {
        kRapidityIntervalMin,   // kCasc_y
        kRapidityIntervalMax,   // == 0
        kTpcDedxPidSigma,       // kCasc_NSigPID
        kDeviationPropLifetime, // kCasc_PropLifetXi
        kLeastTpcClusters,      // kCasc_LeastTPCcls,    // --> DON'T USE
        kCompetingCascRejectOm, // kCasc_CompetingXiMass
        kLeastCRows,            // not used --> USE -> 70 (Emily), 80 (Marek's suggestion - from Michal)
        kLeastCRowsOvF,         // not used --> USE -> (0.8)
        kTrackLengthCut,        // not used --> USE
        kEtaDaughter,           // kCasc_etaDaugh
        kBacBarCosPa,           // kCasc_BacBarCosPA :: not used --> USE
        kV0InvMassWindow,       // kCasc_InvMassLam,     // set to 0.008 instead of 0.005 (lambda == v0)
        kNumCascEvCuts
    };

    enum cascTopoCuts
    {
        kCascTransDecayRadius, // kCasc_CascRad(XiMid)
        kV0TransDecayRadius,   // kCasc_V0RadXi(Mid)
        kDcaBachToPv,          // kCasc_DcaBachToPV
        kDcaV0ToPv,            // kCasc_DcaV0ToPV
        kDcaMesV0ToPv,         // kCasc_DcaMesToPV
        kDcaBarV0ToPv,         // kCasc_DcaBarToPV
        kDcaV0Daughters,       // kCasc_DcaV0Daught(Mid), // DCA V0 daughters (sigma)
        kDcaBachToV0,          // kCasc_DcaCascDaught(Mid)
        kCascCosPa,            // kCasc_CascCosPA(Low)
        kV0CosPa,              // kCasc_V0CosPAXi(Mid)
        kNumCascTopoCuts
    };
    enum cutVars
    {
        kVeryLoose,
        kLoose,
        kTight,
        kVeryTight,
        kNumCutVars
    };

    TString cutNamesEv[kNumCascEvCuts] = {"RapidityIntervalMin",
                                          "RapidityIntervalMax",
                                          "TpcDedxPidSigma",
                                          "DeviationPropLifetime",
                                          "LeastTpcClusters",
                                          "CompetingCascRejectOm",
                                          "LeastCRows",
                                          "LeastCRowsOvF",
                                          "TrackLengthCut",
                                          "EtaDaughter",
                                          "BacBarCosPa",
                                          "V0InvMassWindow"};

    TString cutNamesTopo[kNumCascTopoCuts] = {"CascTransDecayRadius",
                                              "V0TransDecayRadius",
                                              "DcaBachToPv",
                                              "DcaV0ToPv",
                                              "DcaMesV0ToPv",
                                              "DcaBarV0ToPv",
                                              "DcaV0Daughters",
                                              "DcaBachToV0",
                                              "CascCosPa",
                                              "V0CosPa"};

    TString cutVarNames[kNumCutVars] = {"VeryLoose",
                                        "Loose",
                                        "Tight",
                                        "VeryTight"};

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    // cut values setter
    void SetDefOnly(bool);
    void SetEvCutVal(bool, bool, int, double);
    void SetDefCutValue(bool, int, double, int, int);
    void SetVarCutValue(bool, int, int, double, int, int);

    void SetParametricBacBarCosPA(int, float *, float *, int);

    // TrackLength Cut setters
    void SetDeadZoneWidthGeoCut(float DeadZoneWidth) { fDeadZoneWidth_GeoCut = DeadZoneWidth; };
    void SetNcrNclLengthGeoCut(float NcrNclLength) { fNcrNclLength_GeoCut = NcrNclLength; };
    void SetTPCsignalNCut(int TPCsignalNCut) { fTPCsignalNCut = TPCsignalNCut; };
    Float_t GetLengthInActiveZone(AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b);

    // binning setters
    void SetCentbinning(int, int, double *);
    void SetMassbinning(int, int, double, double);
    void SetPtbinning(int, int, double *);

    // set and get which particle specie are analysed
    void SetParticleAnalysisStatus(bool, bool);
    bool GetParticleAnalysisStatus(int);

    // MC-related setters and getters
    void SetIsMC(bool IsMC) { fisMC = IsMC; };
    void SetIsMCassoc(bool IsMCassoc) { fisMCassoc = IsMCassoc; };

    // pile-up rejection setter
    void SetRejectPileUpEvts(bool RejectPileupEvts, int PileupCut = 1)
    {
        if (RejectPileupEvts == kTRUE)
            fPileupCut = PileupCut;
    };

private:
    THistManager *fHistos_eve;   //!
    THistManager *fHistos_XiMin; //!
    THistManager *fHistos_XiPlu; //!
    THistManager *fHistos_OmMin; //!
    THistManager *fHistos_OmPlu; //!

    // objects retreived from input handler
    AliPIDResponse *fPIDResponse; //!
    UInt_t fTriggerMask;          //!

    // AliEventCuts object
    AliEventCuts fEventCuts; //

    // pile-up rejection flag
    int fPileupCut; //

    // MC-realted variables
    bool fisMC;      //
    bool fisMCassoc; //

    // Default cut configuration
    bool fDefOnly; //
    // double fCasc_Cuts[kCasccutsnum]; //
    double fCascTopoCuts[kNumCascTopoCuts]; //

    // particles to be analysed
    bool fParticleAnalysisStatus[ksignednumpart]; //

    // geometrical cut usage
    AliESDtrackCuts fESDTrackCuts; //

    // variables for Cascade analysis
    double fCasc_DcaCascDaught; //!
    double fCasc_CascCosPA;     //!
    double fCasc_CascRad;       //!
    double fCasc_etaPos;        //!
    double fCasc_etaNeg;        //!
    double fCasc_etaBac;        //!
    double fCasc_kinkidx;       //!
    double fCasc_NSigPosProton; //!
    double fCasc_NSigPosPion;   //!
    double fCasc_NSigNegProton; //!
    double fCasc_NSigNegPion;   //!
    double fCasc_NSigBacPion;   //!
    double fCasc_NSigBacKaon;   //!
    double fCasc_LeastCRows;    //!
    double fCasc_LeastCRowsOvF; //!
    double fCasc_LeastTPCcls;   //!
    int fCasc_TrackLengthCut;   //!
    double fCasc_MaxChi2perCls; //!

    double fCasc_InvMassLam;        //!
    double fCasc_DcaV0Daught;       //!
    double fCasc_V0CosPA;           //!
    double fCasc_DcaV0ToPV;         //!
    double fCasc_DcaBachToPV;       //!
    double fCasc_ITSTOFtracks;      //!
    double fCasc_yXi;               //!
    double fCasc_yOm;               //!
    int fCasc_charge;               //!
    double fCasc_Pt;                //!
    double fCasc_DistOverTotP;      //!
    double fCasc_InvMassXiMin;      //!
    double fCasc_InvMassXiPlu;      //!
    double fCasc_InvMassOmMin;      //!
    double fCasc_InvMassOmPlu;      //!
    double fCasc_V0Rad;             //!
    double fCasc_DcaPosToPV;        //!
    double fCasc_DcaNegToPV;        //!
    ULong64_t fCasc_NegTrackStatus; //!
    ULong64_t fCasc_PosTrackStatus; //!
    ULong64_t fCasc_BacTrackStatus; //!
    double fCasc_BacBarCosPA;       //!

    bool fisParametricBacBarCosPA; //
    TH1F *fHist_PtBacBarCosPA;     //
    int fCentLimit_BacBarCosPA;    //

    bool fisParametricTrackLengthCut; //
    TH1I *fHist_CentTrackLengthCut;   //
    float fDeadZoneWidth_GeoCut;      //
    float fNcrNclLength_GeoCut;       //
    int fTPCsignalNCut;               //
    // double fCasc_TrackLength;

    // cut values to be set
    // double cutval_Casc[kCasccutsnum];     //
    double cutValTopo[knumpart][kNumCascTopoCuts][kNumPtInterval];
    double cutValEv[kNumCascEvCuts];
    double def_cutValTopo[knumpart][kNumCascTopoCuts][kNumPtInterval];
    double def_cutValEv[kNumCascEvCuts];
    double var_cutValTopo[knumpart][kNumCascTopoCuts][kNumPtInterval][kNumCutVars];
    double var_cutValEv[kNumCascEvCuts][kNumCutVars];

    int nvarcut_Ev[kNumCascEvCuts];                           //
    double varlowcut_Ev[kNumCascEvCuts];                      //
    double varhighcut_Ev[kNumCascEvCuts];                     //
    int nvarcut_Topo[kNumCascTopoCuts];                           //
    double varlowcut_Topo[kNumCascTopoCuts];                      //
    double varhighcut_Topo[kNumCascTopoCuts];                     //


    // variables to handle binning
    int fncentbins[knumpart];            //
    double fcentbinning[knumpart][50];   //
    int fnmassbins[knumpart];            //
    double fmassbinning[knumpart][1000]; //
    int fnptbins[knumpart];              //
    double fptbinning[knumpart][600];    //

    /// variables to set boundaries for low, mid , high pt ranges
    double ptXiBoundary_LowMid;
    double ptXiBoundary_MidHigh;
    double ptOmBoundary_LowMid;
    double ptOmBoundary_MidHigh;

    // functions to allow flushing part of code out of UserExec
    bool ApplyCuts(int);
    void DataPosting();
    void FillHistCutVariations(double, bool, bool *);
    // functions to allow the correct streaming of the cut variation
    void SetDefCutVals();
    void SetDefCuts();
    void SetCutVariation(bool, int, int, double, double);
    void RandomiseCuts(Int_t);
    // void SetCutVariation(int, int, double, double, double, double);
    void SetDefCutVariations();
    Double_t SetCutValue(bool, int, double, int, int);

    AliAnalysisTaskStrangeCascadesRun2(const AliAnalysisTaskStrangeCascadesRun2 &);            // not implemented
    AliAnalysisTaskStrangeCascadesRun2 &operator=(const AliAnalysisTaskStrangeCascadesRun2 &); // not implemented

    ClassDef(AliAnalysisTaskStrangeCascadesRun2, 1);
};

#endif
