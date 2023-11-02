#ifndef AliCascadeAnalysisMC_Ishaan_H
#define AliCascadeAnalysisMC_Ishaan_H

#include "TString.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"

#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"

class AliCascadeAnalysisMC_Ishaan : public AliAnalysisTaskSE
{
public:
    AliCascadeAnalysisMC_Ishaan();
    AliCascadeAnalysisMC_Ishaan(const char *name, TString lExtraOptions = "");
    virtual ~AliCascadeAnalysisMC_Ishaan();

    // enum and names.
    enum cutnumb_Casc
    {
        kCasc_DcaCascDaughtMid,
        kCasc_CascCosPALow,
        kCasc_CascRadXi,
        kCasc_NSigPID,
        kCasc_LeastCRows,     // not used --> USE -> 70 (emily), 80 (marek's suggestion - from Michal)
        kCasc_LeastCRowsOvF,  // not used --> USE -> (0.8)
        kCasc_LeastTPCcls,    // --> DON'T USE
        kCasc_TrackLengthCut, // not used --> USE
        kCasc_MaxChi2perCls,  // not used
        kCasc_InvMassLam,     // set to 0.008 instead of 0.005 (lambda == v0)
        kCasc_DcaV0DaughtMid,
        kCasc_V0CosPAXi,
        kCasc_DcaV0ToPV,
        kCasc_DcaBachToPV,
        kCasc_ITSTOFtracks, // not used
        kCasc_y,
        kCasc_etaDaugh,
        kCasc_PropLifetXi,
        kCasc_PropLifetOm,
        kCasc_V0RadXi,
        kCasc_DcaMesToPV,
        kCasc_DcaBarToPV,
        kCasc_BacBarCosPA, // not used --> USE
        kCasc_CascRadOm,
        kCasc_V0RadOm,
        kCasc_DcaCascDaughtLow,
        kCasc_DcaCascDaughtHigh,
        kCasc_CascCosPAMid,
        kCasc_CascCosPAHigh,
        kCasc_DcaV0DaughtLow,
        kCasc_DcaV0DaughtHigh,
        kCasc_V0CosPAOmLow,
        kCasc_V0CosPAOmMid,
        kCasc_V0CosPAOmHigh,
        kCasc_CompetingXiMass,
        kCasccutsnum
    }; // kCasc_etaPos, kCasc_etaNeg, kCasc_etaBac, kCasc_kinkidx,
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

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    // cut values setter
    void SetDefOnly(bool);
    void SetCutVal(bool, bool, int, double);
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
    bool fDefOnly;                   //
    double fCasc_Cuts[kCasccutsnum]; //

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
    double fCasc_TrackLength;

    // cut values to be set
    double cutval_Casc[kCasccutsnum];     //
    int nvarcut_Casc[kCasccutsnum];       //
    double varlowcut_Casc[kCasccutsnum];  //
    double varhighcut_Casc[kCasccutsnum]; //

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
    // void FillHistCutVariations(bool, double, bool, bool *, double);
    // functions to allow the correct streaming of the cut variation
    void SetDefCutVals();
    void SetCutVariation(bool, int, int, double, double);
    // void SetDefCutVariations();

    AliCascadeAnalysisMC_Ishaan(const AliCascadeAnalysisMC_Ishaan &);            // not implemented
    AliCascadeAnalysisMC_Ishaan &operator=(const AliCascadeAnalysisMC_Ishaan &); // not implemented

    ClassDef(AliCascadeAnalysisMC_Ishaan, 6);
    // version 6: introduced variations for # of TPC PID clusters and # of ITS-TOF tracks
};

#endif
