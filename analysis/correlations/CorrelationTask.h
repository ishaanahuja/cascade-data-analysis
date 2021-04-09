#ifndef CorrelationTask_H
#define CorrelationTask_H

class TH1D;
class TH2F;
class THnSparse;
class AliEventPoolManager;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include "THnSparse.h"
class CorrelationTask : public AliAnalysisTaskSE

{
public:
    CorrelationTask(const char *name = "default");

    virtual ~CorrelationTask();

    virtual void UserCreateOutputObjects();

    virtual void UserExec(Option_t *option);

    virtual void Terminate(Option_t *option);

    void SetEventMixing(Bool_t EventMixing = kTRUE) { fEventMixing = EventMixing; }
    Bool_t GetEventMixing() { return fEventMixing; }

private:
    Bool_t fEventMixing;           // enable event mixing (default: ON)
    Int_t fMixingTracks;           // size of track buffer for event mixing
    AliEventPoolManager *fPoolMgr; //! event pool manager

    AliAODEvent *fAOD; //! input event

    TList *fOutputList; //! output list

    TH1D *fHistdEta, *fHistEtaTrig, *fHistEtaAssoc;
    TH1D *fHistdPhi, *fHistPhiTrig, *fHistPhiAssoc;
    TH2D *fHistdPhidEta;
    THnSparseD *fHistMixC1, *fHistMixC2;                 // dPhi vs. dEta, mixed events
    CorrelationTask(const CorrelationTask &);            // not implemented
    CorrelationTask &operator=(const CorrelationTask &); // not implemented

    ClassDef(CorrelationTask, 1);
};

/*  AliMixBasicParticle class contains only quantities 
 *	required for the analysis in order to reduce memory consumption for event mixing.
 */
class AliMixBasicParticle : public AliVParticle
{
  public:
    AliMixBasicParticle(Float_t eta, Float_t phi, Float_t pt/*, Short_t candidate*/)
      : fEta(eta), fPhi(phi), fpT(pt)/*, fCandidate(candidate)*/
    {
    }
    virtual ~AliMixBasicParticle() {}

    // kinematics
    virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pt() const { return fpT; }
    virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Phi()        const { return fPhi; }
    virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }


    virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Eta()        const { return fEta; }
    virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }

    virtual Short_t Charge()      const { AliFatal("Not implemented"); return 0; }
    virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
    // PID
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }

    // virtual Short_t WhichCandidate()      const { return fCandidate; }

  private:
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT
    // Short_t fCandidate;   // V0 candidate: 1 - K0sig, 2 - Lamsig, 3 - Alamsig, 4 - K0bg, 5 - Lambg, 6 - Alambg

    ClassDef( AliMixBasicParticle, 1); // class required for event mixing
};

#endif