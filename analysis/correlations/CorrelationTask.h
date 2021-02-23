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
    CorrelationTask();

    CorrelationTask(const char *name);

    virtual ~CorrelationTask();

    virtual void UserCreateOutputObjects();

    virtual void UserExec(Option_t *option);

    virtual void Terminate(Option_t *option);

private:
    Bool_t fFillMixed;             // enable event mixing (default: ON)
    Int_t fMixingTracks;           // size of track buffer for event mixing
    AliEventPoolManager *fPoolMgr; //! event pool manager

    AliAODEvent *fAOD; //! input event

    TList *fOutputList; //! output list

    TH1D *fHistdEta, *fHistEtaTrig, *fHistEtaAssoc;
    TH1D *fHistdPhi, *fHistPhiTrig, *fHistPhiAssoc;
    TH2D *fHistdPhidEta;
    THnSparseD *fHistMixC1, *fHistMixC2;                              // dPhi vs. dEta, mixed events
    CorrelationTask(const CorrelationTask &);            // not implemented
    CorrelationTask &operator=(const CorrelationTask &); // not implemented

    ClassDef(CorrelationTask, 1);
};

#endif