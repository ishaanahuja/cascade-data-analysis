#ifndef CorrelationTask_H
#define CorrelationTask_H

#include <AliAnalysisTaskSE.h>

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
    AliAODEvent *fAOD; //! input event

    TList *fOutputList; //! output list

    TH1D *fHistdEta, *fHistEtaTrig, *fHistEtaAssoc;
    TH1D *fHistdPhi, *fHistPhiTrig, *fHistPhiAssoc;
    TH2F *fHistdEtadPhi;
    CorrelationTask(const CorrelationTask &);            // not implemented
    CorrelationTask &operator=(const CorrelationTask &); // not implemented

    ClassDef(CorrelationTask, 1);
};

#endif