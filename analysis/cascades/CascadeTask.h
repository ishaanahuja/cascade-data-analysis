#ifndef CascadeTask_H
#define CascadeTask_H

class TH1D;


#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class CascadeTask : public AliAnalysisTaskSE

{
public:
    CascadeTask(const char *name = "default");

    virtual ~CascadeTask();

    virtual void UserCreateOutputObjects();

    virtual void UserExec(Option_t *option);

    virtual void Terminate(Option_t *option);

  
private:
  
    AliAODEvent *fAOD; //! input event

    TList *fOutputList; //! output list

    TH1D *fHistXiPt3, *fHistXiPt4, *fHistXiPt5, *fHistXiPt6, *fHistXiPt7;
    TH1D *fHistOmegaPt3, *fHistOmegaPt4, *fHistOmegaPt5, *fHistOmegaPt6, *fHistOmegaPt7;
    TH1D *fHistEtaXi, *fHistEtaOmega;

    CascadeTask(const CascadeTask &);            // not implemented
    CascadeTask &operator=(const CascadeTask &); // not implemented

    ClassDef(CascadeTask, 1);
};

#endif