#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h

#include <GaudiKernel/ToolHandle.h>
#include "G4VUserPrimaryGeneratorAction.hh"
#include <DetSimInterface/IG4PrimaryCnvTool.h>

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(ToolHandle<IG4PrimaryCnvTool>& cnvtool);
  ~PrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);

private:
  ToolHandle<IG4PrimaryCnvTool> tool;
};

#endif

