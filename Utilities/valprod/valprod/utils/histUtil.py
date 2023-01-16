import ROOT

def writeMHist(hist, rootfile):
  hist.SetMarkerColor(ROOT.kRed)
  hist.SetLineColor(ROOT.kRed)
  hist.SetLineWidth(2)
  hist.SetStats(0)

  ROOT.gROOT.SetBatch()
  canv = ROOT.TCanvas("c_"+hist.GetName())
  hist.Draw()

  rootfile.cd()
  canv.Write()

def writeMHistPng(hist, file):
  hist.SetMarkerColor(ROOT.kRed)
  hist.SetLineColor(ROOT.kRed)
  hist.SetLineWidth(2)
  hist.SetStats(0)

  ROOT.gROOT.SetBatch()
  canv = ROOT.TCanvas("c_"+hist.GetName())
  hist.Draw()
  canv.SaveAs(file)
