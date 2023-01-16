import os
import ROOT

class PlotTester:

  TH1Types = ['TH1D', 'TH1F']
  TH2Types = ['TH2D', 'TH2F']
  GraphTypes = ['TGraph']
  testMeth = ['Chi2', 'Kolmogorov']

  def __init__(self, cfg, plotref, datadir = None):
    self.cfg = cfg
    self.cut = cfg.getAttr('histTestCut')
    self.testName = self.cfg.getAttr('histTestMeth')
    assert self.testName in PlotTester.testMeth, 'ERROR: unknown test method ' + self.testName 
    self.refFileNames = (type(plotref) == list and (plotref,) or ([plotref],))[0]
    self.datadir = datadir
    #print('@PlotTester: refFileNames=%s'%self.refFileNames)
    self.fileNames = []
    for fn in self.refFileNames:
      # We assume that the file has the same name with the reference file
      if self.datadir:
        self.fileNames.append(datadir+'/'+os.path.basename(fn))
      else:
        self.fileNames.append(os.path.basename(fn))
    self.jointFile = None
    self.curPath = ''
    self.tempPlots = []
    self.okay = True
    self.testResult = {}

  def _extractPlot(self, file, refFile):
    # Get all plots in the file recursively, and find the target ones
    self._recursiveGet(file)
    for v in self.tempPlots:
      v.append(refFile.Get(v[0]))

  def _recursiveGet(self, dir):
    for k in dir.GetListOfKeys():
      clName = k.GetClassName()
      objName = k.GetName()
      if clName == 'TDirectoryFile':
        tempPath = self.curPath
        self.curPath += objName + '/'
        self._recursiveGet(k.ReadObj())
        self.curPath = tempPath
      elif clName in PlotTester.TH1Types + PlotTester.TH2Types + PlotTester.GraphTypes:
        self.tempPlots.append([self.curPath + objName, k.ReadObj()])

  def _drawTH1CMP(self, path, plot, refPlot, testResult):
    canvas = self._mkCanvas("compare_"+plot.GetName())
    plot.SetLineColor(ROOT.kRed)
    plot.SetStats(False)
    refPlot.SetStats(False)
    plot.Draw()
    refPlot.Draw("same")
    leg = ROOT.TLegend(0.13, 0.77, 0.3, 0.87)
    leg.AddEntry(plot, "new_" + plot.GetName(), "lp")
    leg.AddEntry(refPlot, "old_" + plot.GetName(), "lp")
    leg.Draw("same")
    if testResult:
      text = ROOT.TPaveText(0.67, 0.74, 0.83, 0.87, "NB NDC")
      text.SetFillColor(0)
      text.SetFillStyle(0)
      text.SetBorderSize(0)
      text.AddText(self.testName + ': ' + str(testResult))
      text.SetTextSize(0.04)
      if testResult < self.cut:
        text.AddText('INCONSISTENT')
        text.SetTextColor(ROOT.kRed)
      else:
        text.AddText('CONSISTENT')
        text.SetTextColor(ROOT.kGreen)
      text.Draw()
    canvas.Write()

  def _drawTH2CMP(self, path, plot, refPlot, testResult):
    canvas = self._mkCanvas("compare_"+plot.GetName())
    canvas.Divide(2,1)
    canvas.GetPad(1).cd()
    plot.SetTitle(plot.GetTitle() + ' (new)')
    plot.Draw('Colz')
    if testResult:
      text = ROOT.TPaveText(0.13, 0.77, 0.3, 0.87, "NB NDC")
      text.SetFillColor(0)
      text.SetFillStyle(0)
      text.SetBorderSize(0)
      text.AddText(self.testName + ': ' + str(testResult))
      text.SetTextSize(0.04)
      if testResult < self.cut:
        text.AddText('INCONSISTENT')
        text.SetTextColor(ROOT.kRed)
      else:
        text.AddText('CONSISTENT')
        text.SetTextColor(ROOT.kGreen)
      text.Draw()    
    canvas.GetPad(2).cd()
    refPlot.SetTitle(refPlot.GetTitle() + ' (old)')
    refPlot.Draw('Colz')
    canvas.Write()

  def _drawGraphCMP(self, path, plot, refPlot):
    canvas = self._mkCanvas("compare_"+plot.GetTitle())
    tmg = ROOT.TMultiGraph()
    plot.SetLineColor(ROOT.kRed)
    refPlot.SetLineColor(ROOT.kBlue)
    tmg.Add(plot)
    tmg.Add(refPlot)
    tmg.Draw('APL')
    leg = ROOT.TLegend(0.13, 0.77, 0.3, 0.87)
    leg.AddEntry(plot, "new_" + plot.GetTitle(), "lp")
    leg.AddEntry(refPlot, "old_" + plot.GetTitle(), "lp")
    leg.Draw("same")
    canvas.Write()
    '''
    canvas = self._mkCanvas("compare_"+plot.GetTitle())
    canvas.Divide(2,1)
    canvas.GetPad(1).cd()
    plot.SetTitle(plot.GetTitle() + ' (new)')
    plot.Draw()
    canvas.GetPad(2).cd()
    refPlot.SetTitle(refPlot.GetTitle() + ' (old)')
    refPlot.Draw()
    canvas.Write()
    '''

  def _mkCanvas(self, name):
    if not self.jointFile:
      jointName = self.cfg.getAttr('cmpOutput')
      self.jointFile = ROOT.TFile(jointName, 'recreate')
    self.jointFile.cd()
    #makeROOTDir(self.jointFile, v[0])
    canvas = ROOT.TCanvas(name)
    return canvas

  def _test(self, h1, h2):
    if self.testName == 'Kolmogorov':
      return h1.KolmogorovTest(h2)
    elif self.testName == 'Chi2':
      return h1.Chi2Test(h2)

  def run(self):
    print("@PlotTester: self.fileNames=%s"%self.fileNames)
    print("@PlotTester: self.refFileNames=%s"%self.refFileNames)
    for i in range(len(self.refFileNames)):
      assert os.path.exists(self.fileNames[i]), '%s does not exsist!' % self.fileNames[i]
      assert os.path.exists(self.refFileNames[i]), '%s does not exsist!' % self.refFileNames[i]
      file = ROOT.TFile.Open(self.fileNames[i])
      refFile = ROOT.TFile.Open(self.refFileNames[i])
      self._extractPlot(file, refFile)
      # plots: [path, plot, plotref]
      for plot in self.tempPlots:
        clName = plot[1].ClassName()
        if clName in PlotTester.TH1Types + PlotTester.TH2Types:
          value = self._test(plot[1], plot[2])
          self.okay = self.okay and value > self.cut
          self.testResult[plot[0]] = value
          if clName in PlotTester.TH1Types:
            self._drawTH1CMP(plot[0], plot[1], plot[2], value)
          elif clName in PlotTester.TH2Types:
            self._drawTH2CMP(plot[0], plot[1], plot[2], value)
        elif clName in PlotTester.GraphTypes:
          self._drawGraphCMP(plot[0], plot[1], plot[2])
      file.Close()
      refFile.Close()
      self.tempPlots = []
      if self.jointFile:
        self.jointFile.Close()

    return self.okay, self.testResult
