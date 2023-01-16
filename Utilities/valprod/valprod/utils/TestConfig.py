class TestConfig:

  defaults = {
               'verbose' : False,
               'genLog' : True,
               'logFileName' : None,
               'maxTime' : None,
               'walltime' : None,
               'maxVIR' : None,
               'parser' : False,
               'profile' : False,
               'timeInterval' : 5,
               'monitorBackend': None,  # ps or prmon
               'plottingBackend': 'root', # root or matplotlib
               'perf': False,
               'profileIO': False,
               'fatalPattern' : None,
               'plotRef' : None,
               'cmpOutput' : 'plotcmp.root',
               'histTestMeth' : 'Kolmogorov',
               'histTestCut' : 0.9,
               'logName' : None
             }

  def __init__(self):
    self.config = dict([(k,v) for (k,v) in TestConfig.defaults.items()])

  def update(self, **kwa):
    self.config.update(kwa)

  def setAttr(self, name, value):
    assert type(name) == str, "ERROR: attribute must be of String type!"
    self.config[name] = value

  def getAttr(self, name):
    if name in self.config:
      return self.config[name]
    return None

  def getConfig(self):
    return self.config

globalConfig = TestConfig()
