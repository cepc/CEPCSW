import re, os, copy

class Parser:

  fatalPattern = [
                    # C
                    '.*Segmentation fault',
                    # SNiPER
                    '.*ERROR:',
                    '.*FATAL:',
                    # Python
                    '.*IOError',
                    '.*ImportError',
                    '.*TypeError',
                    '.*MemoryError',
                    '.*SyntaxError',
                    '.*NameError',
                    '.*RuntimeError',
                    # Other
                    '.*\*\*\* Break \*\*\* segmentation violation',
                    '.*Warning in <TApplication::GetOptions>: macro .* not found',
                 ]

  successPattern = []

  def __init__(self, cfg):
    self.fatal = []
    self.success = {}
    fs = cfg.getAttr('fatalPattern') or Parser.fatalPattern
    for fatal in fs:
      self.fatal.append(re.compile(fatal))
    sp = cfg.getAttr('successPattern') or Parser.successPattern
    for s in sp:
      self.success[s] = re.compile(s)

  def parse(self, word):
    if not word:
      return True, None
    wordl = word.split('\n')
    for pattern in self.fatal:
      for w in wordl:
        match = pattern.match(word)
        if match:
          return False, "Fatal line: %s" %w
    return True, None

  def parseFile(self, file):

    if not os.path.exists(file):
      return False, "Can not find log file: %s" % file

    s = copy.copy(self.success)
    f = open(file, 'r')
    for line in f:
      w = line.strip('\n')
      for pattern in self.fatal:
        match = pattern.match(w)
        if match:
          f.close()
          return False, "Fatal line: %s detected in log file: %s" % (w,file)
      for ps, pattern in s.items():
        match = pattern.match(w)
        if match:
          del s[ps]
    f.close()
    if len(s):
      return False, "Key words not found in %s:\n %s" % (file,'\n'.join(s.keys()))
    return True, None
