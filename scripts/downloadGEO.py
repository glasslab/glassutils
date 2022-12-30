import GEOparse, sys, os, re

class suppress_output:
    def __init__(self, suppress_stdout=False, suppress_stderr=False):
        self.suppress_stdout = suppress_stdout
        self.suppress_stderr = suppress_stderr
        self._stdout = None
        self._stderr = None

    def __enter__(self):
        devnull = open(os.devnull, "w")
        if self.suppress_stdout:
            self._stdout = sys.stdout
            sys.stdout = devnull

        if self.suppress_stderr:
            self._stderr = sys.stderr
            sys.stderr = devnull

    def __exit__(self, *args):
        if self.suppress_stdout:
            sys.stdout = self._stdout
        if self.suppress_stderr:
            sys.stderr = self._stderr

def MsgHelp():
  print("Download public data from GEO")
  print("Please provide accession number")
  exit()
  
def MsgError(strMsg=""):
  print(strMsg)
  exit()

def dwOne(gsm,uID,tryN):
  name_regex = r"[\s\*\?\(\),\.;]"
  strDIR = "%s_%s_%s"%("Supp",gsm.get_accession(),re.sub(name_regex, "_", gsm.metadata["title"][0]))
  extF = os.listdir(strDIR)
  if any(s.endswith('.gz') for s in extF) and not any(s.endswith('.sra') for s in extF):
    print("\tSkip! *.gz existed in %s"%strDIR)
    return()
  with suppress_output(suppress_stdout=True, suppress_stderr=True):
      try:
        gsm.download_SRA('%s@health.ucsd.edu'%uID)
      except Exception as e:
        print(e)
        if tryN<3:
          print("\tTry again! Max 3 times")
          dwOne(gsm,uID,tryN+1)
      else:
        print("Finished %s\n"%gsm.get_accession())

def downloadGEO(strGEO,uID):
  gse = GEOparse.get_GEO(geo=strGEO, destdir="./")
  for gsm_id,gsm in gse.gsms.items():
    print("\n*****\ndownloading %s ..."%gsm_id)
    dwOne(gsm,uID,0)

def main():
  if len(sys.argv)<2:
    MsgHelp()
  
  prjID = sys.argv[1]
  uID = sys.argv[2]
  print("Processing %s"%prjID)
  print("%s"%uID)
  if prjID.startswith("GSE"):
    downloadGEO(prjID,uID)
  else:
    MsgError("Currently only support GSE accessions.")

if __name__ == "__main__":
  main()

