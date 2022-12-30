import GEOparse, sys, os

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

def downloadGEO(strGEO,uID):
  gse = GEOparse.get_GEO(geo=strGEO, destdir="./")
  for gsm_id,gsm in gse.gsms.items():
    print("downloading %s ..."%gsm_id)
    with suppress_output(suppress_stdout=True, suppress_stderr=True):
      gsm.download_SRA('%s@health.ucsd.edu'%uID)
    print("Finished %s"%gsm_id)

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

