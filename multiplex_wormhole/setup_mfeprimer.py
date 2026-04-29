#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: SET UP MFEPRIMER
Purpose: Auto-downloads MFEprimer dependency & stores with multiplex_wormhole src code
Requires internet for download!

Created on Mon Apr 27 09:55:08 2026
@author: maggiehallerud
"""

# load dependency modules
import sys
import os
import urllib
import gzip
import shutil
import subprocess
pwd = os.path.dirname(__file__)



def main():
    """
    ---
    Auto-downloads MFEprimer binaries and stores within multiplex_wormhole src code
    *NOTE*: Requires internet connection.
    """
    # define filepath based on platform
    pf = sys.platform
    if pf == "darwin":
        mfeprimer = "https://github.com/quwubin/MFEprimer-3.0/releases/download/v3.2.7/mfeprimer-3.2.7-darwin-10.6-amd64.gz"
    elif pf in ['linux', 'linux2']:
        mfeprimer = "https://github.com/quwubin/MFEprimer-3.0/releases/download/v3.2.7/mfeprimer-3.2.7-linux-amd64.gz"
    elif pf in ['win32','cygwin','msys']:
        mfeprimer = "https://github.com/quwubin/MFEprimer-3.0/releases/download/v3.2.7/mfeprimer-3.2.7-windows-4.0-amd64.exe.gz"
    else:
        print("An MFEprimer package distribution could not be found for platform: "+pf)
        
    # attempt to download
    print("Downloading MFEprimer......")
    outpath = os.path.join(pwd, os.path.basename(mfeprimer))
    try:
        urllib.request.urlretrieve(mfeprimer, outpath)
    except Exception:
        print("Failed to download - try manual download from "+mfeprimer)
        print("Save to "+pwd)
    
    # unzip file
    if outpath.endswith(".gz"):
        print("Unzipping downloaded file......")
        try:
            with gzip.open(outpath, "rb") as f_in:
                with open(outpath.replace(".gz",""), "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            # delete gz
            try:
                os.remove(outpath)
            except Exception:
                pass
        except Exception:
            print("MFEprimer successfully downloaded to "+outpath+
                  " but could not be unzipped. To proceed, please manually unzip the package"+
                  " so that binaries are accessible to multiplex_wormhole.")
    
    # change permissions to allow usage
    f_out = outpath.replace(".gz","")
    if os.path.exists(f_out):
        print("Allowing usage as executable....")
        try:
            subprocess.call("chmod +x "+f_out, shell=True)
        except Exception:
            print("Permissions could not be changed for MFEprimer- try running:")
            print("    chmod +x "+f_out)
            print("in terminal/command line, or change permissions manually.")



if __name__ == "__main__":
    main()
