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
import glob
import urllib
import gzip
import shutil
import subprocess



def main():
    """
    Checks for MFEprimer binary, downloads & configures if can't be found.
    ---
    Returns path to MFEprimer
    """
    outdir = os.path.dirname(os.path.dirname(__file__))
    path = glob.glob(outdir+"/*mfeprimer*")
    if len(path)==0:
        path = install_mfeprimer()
    return path[0]



def install_mfeprimer(outdir):
    """
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
    outpath = os.path.join(outdir, os.path.basename(mfeprimer))
    try:
        urllib.request.urlretrieve(mfeprimer, outpath)
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
                
                # change permissions to allow usage
                f_out = outpath.replace(".gz","")
                if os.path.exists(f_out):
                    print("Allowing usage as executable....")
                    try:
                        subprocess.call("chmod +x "+f_out, shell=True)
                        return [f_out]
                    except Exception:
                        print("Permissions could not be changed for MFEprimer- try running:")
                        print("    chmod +x "+f_out)
                        print("in terminal/command line, or change permissions manually.")
            
            except Exception:
                print("MFEprimer successfully downloaded to "+outpath+
                      " but could not be unzipped. To proceed, please manually unzip the package"+
                      " and ensure that permissions allow execution (chmod + x <file>) "+
                      "so that binaries are accessible to multiplex_wormhole.")
            
    except Exception:
        print("Failed to download - try manual download from "+mfeprimer)
        print("Save to "+outdir, " then unzip and ensure that permissions allow execution (chmod +x <file>)")



if __name__ == "__main__":
    main()
