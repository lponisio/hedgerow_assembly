"""runNetworkChangePoint.py - module to run changepoint detection on a sequence of networks
    Copyright (C) 2014 Leto Peel

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA"""



if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run network changepoint detection on a sequence of networks")
    
    parser.add_argument("nodenamesfile", help="Input node names file e.g. names.lut")
    parser.add_argument("windowsize", help="Length of sliding window")
    parser.add_argument("networkfilesequence", help='Input sequence of network files: e.g. "network1.pairs network2.pairs network3.pairs"')
    parser.add_argument("-p","--path", help="Path to files (if not in current directory)")
    args = parser.parse_args()




    import changepointDetection

    if args.path:
        path=args.path
    else:
        path="."

    networkFileSequence=args.networkfilesequence.split()
    
    cpDetector = changepointDetection.anomalyDetection()
    window = int(args.windowsize)
    cpDetector.detectAnomaliesInSequence(networkFileSequence,args.nodenamesfile,window,path)
