#!/usr/bin/env python

import sys, ROOT, math

beam = ROOT.TLorentzVector(0,0,10.6041,10.6041)
targ = ROOT.TLorentzVector(0,0,0,0.938)

with open(sys.argv[1]) as ff:
  lines=[next(ff)]
  for line in ff:
    vvs = line.split()
    if len(vvs)==14:
      px,py,pz,ee,mm = [float(v) for v in vvs[6:11]]
      pid=int(vvs[3])
      if pid==11:
        ele = ROOT.TLorentzVector(px,py,pz,ee)
        q2 = (beam-ele).M2()
      elif pid==211:
        pip = ROOT.TLorentzVector(px,py,pz,ee)
      elif pid==2112:
        neu = ROOT.TLorentzVector(px,py,pz,ee)
        tt = -(neu-targ).M2()

      lines.append(line)
    elif lines:
      if tt<2.5 and math.degrees(pip.Theta())<50:
        for ll in lines:
          print(ll.rstrip())
      lines = [line]

