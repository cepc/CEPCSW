#!/usr/bin/env python
#Author: Zhan Li <lizhan@ihep.ac.cn>
#Created [2024-03-07 Thu 14:53]

import ROOT
import math

#read
#path1="/scratchfs/atlas/lizhan/cepc/CEPCSW/test-onlyVXD-single.root"
#path1="/scratchfs/atlas/lizhan/cepc/CEPCSW/test-detsim10.root"

def GetVXDHitMapFromFile(infile_name,out_file_name):
    file = ROOT.TFile(infile_name)
    tree = file.Get("events")
    posx_VXD=[]
    posy_VXD=[]
    posz_VXD=[]
    posr_VXD=[]
    for entry in tree:
        for elem in entry.VXDCollection:
            x=elem.position.x
            y=elem.position.y
            z=elem.position.z
            r = math.sqrt(x*x + y*y)
            posx_VXD.append(x)
            posy_VXD.append(y)
            posz_VXD.append(z)
            posr_VXD.append(r)
    file.Close()
    VXDHitMap = ROOT.TH2F("VXDHitMap","Hit Map of VXD",40,-125,125,20,0,80)
    for i in range(0,len(posz_VXD)):
        VXDHitMap.Fill(posz_VXD[i],posr_VXD[i])
    ofile=ROOT.TFile(out_file_name,"Recreate")
    VXDHitMap.Write()
    ofile.Close()
    
    
def GetSITHitMapFromFile(infile_path,out_file_name):    
    file = ROOT.TFile(infile_path)
    tree = file.Get("events")
    posx_SIT=[]
    posy_SIT=[]
    posz_SIT=[]
    posr_SIT=[]
    for entry in tree:
        for elem in entry.SITCollection:
            x=elem.position.x
            y=elem.position.y
            z=elem.position.z
            r = math.sqrt(x*x + y*y)
            posx_SIT.append(x)
            posy_SIT.append(y)
            posz_SIT.append(z)
            posr_SIT.append(r)
    file.Close()
    SITHitMap = ROOT.TH2F("SITHitMap","Hit Map of SIT",200,-2400,2400,20,120,320)
    for i in range(0,len(posz_SIT)):
        SITHitMap.Fill(posz_SIT[i],posr_SIT[i])
    ofile=ROOT.TFile(out_file_name,"Recreate")
    SITHitMap.Write()
    ofile.Close()

print("test")

def GetSETHitMapFromFile(infile_path,out_file_name):    
    file = ROOT.TFile(infile_path)
    tree = file.Get("events")
    posx_SET=[]
    posy_SET=[]
    posz_SET=[]
    posr_SET=[]
    for entry in tree:
        for elem in entry.SETCollection:
            x=elem.position.x
            y=elem.position.y
            z=elem.position.z
            r = math.sqrt(x*x + y*y)
            posx_SET.append(x)
            posy_SET.append(y)
            posz_SET.append(z)
            posr_SET.append(r)
    file.Close()
    SETHitMap = ROOT.TH2F("SETHitMap","Hit Map of SET",200,-2400,2400,20,1700,2000)
    for i in range(0,len(posz_SET)):
        SETHitMap.Fill(posz_SET[i],posr_SET[i])
    ofile=ROOT.TFile(out_file_name,"Recreate")
    SETHitMap.Write()
    ofile.Close()



#draw
def DrawHitMap(in_file_name,TH2name,out_file_name):
    file = ROOT.TFile(in_file_name)
    TH2_HitMap=file.Get(TH2name)
    c = ROOT.TCanvas("c", "c", 800, 800)
    ROOT.gStyle.SetOptStat(0)
    TH2_HitMap.GetXaxis().SetTitle("Z [mm]")
    TH2_HitMap.GetXaxis().SetTitleOffset(0.77)
    TH2_HitMap.GetXaxis().SetTitleSize(0.05)
    TH2_HitMap.GetYaxis().SetTitle("R [mm]")
    TH2_HitMap.GetYaxis().SetTitleOffset(0.8)
    TH2_HitMap.GetYaxis().SetTitleSize(0.045)
    TH2_HitMap.Draw("COL Z CJUST SAME")

    c.SaveAs(out_file_name+".root")
    c.SaveAs(out_file_name+".pdf")

path1="/scratchfs/atlas/lizhan/cepc/CEPCSW/test-SIT.root"

VXDOutFile="/scratchfs/atlas/lizhan/cepc/CEPCSW/VXDHM.root"
SITOutFile="/scratchfs/atlas/lizhan/cepc/CEPCSW/SITHM.root"
SETOutFile="/scratchfs/atlas/lizhan/cepc/CEPCSW/SETHM.root"
GetVXDHitMapFromFile(path1,VXDOutFile)
GetSITHitMapFromFile(path1,SITOutFile)
GetSETHitMapFromFile(path1,SETOutFile)

DrawHitMap(VXDOutFile,"VXDHitMap","Full_VXDHM_1evt")
DrawHitMap(SITOutFile,"SITHitMap","Full_SITHM_1evt")
DrawHitMap(SETOutFile,"SETHitMap","Full_SETHM_1evt")