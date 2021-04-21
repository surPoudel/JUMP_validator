#!/usr/bin/env python
# coding: utf-8

# In[27]:


import sys
import pandas as pd
import glob
import os
import re
from collections import Counter
import configparser
config = configparser.ConfigParser()

import pyteomics as pyteomics
import pandas as pd
from pyteomics import mzxml, auxiliary, pepxml

from pyteomics import mass
# from ipynb.fs.full.dunpFunctions_JUMP_l import *
import numpy as np
import scipy.stats as ss
import matplotlib
from operator import attrgetter
from scipy.stats import hypergeom

import matplotlib.pyplot as plt
import seaborn as sns
import time
from os.path import dirname
import itertools   
# In[2]:

'''
These are jump validator specific functions 
'''


def makeJumpPlainPeptideAndModCols(row,sta_AA,jump_mod_dict):
  
  peptide = row.Peptides
  pepNoFlank = peptide.split(".")[1]
  pepNoFlankList = list(pepNoFlank)
  
  plain_peptide_list =[]
  for aa in pepNoFlankList:
    if aa in pyteomics.parser.std_amino_acids:
      plain_peptide_list.append(aa)
  plain_peptide = "".join(plain_peptide_list)
  modifications = ",".join(DynModPositionPeptide(pepNoFlank,sta_AA,jump_mod_dict))
  return pd.Series([plain_peptide,modifications])


def prepareInput(inputFile):
    df = pd.read_excel(inputFile)
    df["spectrum"] = df.Exp+"."+df.scan.astype("str")+"."+df.charge.astype("str")
    return df

def splitTwoCols(df):
    parseCols = ['Protein Accession', 'Site_LScore']
    modDict = {}
    for col in list(df.columns):
        modDict[col]=[]

    for row in df.itertuples():
        if "," in row[-2]:
            mod_list = row[-3].split(",")
            score_list = row[-2].split(",")
            for index, val in enumerate(mod_list):
                #c = [x for x in a if x not in b]
                new_list = [x for x in list(df.columns) if x not in parseCols]
                for key in new_list:
                    modDict[key].append(row[list(df.columns).index(key)+1])

                modDict[parseCols[0]].append(val)
                modDict[parseCols[1]].append(score_list[index])
    return modDict


def DynModPositionPeptide(PepSeq,sta_AA,jump_mod_dict):
    mods = []
    mod_pos = []
    #   aa_site = []
    aalist=list(PepSeq)
    #   print (aalist)
    for i,aa in enumerate(aalist):
        if aa in jump_mod_dict.keys():
            substract_val = len(mod_pos)
            pos = i-substract_val
            aaRes = PepSeq[i-1]
            mod = str(pos)+"_"+"V_"+str(jump_mod_dict[aa])
            mods.append(mod)
            mod_pos.append(pos)

        if aa in sta_AA.keys():
            mod = str(i+1-len(mod_pos))+"_"+"S_"+str(sta_AA[aa])
            mods.append(mod)
        if "n" in sta_AA.keys():
            nterm = "1_S_"+str(sta_AA["n"])+"_n"
            mods.append(nterm)

    return mods

def mkdir(dir1):
  cmd = "mkdir "+dir1
  try:
    os.system(cmd)
  except:
    "Dynamic Modification Directory exits!!"



def createSymbolDict(posMassDict):
  seqSymbol = {}
  for vals in posMassDict.values():
    if vals==229.162931:
      seqSymbol["(229.162931)"]  = "%"
    elif vals==15.99492:
      seqSymbol["(15.99492)"]  = "@"
    elif vals == 57.021464:
      seqSymbol["(57.021464)"]  = "$"
    else:
      newVal = "("+str(vals)+")"
      seqSymbol[newVal]  = "*"
  return seqSymbol



def createSymbolDict2(posMassDict, jump_mod_dict, sta_AA):
  seqSymbol = {}
  for keys, vals in posMassDict.items():
    for keyJ, valsJ in jump_mod_dict.items():
        if vals == valsJ:
            seqKey = "("+str(vals)+")"
            seqVal = keyJ
            seqSymbol[seqKey] = seqVal
    
    for keySta, valsSta in sta_AA.items():
        if vals == valsSta:
            seqKey = "("+str(vals)+")"
            seqVal = "+"+str(int(valsSta))
            seqSymbol[seqKey] = seqVal
    
  return seqSymbol


def rightInterval(row):
  val = int(row.right)
  if val == 10000:
    new_val = "<10e4"
  if val == 100000:
    new_val = "10e4-10e5"
  if val == 1000000:
    new_val = "10e5-10e6"
  if val == 10000000:
    new_val = "10e6-10e7"
  if val == 100000000:
    new_val = "10e7-10e8"
  if val == 1000000000:
    new_val = "10e8-10e9"
  if val == 10000000000:
    new_val = "10e9-10e10"
  return new_val

def paramsOptToDict(parameter):
    dict_parameter = {} #dictionary to create using paramter option 
    listParameter = parameter.split(",")
    for x in listParameter:
        key_val = x.split(":")
        dict_parameter[key_val[0]]=key_val[1]
    return dict_parameter

#this function helps to quickly concatenate the txt files converting them to dataframe

def txtToConcatDF(pepVariantsScoreTxt):
  start = time.time()
  super_x =[] #dataframe are updated in this list
  df = pd.read_csv(pepVariantsScoreTxt[0], delimiter="\t")
  super_x.append(df)
  
  for x in range(1,len(pepVariantsScoreTxt)):
    df2 = pd.read_csv(pepVariantsScoreTxt[x], delimiter="\t")

    super_x.append(df2)
    
  finalDF = pd.concat(super_x, axis=0) #concat the list of dataframe along the axis=0
  end = time.time()
  print ("total time required to concatenate all mzXML files = ", end-start,"seconds\n")
  return finalDF


def modsForReport(mods, peptide):
  newMods = []
  modSplit = mods.split(",")
  for x in modSplit:
    xSplit = x.split("_")
    pos = xSplit[0]
    modKind = xSplit[1]
    aa = peptide[int(pos)-1]
    if modKind == "S":
      modType = "St"
    else:
      modType = "Dy"
    finalMods = str(pos)+"_"+aa+"_"+str(xSplit[2])+"_"+modType
    newMods.append(finalMods)
  modsPrint = ",".join(newMods)
  return modsPrint


def reformat_dataframe2(df, posMassDict,jump_mod_dict, sta_AA, match_list):
  matched_list_array = list(np.array(match_list).round(4))
  dfNew = df.copy() #duplicate dataframe to work

  df2 = dfNew.rename(columns={"Peptide_Mod_Seq":"Seq"})
  
  seqSymbol = createSymbolDict2(posMassDict,jump_mod_dict, sta_AA)

  matched_list_array = list(np.array(match_list).round(4))

  for name, values in df2.iteritems():
    for val in range(0, df2.shape[0]):
  #     print('{name}: {value}'.format(name=name, value=values[val]))
      value = values[val]
      
      if name == "Seq":
        for keys in seqSymbol.keys():
          if keys in value:
            new_val = value[0]
            df2 = df2.replace(value,value[0]+seqSymbol[keys]) 

  return df2, seqSymbol


def excelWriter2(writer, df,  worksheetName, figure1, figure2,spectrumSplit,xcorr, prob,lscoreSite, massSeriesLength, matched, seqSymbol,jump_mod_dict, sta_AA, peptide, N=34, updater=0):
  text1 = "Spectrum = "+spectrumSplit[0]+"; Plain Peptide Sequence = "+spectrumSplit[1]+"; Modifications on the Peptide = "+modsForReport(spectrumSplit[2], peptide)
  text2 = "Xcorr = "+xcorr.astype("str")+"; Localization Probablity Score = "+ str(prob)+"; Localization Site:Score = "+ lscoreSite
  text3 = "Figure a) Intensity Plot"
  text4 = "Figure b) Mass deviation of matched ions (ppm tolerance) Plot"
  text5 = "Figure c) Ion series b and y ions with different charge states Plot" 
  text6 = "Matched ion = "+matched
  text7 = "Symbols and their meaning"
  text8 = ""
  for symbols in seqSymbol.keys(): 
    text8+=seqSymbol[symbols]+"="+symbols+"; "
    
  
  df.to_excel(writer, sheet_name=worksheetName, startrow=N, index=None)
  workbook  = writer.book
  worksheet = writer.sheets['Sheet1']
  format1 = workbook.add_format({'num_format': '0.00'})
  worksheet.set_column('C:C', None, format1)  # Adds formatting to column C
    
    

  length = 34+massSeriesLength+1
  # Get the xlsxwriter workbook and worksheet objects.
  workbook  = writer.book
  worksheet = writer.sheets[worksheetName]
    
  worksheet.write(updater, 0, text1)
  worksheet.write(updater+1, 0, text2)
  worksheet.write(updater+2, 0, text3)

  location = "A"+str(updater+3)
  worksheet.insert_image(location,figure1)
  
  worksheet.write(updater+16, 0, text4)  
    
  location2 = "A"+str(updater+18)
  worksheet.insert_image(location2,figure2)
 
  worksheet.write(updater+33, 0, text5)  

  worksheet.write(updater+length+2, 0, text6)
  worksheet.write(updater+length+3, 0, text7)
  worksheet.write(updater+length+4, 0, text8)

  updater = updater+length+6
  N = updater + 35

  return writer, updater, N

#this function is different in consensus library as we cannot sum modifications there
#as the PTMs are used to extract unimod information
def spectrumToDict(spectrum):
    dict1 = {}
    spectrumCommaSplit = spectrum.split(",")
    for x in spectrumCommaSplit:
        y=x.split("_")
        if y[0] not in dict1.keys():
            dict1[y[0]] = [float(y[2])]
        else:
            dict1[y[0]].append(float(y[2]))
    dict2 = {}

    for key in dict1.keys():
        value = np.sum(list(set(dict1[key])))
        dict2[key] = value
    return dict2

def displaySeqTable(row,massPosDict): #massPosDict example is massPosList[0]
  aa = row.Seq
  pos = row["b-series"]
  if str(pos) in massPosDict.keys():
    massShift = massPosDict[str(pos)]
    aa2 = aa+"("+str(massShift)+")"
  else:
    aa2 = aa
  return aa2
  
#These functions are used for jump -l and jump -validator

#This is time consuming so now each mzxml file will be individually parsed and look for id in txt file so will be submitted
#as single job per mzXML file and later converted to txt file
def mzmlDF_SingleFile(mzxmlFile, spec_list):
  start = time.time()
  super_x =[]
  x1 = pyteomics.mzxml.read(mzxmlFile,iterative=True)  #reading mzXML file using pyteomics iterative (bool, optional) – Defines whether iterative parsing should be used. It helps reduce memory usage at almost the same parsing speed. Default is True.
  df = pd.DataFrame([x for x in x1])  #dataframe of the mzXML file
  ms2 = df.copy().loc[df.msLevel==2]     #ms2 level scans
  run = mzxmlFile.split("/")[-1].split(".mzXML")[0]
  ms2["expScan"] = run+"."+ms2.num.astype("str")
  ms2F = ms2.copy().loc[ms2.expScan.isin(spec_list)]
  ms2F.rename(columns={"num":"scan","m/z array":"m/z","intensity array":"intensity"}, inplace=True)
  reqCols = ["expScan","scan","m/z","intensity"]
  ms2F2 = ms2F[reqCols]
  return ms2F2


#each mzXML or ms2 file are now converted to dataframe based on teh ID in the text file and converted to single txt file for each raw file based on the id
#later these new txt file are used to merge together the txt files into dataframe
def peptideVariantsScoringEachFile(mzxmlFile, df_phoSym,spec_list,dyn_AA, tol, psmTxtFile, ms2_fileType, ion_types, ionLoss):
    if ms2_fileType.upper() == "MZXML":
        testID_mzXML = mzmlDF_SingleFile(mzxmlFile, list(set(df_phoSym["expScan"])))
        
        new_list = list(set(testID_mzXML.expScan))
        
        #extract these scans from df_pho (ID.txt file)
        df_pho = df_phoSym.copy().loc[df_phoSym.expScan.isin(new_list)]
        
        mz_spectrumDict = dict(zip(testID_mzXML.expScan, testID_mzXML["m/z"]))
        int_spectrumDict = dict(zip(testID_mzXML.expScan, testID_mzXML["intensity"]))

        df_pho["m/z"]=df_pho["expScan"].map(mz_spectrumDict)
        df_pho["intensity"]= df_pho["expScan"].map(int_spectrumDict)
        ms2DF = df_pho[["spectrum","exp", "scan", "charge","expScan", "m/z", "intensity"]]
    else:
#         ms2DF=msToDF(mzXML_list, list(set(df_pho["spectrum"])))
        df = ms2ToDf(mzxmlFile)
        ms2DF = df.copy().loc[df.spectrum.isin(spec_list)]
        run = mzxmlFile.split("/")[-1].split(".ms2")[0]
        ms2DF["expScan"] = run+"."+ms2DF.scan.astype("str")
        new_list = list(set(ms2DF.spectrum))
        
        #extract these scans from df_pho (ID.txt file)
        df_pho = df_phoSym.loc[df_phoSym.spectrum.isin(new_list)]
        
        
    pval_list = []
    revPval_list = []
    seq_probability_list = []
    spectrum_list = []
    peptide_list = []
    all_mods_list =[]
    # original_Peptide_mod = []

    all_scans_ms2_dict = {}

    all_scanDF = {}
    counter = 0
    mz_cols = list(df_pho.columns)
    for row in df_pho.itertuples():
        scan = str(row[mz_cols.index("scan")+1])
        peptide = row[mz_cols.index("plain_peptide")+1]
        #original peptide modification

    #     ori_peptide = row[mz_cols.index("Peptides")+1]

        spectrum = row[mz_cols.index("spectrum")+1]
        maxCharge = int(row[mz_cols.index("charge")+1])
        mods = row[mz_cols.index("modifications")+1]
        mz = ms2DF.loc[ms2DF.scan == scan]['m/z']
        newDF = ms2DF.loc[ms2DF["scan"] == scan]
        # precMZ = list(newDF["prec_MZ"])[0]

        intensity = ms2DF.loc[ms2DF.scan == scan]['intensity']
        exp_mz_list = list(mz.tolist()[0])
        intensity_list = list(intensity.tolist()[0])
        #print ("Peptide = ", peptide)
        #print ("Mods = ", mods)
        #print ("dyn_AA = ", dyn_AA)
        #print ("Spectrum = ",spectrum)

        modsVariantsAll, massPosList = statModDynModAllVariants(peptide, mods, dyn_AA)
        #print ("modVariantsAll = ", modsVariantsAll)
        #print ("massPosList = ", massPosList)


        one_scan_revPval_list = []
        pval_Rev_sum = 0.0
        counter+=1
        #Turn this off if you have already created ms2 files or it is not required
        #createMS2FileFromMzXML(scan, precMZ, maxCharge, exp_mz_listC, intensity_list, MH, dir2+"/"+spectrum+".ms2")
        print ("The spectrum analyzed = ", counter) 
        for i,dicts in enumerate(massPosList):
            
    #         original_Peptide_mod.append(ori_peptide)

            df_MOD = ionSeriesIonLossSpeRes(peptide,maxcharge=maxCharge,massPosDict=dicts,useMod ="Yes")
            #         if spectrum == "w010.23411.3":
            #           checkDF = df_MOD

            pvalue, true_ions_list, matched_int_list, true_ions_ppm_list = hypergeometricTest(exp_mz_list,intensity_list, df_MOD, tol=tol, ion_types = ion_types, ionLoss = ionLoss)
            keySpectrum = spectrum+"\t"+peptide+"\t"+modsVariantsAll[i]

            if len(true_ions_list) >=1:
                peptide_list.append(peptide)
                spectrum_list.append(spectrum)
                all_mods_list.append(modsVariantsAll[i])                

                pvalReverse = 1/pvalue
                pval_Rev_sum+=pvalReverse

                pval_list.append(pvalue)
                revPval_list.append(pvalReverse)
                one_scan_revPval_list.append(pvalReverse)
                #keySpectrum = spectrum+"\t"+peptide+"\t"+modsVariantsAll[i]

            else:
                print ("No true matches for ", keySpectrum)
        for valRev in one_scan_revPval_list:
            seqPr = valRev/pval_Rev_sum*100.00
            seq_probability_list.append(seqPr)



    all_scanDF = {"spectrum":spectrum_list,"plain_peptide":peptide_list,"modifications":all_mods_list,"pvalue":pval_list,"1/p-value":revPval_list,"SequenceProbablity":seq_probability_list}  
    ptmsite = pd.DataFrame(all_scanDF)
    ptmsite.to_csv(psmTxtFile, sep="\t", index=None) #generates a single file for each dataframe for all the ids that are present in 
    #raw file and corresponding scans in the ID.txt
    


def peptideVariantsScoring(df_pho,ms2DF,dyn_AA, tol):
    
    pval_list = []
    revPval_list = []
    seq_probability_list = []
    spectrum_list = []
    peptide_list = []
    all_mods_list =[]
    # original_Peptide_mod = []

    all_scans_ms2_dict = {}

    all_scanDF = {}
    counter = 0
    mz_cols = list(df_pho.columns)
    for row in df_pho.itertuples():
        scan = str(row[mz_cols.index("scan")+1])
        peptide = row[mz_cols.index("plain_peptide")+1]
        #original peptide modification

    #     ori_peptide = row[mz_cols.index("Peptides")+1]

        spectrum = row[mz_cols.index("spectrum")+1]
        maxCharge = int(row[mz_cols.index("charge")+1])
        mods = row[mz_cols.index("modifications")+1]
        mz = ms2DF.loc[ms2DF.scan == scan]['m/z']
        newDF = ms2DF.loc[ms2DF["scan"] == scan]
        # precMZ = list(newDF["prec_MZ"])[0]

        intensity = ms2DF.loc[ms2DF.scan == scan]['intensity']
        exp_mz_list = list(mz.tolist()[0])
        intensity_list = list(intensity.tolist()[0])
        #print ("Peptide = ", peptide)
        #print ("Mods = ", mods)
        #print ("dyn_AA = ", dyn_AA)
        #print ("Spectrum = ",spectrum)

        modsVariantsAll, massPosList = statModDynModAllVariants(peptide, mods, dyn_AA)
        #print ("modVariantsAll = ", modsVariantsAll)
        #print ("massPosList = ", massPosList)


        one_scan_revPval_list = []
        pval_Rev_sum = 0.0
        counter+=1
        #Turn this off if you have already created ms2 files or it is not required
        #createMS2FileFromMzXML(scan, precMZ, maxCharge, exp_mz_listC, intensity_list, MH, dir2+"/"+spectrum+".ms2")
        print ("The spectrum analyzed = ", counter) 
        for i,dicts in enumerate(massPosList):
            
    #         original_Peptide_mod.append(ori_peptide)

            df_MOD = ionSeriesIonLossSpeRes(peptide,maxcharge=maxCharge,massPosDict=dicts,useMod ="Yes")
            #         if spectrum == "w010.23411.3":
            #           checkDF = df_MOD

            pvalue, true_ions_list, matched_int_list, true_ions_ppm_list = hypergeometricTest(exp_mz_list,intensity_list, df_MOD, tol=tol, ion_types = ["b","y"], ionLoss = ["H2O","NH3"])
            keySpectrum = spectrum+"\t"+peptide+"\t"+modsVariantsAll[i]

            if len(true_ions_list) >=1:
                peptide_list.append(peptide)
                spectrum_list.append(spectrum)
                all_mods_list.append(modsVariantsAll[i])                

                pvalReverse = 1/pvalue
                pval_Rev_sum+=pvalReverse

                pval_list.append(pvalue)
                revPval_list.append(pvalReverse)
                one_scan_revPval_list.append(pvalReverse)
                #keySpectrum = spectrum+"\t"+peptide+"\t"+modsVariantsAll[i]

            else:
                print ("No true matches for ", keySpectrum)
        for valRev in one_scan_revPval_list:
            seqPr = valRev/pval_Rev_sum*100.00
            seq_probability_list.append(seqPr)



    all_scanDF = {"spectrum":spectrum_list,"plain_peptide":peptide_list,"modifications":all_mods_list,"pvalue":pval_list,"1/p-value":revPval_list,"SequenceProbablity":seq_probability_list}  
    ptmsite = pd.DataFrame(all_scanDF)
    return ptmsite

def valAddKey(dict1, key, val):
  if key not in dict1.keys():
    dict1[key] = [val]
  else:
    dict1[key].append(val)
  return dict1


# In[6]:


def getPhoPosition(massPosDict):
    phoPos = []
    for key, value in massPosDict.items():
        if "79.96" in str(value):
            phoPos.append(int(key))
    return sorted(phoPos)
            

def ionSeriesIonLossSpeRes(peptide,massPosDict, maxcharge=1,useMod ="Yes"):
  #first check if the mods has phosphorylation if yes get the maximum and minimum position
  phoPos = getPhoPosition(massPosDict)
  minPos = 1  #assigns minPosition which is later updated
  maxPos = len(peptide) #assigns maximum position which is later updated
  if len(phoPos) >= 1:
    minPos = np.min(phoPos) #computes minimum phospho position so tha b and y ions can be computed according for phospho loss
    maxPos = np.max(phoPos) #computes maximum phospho position so that b and y ions can be computed according ly

  if useMod != "Yes": #this checks whether the modificaiton is searched or not, generally this is always Yes for localization 
    massPosDict = {0:0.0} #if modificaiton is no than the dictionary has no information
  h2o = mass.calculate_mass(formula='H2O') #h2o mono mass
  co = mass.calculate_mass(formula='CO') #co mono mass
  nh3 = mass.calculate_mass(formula='NH3') #nh3 mono mass
  xmassMore = co-mass.calculate_mass(formula='H2') #hydrogen mono mass
  proton = mass.calculate_mass(formula='H+') #proron mono mass
  hydrogenAtom = mass.calculate_mass(formula='H') #hydrogen atom mono mass
  hpo3 = mass.calculate_mass(formula='HPO3') #computes hpo3 mono mass
  h3po4 = mass.calculate_mass(formula='H3PO4') #computes h3po4 mono mass

  all_ions_dict = {}#iniitates a dictionary for theoretical ions which is later conveted to dataframe
#possible ion loss or no ion loss
  ionloss = {"":0,"-H2O":h2o,"-HPO3":hpo3,"-H3PO4":h3po4,"-NH3":nh3}
  
#this section computes a,b,c ions  
  for i in range(1, len(peptide)+1):
    addModMass = 0.0
    for massKey in massPosDict.keys():
      if int(massKey) <= i:
        addModMass += float(massPosDict[massKey])
        
    valAddKey(all_ions_dict,"Seq",peptide[i-1])
#     print (peptide[0:i-1])
    for losses in ionloss.keys():
      if losses == "-H2O":
#water loss if the amino acid sequence is STED
        fate = checkAA(peptide[0:i],["S","T","E","D"])
        for charge in range(1, maxcharge+1):
          if (fate == "True") and (i < len(peptide)):
            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
          else:
            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
      if losses == "-NH3":
#ammonia loss if the amono acid sequence is RKQN
        fate = checkAA(peptide[0:i],["R","K","Q","N"])
        for charge in range(1, maxcharge+1):
          if (fate == "True") and (i < len(peptide)):
            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
            valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
          else:
            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
                      
      if (losses == "-H3PO4") or (losses == "-HPO3"):
        for charge in range(1, maxcharge+1):
          if "79.96" in str(massPosDict.values()):
#             fate = checkAA(peptide[0:i],["S","T","Y"])
#             if fate == "True":
#phosoho loss if massPosDict have phopsho loss in it. If yes based on minPos and maxPos ions are calculated
            if (i >= minPos) and (i < len(peptide)):
              valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
              valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
            else:
              valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
              valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
                      
      if losses == "":
        for charge in range(1, maxcharge+1):
          if i < len(peptide):
            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
            valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
            valAddKey(all_ions_dict,"c"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge)
            valAddKey(all_ions_dict,"c(-1)"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge-hydrogenAtom)
            valAddKey(all_ions_dict,"c1"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge+hydrogenAtom)
            valAddKey(all_ions_dict,"c2"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge+(2*hydrogenAtom))

          else:
            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"c"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"c(-1)"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"c1"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"c2"+losses+"+"+str(charge),"---")

#this section computes x,y,z ions
  for i in range(0, len(peptide)):
    
    addModMass = 0.0
    for massKey in massPosDict.keys():
      if int(massKey) > i:
        addModMass += float(massPosDict[massKey])
    for losses in ionloss.keys():  
      if losses == "-H2O":
        fate = checkAA(peptide[i:len(peptide)],["S","T","E","D"])
        for charge in range(1, maxcharge+1):
          if (fate == "True") and (i >=1):
            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
          else:
            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
      
      if losses == "-NH3":
        fate = checkAA(peptide[i:len(peptide)],["R","K","Q","N"])
        for charge in range(1, maxcharge+1):
          if (fate == "True") and (i >=1):
            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
          else:
            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
      
      
      if (losses == "-H3PO4") or (losses == "-HPO3"):
#         fate = checkAA(peptide[i:len(peptide)],["S","T","Y"])
        for charge in range(1, maxcharge+1):
          if "79.96" in str(massPosDict.values()):
            if (i < maxPos) and (i >=1):
              valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
            else:
              valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
      if losses == "":
        for charge in range(1, maxcharge+1):
          if (i >=1):
            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
            valAddKey(all_ions_dict,"x"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)+xmassMore-ionloss[losses]+addModMass)/charge)

            valAddKey(all_ions_dict,"z"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+hydrogenAtom)
            valAddKey(all_ions_dict,"z1"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+(2*hydrogenAtom))
            valAddKey(all_ions_dict,"z2"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+(3*hydrogenAtom))
            valAddKey(all_ions_dict,"z3"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+(4*hydrogenAtom))
          else:
            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"x"+losses+"+"+str(charge),"---")

            valAddKey(all_ions_dict,"z"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"z1"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"z2"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"z3"+losses+"+"+str(charge),"---")
         
                

  df = pd.DataFrame(all_ions_dict)
  return df
#for z ion calculation please see the notes:
#https://prospector.ucsf.edu/prospector/html/misc/faq.htm#z_ions
# In[5]:
def ionMatches(exp_mz_list, intensity_list, theor_fragments_df, ion_types = ["b","y"], ionLoss = ["NH3","H2O"], tol=10):
    
  cols = theor_fragments_df.columns
  ionsCols = []
  for ions in ion_types:
    for ionC in cols:
      if ions+"+" in ionC:
        ionsCols.append(ionC)
      for ionLosses in ionLoss:
        if ions+"-"+ionLosses in ionC:
          ionsCols.append(ionC)   

  checkDF = theor_fragments_df[ionsCols]
  newDF = checkDF.replace("---",np.nan)
  arr = newDF.to_numpy().flatten()
  theor_ions_list = list(arr[~np.isnan(arr)])
    
  true_ions_ppm_list = []  
  true_ion_int_list = []
  true_theretical_ions_list = []
  for i,exp_ion in enumerate(exp_mz_list):
    for theor_ion in theor_ions_list:
      ppm = ppmCalc(theor_ion,exp_ion, tol=tol)
      
      if (ppm <= float(tol)) and (ppm >= (-1*float(tol))):
        true_theretical_ions_list.append(theor_ion)
        true_ion_int_list.append(intensity_list[i])
        true_ions_ppm_list.append(ppm)
        
  return true_theretical_ions_list,true_ion_int_list,true_ions_ppm_list


def checkAA(pepSeq, check_list):
  update_val = "False"
  for aa in list(pepSeq):
    if aa in check_list:
      update_val = "True"
  return update_val





def characterToIndices(string, char):
  string_list = list(string)
  indices = [i for i, x in enumerate(string_list) if x == char]
  return indices


# In[2]:


def multipleDynModToIndices(peptide, dynAA):
  #This portion was used when Yuxin made a program to replace n-term with J
  # if peptide.startswith("K"):
  #   peptide = peptide
  # else:
  #   peptide = "J"+peptide[1:]
  finalIndices = []
  for val in dynAA:
    listAA = list(val)
    allIndices = []
    for aa in listAA:
      charList = characterToIndices(peptide, aa)
      allIndices+=charList
    finalIndices.append(allIndices)
  return finalIndices





def statModDynModAllVariants(peptide, mods,dyn_AA): #mods is modification column from comet, peptide is plain peptide
  dynAA = []          #list that updates dynamic amino acids in the sequential order
  dynDelMass = []     #list that updates delta masses of dynamic modification  
  modsSplit = mods.split(",")              
  all_static = {}  #all static modification 
  all_static_nterm = {} #all static n term mods
  dynamicNterm = {}   #stores dynamic modification Nterm
  all_dyn_mod_combinations = []   #list that updates all dynamic modifications combinations

  countVariableMods = 0  #does not include n-term dynamic here
  for x in modsSplit:
    new_mod = ""
    new_x = x.split("_")
    changedMods = ""
    if (new_x[1] == "S") and (new_x[-1] == "n"):   #checks for static modification such as cysteine 
      all_static_nterm[new_x[0]] = new_x[2]
    elif new_x[1] == "S":   #checks for static modification such as cysteine 
      all_static[new_x[0]] = new_x[2]
    #elif len(new_x) == 4:   #checks for dynamic modification with n-term #This does not work for multiple amino acid for dynamic modification such as STY
    elif new_x[-1] == "n":
      dynamicNterm[new_x[0]] = new_x[2]   
    else:
      countVariableMods+=1   
      delMass = new_x[2]    #looks for all other dynamic modififcations 
      dynDelMass.append(delMass)   #appends delta mass to list
      mod_aa = peptide[int(new_x[0])-1] #since new_x[0] is real position and not index we need to subtract 1 to make index
    
      for val in dyn_AA:
        if mod_aa in val:
          dynAA.append(val) #appends dynamic amino acid to list
    
    all_dyn_mod_combinations= multipleDynModToIndices(peptide, dynAA)     #calls the function characterToIndices to generate combinations for all positions and makes list of list
    
  ##############
  
  #Creates combination of different list using itertolls product from for all combinations
  out_list = []
  combination_list = list(itertools.product(*all_dyn_mod_combinations)) #the list of list to check for all combinations

  for value in combination_list:
    if len(set(value)) == len(value):  #removes repetition of same position in the different combination
      if sorted(value) not in out_list:
        out_list.append(sorted(value))
#      out_list.append(value)
#   print(out_list)
#   print (dynAA)
  modsVariantsAll = []
  massPosList = []
  
  for element in out_list:
    massPosDict = {}
    modsVariants = []
    if len(dynamicNterm) >=1:
      for keyD in dynamicNterm.keys():
        modsVariants.append(keyD+"_V_"+dynamicNterm[keyD]+"_n")
        valAddKey(massPosDict, str(keyD), dynamicNterm[keyD])
#         massPosDict[keyD] = dynamicNterm[keyD]

    if len(all_static_nterm) >=1:
        for keySn in all_static_nterm.keys():
            modsVariants.append(keySn+"_S_"+all_static_nterm[keySn]+"_n")
            valAddKey(massPosDict, str(keySn), all_static_nterm[keySn])

    if len(all_static) >=1:
      for keyS in all_static.keys():
        modsVariants.append(keyS+"_S_"+all_static[keyS])
        valAddKey(massPosDict, str(keyS), all_static[keyS])
#         massPosDict[keyS] = all_static[keyS]
    for i,eachVal in enumerate(element):
      aaMod = peptide[eachVal]
#       aaMod = dynAA[i]
      delM = dynDelMass[i]
      modsVariants.append(str(eachVal+1)+"_V_"+delM+"_"+aaMod)
      valAddKey(massPosDict, str(eachVal+1), delM)
#       massPosDict[str(eachVal+1)] = delM
   #add the mods if they are in same residue and re-create the dictionary
    massPosDict2 = {}
#     print (modsVariants)
#     print (massPosDict)
    for posKey in massPosDict.keys():
      new_val = 0.0
      for eachVal in massPosDict[posKey]:
        new_val += float(eachVal)
      massPosDict2[posKey] = new_val

    massPosList.append(massPosDict2)
    modsVariantsAll.append(",".join(modsVariants))
  return modsVariantsAll, massPosList
  
  


def return_skiprows(file, delimiter, peptide):
  with open(file, "r") as f:
    skiprows = 0
    for line in f:
      if peptide+delimiter in line:
        break
      else:
        skiprows+=1
  return skiprows


# In[10]:



def createOutfile(row, df):
  columns = list(df.columns)
  if "Outfile" in columns:
    outfile_split = row.Outfile.split(".")
    run = outfile_split[-5].split("/")[-1]
    scan = int(outfile_split[-4])
    charge = outfile_split[-2]
    spectrum = run+"."+str(scan)+"."+str(charge)
  else:
    run = row["Run#"]
    scan = row["Scan#"]
    charge = row["z"]
    spectrum = run+"."+str(scan)+"."+str(charge)
#   print (spectrum)
  return spectrum


# In[13]:



def tidy_split(df, column, sep='|', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df


# In[ ]:
#this function extracts dynamic and static modifications using pepxml files and stores them as the dictionary
def getDynStatModsInfoPepXml(pepxml):
    f = open(pepxml,"r") #read the file
    line = f.readline()
    var_AA_mass = {} 
    var_AA_symbol = {} #symbol mass dictionary
    stat_AA_mass = {}
    while "<spectrum_query" not in line: #end reading the file if this is sen
    # look for aminocaid modification keyword so that the modification infomration can be parsed
        if "<aminoacid_modification" in line: #found modification information
            if "symbol=" in line.strip(): #only dynamic modification has symbols as static are fixed

                #check the patter
                pattern = '<aminoacid_modification aminoacid="(\w)" massdiff="(\d+\.\d+)" mass="\d+\.\d+" variable="(\w)" symbol="(.+?)"/>'

                #match for the patter in the line and store them as the tuples there are 4 matches () this is changed to list with list function
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3,4))

                modAA = mods_Description[0] #first element in list is amino acid 
                varMass = float(mods_Description[1]) #second element in list is mass and it is converted to float
                variable = mods_Description[2] #third element in the list is determination of variable of static modification "Y" is variable and "N" is static
                symbol = mods_Description[3] #Symbol. this is used in the dictionary
                valAddKey(var_AA_mass,  varMass, modAA)
    #             valAddKey(var_AA_symbol, symbol, varMass)
                var_AA_symbol[symbol] = varMass #symbol as key and mass as values in the dictionary
                line = f.readline()
            else:
                #this is for static modification so there is no symbol hence we have only three values in the list
                pattern = '<aminoacid_modification aminoacid="(\w)" massdiff="(\d+\.\d+)" mass="\d+\.\d+" variable="(\w)"/>'
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3))
                modAA = mods_Description[0]
                varMass = float(mods_Description[1])
                variable = mods_Description[2]
                stat_AA_mass[modAA] = varMass
    #             valAddKey(stat_AA_mass, modAA, varMass)
                line = f.readline()

        elif "<terminal_modification terminus" in line: #This is for terminal modification such as N-term or C-term
            if "symbol=" in line.strip():
                pattern = '<terminal_modification terminus="(\w)" massdiff="(\d+\.\d+)" mass="\d+\.\d+" variable="(\w)".+symbol="(.+?)"/>'
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3,4))

                modAA = mods_Description[0].lower()
                varMass = float(mods_Description[1])
                variable = mods_Description[2]
                symbol = mods_Description[3]
                valAddKey(var_AA_mass,  varMass, modAA)
    #             valAddKey(var_AA_symbol, symbol, varMass)
                var_AA_symbol[symbol] = varMass
                line = f.readline()
            else:
                pattern = '<terminal_modification terminus="(\w)" massdiff="(\d+\.\d+)" mass="\d+\.\d+" variable="(\w)".+/>'
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3))
                modAA = mods_Description[0].lower()
                varMass = float(mods_Description[1])
                variable = mods_Description[2]
                stat_AA_mass[modAA] = varMass
    #             valAddKey(stat_AA_mass, modAA, varMass)
                line = f.readline()
        else:
            line = f.readline()
    return var_AA_mass,var_AA_symbol, stat_AA_mass



#This fucntion computes modifications with all the static and dynamic represenatation
def computeModifications(row, jump_mod_dict, sta_AA):
    mod_peptide = row.Peptides
    peptideNoFlank = mod_peptide.split(".")[1]
    pepNoFlankList = list(peptideNoFlank)

    plain_peptide_list =[]
    for aa in pepNoFlankList:
        if aa in pyteomics.parser.std_amino_acids:
            plain_peptide_list.append(aa)
    plain_peptide = "".join(plain_peptide_list)
    
    dynAA_count = 0
    for key in jump_mod_dict.keys():
        if key in mod_peptide:
            dynAA_count +=1

    mods = []
    if "n" in peptideNoFlank:
        dynAA_count-=1
        mods.append("1_V_"+str(jump_mod_dict["n"])+"_n")
        peptideNoFlank = peptideNoFlank[1:]
    #find index of other dynamic modification symbols

    mod_symbol_index = []
    for key in jump_mod_dict.keys():
        for aa_index, aa in enumerate(peptideNoFlank):
            if aa == key:
                mod_symbol_index.append(aa_index)
    no_remaining_dyn_mods = len(mod_symbol_index)

    iterate = 1
    for mod_no in range(0,no_remaining_dyn_mods):
        aa = peptideNoFlank[mod_symbol_index[mod_no]-1]
        position = mod_symbol_index[mod_no]-iterate+1 #position is always greater than index
        mod_mass = jump_mod_dict[peptideNoFlank[mod_symbol_index[mod_no]]]

        iterate+=1
        new_mod = str(position)+"_V_"+str(mod_mass)
        mods.append(new_mod)
    #static modification

    for sta_key in sta_AA.keys():
        if sta_key == "n":
            sta_mod = "1_S_"+str(sta_AA["n"])+"_n"
            mods.append(sta_mod)
        else:
            plain_pep_list = list(plain_peptide)
            for index, aa in enumerate(plain_pep_list):
                if aa == sta_key:
                    sta_mod = str(index+1)+"_S_"+str(sta_AA[aa])
                    mods.append(sta_mod)
    modifications = ",".join(mods)              
    return pd.Series([plain_peptide,modifications])


#if the input is ms2File it creates a dataframe with columns = scan, charge, m/z, intensity and so on
def ms2ToDf(ms2File):
    g = open(ms2File,"r")
    line = g.readline()

    scan_list = []
    charge_list = []
    MH_list = []
    precursorIon_list = []
    ms2_mz = []
    ms2_int = []
    while line:
        specMZ = []
        specInt = []
        while "S\t" not in line:
            line = g.readline()
        if "S\t" in line: #spectrum found
            temp_line = line.strip().split("\t")
            scan = temp_line[1]
            scan_list.append(scan)
            precursorIon = temp_line[-1]
            precursorIon_list.append(precursorIon)
            nextline = g.readline()

            temp_line2 = nextline.strip().split("\t")
            charge = int(temp_line2[1])
            charge_list.append(charge)
            MH = temp_line2[-1]
            MH_list.append(MH)
            line = g.readline()
    #         print (line)
        while ("S\t" not in line) and (line != ""):

            temp_line = line.strip().split("\t")
            specMZ.append(float(temp_line[0]))
            specInt.append(float(temp_line[1]))
            line = g.readline()
        ms2_mz.append(specMZ)
        ms2_int.append(specInt)

    dict1 = {"scan":scan_list,"charge":charge_list,"[M+H]+":MH_list,
            "prec_MZ":precursorIon_list,"m/z":ms2_mz,"intensity":ms2_int}

    ms2Df = pd.DataFrame.from_dict(dict1)
    exp = ms2File.split("/")[-1][0:-4]
    ms2Df["spectrum"] = exp+"."+ms2Df.scan.astype("str")+"."+ms2Df.charge.astype("str")
    return ms2Df


#converts mzXML_list to the dataframe using by selecting spectrum list 
def mzmlDF(mzXML_list, spec_list):
  start = time.time()
  super_x =[]
  x1 = pyteomics.mzxml.read(mzXML_list[0],iterative=True)  #reading mzXML file using pyteomics iterative (bool, optional) – Defines whether iterative parsing should be used. It helps reduce memory usage at almost the same parsing speed. Default is True.
  df = pd.DataFrame([x for x in x1])  #dataframe of the mzXML file
  ms2 = df.copy().loc[df.msLevel==2]     #ms2 level scans
  run = mzXML_list[0].split("/")[-1].split(".mzXML")[0]
  ms2["expScan"] = run+"."+ms2.num.astype("str")
  ms2F = ms2.copy().loc[ms2.expScan.isin(spec_list)]
  ms2F.rename(columns={"num":"scan","m/z array":"m/z","intensity array":"intensity"}, inplace=True)
  reqCols = ["expScan","scan","m/z","intensity"]
  ms2F2 = ms2F[reqCols]
  super_x.append(ms2F2)
  n=1
  print ("The mzxml file are now concatenated")
  print ("File ",n," is concatenated")
  for x in range(1,len(mzXML_list)):
    x2 = pyteomics.mzxml.read(mzXML_list[x],iterative=True)  #reading mzXML file using pyteomics
    df2 = pd.DataFrame([x for x in x2])  #dataframe of the mzXML file
    ms2_1 = df2.copy().loc[df2.msLevel==2]     #ms2 level scans
    run2 = mzXML_list[x].split("/")[-1].split(".mzXML")[0]
    ms2_1["expScan"] = run2+"."+ms2_1.num.astype("str")
    ms2F_1 = ms2_1.copy().loc[ms2_1.expScan.isin(spec_list)]
    ms2F_1.rename(columns={"num":"scan","m/z array":"m/z","intensity array":"intensity"}, inplace=True)
    ms2F_2 = ms2F_1[reqCols]
    #frames = [ms2F2, ms2F_2]
    #ms2F2 = pd.concat(frames)
    super_x.append(ms2F_2)
    n+=1
    print ("Files ",n," is concatenated")
  print ("concatenating ",n,"  files. This will take time depending on number of mzXML file\n")
  finalDF = pd.concat(super_x, axis=0)
  end = time.time()
  print ("total time required to concatenate all mzXML files = ", end-start,"seconds\n")
  return finalDF


#converst ms2 files from folders and subsets using spectrum list
def msToDF(allms2Folder, spec_list):
  df = ms2ToDf(allms2Folder[0])
  df_pre3 = df.copy().loc[df.spectrum.isin(spec_list)]
  run = allms2Folder[0].split("/")[-1].split(".ms2")[0]
  df_pre3["expScan"] = run+"."+df_pre3.scan.astype("str")
    
  n=1
  print ("The ms2 file are now concatenated")
  print ("File ",n," Fraction name = ", allms2Folder[0], " is concatenated")

  
    
  for x in range(1,len(allms2Folder)):
    df_2 = ms2ToDf(allms2Folder[x])
    df_pre33 = df_2.copy().loc[df_2.spectrum.isin(spec_list)]
    run2 = allms2Folder[0].split("/")[-1].split(".ms2")[0]
    df_pre33["expScan"] = run2+"."+df_pre33.scan.astype("str")
    frames = [df_pre3, df_pre33]
    df_pre3 = pd.concat(frames)
    n+=1
    print ("File ",n," Fraction name = ", allms2Folder[x], " is concatenated")


  return df_pre3


#this program creates the symbols to filter and make sure that wildcard characters escape
def symbolsToFilter(jump_mod_dict, delta_mass_to_localize): #delta_mass_to_localize = parameter that users defines to localize the delta mass for example for pho search STY = 79.966331 
    wildcharacters = [".","^","$","*","+","?"] # we need to escape these symbols if we want to filter the peptides with these symbols
    symbol_to_localize = []
    #the dynamic modifications will be selected based on the dictionary
    for key,values in jump_mod_dict.items():
        for mass in delta_mass_to_localize:
            if values == float(mass): #float the values as the parameter will be string
                symbol_to_localize.append(key)
    # make a string to filter df_pho (ID.txt) matrix using symbols. We have to escape wild characters. THis is also important to decrease the time for JUMP -localization    
    filter_symbols = []
    for sym in symbol_to_localize:
        if sym in wildcharacters:
            filter_symbols.append("\\"+sym)
        else:
            filter_symbols.append(sym)
    final_symbols_filter = "|".join(filter_symbols) #concatenates symbols as one string using | symbols which is used as "OR" in python
    return final_symbols_filter

def ppmCalc(a, b, tol=10):
  a = float(a)
  b = float(b)
#   massError = abs(a-b)*1e6/a
  massError = (b-a)*1e6/a  #calculates ppm error of the a(theoretical) from b(observed)
  return float(massError)


#https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html
# pmf(k, M, n, N, loc=0)  ---- this is what the website describes
#pmf(k, N, K, n, loc=0)  ---- this is what we have used https://en.wikipedia.org/wiki/Hypergeometric_distribution
def hyperTestPval(k, N, K, n):
  return hypergeom.pmf(k, N, K, n)

#calcualation of scores
def hypergeometricTest(exp_mz_list, intensity_list, theor_fragments_df, tol=10, ion_types = ["b","y"], ionLoss = ["NH3","H2O"]):
  cols = theor_fragments_df.columns
  ionsCols = []
  for ions in ion_types:
    for ionC in cols:
      if ions+"+" in ionC:
        ionsCols.append(ionC)
      for ionLosses in ionLoss:
        if ions+"-"+ionLosses in ionC:
          ionsCols.append(ionC) 
    

  checkDF = theor_fragments_df[ionsCols]
  newDF = checkDF.replace("---",np.nan)
  arr = newDF.to_numpy().flatten()
  theor_ions_list = list(arr[~np.isnan(arr)])
  
  mmax = np.max(exp_mz_list)
  mmin = np.min(exp_mz_list)
  
#   print ("mmax = ",mmax)
#   print ("mmin = ",mmin)

  Nfloat = (mmax-mmin)/((tol*(mmin+mmax)/1000000)*2)
    
#   print ("Nfloat = ",Nfloat)

  N = int(Nfloat.round())
  K = len(exp_mz_list) 
  n = len(theor_ions_list)


  true_ions_ppm_list = []
  true_ion_int_list = []
  true_ions_list = []
  true_theretical_ions_list = []
  for i,exp_ion in enumerate(exp_mz_list):
    for theor_ion in theor_ions_list:
      ppm = ppmCalc(theor_ion,exp_ion, tol=tol)
      
      if (ppm <= float(tol)) and (ppm >= (-1*float(tol))) :
        true_ions_ppm_list.append(ppm)
        true_ions_list.append(exp_ion)
        true_ion_int_list.append(intensity_list[i])
        true_theretical_ions_list.append(theor_ion)
#         print ("True ions, Theor ion and PPM,", exp_ion, theor_ion, ppm)
  k = len(true_ions_list)
    
  
  pval = hyperTestPval(k, N, K, n)
  if pval == 0:
    pval = 1e-40
  return pval,true_ions_list, true_ion_int_list,true_ions_ppm_list

#calculates psm level scores
def psmScore1(df, spec):
    #extracts one spectrum information on modification, peptide, pvalue, and sequence probabilities
    psmAll = df.copy().loc[df.spectrum == spec]
    
    #conts the number of dynamic modifications in the peptide by splitting "," as "," separates different modifications
    psmAll["countMods"] = psmAll.modifications.apply(lambda x: len(x.split(",")))
    
    #makes a new row for each modification for example if there are two modifications, two new rows will be created with all same information except the modification 
    psmAll = tidy_split(psmAll,"modifications",",",keep=False)
    
    #group the dataframe with "spectrum","plain_peptide","modifications","countMods" and sums the sequence probability. 
    #If the same modification is observed many times that will be summed and kept as unique probability
    psmAll_grp = psmAll.groupby(["spectrum","plain_peptide","modifications","countMods"])["SequenceProbablity"].sum().reset_index()
    
    #gives the required number of modifications in the peptide. This would be same for all rows so we just take one element
    modRequired = list(set(psmAll_grp.countMods))[0]
    
    #modRequired is used to pick up highest scoring modification of that number for example if modRequired = 2, highest 2 probab score are picked up. Low score is the second highest in this case, if 3 mods then it is 3rd highest
    #new dataframe with rows == number of modification (modRequired) is selected. 
    psmAll_grpsort = psmAll_grp.nlargest(modRequired,['SequenceProbablity'])
    
    #Next we sort the matrix by the probability
    psmAll_grpsort_select = psmAll_grpsort.sort_values(by=["SequenceProbablity"], ascending=False)
    
    #Get the scores of probability to 2 decimal places
    psmAll_grpsort_select["SequenceProbablity"] = psmAll_grpsort_select["SequenceProbablity"].apply(lambda x: '{0:.2f}'.format(x))
    
    #two rows of the peptide are concatenated into one row with the scores separated by ",". This is achieved by grouping with "spectrum","plain_peptide","countMods" columns and joining remaining columns by "," applying tolist
    psmAll_Final = psmAll_grpsort_select.groupby(["spectrum","plain_peptide","countMods"]).agg(lambda x: ",".join(x.tolist())).reset_index()  
    
    #same number of rows as the input spectrum
    return psmAll_Final

#This function helps to map the peptide to the protein level and get the sites information
def proteinLevelSites(row):
    prot_level_mod = []
    
    #First modificaition information is retrieved because that has sites for modified amino acids site in the peptide
    jump_l_mods = row.modifications
    
    #Next the position of peptide mapped is retrieved from ID.txt
    pos_jumpf = row.pos
    
    #position start is extracted using regex findall function with 0 index being the start (beginning number)
    pos_start = int(re.findall("\d+", pos_jumpf)[0])
    
    #the modification are splited by comma
    commaSplit = jump_l_mods.split(",")

    #this loops within different modication that are in commaSplit and the new protein positions replaces old number and the modification nomenclature is recreated
    for mod in commaSplit:
        #splits by "_"
        modSplit = mod.split("_")
        #the new value is defined as the previous number -1 + new protein mapped number
        newVal = str(int(modSplit[0])-1+pos_start)+"_"+"_".join(modSplit[1:])
        #list containing protein level modification is created
        prot_level_mod.append(newVal)
    
    return ",".join(prot_level_mod) #the protein modfiications is joined by "," to make the similar modificaiotn as peptide


def dfColsCommaToList(df, column, delimiter=","):
    finalList = []
    colToList = list(df[column])
    for val in colToList:
        valSplit = val.split(delimiter)
        finalList+=valSplit
    return finalList


#each row of teh dataframe protAll is applied with this function
def protLevelScoreMap(row, dict1): #dict1 is the site dictionary with key as sites (protein) and values as the maximum protein level probability score
    #print ("modification dictionary = \n", dict1)
    total_mods = int(row.countMods) #checks the total number of modifications
    #print ("total modification = ", total_mods)
    protMatchedScores = [] ##list of protein level probability scores for modification sites protein level for particular protein row
    protMatchedMods = [] #list of protein level modifications or sites for particular protein row
    protSite = row.prot_level_site #these are protein level sites
    #print ("protein sites = ", protSite)
    protSiteSplit = protSite.split(",") #prot level sites splitted
    for val in protSiteSplit:
        protMatchedMods.append(val) #protein level mods or sites
        protMatchedScores.append(dict1[val]) #maximum score of sequence probability
    
#     print ("protMatchedScores\n",protMatchedScores)
#     print ("protMatchedMods\n",protMatchedMods)
    score = np.array(protMatchedScores) #list of scores converted to numpy array
    sites = np.array(protMatchedMods) #list of sites converted to numpy array
    
#     print ("score = ", score)
#     print ("sites = ", sites)
    
    inds = score.argsort() #index of sorted scores [10,2,30,4] will give array([1, 3, 0, 2]) which means 0 index is highest score and is at third position and so on
#     print ("Inds ", inds)
    sortedSites = sites[inds] #corresponding sorting of sites by ascending order
#     print ("sorted sites ",sortedSites)
    
    sortedScore = score[inds] #this will sort scores based on the number that it returns as inds. Above example will be array([ 2,  4, 10, 30]) sorted by ascending order 
    sortedScore_rev = sortedScore[::-1] # above result will be  array([30, 10,  4,  2]) .. reverse sorting of ascending order
#     print ("sortedScore_rev ",sortedScore_rev)
    
    sortedSites_rev = sortedSites[::-1] #reverse sorting of sites too
#     print ("sortedSites_rev ",sortedSites_rev)
    
    lowestScore = sortedScore_rev[0:total_mods][-1] #if total mods = 2, the value will be 10. That means if 2 values are allowed second highest is 10 
#     print ("lowestScore ",lowestScore)
    
    new_scores = [x for x in sortedScore_rev if x >= lowestScore] #list comprehesion that determines all the scores above lowestScore in this case [30, 10] based on number of mods. 
    #Since there are 2 mods we can only accept 2 scores and highest 2 scores are [30, 10]
    
    #new sites are determined by the number of scores selected for example 2 in this case. So, top 2 sites from sortedSites_rev are selected
    new_sites = sortedSites_rev[0:len(new_scores)]
    
    #The selected new sites are joined by ","
    NSites  = ",".join(new_sites)
    #Selcted new scores are joined by ","
    NScores = ",".join(str(n) for n in new_scores)
#     print ("NSites")
#     print ("NScores")
    return pd.Series([NSites, NScores])


def publicationTableSiteScoreColumns(row, delta_mass_to_localize): #delta_mass_to_localize = use user defined masses to localize and give sites information and scores
    prot_site = row.prot_level_site.split(",") #protein mod site is separated by ","
    site_scores = row.Protein_level_site_scores.split(",") #Protein localization scores is separated by ","
    keep_sites = [] #updates the sites to keep
    keep_scores = [] #updates the scores to keep
    for index, val in enumerate(prot_site): 
        val_split = val.split("_") #splits the mods by "_"
        if val_split[1] == "V": #check these sites
            if val_split[2] in delta_mass_to_localize: #keep these sites
                site = val_split[-1]+val_split[0] #creates the site informaation for example S170 = serine is modified at 170 position
                keep_sites.append(site) #append those sites to the list
                keep_scores.append(site_scores[index]) #add the corresponding score to the score list
    #Selcted new sites are joined by ","
    NSites  = ",".join(keep_sites)
    #Selcted new scores are joined by ","
    NScores  = ",".join(keep_scores)
    
    return pd.Series([NSites, NScores])     
    


#each protein accession are checked for this section "locProtein = list(set(psmLocDF.Protein))" gives the unique list of protein accession
def proteinScore(df, protein, delta_mass_to_localize): #protein is one accession from locProtein list
    protAll = df.copy().loc[df.Protein == protein] #creates a dataframe that has only one protein (subset). The idea is to just look at one protein
#     print (protAll)
    allprobability = list(protAll.SequenceProbablity) #Gives the list of sequence probability for the peptide that maps to this particular protein
    
    #this dfColsCommaToList converts entire column to a list. If there is "," separated values as in some sequence probability
    #it opens that up and make a new linear list with all scores
    
    seqProbab = dfColsCommaToList(protAll, "SequenceProbablity", delimiter=",")
    
    #does same for the protein level sites too ex.418_S_229.162931,419_V_79.966331_T (separated by ",")
    modSites = dfColsCommaToList(protAll, "prot_level_site", delimiter=",")
    
    #modSitesProbDict is a dictionary that will be updated with sites as key and if same sites present the seqProbab variable will be added as a list of values
    
    modSitesProbDict = {}

    for i, sites in enumerate(modSites):
        if sites not in modSitesProbDict:
            modSitesProbDict[sites] = [float(seqProbab[i])] #sequence prob values added as the list
        else:
            modSitesProbDict[sites].append(float(seqProbab[i])) #if new key (sites) is seen a new list is generated with values 

    maxProbDict = {} #this is another dictionary that store maximum score from the modSitesProbDict
    for key in modSitesProbDict.keys(): #key is the site (ex. 418_S_229.162931)
        maxProbDict[key] = max(modSitesProbDict[key]) #maximum score for each sites is retained as values
        
        #in this step the new function protLevelScoreMap is applied. see this fucntion for more details but it takes the max prob dictionary
#         print (maxProbDict)
    protAll[["Protein_final_mods","Protein_level_site_scores"]]=protAll.apply(protLevelScoreMap, dict1=maxProbDict, axis=1) #This will return protein level sites and prtein level scores
    
    protAll[["Mod sites","Lscore"]]=protAll.apply(publicationTableSiteScoreColumns, delta_mass_to_localize=delta_mass_to_localize, axis=1) #This will return protein level sites and prtein level scores to keep
    
    return protAll


def useValuesGetKey(dict1, string):
    result = ""
    for key,values in dict1.items():
        if values == string:
            result = key
    return result

def jump_l_peptide(row, jump_mod_dict):
    #Add modifications symbols in the peptide sequence with those modifications that have highest scores.
    mods = row.modifications  #stores the new modification with best scores
    plain_peptide = row.plain_peptide #plain peptide for counting the position of mods
    
    prev_mod_peptide = row.Peptides_beforeJUMP_l  #previous peptide to obtain flanking aa information
    prev_mod_peptide_split = prev_mod_peptide.split(".") #split peptide for flanking aa
    
    left_flank =prev_mod_peptide_split[0] #left flank
    right_flank = prev_mod_peptide_split[-1] #right flank
    
    modsSplit = mods.split(",") #this splits each modification so that each mods symbols can be added back to the peptide             
    new_mod_peptide = [left_flank+"."] #the new peptide will be in the forma of list. The list contains left flank amino acid followed by "."
    #this new_mod_peptide will now be updated with new amino acids and symbols and at the end converted into one peptide
    for index_aa, aa in enumerate(plain_peptide): #this looks at each index poistion and amino acid with enumerate function
        for x in modsSplit: #considers each modifications
            new_x = x.split("_") #splits the modifications
    #         print (new_x)
            if str(index_aa+1) == new_x[0]: #index starts from 0 so we add +1 to make it same as the position of modification
                if (new_x[1] == "V") and (new_x[-1] == "n"):   #Found a Nterm dynamic modification
                    newSymbol = useValuesGetKey(jump_mod_dict, float(new_x[2])) #symbol is retrieved from jump_mod_dict key if the value matches the floating number in modificaiotn
    #                 print ("n"+aa+newSymbol)
                    new_mod_peptide.append("n"+aa+newSymbol) #the new modified peptide is appended

                elif (new_x[1] == "V"):   #checks for dynamic modification besides Nterm
                    newSymbol = useValuesGetKey(jump_mod_dict, float(new_x[2])) #another symbol of modification that is not Nterm
                    new_mod_peptide.append(aa+newSymbol) #appends the amino acid with new symbol

                else: #there are no more dynamic modification. So just append amino acid 
                    new_mod_peptide.append(aa) 

                    break #very important, if the criteria matches, we do not search other modificaition, the loop breaks
        else: #interesting for else concept. This is if there is no matching numbers in modifications and we put the aa as it is
            new_mod_peptide.append(aa)
    new_mod_peptide.append("."+right_flank)
    return "".join(new_mod_peptide)


#This is use to grab the non structured lines above the dataframe
#file = ID.txt, delimiter = ;, peptide = Peptide column, outFile = new ID.txt
def grabHeader(file, delimiter, peptide,outFile):
  
  with open(file, "r") as f, open(outFile, "w") as g:
    skiprows = 0
    for line in f:
      if peptide+delimiter in line:
        break
      else:
        g.write(line.strip())
        skiprows+=1
        if skiprows >= 1:
            g.write("\n")

#delta_mass_to_localize -- count the dynamic modifications used by the users
#file = ew publication tables for jump localization, delimiter = "\t", peptide = Peptides column
#This function takes the jump_mod_dict {'*': 79.966331, '#': 15.994915} of aa and use the symbols and delta mass to find the modification to count the modified peptides and sites
#This function also uses jump_modAA_dict that has {79.966331: ['S', 'T', 'Y'], 15.994915: ['M']} delta mass and list of aa
def writePublicationTableHeaderPept(file, peptide_df, jump_mod_dict, jump_modAA_dict): #peptide_df after jump_l 
  with open(file, "w") as g: #new txt (peptides) files
    g.write("Unique modified peptides identified by mass spectrometry\n")
    g.write("n = "+str(peptide_df.shape[0])+" modified peptides\n")
    number_peptides = peptide_df.shape[0] #this gives the number of rows of peptides. All of these peptides are modified peptides
    mod_peptide_aa = [] #this list updates the amino acids that are modified and the total number of modified sites
    all_mod_sites = 0 #gives the sum of all modified sites for example if S= 10, T= 3, Y=1 sum = 14
    for key,value in jump_mod_dict.items(): #iterate over key and values of jump_mod_dict
        aa_list = jump_modAA_dict[jump_mod_dict[key]] #get the amino acid list from jump_modAA_dict using the key as the value of jump_mod_dict[key]
        symbol=symbolsToFilter(jump_mod_dict, [value]) #gets the symbol to search and also make sure if a symbol needs to be escaped as being the wildcharacter for regex
        for aa in aa_list: #checks each amino acid in dynamic modification to count the number of sites
            newDF= peptide_df.copy().loc[peptide_df.Peptides.str.contains(aa+symbol)] #creates a dataframe with that amino acid and modification
            newDF["Count"] = newDF.Peptides.str.count(aa+symbol) #adds a new column to the dataframe that counts the modifications in the peptide (if there are two modification count = 2)
#             print(aa,np.sum(newDF.Count))
            all_mod_sites+=np.sum(newDF.Count) #sums the total count of modification across the dataframe that sees the amino acid modificaiotn
            mod_peptide_aa.append(str(np.sum(newDF.Count))+" "+aa)
#     print ("n = "+str(all_mod_sites)+" modified sites ("+",".join(mod_peptide_aa)+")")
    g.write("n = "+str(all_mod_sites)+" modified sites ("+",".join(mod_peptide_aa)+")\n")



#delta_mass_to_localize -- count the dynamic modifications used by the users
#file = ew publication tables for jump localization, delimiter = "\t", peptide = Peptides column
#This function takes the jump_mod_dict {'*': 79.966331, '#': 15.994915} of aa and use the symbols and delta mass to find the modification to count the modified peptides and sites
#This function also uses jump_modAA_dict that has {79.966331: ['S', 'T', 'Y'], 15.994915: ['M']} delta mass and list of aa
def writePublicationTableHeaderProt(file, protein_df, jump_mod_dict, jump_modAA_dict): #peptide_df after jump_l 
  with open(file, "w") as g: #new txt (peptides) files
    g.write("Unique proteins identified by mass spectrometry (n = "+str(protein_df.shape[0])+")\n")

        
