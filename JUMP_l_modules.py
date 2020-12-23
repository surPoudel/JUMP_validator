#!/usr/bin/env python
# coding: utf-8

# In[27]:


import pandas as pd
import os, sys, glob
import re
from pyteomics import mgf, pepxml, mass
import os
from urllib.request import urlretrieve
import pylab
from pyteomics import mass
import itertools
from operator import attrgetter
import numpy as np

# In[2]:


def valAddKey(dict1, key, val):
  if key not in dict1.keys():
    dict1[key] = [val]
  else:
    dict1[key].append(val)
  return dict1


# In[6]:


def ionSeriesIonLossSpeRes(peptide,massPosDict, maxcharge=1,useMod ="Yes"):
  if useMod != "Yes":
    massPosDict = {0:0.0}
  h2o = mass.calculate_mass(formula='H2O')
  co = mass.calculate_mass(formula='CO')
  nh3 = mass.calculate_mass(formula='NH3')
  xmassMore = co-mass.calculate_mass(formula='H2')
  proton = mass.calculate_mass(formula='H+')
  
  hpo3 = mass.calculate_mass(formula='HPO3')
  h3po4 = mass.calculate_mass(formula='H3PO4')

  all_ions_dict = {}

  ionloss = {"":0,"-H2O":h2o,"-HPO3":hpo3,"-H3PO4":h3po4,"-NH3":nh3}
  
  
  for i in range(1, len(peptide)+1):
    addModMass = 0.0
    for massKey in massPosDict.keys():
      if int(massKey) <= i:
        addModMass += float(massPosDict[massKey])
        
    valAddKey(all_ions_dict,"Seq",peptide[i-1])
#     print (peptide[0:i-1])
    for losses in ionloss.keys():
      if losses == "-H2O":
        fate = checkAA(peptide[0:i],["S","T","E","D"])
        for charge in range(1, maxcharge+1):
          if fate == "True":
            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
          else:
            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
      if losses == "-NH3":
        fate = checkAA(peptide[0:i],["R","K","Q","N"])
        for charge in range(1, maxcharge+1):
          if fate == "True":
            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
            valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
          else:
            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
            valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
                      
      if (losses == "-H3PO4") or (losses == "-HPO3"):
        for charge in range(1, maxcharge+1):
          if "79.96" in str(massPosDict.values()):
            fate = checkAA(peptide[0:i],["S","T","Y"])
            if fate == "True":
              valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
              valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
            else:
              valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
              valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
                      
      if losses == "":
        for charge in range(1, maxcharge+1):
          valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
          valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
          valAddKey(all_ions_dict,"c"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge)
         

  for i in range(0, len(peptide)):
    
    addModMass = 0.0
    for massKey in massPosDict.keys():
      if int(massKey) > i:
        addModMass += float(massPosDict[massKey])
    for losses in ionloss.keys():  
      if losses == "-H2O":
        fate = checkAA(peptide[i:len(peptide)],["S","T","E","D"])
        for charge in range(1, maxcharge+1):
          if fate == "True":
            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
          else:
            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
      
      if losses == "-NH3":
        fate = checkAA(peptide[i:len(peptide)],["R","K","Q","N"])
        for charge in range(1, maxcharge+1):
          if fate == "True":
            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
          else:
            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
      
      
      if (losses == "-H3PO4") or (losses == "-HPO3"):
        fate = checkAA(peptide[i:len(peptide)],["S","T","Y"])
        for charge in range(1, maxcharge+1):
          if "79.96" in str(massPosDict.values()):
            if fate == "True":
              valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
            else:
              valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
      if losses == "":
        for charge in range(1, maxcharge+1):
          valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
          valAddKey(all_ions_dict,"x"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)+xmassMore-ionloss[losses]+addModMass)/charge)
          valAddKey(all_ions_dict,"z"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge)
      
    
                

  df = pd.DataFrame(all_ions_dict)
  return df


# In[5]:


def checkAA(pepSeq, check_list):
  update_val = "False"
  for aa in list(pepSeq):
    if aa in check_list:
      update_val = "True"
  return update_val


# In[3]:


def ionSeriesIonLossAllResidue(peptide,massPosDict, maxcharge=1,useMod ="Yes"):
  if useMod != "Yes":
    massPosDict = {0:0.0}
  h2o = mass.calculate_mass(formula='H2O')
  co = mass.calculate_mass(formula='CO')
  nh3 = mass.calculate_mass(formula='NH3')
  xmassMore = co-mass.calculate_mass(formula='H2')
  proton = mass.calculate_mass(formula='H+')
  
  hpo3 = mass.calculate_mass(formula='HPO3')
  h3po4 = mass.calculate_mass(formula='H3PO4')

  all_ions_dict = {}

  ionloss = {"":0,"-H2O":h2o,"-HPO3":hpo3,"-H3PO4":h3po4,"-NH3":nh3}
  
  
  for i in range(1, len(peptide)+1):
    addModMass = 0.0
    for massKey in massPosDict.keys():
      if int(massKey) <= i:
        addModMass += float(massPosDict[massKey])
    valAddKey(all_ions_dict,"Seq",peptide[i-1])
    for losses in ionloss.keys():
      for charge in range(1, maxcharge+1):
      
        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
        valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
        valAddKey(all_ions_dict,"c"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge)
      


  for i in range(0, len(peptide)):
    addModMass = 0.0
    for massKey in massPosDict.keys():
      if int(massKey) > i:
        addModMass += float(massPosDict[massKey])
    for losses in ionloss.keys():  
      for charge in range(1, maxcharge+1):
        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
        valAddKey(all_ions_dict,"x"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)+xmassMore-ionloss[losses]+addModMass)/charge)
        valAddKey(all_ions_dict,"z"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge)

  df = pd.DataFrame(all_ions_dict)
  return df


# In[6]:


def characterToIndices(string, char):
  string_list = list(string)
  indices = [i for i, x in enumerate(string_list) if x == char]
  return indices


# In[2]:


def multipleDynModToIndices(peptide, dynAA):
  if peptide.startswith("K"):
    peptide = peptide
  else:
    peptide = "J"+peptide[1:]
  finalIndices = []
  for val in dynAA:
    listAA = list(val)
    allIndices = []
    for aa in listAA:
      charList = characterToIndices(peptide, aa)
      allIndices+=charList
    finalIndices.append(allIndices)
  return finalIndices


# In[1]:


def multipleCharacterToIndices(peptide, dyn_AA):
  
  listAA = list(dyn_AA)
  allIndices = []
  for aa in listAA:
    charList = characterToIndices(peptide, aa)
    allIndices+=charList
  return allIndices


# In[9]:


def statModDynModAllVariants(peptide, mods,dyn_AA): #mods is modification column from comet, peptide is plain peptide
  dynAA = []          #list that updates dynamic amino acids in the sequential order
  dynDelMass = []     #list that updates delta masses of dynamic modification  
  modsSplit = mods.split(",")              
  all_static = {}  #all static modification or n-term if dynamic 
  dynamicNterm = {}   #stores dynamic modification Nterm
  all_dyn_mod_combinations = []   #list that updates all dynamic modifications combinations

  countVariableMods = 0  #does not include n-term dynamic here
  for x in modsSplit:
    new_mod = ""
    new_x = x.split("_")
    changedMods = ""
    if new_x[1] == "S":   #checks for static modification such as cysteine 
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
  
  


# In[8]:


def makeMSViewer(row):
  new_mods_list = []
  mods = row.modifications
  if mods != "-":
    modsSplit = mods.split(",")
    for x in modsSplit:
      new_mod = ""
      new_x = x.split("_")
      if (len(new_x) == 4) and (new_x[-1] == "n"):
        new_mod = new_x[2]+"@N-term"
      else:
        new_mod = new_x[2]+"@"+new_x[0]
      new_mods_list.append(new_mod)
  else:
    new_mods_list = ["-"]
  return ";".join(new_mods_list)


# In[9]:


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


def correctNterm(row):
  pep = row.Peptides
  if ("n#" in pep) or ("n$" in pep):
    pep1 = pep[0:3]
    pep2 = pep[3:]
    final_pep = pep1[::-1][0:-1]+pep2
  else:
    final_pep = pep
  return final_pep

def makeCometDf(allPepToTsvFolder, spec_list, exp):
  comet = allPepToTsvFolder+"/"+exp+"/"+exp+".1.txt"
#   comet = glob.glob(allPepToTsvFolder+"/*.1.txt")
  skiprows=return_skiprows(comet, "\t", "modified_peptide")
  df = pd.read_csv(comet, delimiter = "\t", skiprows=skiprows)
  head = list(df.columns)
  df_pre = df[head[0:-1]]
  df_pre2 = df_pre.reset_index()
  df_pre2.columns = head
  df_pre2["Run"] = comet.split("/")[-1].split(".")[0]
  df_pre2["spectrum"]=df_pre2.Run+"."+df_pre2.scan.astype("str")+"."+df_pre2.charge.astype("str")
  df_pre3 = df_pre2.loc[df_pre2.spectrum.isin(spec_list)]
  dfRank1 = df_pre3.loc[df_pre3["num"] == 1]
  return dfRank1
# In[17]:


def merge_comet_searches(allPepToTsvFolder, spec_list):
  comet = glob.glob(allPepToTsvFolder+"/*/*.1.txt")
#   comet = glob.glob(allPepToTsvFolder+"/*.1.txt")
  skiprows=return_skiprows(comet[0], "\t", "modified_peptide")
  df = pd.read_csv(comet[0], delimiter = "\t", skiprows=skiprows)
  head = list(df.columns)
  df_pre = df[head[0:-1]]
  df_pre2 = df_pre.reset_index()
  df_pre2.columns = head
  df_pre2["Run"] = comet[0].split("/")[-1].split(".")[0]
  df_pre2["spectrum"]=df_pre2.Run+"."+df_pre2.scan.astype("str")+"."+df_pre2.charge.astype("str")
  df_pre3 = df_pre2.loc[df_pre2.spectrum.isin(spec_list)]
  
  n=1
  print ("The comet txt files are now concatenated")
  print ("File ",n," Fraction name = ", comet[0], " is concatenated")

  
    
  for x in range(1,len(comet)):
    df_2 = pd.read_csv(comet[x], delimiter = "\t", skiprows=skiprows)
    df_pre_2 = df_2[head[0:-1]]
    df_pre_22 = df_pre_2.reset_index()
    df_pre_22.columns = head

    df_pre_22["Run"] = comet[x].split("/")[-1].split(".")[0]
    df_pre_22["spectrum"]=df_pre_22.Run+"."+df_pre_22.scan.astype("str")+"."+df_pre_22.charge.astype("str")
    df_pre33 = df_pre_22.loc[df_pre_22.spectrum.isin(spec_list)]
    frames = [df_pre3, df_pre33]
    df_pre3 = pd.concat(frames)
    n+=1
    print ("File ",n," Fraction name = ", comet[x], " is concatenated")

  
  dfRank1 = df_pre3.loc[df_pre3["num"] == 1]
  return dfRank1


# In[14]:


def merge_comet_searchesMapJUMPf(allPepToTsvFolder, dyn_mod, spec_list):
  comet = glob.glob(allPepToTsvFolder+"/*/*.1.txt")
#   comet = glob.glob(allPepToTsvFolder+"/*.1.txt")
  skiprows=return_skiprows(comet[0], "\t", "modified_peptide")
  df = pd.read_csv(comet[0], delimiter = "\t", skiprows=skiprows)
  head = list(df.columns)
  df_pre = df[head[0:-1]]
  df_pre2 = df_pre.reset_index()
  df_pre2.columns = head
  df_pre2["Run"] = comet[0].split("/")[-1].split(".")[0]
  df_pre2["spectrum"]=df_pre2.Run+"."+df_pre2.scan.astype("str")+"."+df_pre2.charge.astype("str")
  df_pre3 = df_pre2.loc[df_pre2.spectrum.isin(spec_list)]
  
  n=1
  print ("The comet txt files are now concatenated")
  print ("File ",n," Fraction name = ", comet[0], " is concatenated")

  
    
  for x in range(1,len(comet)):
    df_2 = pd.read_csv(comet[x], delimiter = "\t", skiprows=skiprows)
    df_pre_2 = df_2[head[0:-1]]
    df_pre_22 = df_pre_2.reset_index()
    df_pre_22.columns = head

    df_pre_22["Run"] = comet[x].split("/")[-1].split(".")[0]
    df_pre_22["spectrum"]=df_pre_22.Run+"."+df_pre_22.scan.astype("str")+"."+df_pre_22.charge.astype("str")
    df_pre33 = df_pre_22.loc[df_pre_22.spectrum.isin(spec_list)]
    frames = [df_pre3, df_pre33]
    df_pre3 = pd.concat(frames)
    n+=1
    print ("File ",n," Fraction name = ", comet[x], " is concatenated")

  


# In[12]:


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


def from_cometPep_cleanPep(df, peptide, cometDF, filename):
  columns = list(df.columns)[1:]  #takes peptide column off
#   df_1 =df.loc[df[peptide].str.contains("\*|\$|%")]
  df_1 = df
  print ("length of the dataframe =", df.shape[0]) 

  cometMergedJump = pd.merge(df_1,cometDF,how="inner",on="spectrum")
  cometMergedJump["New_MODS"] = cometMergedJump.apply(makeMSViewer, axis=1)
  msviewer_cols = ['plain_peptide','New_MODS','modified_peptide','spectrum','charge','PeptideNtermChanged','Protein','scan', 'exp_neutral_mass', 'calc_neutral_mass', 'e-value', 'xcorr',
       'delta_cn', 'sp_score', 'ions_matched', 'ions_total']

  filename1 = filename+".txt"
  filename2 = filename+".xlsx"
  filename3 = filename+"_MSviewerInput.txt"

  cometMergedJump.to_csv(filename1, sep="\t", index=None)
  cometMergedJump.to_excel(filename2, index=None)  
  
  cometMergedJump2 = cometMergedJump[msviewer_cols]
  cometMergedJump2.to_csv(filename3, sep="\t", index=None)


# In[11]:


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


def calibratedMass(mz, massError):
  calibMass = mz + ((massError*mz)/1e6)
  return calibMass

def massCorrectionFunction(exp_list, massError=12):
  massCorrectedList = []
  for mz in exp_list:
    calibMass = calibratedMass(mz, massError)
    massCorrectedList.append(calibMass)
  return massCorrectedList

def ppmCalc(a, b, tol=10):
  a = float(a)
  b = float(b)
#   massError = abs(a-b)*1e6/a
  massError = (b-a)*1e6/a  #calculates ppm error of the a(theoretical) from b(observed)
  return float(massError)

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
