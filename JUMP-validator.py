#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from docx import Document
from docx.shared import Inches
import numpy as np
from ast import literal_eval
import configparser
config = configparser.ConfigParser()
# from ipynb.fs.full.JUMP_l_modules import *
from JUMP_l_modules import *
import matplotlib.pyplot as plt
import seaborn as sns
import os, path
from os.path import dirname
from pandas.plotting import table
import dataframe_image as dfi
from docx.shared import Pt
import pyteomics as pyteomics
from pyteomics import mzxml, auxiliary, pepxml


# In[2]:


params_file = sys.argv[1]
# params_file = "../JUMP_localizationProgram/jump_l.params"
#params_file = "../parameterFile/map_comet_jump_fJUMP.params"
config.read(params_file)


# In[3]:


'''
#this is the input excel file for the each peptide that are to be manually validated. Headers = Exp	scan	charge (Exp = experiment/fraction name, scan = scan number, charge = z)
peptideFile = /Users/spoudel1/Desktop/PTM_study/ChaoPeng/Check_pho_peptides_Chao.xlsx

#this is the input ID.txt file. It can also be input from jump_l program ID.lscore file
jump_f_id = /Users/spoudel1/Desktop/PTM_study/ChaoPeng/solubleTauPTMs/phosphorylation/ID.lscore

#this option says whether the output is from jump_l program or not. Type "yes" if JUMP_l was used for pho site localization else type "no"
jumpl = Yes

#path (folder) for .ms2 if .ms2 was used for searching
ms2_path = /Users/spoudel1/Desktop/PTM_study/ChaoPeng

#path (folder) for .mzXML 
mzXML_path = /Users/spoudel1/Desktop/PTM_study/ChaoPeng

#files used to extract ms2 information. options (ms2 or mzxml). If option is .ms2, the program will only use ms2_path and if option is .mzXML, the program will use mzXML_path

ms2_fileType = ms2

#list of symbols used in the modifications (dynamic modificaitons). This information you can get if from pep.xml files
dynMods = @:15.99492,*:79.966331,#:79.966331,%:79.966331

#list static modification in the parameter file. If you have n-term modification n:229.162932
staticMods = C:57.021464

#parameters for matching ions
#ion_types to analyze a.b.c,x,y,z
ion_types = b,y

#neutral losses to analyze H2O, NH3, H3PO4
ionLoss = NH3,H2O 

#tolerance for ms2 matching = fragment ion tolerance
tol=10
'''


# In[4]:


peptideFile = config["caseEvaluation"]["peptideFile"]
jump_f_id = config["caseEvaluation"]["jump_f_id"]
jumpl = config["caseEvaluation"]["jumpl"]
ms2_path = config["caseEvaluation"]["ms2_path"]
mzXML_path = config["caseEvaluation"]["mzXML_path"]
ms2_fileType = config["caseEvaluation"]["ms2_fileType"]
dynMods = config["caseEvaluation"]["dynMods"]
staticMods = config["caseEvaluation"]["staticMods"]
ion_types = config["caseEvaluation"]["ion_types"]
ionLoss = config["caseEvaluation"]["ionLoss"]
tol = config["caseEvaluation"]["tol"]
out_fol = config["caseEvaluation"]["out_fol"]


# In[5]:


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


# In[6]:


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


# In[7]:


# peptideFile = "../PTM_study/ChaoPeng/solubleTauPTMs/phosphorylation/ID.lscore"
df_pho2 = pd.read_csv(jump_f_id,delimiter=";",skiprows=return_skiprows(jump_f_id,";", "Peptide"))


# In[8]:


# inputCheck = "../PTM_study/ChaoPeng/Check_pho_peptides_Chao.xlsx"
inputCheck =peptideFile


# In[9]:


# jumpl = "Yes"


# In[10]:


def prepareInput(inputFile):
    df = pd.read_excel(inputFile)
    df["spectrum"] = df.Exp+"."+df.scan.astype("str")+"."+df.charge.astype("str")
    return df


# In[11]:


inputDF = prepareInput(inputCheck)


# In[12]:


df_pho2["spectrum"] = df_pho2.apply(createOutfile, df=df_pho2, axis=1)

if jumpl.upper()=="YES":
    df_pho2.rename(columns={"Peptide":"JUMP-f_Peptide","JUMPl_site":"Peptides"}, inplace=True)
    scoreDict = dict(zip(df_pho2.spectrum,df_pho2.JUMPl_score))
else:
    df_pho2.rename(columns={"Peptide":"Peptides"}, inplace=True)
    scoreDict = {}
    
df_pho2[["exp","scan","charge"]] = df_pho2["spectrum"].str.split(".",expand=True)


# In[13]:


df_pho2.drop_duplicates(subset="spectrum", inplace=True, keep="first")


# In[14]:


df_pho3 = df_pho2.loc[df_pho2.spectrum.isin(list(inputDF.spectrum))]
df_pho = df_pho3.copy()
df_pho["expScan"] = df_pho.exp+"."+df_pho.scan.astype("str")


# In[15]:


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


# In[16]:


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


# In[17]:


if ms2_fileType.upper() == "MZXML":
    mzXML_list = glob.glob(mzXML_path+"/*.mzXML")
else:
    mzXML_list = glob.glob(ms2_path+"/*.ms2")


# In[18]:


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


# In[19]:


def mzmlDF(mzXML_list, spec_list):
  x1 = pyteomics.mzxml.read(mzXML_list[0])  #reading mzXML file using pyteomics
  df = pd.DataFrame([x for x in x1])  #dataframe of the mzXML file
  ms2 = df.copy().loc[df.msLevel==2]     #ms2 level scans
  run = mzXML_list[0].split("/")[-1].split(".mzXML")[0]
  ms2["expScan"] = run+"."+ms2.num.astype("str")
  ms2F = ms2.copy().loc[ms2.expScan.isin(spec_list)]
  ms2F.rename(columns={"num":"scan","m/z array":"m/z","intensity array":"intensity"}, inplace=True)
  reqCols = ["expScan","scan","m/z","intensity"]
  ms2F2 = ms2F[reqCols]
  for x in range(1,len(mzXML_list)):
    x2 = pyteomics.mzxml.read(mzXML_list[x])  #reading mzXML file using pyteomics
    df2 = pd.DataFrame([x for x in x2])  #dataframe of the mzXML file
    ms2_1 = df2.copy().loc[df2.msLevel==2]     #ms2 level scans
    run2 = mzXML_list[x].split("/")[-1].split(".mzXML")[0]
    ms2_1["expScan"] = run2+"."+ms2_1.num.astype("str")
    ms2F_1 = ms2_1.copy().loc[ms2_1.expScan.isin(spec_list)]
    ms2F_1.rename(columns={"num":"scan","m/z array":"m/z","intensity array":"intensity"}, inplace=True)
    ms2F_2 = ms2F_1[reqCols]
    frames = [ms2F2, ms2F_2]
    ms2F2 = pd.concat(frames)
  return ms2F2


# In[20]:


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


# In[21]:


if ms2_fileType.upper() == "MZXML":
    testID_mzXML = mzmlDF(mzXML_list, list(set(df_pho["expScan"])))
    mz_spectrumDict = dict(zip(testID_mzXML.expScan, testID_mzXML["m/z"]))
    int_spectrumDict = dict(zip(testID_mzXML.expScan, testID_mzXML["intensity"]))

    df_pho["m/z"]=df_pho["expScan"].map(mz_spectrumDict)
    df_pho["intensity"]= df_pho["expScan"].map(int_spectrumDict)
    ms2DF = df_pho[["spectrum","exp", "scan", "charge","expScan", "m/z", "intensity","JUMPl_score"]]
else:
    ms2DF=msToDF(mzXML_list, list(inputDF.spectrum))
    ms2DF["JUMPl_score"] = ms2DF.spectrum.map(scoreDict)


# In[22]:


def paramsOptToDict(parameter):
    dict_parameter = {} #dictionary to create using paramter option 
    listParameter = parameter.split(",")
    for x in listParameter:
        key_val = x.split(":")
        dict_parameter[key_val[0]]=key_val[1]
    return dict_parameter




jump_mod_dict = paramsOptToDict(dynMods)
sta_AA = paramsOptToDict(staticMods)


# In[27]:


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
            mod = str(pos)+"_"+"V_"+jump_mod_dict[aa]
            mods.append(mod)
            mod_pos.append(pos)

        if aa in sta_AA.keys():
            mod = str(i+1-len(mod_pos))+"_"+"S_"+sta_AA[aa]
            mods.append(mod)
        if "n" in sta_AA.keys():
            nterm = "1_S_"+sta_AA["n"]+"_n"
            mods.append(nterm)

    return mods


# In[28]:


df_pho[["plain_peptide","modifications"]] = df_pho.apply(makeJumpPlainPeptideAndModCols, sta_AA=sta_AA,jump_mod_dict=jump_mod_dict, axis=1)


# In[29]:


if "SequenceProbablity" not in df_pho.columns:
    df_pho["SequenceProbablity"] = "NA"


# In[30]:


reqdCols = ['spectrum', 'plain_peptide', 'modifications',
       'SequenceProbablity', 'XCorr','JUMPl_score']


# In[31]:


df_pho2 = df_pho[reqdCols]
inputFileDf = df_pho2.copy()
# inputFileDf = df_pho2.loc[df_pho2.spectrum.str.contains("q190622_VL1377_A")]


# In[32]:


inputFileDf['combined'] = inputFileDf.apply(lambda row: '\t'.join(row.values[0:3].astype(str)), axis=1)


# In[33]:


combinedDict = dict(zip(inputFileDf.spectrum,inputFileDf.combined))
plainPepDict = dict(zip(inputFileDf.spectrum,inputFileDf.plain_peptide))
modSiteDict = dict(zip(inputFileDf.spectrum,inputFileDf.modifications))



ms2DF["spectrum_peptide_mod"] = ms2DF.spectrum.map(combinedDict)
ms2DF["plain_peptide"] = ms2DF.spectrum.map(plainPepDict)
ms2DF["mod_site"] = ms2DF.spectrum.map(modSiteDict)


# In[36]:


ms2DF2=ms2DF.dropna(subset=["spectrum_peptide_mod"]).reset_index()


# In[37]:


newDF = ms2DF2.rename(columns={"m/z":"exp_mz_list","intensity":"intensity_list"})


# In[38]:


file_path_dyn = os.getcwd()


# In[39]:


def mkdir(dir1):
  cmd = "mkdir "+dir1
  try:
    os.system(cmd)
  except:
    "Dynamic Modification Directory exits!!"


# In[40]:


def modsForReport(mods):
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


# In[41]:


def reformat_dataframe2(df, posMassDict, match_list):
  matched_list_array = list(np.array(match_list).round(4))
  dfNew = df.copy() #duplicate dataframe to work

  df2 = dfNew.rename(columns={"Peptide_Mod_Seq":"Seq"})
  
  seqSymbol = createSymbolDict(posMassDict)

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


# In[42]:


def excelWriter2(df,  worksheetName, figure1, figure2,spectrumSplit,xcorr, prob,lscoreSite, massSeriesLength, matched, seqSymbol, N=34, updater=0):
  text1 = "Spectrum = "+spectrumSplit[0]+"; Plain Peptide Sequence = "+spectrumSplit[1]+"; Modifications on the Peptide = "+modsForReport(spectrumSplit[2])
  text2 = "Xcorr = "+str(xcorr)+"; Localization Probablity Score = "+ str(prob)+"; Localization Site:Score = "+ lscoreSite
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

  


# In[43]:


def spectrumToDict(spectrum):
  dict1 = {}
  spectrumCommaSplit = spectrum.split(",")
  for x in spectrumCommaSplit:
    y=x.split("_")
    dict1[y[0]] = float(y[2])
  return dict1


# In[44]:


def displaySeqTable(row,massPosDict): #massPosDict example is massPosList[0]
  aa = row.Seq
  pos = row["b-series"]
  if str(pos) in massPosDict.keys():
    massShift = massPosDict[str(pos)]
    aa2 = aa+"("+str(massShift)+")"
  else:
    aa2 = aa
  return aa2


# In[45]:


newDir = out_fol+"/ManualValidation_dta_ms2_Deiostope_test"
mkdir(newDir)
workbookName = newDir+"/Report_of_manual_validation.xlsx"
writer = pd.ExcelWriter(workbookName, engine='xlsxwriter')

N=34
updater=0

for specs in list(inputFileDf['combined']):
#   print (specs)
  xcorr= "missing"
  prob = "missing"
  lscoreSite = "missing"
  spectrumSplit = specs.split("\t")
  filename = newDir+"/"+"__".join(spectrumSplit)
  if "XCorr" in inputFileDf.columns:
    xcorr = inputFileDf.loc[inputFileDf['combined'] == specs].XCorr.values[0]
  if "SequenceProbablity" in inputFileDf.columns:
    prob = inputFileDf.loc[inputFileDf['combined'] == specs].SequenceProbablity.values[0]
  if "JUMPl_score" in inputFileDf.columns:
    lscoreSite = inputFileDf.loc[inputFileDf['combined'] == specs].JUMPl_score.values[0]
  spectrum_DF = newDF.loc[newDF.spectrum_peptide_mod == specs] 
#   print (spectrum_DF)
  peptide = spectrum_DF.loc[spectrum_DF.spectrum_peptide_mod == specs].plain_peptide.values[0]
  modsOri = spectrum_DF.loc[spectrum_DF.spectrum_peptide_mod == specs].mod_site.values[0]
  maxCharge = int(spectrum_DF.loc[spectrum_DF.spectrum_peptide_mod == specs].charge.values[0])
  mods = modsForReport(modsOri)
#   modsVariantsAll, massPosList = statModDynModAllVariants(peptide, mods, dyn_AA)
  massPosDict1 = spectrumToDict(spectrumSplit[-1])

#   print (massPosDict1)

  ion_types = ["b","y"]
  ionLoss = ["H2O","NH3","H3PO4"]
  
  df_pep = ionSeriesIonLossSpeRes(peptide,maxcharge=maxCharge,massPosDict=massPosDict1,useMod ="Yes")
  

  reqdCols = ["Seq"]
  for x in df_pep.columns:
    for y in ion_types:
      if y+"+" in x:
        reqdCols.append(x)
      for z in ionLoss:
        if y+"-"+z in x:
          reqdCols.append(x)
        
  df_pep2 = df_pep.copy()[reqdCols]   
  df_pep2["b-series"] = [*range(1, len(df_pep2)+1, 1)] 
  df_pep2["y-series"] = [*range(len(df_pep2),0, -1)]
  
  exp_mz_list = spectrum_DF.exp_mz_list.values[0]
  exp_int_list = spectrum_DF.intensity_list.values[0]
  
  match_list, match_int_list, true_ions_ppm_list=ionMatches(exp_mz_list, exp_int_list, df_pep, ion_types =ion_types, ionLoss=ionLoss, tol=10)
  
  
  
  
    
  df_pep2["Peptide_Mod_Seq"] = df_pep2.apply(displaySeqTable, massPosDict=massPosDict1, axis=1)

  displayTableCols = []
  for cols in df_pep2.columns:
    if "b" in cols:
      if "series" not in cols:
        displayTableCols.append(cols)
  displayTableCols=displayTableCols[::-1]
  displayTableCols.append("b-series")
  displayTableCols.append("Peptide_Mod_Seq")
  displayTableCols.append("y-series")
  for cols in df_pep2.columns:
    if "y" in cols:
      if "series" not in cols:
        displayTableCols.append(cols)

  df=df_pep2.copy()[displayTableCols]

  seriesName = []
  
  for vals in match_list:
    val = float(vals)
    row_column_pair = df_pep2[df_pep2.isin([val])].stack().index[0]
#     print (row_column_pair)
    row = row_column_pair[0]
    column = row_column_pair[1]
  #   print (row, column) 
    if "-" in column:
      columnSplit = column.split("-")
      if ("x" in column) or ("y" in column) or ("z" in column):
        ionNu = columnSplit[0]+str(df_pep2["y-series"][row])+"-"+columnSplit[1]
      if ("a" in column) or ("b" in column) or ("c" in column):
        ionNu = columnSplit[0]+str(df_pep2["b-series"][row])+"-"+columnSplit[1]
    else:
      columnSplit = column.split("+")
      if ("x" in column) or ("y" in column) or ("z" in column):
        ionNu = columnSplit[0]+str(df_pep2["y-series"][row])+"+"+columnSplit[1]
      if ("a" in column) or ("b" in column) or ("c" in column):
        ionNu = columnSplit[0]+str(df_pep2["b-series"][row])+"+"+columnSplit[1]
    seriesName.append(ionNu)
  seriesName = list(map(lambda x: x.replace('+1',''),seriesName))
#   print (len(list(spectrum_DF.matched_ions_list)[0]))
  plt.rcParams.update({'font.size': 8,'figure.max_open_warning': 0})
  plt.figure(figsize=(7,2.5))
  #plt.title(" ".join(specs.split("\t")))
  plt.xlabel('m/z')
  plt.ylabel('Intensity, absolute. units')
  plt.bar(list(spectrum_DF.exp_mz_list)[0], list(spectrum_DF.intensity_list)[0], width=0.1, linewidth=0.5,
        edgecolor='black')

  plt.bar(match_list, match_int_list, width=0.1, linewidth=0.5,
        edgecolor='red', alpha=0.7)

  plt.xlabel("m/z", color='black')

  plt.ylabel('Absolute Intensity', color='black')
  plt.tick_params(color='black')
  plt.xticks(color="black")
  plt.yticks(color="black")

  for i, txt in enumerate(seriesName):
    for label in ["b","y"]:
      if label in txt:
        if "-" not in txt:
          plt.annotate(txt, (match_list[i],match_int_list[i]), ha="center")


  
  #figurename = filename+"__intensityPlot.pdf"
  figurename1 = filename+"__intensityPlot.png"
  #plt.savefig(figurename, bbox_inches="tight", dpi=600 )
  plt.savefig(figurename1, bbox_inches="tight", dpi=600 )


  ionDFmatch = pd.DataFrame(columns=["matched_ions","matched_intensity","matched_ppm"])
  
  
  ionDFmatch["matched_ions"] = match_list
  ionDFmatch["matched_intensity"] = match_int_list
  ionDFmatch["matched_ppm"] = true_ions_ppm_list 
  bins = [0,10000,100000,1000000, 10000000,100000000,1000000000,10000000000]
  ionDFmatch['intRange'] = pd.cut(ionDFmatch.matched_intensity, bins=bins, include_lowest=True)
  ionDFmatch['right'] = ionDFmatch['intRange'].map(attrgetter('right'))
  ionDFmatch["intensity"] = ionDFmatch.apply(rightInterval,axis=1)


  fig, ax = plt.subplots(figsize=(6,2.5))
  plt.style.use("ggplot")


  plt.rcParams['axes.edgecolor'] = "#010101"
  plt.rcParams['axes.facecolor'] = '#FFFFFF'


  plt.rcParams.update({'font.size': 8,'figure.max_open_warning': 0})
  ax = plt.axes(facecolor="w")

  ionDFmatch2 =ionDFmatch.sort_values(by="right", ascending=True)
  sns.scatterplot(data=ionDFmatch2, x="matched_ions", y="matched_ppm",s=25,hue="intensity")
  
  plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  plt.hlines(0, 50, 2200, color = "green", linestyles="dotted",linewidth=0.5)
  ax.set_ylabel("mass error (ppm)", color="black")
  ax.set_xlabel("all matched ions", color="black")
  ax.set_xlim(0,max(match_list)+100)
  ax.tick_params(axis='x', colors='black')
  ax.tick_params(axis='y', colors='black')
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  for i, txt in enumerate(seriesName):
    for label in ["b","y"]:
      if label in txt:
        if "-" not in txt:
          plt.annotate(txt, (list(ionDFmatch.matched_ions)[i],list(ionDFmatch.matched_ppm)[i]+0.3), ha="center")

          name = "_".join(specs.split("\t"))
  
  #figurename = filename+"__ToleranceMap.pdf"
  figurename2 = filename+"__ToleranceMap.png"
    
  #plt.savefig(figurename, bbox_inches="tight", dpi=600 )
  plt.savefig(figurename2, bbox_inches="tight", dpi=600 )
  
  matched = str(len(ionDFmatch))+"/"+str(len(list(spectrum_DF.exp_mz_list)[0]))

  figurename3 = filename+"__MasterSeries.png"
  

  finalDF, seqSymbol = reformat_dataframe2(df, massPosDict1, match_list)
  finalDF2 = finalDF.style.applymap(lambda x: 'color: red' if x in match_list else 'color: black')

  massSeriesLength = df.shape[0]
  

  writer,updater, N = excelWriter2(finalDF2, "Sheet1", figurename1, figurename2, spectrumSplit,xcorr, prob,lscoreSite, massSeriesLength,matched, seqSymbol,N, updater)
writer.save()


