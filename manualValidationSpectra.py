#!/usr/bin/env python
# coding: utf-8

# In[198]:


import pandas as pd
from docx import Document
from docx.shared import Inches
import numpy as np
from ast import literal_eval
import configparser
config = configparser.ConfigParser()
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
from os.path import dirname
from pandas.plotting import table
import dataframe_image as dfi
from docx.shared import Pt
sys.path.append("/home/spoudel1/bin/python/JUMP_l/")
from JUMP_l_modules import *

# In[169]:

params_file = sys.argv[1]
#params_file = "JUMPl_phosphoValidate/map_comet_jump_f.params"
config.read(params_file)


# In[170]:


file_path_dyn = config["caseEvaluation"]["file_path_dyn"]
dyn_AA = config["caseEvaluation"]["dyn_AA"]
dyn_mod = config["caseEvaluation"]["dyn_mod"]
jump_f_dyn = config["caseEvaluation"]["jump_f_dyn"]
ptm = config["caseEvaluation"]["ptm"]
ion_series = config["caseEvaluation"]["ion_series"]
ion_losses_MS2 = config["caseEvaluation"]["ion_losses_MS2"]
tol_max =  config["caseEvaluation"]["tol_max"]

ion_types = ion_series.split(",")
ionLoss = ion_losses_MS2.split(",")
#print (ionLoss)
# In[171]:

#inputFile = dirname(file_path_dyn)+"/SampleInputFileManualValidation_v001.xlsx"
inputFile = dirname(file_path_dyn)+"/PhoPeptidesSelected_v001.xlsx"
inputFileDf = pd.read_excel(inputFile)
inputFileDf['combined'] = inputFileDf[list(inputFileDf.columns)].apply(lambda row: '\t'.join(row.values[0:3].astype(str)), axis=1)
inputFileDf['Run'] = inputFileDf.Spectrum.str.split(".",expand=True)[0]

all_exp = list(set(inputFileDf['Run']))

interest_files = []
for run in all_exp:
  file1 = "ManualValidationFiles/"+run+"SaveDfToFile.txt"
  interest_files.append(file1)

print ("Interested Files are ",interest_files)

def merge_df_vertical(file_list):
  df = pd.read_csv(file_list[0], delimiter="\t")
  for x in range(1,len(file_list)):
    dfnew = pd.read_csv(file_list[x], delimiter="\t")
    frames = [df, dfnew]
    df = pd.concat(frames)
  return df

#newDF = pd.read_csv("SaveDfToFile.txt", delimiter="\t")
newDF2 =  merge_df_vertical(interest_files)

print ("The interested files are merged to create a new dataframe of shape ",newDF2.shape)

neededSpecsPeptide = dict(zip(inputFileDf.combined, inputFileDf.newPeptide))
#neededSpecs = list(inputFileDf['combined'])

# In[172]:
newDF = newDF2.loc[newDF2.spectrum_peptide_mod.isin(neededSpecsPeptide.keys())]

cols = ['exp_mz_list', 'intensity_list',
       'matched_ions_list', 'matched_int_list', 'true_ions_ppm_list',
       'theoretical_ions_matched']
for i,val in enumerate(cols):
  newDF.loc[:,val] = newDF.loc[:,val].apply(lambda x: literal_eval(x))
  print ("The string to list conversion of ", val ," column is completed")
  print (i," columns is/are done")


# In[173]:


def displaySeqTable(row,massPosDict): #massPosDict example is massPosList[0]
  aa = row.Seq
  pos = row["b-series"]
  if str(pos) in massPosDict.keys():
    massShift = massPosDict[str(pos)]
    aa2 = aa+"("+str(massShift)+")"
  else:
    aa2 = aa
  return aa2


# In[174]:


def render_mpl_table(data, col_width=1, row_height=0.3, font_size=7,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')
    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in mpl_table._cells.items():
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    return ax.get_figure(), ax


# In[175]:


def replaceSomeVals(row, seqSymbol):
  seq = row.Seq
  for key in seqSymbol.keys():
    if key in seq:
      seq = seq.replace(key, seqSymbol[key])
  return seq


# In[176]:

def createSymbolDict(posMassDict):
  seqSymbol = {}
  for vals in posMassDict.values():
    if vals==229.162931:
      seqSymbol["(229.162931)"]  = "%"
    elif vals==229.162932:
      seqSymbol["(229.162932)"]  = "%"
    elif vals==15.994915:
      seqSymbol["(15.994915)"]  = "@"
    elif vals == 57.021464:
      seqSymbol["(57.021464)"]  = "$"
    else:
      newVal = "("+str(vals)+")"
      seqSymbol[newVal]  = "#"
  return seqSymbol


# In[177]:


def reformat_dataframe(df, posMassDict, match_list):
  
  dfNew = df #duplicate dataframe to work

  df2 = dfNew.rename(columns={"Peptide_Mod_Seq":"Seq"})
  df3 = df2.round(4)
  seqSymbol = createSymbolDict(posMassDict)

  matched_list_array = list(np.array(match_list).round(4))
  for name, values in df3.iteritems():
    for val in range(0, df3.shape[0]):
  #     print('{name}: {value}'.format(name=name, value=values[val]))
      value = values[val]
      if value in matched_list_array:
        df3 = df3.replace(value, str(value)+"*")
      if name == "Seq":
        for keys in seqSymbol.keys():
          if keys in value:
            new_val = value[0]
            df3 = df3.replace(value,value[0]+seqSymbol[keys]) 
      
  return df3, seqSymbol
  


# In[178]:


def mkdir(dir1):
  cmd = "mkdir "+dir1
  try:
    os.system(cmd)
  except:
    "Dynamic Modification Directory exits!!"


# In[179]:


def spectrumToDict(spectrum):
  dict1 = {}
  spectrumCommaSplit = spectrum.split(",")
  for x in spectrumCommaSplit:
    y=x.split("_")
    dict1[y[0]] = float(y[2])
  return dict1
    


# In[180]:


#inputFile = dirname(file_path_dyn)+"/SampleInputFileManualValidation_v001.xlsx"
#inputFileDf = pd.read_excel(inputFile)
#inputFileDf['combined'] = inputFileDf[list(inputFileDf.columns)].apply(lambda row: '\t'.join(row.values[0:3].astype(str)), axis=1)


# In[181]:


def reformat_dataframe2(df, posMassDict, match_list):
  matched_list_array = list(np.array(match_list).round(4))
  dfNew = df #duplicate dataframe to work

  df2 = dfNew.rename(columns={"Peptide_Mod_Seq":"Seq"})
  df3 = df2.round(4)
  seqSymbol = createSymbolDict(posMassDict)

  matched_list_array = list(np.array(match_list).round(4))
  for name, values in df3.iteritems():
    for val in range(0, df3.shape[0]):
  #     print('{name}: {value}'.format(name=name, value=values[val]))
      value = values[val]
      
      if name == "Seq":
        for keys in seqSymbol.keys():
          if keys in value:
            new_val = value[0]
            df3 = df3.replace(value,value[0]+seqSymbol[keys]) 
            
  df4=df3.style.applymap(color, match_list = matched_list_array)

  return df4, seqSymbol


# In[182]:


def color(val, match_list):
    matched_list_array = list(np.array(match_list).round(4))
    if val in matched_list_array:
        color = 'red'
    else:
        color = 'black'
#     return 'background-color: %s' % color
    return 'color: %s' % color


# In[201]:
def excelWriter2(df,  worksheetName, figure1, figure2,spectrumSplit,newPeptide,xcorr, prob, massSeriesLength, matched, seqSymbol, N=34, updater=0):
  text1 = "Spectrum = "+spectrumSplit[0]+"; Peptide Sequence = "+newPeptide+"; Mods = "+modsForReport(spectrumSplit[2])
  text2 = "Xcorr = "+str(xcorr)+"; Localization Score = "+ str(prob)
  text3 = "Figure a) Intensity Plot"
  text4 = "Figure b) Mass deviation of matched ions"
  text5 = "Figure c) Product Ion series" 
  text6 = "Matched ion = "+matched
  text7 = "Symbols and their meaning"
  text8 = ""
  for masses in seqSymbol.keys(): 
    text8+=seqSymbol[masses]+"="+masses[1:-1]+"; "
    
  
  df.to_excel(writer, sheet_name=worksheetName, startrow=N, index=None)
  length = 40+massSeriesLength+1
  # Get the xlsxwriter workbook and worksheet objects.
  workbook  = writer.book
  worksheet = writer.sheets[worksheetName]
    
  worksheet.write(updater, 0, text1)
  worksheet.write(updater+1, 0, text2)
  worksheet.write(updater+2, 0, text3)

  location = "A"+str(updater+5)
  worksheet.insert_image(location,figure1)
  
  worksheet.write(updater+20, 0, text4)  
    
  location2 = "A"+str(updater+23)
  worksheet.insert_image(location2,figure2)
 
  worksheet.write(updater+38, 0, text5)  

  worksheet.write(updater+length+2, 0, text6)
  worksheet.write(updater+length+3, 0, text7)
  worksheet.write(updater+length+4, 0, text8)

  updater = updater+length+6
  N = updater + 41

  return writer, updater, N

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


newDir = file_path_dyn+"/ManualValidationFiguresAndDocuments"
mkdir(newDir)
workbookName = file_path_dyn+"/Report_of_manual_validation.xlsx"
writer = pd.ExcelWriter(workbookName, engine='xlsxwriter')

N=40
updater=0
for specs in neededSpecsPeptide.keys():
  #print (specs)
#for specs in list(inputFileDf['combined']):
  
  xcorr= "missing"
  prob = "missing"
  spectrumSplit = specs.split("\t")
  filename = newDir+"/"+"__".join(spectrumSplit)
  if "Xcorr" in inputFileDf.columns:
    xcorr = inputFileDf.loc[inputFileDf['combined'] == specs].Xcorr.values[0]
  if "SequenceProbablity" in inputFileDf.columns:
    prob = inputFileDf.loc[inputFileDf['combined'] == specs].SequenceProbablity.values[0]
  spectrum_DF = newDF.loc[newDF.spectrum_peptide_mod == specs] 
#   print (spectrum_DF)
  peptide = spectrum_DF.loc[spectrum_DF.spectrum_peptide_mod == specs].plain_peptide.values[0]
  modsOri = spectrum_DF.loc[spectrum_DF.spectrum_peptide_mod == specs].mod_site.values[0]
  maxCharge = spectrum_DF.loc[spectrum_DF.spectrum_peptide_mod == specs].charge.values[0]
  mods = modsForReport(modsOri)
#   modsVariantsAll, massPosList = statModDynModAllVariants(peptide, mods, dyn_AA)
  massPosDict1 = spectrumToDict(spectrumSplit[-1])
  newPeptide = neededSpecsPeptide[specs]
#   print (massPosDict1)

  #ion_types = ["b","y"]
  #ionLoss = ["H2O","NH3","H3PO4"]
  
  df_pep = ionSeriesIonLossSpeRes(peptide,maxcharge=maxCharge,massPosDict=massPosDict1,useMod ="Yes")
  

  reqdCols = ["Seq"]
  for x in df_pep.columns:
    for y in ion_types:
      y=y.strip()
      if y+"+" in x:
        reqdCols.append(x)
      for z in ionLoss:
        z=z.strip()
        if y+"-"+z in x:
          reqdCols.append(x)
        
  df_pep2 = df_pep[reqdCols]   
  df_pep2["b-series"] = [*range(1, len(df_pep2)+1, 1)] 
  df_pep2["y-series"] = [*range(len(df_pep2),0, -1)]
  
  exp_mz_list = spectrum_DF.exp_mz_list.values[0]
  exp_int_list = spectrum_DF.intensity_list.values[0]
  
  match_list, match_int_list, true_ions_ppm_list=ionMatches(exp_mz_list, exp_int_list, df_pep, ion_types =ion_types, ionLoss=ionLoss, tol=3)
  
  
  
    
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

  df=df_pep2[displayTableCols]

  seriesName = []
  
  for vals in match_list:
    val = float(vals)
    row_column_pair = df_pep2[df_pep2.isin([val])].stack().index[0]
  #   print (row_column_pair)
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
  #print (len(list(spectrum_DF.matched_ions_list)[0]))
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
  ax.set_xlim(0,max(list(spectrum_DF.matched_ions_list)[0])+100)
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
  
  matched = str(len(list(spectrum_DF.matched_ions_list)[0]))+"/"+str(len(list(spectrum_DF.exp_mz_list)[0]))

  figurename3 = filename+"__MasterSeries.png"
  

  finalDF, seqSymbol = reformat_dataframe2(df, massPosDict1, match_list)
   
  massSeriesLength = df.shape[0]
  
  
  writer,updater, N = excelWriter2(finalDF, "Sheet1", figurename1, figurename2, spectrumSplit, newPeptide, xcorr, prob, massSeriesLength,matched, seqSymbol,N, updater)
writer.save()
