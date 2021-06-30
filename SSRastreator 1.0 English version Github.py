# -*- coding: utf-8 -*-
"""
Created on June 2021

SHORT SEQUENCE RASTREATOR v.1.0 

@author: Francisco Javier Bielsa González || Graduate in Biotechnology, 
         Master in Molecular and Cell Biology.

SSRastreator works with Allele Data from SSR sequencing and can be used to study genetic diversity.
SSRastreator also fastens workflow between genetic structure and diversity analysis softwares such as 
STRUCTURE 2.3.4 (Stanford University) or DARwin (CIRAD) preparing input files to both softwares with
your original data sheet.

This code was made during the development of my Master's Thesis and was tested with Genomic Data from Pyrus Communis L.


"""


import pandas as pd
import numpy as np
import itertools


def convert_to_list(x):
  if x == '//' or x == '/' or isinstance(x, float):
    return []
  else:
    return [int(i) for i in list(filter(None, x.split('/')))]


def possible_Alleles(Column):
  total = set()
  for accession in Column:
    total = total | set(accession)
  return sorted(list(total))


def Count_hits(dataframe, new_dataframe, Locus, Alleles):
  for allele in Alleles:
    new_Column = Locus + '.' + str(allele)
    Column = []
    # print(new_Column)
    for accession in dataframe[Locus]:
      Column.append(accession.count(allele))
    new_dataframe[new_Column] = Column
  return new_dataframe


def Relative_frequencies(TableF, Locus, Alleles):
        for allele in Alleles:
            freq = PoolOfAlleles.count(allele)
            relfreq = freq/len(PoolOfAlleles)
            print(Locus,allele,relfreq)
            allele = str(allele)
            RelFreq = pd.Series([Locus, allele, relfreq], index = TableF.columns)
            TableF = TableF.append(RelFreq, ignore_index=True).copy()
        return TableF
    
    
def convert_to_list_het(x):
  if x == '//' or x == '/' or x =='///' or x =='////' or isinstance(x, float):
    return ['nan']
  else:
    return [int(i) for i in list(filter(None, x.split('/')))]


def Observed_Heterozigosity(HObs, HetM):
   Heterozygotes = 0
   Accessions = HetM[Locus].count()
   for accession in HetM[Locus]:
       print(accession)
       if accession == ['nan']:
           continue
       if accession[0] != accession[1]:
           Heterozygotes +=  1
       print(Heterozygotes)
   HO = Heterozygotes/Accessions
   Obser = pd.Series([HO,Locus], index = HObs.columns)
   HObs = HObs.append(Obser, ignore_index = True).copy()
   return(HObs)


def Expected_Heterozigosity(HJ,HetM):
    nume = 0  
    for allele in Alleles:
        print(allele)
        freq = PoolOfAlleles.count(allele)
        relfreq = freq/len(PoolOfAlleles)
        pi2 = pow(relfreq, 2)
        nume = nume + pi2
    Hj = 1 - (nume)
    Expected = pd.Series([Hj,Locus], index = HJ.columns)
    HJ = HJ.append(Expected, ignore_index = True).copy()
    return(HJ)


def Unique_Alleles(UniqueAlleles,Locus,Alleles):
    for allele in Alleles:
        if PoolOfAlleles.count(allele) == 1:
            UniqueAllele = pd.Series([allele, Locus], index = UniqueAlleles.columns)
            UniqueAlleles = UniqueAlleles.append(UniqueAllele, ignore_index = True)
            print(allele, 'is an unique allele for Locus', Locus)
    return(UniqueAlleles)


def Number_of_Alleles(NumberOfAlleles, Locus, Alleles):
    NuAl = pd.Series([len(Alleles),Locus], index = NumberOfAlleles.columns)
    NumberOfAlleles = NumberOfAlleles.append(NuAl, ignore_index = True)
    return (NumberOfAlleles)


"""

Insert path of the file you want to work with - f.e. 'C:/Users/Admin/Documents/SSR/Dataframe1.csv'
This program works with CSV-UTF8 encoded dataframes. Header is expected  in the first row. 
Rows contain sample data, each row represents an accession.
The first columns of your dataframe should contain important information such as name, origin or ploidy.
This program supports diploid and triploid SSR analysis. 

"""

    
path = 'C:/Users/Admin/Documents/SSR/Dataframe1.csv'    


#Creation of working dataframe and another dataframe that will result into DARwin working Matrix.
#loci: Set of all locus. Columns containing information different from genomic data must be included
#      in the second part of the code line 'set(['Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy'])'. 


df=pd.read_csv(path , sep=';',header=0)
loci = set(df.columns.values) - set(['Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy'])
new_df = df[['Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy']].copy()


#Accessions summary


print('Number of Accessions is', len(df),'\nNumber of loci is', len(loci))
Diploids = 0
Triploids = 0
for x in df['Ploidy']:
    if x == 2:
        Diploids +=1
    else:
        Triploids +=1
print('Diploid percentage is', Diploids/len(df))
print('Triploid percentage is', Triploids/len(df))


NDiploids = ('Diploid percentage is', Diploids/len(df))
NTriploids = ('Triploid percentage is', Triploids/len(df))


"""

CREATION OF HITS MATRIX / DARwin MATRIX 

Here we are creating a 'Hits' dataframe (Accession vs number of times an allele appears in a given locus)
This matrix can be used in DARwin software and can be saved to the path we are working on. 
Be careful with DARwin parameters, because the matrix may contain values different from 1 or 0 and 
may give errors in DARwin. This can happen because homozygotes/triploids give 2 hits in some alleles.
In order to make the matrix work in DARwin you should replace all values=2 to values=1.

"""


for Locus in loci:
  #print('Procesando Locus', Locus, '...')
  df[Locus] = df[Locus].apply(lambda x: convert_to_list(x))
  Alleles = possible_Alleles(df[Locus])
  new_df = Count_hits(df, new_df, Locus, Alleles)


ExcelHits = new_df
ExcelHits.to_csv('Darwin.csv', index = False)


"""

CREATION OF ALLELIC FREQUENCY TABLE

Outputs: 'TableF.csv' contains relative frequencies of alleles.
      'UniqueAlleles.csv' containes unique alleles in a given locus.
      'NumberOfAlleles.csv' contains number of alleles per locus.
      'RareAlleles' contains rare alleles found at a given locus.
      
"""


Tf = ['Locus', 'allele', 'Frequency']
TableF = pd.DataFrame([], columns = Tf)
UniqueAlleles = pd.DataFrame([], columns = ['Unique Allele', 'Locus'])
NumberOfAlleles = pd.DataFrame([], columns = ['Number of Alleles','Locus'])


for Locus in loci:
    AList = df[Locus].values.tolist()
    PoolOfAlleles = list(itertools.chain(*AList))
    Alleles = possible_Alleles(df[Locus])
    UniqueAlleles = Unique_Alleles(UniqueAlleles,Locus,Alleles)
    TableF = Relative_frequencies(TableF,Locus, Alleles)
    TableF.to_csv('TableF.csv', index = False)
    UniqueAlleles.to_csv('UniqueAlleles.csv', index = False)
    NumberOfAlleles = Number_of_Alleles(NumberOfAlleles, Locus, Alleles)
    NumberOfAlleles.to_csv('NumberOfAlleles.csv', index = False)
    
    
RareAlleles = TableF[TableF['Frequency'] < 0.05]
RareAlleles = RareAlleles.reset_index(drop=True)
RareAlleles.to_csv('RareAlleles.csv', index = False)


"""

DISCRIMINANT POWER OF EACH Locus - PD = 1- Epi^2 (Genotypes)

Output: 'PD.csv' Contains DP of each locus.

"""


df=pd.read_csv(path , sep=';',header=0)
loci = set(df.columns.values) - set(['Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy'])
DiscriminantPower = pd.DataFrame([])


for Locus in loci:
    Lista = df[Locus].values.tolist()
    nume = 0
    Genotipos = set(Lista)
    for genotipo in Genotipos:
        freq = Lista.count(genotipo)
        RelFreq = freq/len(Lista)
        pg2 = pow(RelFreq, 2)
        nume = nume + pg2
    PD = 1 - (nume)
    serie = pd.Series([Locus,PD])
    serie = pd.DataFrame([serie])
    DiscriminantPower = pd.concat([serie,DiscriminantPower], ignore_index=True)
    
    
DiscriminantPower.columns = ['Locus', 'PD']
DiscriminantPower.to_csv('PD.csv', index = False)    


"""

OBSERVED AND EXPECTED HETEROZYGOSITY

Outputs: 'HObs.csv' and 'HJ.csv'

"""
    

HetM = pd.read_csv(path , sep=';',header=0)
loci = set(HetM.columns.values) - set(['Ord', 'Nclon', 'Origin', 'Ploidy', 'AccessionName'])


#Observed Heterozygosity


Heter = ['Hobs','Locus']
HObs = pd.DataFrame([], columns = Heter)


for Locus in loci:
    HetM[Locus] = HetM[Locus].apply(lambda x: convert_to_list_het(x))
for Locus in loci:
    HObs = Observed_Heterozigosity(HObs, HetM)
print('Observed Heterozigosity is \n', HObs)


HObs.to_csv('HObs.csv', encoding = 'utf-8')


#Expected Heterozygosity


HetM = pd.read_csv(path , sep=';',header=0)
loci = set(HetM.columns.values) - set(['Ord', 'Nclon', 'Origin', 'Ploidy', 'AccessionName'])
Expected = ['Hj', 'Locus']
HJ = pd.DataFrame([], columns = Expected)
for Locus in loci:
    HetM[Locus] = HetM[Locus].apply(lambda x: convert_to_list(x))
    Lista = HetM[Locus].values.tolist()
    PoolOfAlleles = list(itertools.chain(*Lista))
    Alleles = possible_Alleles(HetM[Locus])
    HJ = Expected_Heterozigosity(HJ,HetM)
print('Heterocigosidad Expectedada es \n',  HJ)


HJ.to_csv('HJ.csv', encoding = 'utf-8')


HObsTotal = sum(HObs['Hobs'].values)/14
HJTotal = sum(HJ['Hj'].values)/14


"""

Wright F-statistics

F value || F = 1- (Hobs/Hesp) 
This function considers every individual as a part of an unique population

"""


F = 1 - (HObsTotal/HJTotal)
print('Fit value is', F)
if F > 0.2:
    print('Slight differentiation')
if F > 0.15 and F < 0.2:
    print('Strong differentiation')


"""

#Nule alleles (r=(He-Ho)/(1+He))

"""


DfEst = pd.DataFrame([HJ['Locus'],HJ['Hj'],HObs['Hobs']])
DfEst = DfEst.transpose()
DfEst['Na'] = (DfEst['Hj']-DfEst['Hobs'])/(1+DfEst['Hj'])
DfEst.to_csv('DfEst.csv', encoding = 'utf-8')



"""

Both Fvalue and Nule alleles function require further revision. Use under your own risk.

"""


"""

SUMMARY TABLE

#'Results.csv'.

"""


Summary = pd.DataFrame([])
PD_sorted = DiscriminantPower.sort_values(by = ['Locus'], inplace = False, ascending = True)
PD_sorted = PD_sorted.reset_index()
PD_sorted = PD_sorted.drop('index', axis =1 )
DfEst_sorted = DfEst.sort_values(by = ['Locus'],inplace = False, ascending = True)
NumberOfAlleles_sorted = NumberOfAlleles.sort_values(by = ['Locus'], inplace = False,ascending = True)
NumberOfAlleles_sorted = NumberOfAlleles_sorted.reset_index()
NumberOfAlleles_sorted = NumberOfAlleles_sorted.drop('index', axis = 1)
Summary = DfEst_sorted
Summary = Summary.reset_index()
Summary = Summary.drop('index', axis = 1)
Summary['PD'] = PD_sorted['PD']
Summary['Number of Alleles'] = NumberOfAlleles_sorted['Number of Alleles']
Summary['Ploidy'] = ''
Summary['Ploidy'][0]=NDiploids
Summary['Ploidy'][1]=NTriploids
Summary['Information'] = ''
Summary['Information'][0] = ('Number of Accessions is', len(df),'Number of loci is' ,len(loci))
Summary.to_csv('Results.csv', encoding = 'utf-8')
print(Summary)


""""

DUPLICATES TABLE

Output: 'Duplicates.csv'
In 'col[6:(len(new_df.columns))]' number 6 need to be checked and replaced 
with the last position of information columns that you have in your dataframe

"""


col = new_df.columns
sel = col[6:(len(new_df.columns))]
Duplicates = new_df[new_df.duplicated(sel, keep = False)]
sorti = list(sel)
Duplicates = Duplicates.sort_values(by = sorti)
Duplicates.to_csv('Duplicates.csv', encoding = 'utf-8')


""""

STRUCTURE FILE

Insert file path of your unique genotypes matrix

Change Information columns in '{'Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy'}'
and St = Stru.drop(['Ord', 'Nclon', 'Origin', 'Ploidy'], axis = 1)

Execute and enjoy saving 2 hours of excel work.

Output: 'Structure.csv'

"""


path = 'C:/Users/Admin/Documents/SSR/DataframeUniqueGenotypes.csv'    
Stru = pd.read_csv(path, sep=';', header=0)
loci = set(Stru.columns.values) - {'Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy'}
for Locus in loci:
    Stru[Locus]= Stru[Locus].apply(lambda x: convert_to_list(x))
 
    
#Convert to structure format
 
    
for i in range(len(Stru)):
    for Locus in loci:
        while len(Stru[Locus].loc[i]) < 3:
            Stru[Locus].loc[i].append('')


St = Stru.drop(['Ord', 'Nclon', 'Origin', 'Ploidy'], axis = 1)           
St = St.apply(pd.Series.explode)
St.replace('', np.nan, inplace= True)     
St.replace(np.nan,-9, inplace = True)
St.to_csv('Structure.csv', encoding = 'UTF-8')


""""

SPAGeDi File

Insert your genotype matrix file path. Does not need to contain only unique genotypes.

Change Information columns in 'set(['Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy'])'

Output : 'Spagedi.csv'

""""


path = 'C:/Users/Admin/Documents/SSR/Dataframe1.csv'     
df=pd.read_csv(path, sep=';',header=0)
loci = set(df.columns.values) - set(['Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy'])


for Locus in loci:    
    df = df.join(df[Locus].str.split('/',expand = True).add_prefix(Locus).fillna(np.nan))
    df = df.drop(columns = Locus)


newset = set(df.columns.values) - set(['Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy'])


for col in newset:
    df[col] = df[col].apply(lambda x: '{0:0>3}'.format(x))


for Locus in loci:
    df[Locus] = df[Locus + '0'] + df[Locus + '1'] + df[Locus + '2']
    del (df[Locus + '0'], df[Locus + '1'], df[Locus + '2'])
    
    
Spagedi = df
Spagedi.to_csv('Spagedi.csv', encoding = 'UTF-8')