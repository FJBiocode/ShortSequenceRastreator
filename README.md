# ShortSequenceRastreator
Developed by Francisco Javier Bielsa (Graduated in Biotechnology and Master in Molecular and Cell Biology. University of Zaragoza). 
Greetings! I hope my work helps yours and your SSR studies. 
This code was made during the development of my Master's Thesis and was tested with Genomic Data from Pyrus Communis L. held at CITA of Aragón Germplasm bank. 

ShortSequenceRastreator allows Genetic diversisity analysis and input file preparation for STRUCTURE 2.3.4 (Stanford), DARwin 6.0 (CIRAD) and SPAGeDi (Iowa State Uiversity).

This program works with CSV-UTF8 encoded dataframes. Header is expected  in the first row. 
Rows contain sample data, each row represents an accession.
The first columns of your dataframe should contain important information such as name, origin or ploidy. These should be followed by loci columns containing genomic data.
I.E:

| Accession  | Name  | Origin | Ploidy  | Locus.1  | Locus.2 |
| :------------ |:---------------:| :-----:| :------------ |:---------------:| -----:|
| 100      | M.1 | SPA | 2      | 123/123 | 241/243 |
| 101      | M.2 | SPA | 3      | 123/125/127 | 241/243/247 |
| 102      | M.3 | FRA | 2      | 123/125 | 241/247 |
| 103      | M.4 | SPA | 2      | 121/123 | 241/243 |
 
This program supports diploid and triploid SSR analysis. 
Insert path of the file you want to work with - f.e. 'C:/Users/Admin/Documents/SSR/Dataframe1.csv'

Before executing the code set of columns containing information should be checked and changed:  set(['Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy']).

FEATURES:

ACCESSIONS SUMMARY: Number of accessions, Number of loci, Percentage of diploid and triploids.

CREATION OF HITS MATRIX / DARwin MATRIX : Here we are creating a 'Hits' dataframe (Accession vs number of times an allele appears in a given locus).
                                          This matrix can be used in DARwin software and can be saved to the path we are working on. 
                                          Be careful with DARwin parameters, because the matrix may contain values different from 1 or 0 and 
                                          may give errors in DARwin. This can happen because homozygotes/triploids give 2 hits in some alleles.
                                          In order to make the matrix work in DARwin you should replace all values=2 to values=1.

ALLELE SUMMARY: Relative frequencies, rare alleles and unique alleles.

DISCRIMINANT POWER

OBSERVED AND EXPECTED HETEROZYGOSITY

F-STAT AND NULE ALELES 

DUPLICATE IDENTIFICATION : Get a dataframe with your duplicates! WARNING: In 'col[6:(len(new_df.columns))]' number 6 need to be checked and replaced 
with the last position of information columns that you have in your dataframe.

TRIALELIC LOCI SEARCH

CRETION OF STRUCTURE INPUT FILE: Enjoy saving 2 hours of excel work! WARNING: Insert file path of your unique genotypes matrix. Change Information columns in '{'Ord', 'Nclon', 'AccessionName', 'Origin', 'Ploidy'}' and St = Stru.drop(['Ord', 'Nclon', 'Origin', 'Ploidy'], axis = 1).


CREATION OF SPAGeDI INPUT FILE: Feel free to dance to the rythm of https://www.youtube.com/watch?v=Z3w5gVM_4y8 while python does your work. 
