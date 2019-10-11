## A novel and accurate tool for predicting lncRNAs based on multiple features. 
PreLnc is mainly used to distinguish long non-coding transcripts from protein-coding transcripts, supporting multiple species of plants and animals.  

---
### Start
Clone this repository to use this tool:  
`git clone https://github.com/LeiCao97/PreLnc.git`  

---
### Requirements  
1. Python(>2.7) 
    * [Biopython](https://biopython.org)   
    * [StatsModel](http://www.statsmodels.org/stable/index.html)  
2. txCdsPredict  

--- 
### Usage
PreLnc currently provides models and documents for 6 species, as well as support for building models of the species you need. 
#### Introduction  
1. Animals--Human, Mouse, Cow
	- prelnc.py : predicting long noncoding transcripts through 21 features.  
	- createmodel.py : create prediction model. 
2. Plants--Arabidopsis thaliana, Oryza sativa, Zea mays  
	- prelnc\_plant.py : predicting long noncoding transcripts through 10 features.
	- createmodel\_plant.py :  create prediction model.

#### Example  
Predicting lncRNAs of plants and animals requires different script files (see **Introduction**).  
***For human and animals:***  
`prelnc.py -i transcripts.fa -m human.pkl -r human_hexamer.csv -o prediction.txt`

If you need to model another species, you need to create an hexamer file through a script in CPAT called "make\_hexamer\_tab.py".  
1. prepare the CDS sequence and long non-coding transcripts. All files are in fasta format.   
`make_hexamer_tab.py -c Human_coding_transcripts_CDS.fa   -n Human_noncoding_transcripts_RNA.fa >Human_Hexamer.tsv`
2. Input high-confident long non-coding transcripts and protein-coding transcripts to build the model (pkl format).
`createmodel.py -l human_lnc.fa -p human_pc.fa -r human_hexamer.csv -o human_model.pkl`

***For plants:***
Use the same script as the animals, the script files used are 'prelnc\_plant.py' and 'createmodel\_plant.py'.

--- 
### Output File 
The main output file contains the feature information and prediction results of the transcript.  
***Human and animals (columns describe):***  
1. Sequence length  
2. GC content  
3. Standard deviation of stop codon counts    
4. CDS percentage   
5. CDS score of txCdsPredict prediction  
6. Fickett score  
7. Isoelectric point    
8. Open reading frame integrity  
9. Hexamer score  
10. ACG  
11. AGC  
12. CAG  
13. CAT    
14. CGG   
15. CGT  
16. GAC  
17. GAG    
18. GGC  
19. GGG  
20. TAG  
21. TCA   
22. The possibility of lncRNA:1-lncRNA 

***Plants (columns describe):***  
1. Sequence length  
2. GC content  
3. Standard deviation of stop codon counts    
4. CDS percentage   
5. CDS score of txCdsPredict prediction  
6. Fickett score     
7. Open reading frame integrity  
8. Hexamer score  
9. CTA  
10. TGG   
11. The possibility of lncRNA:1-lncRNA

---
### Statement: 
The feature extraction method in the scripts mainly comes from from CPAT[1], CPC2[2] and LncRscan-SVM[3].

---
### References:  
[1] Wang L , Park H J , Dasari S , et al. CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model[J]. Nucleic Acids Research, 2013, 41(6):e74-e74.  
[2] Kang Y J , Yang D C , Kong L , et al. CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features[J]. Nucleic Acids Research, 2017.  
[3] Sun L , Liu H , Zhang L , et al. lncRScan-SVM: A Tool for Predicting Long Non-Coding RNAs Using Support Vector Machine[J]. PLOS ONE, 2015, 10.





