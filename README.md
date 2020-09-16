## An accurate tool for predicting lncRNAs based on multiple features. [[Paper]](https://www.mdpi.com/2073-4425/11/9/981/htm)
PreLnc is mainly used to distinguish long non-coding transcripts from protein-coding transcripts, supporting multiple species of plants and animals.   

---
### Requirements (Ubuntu 16.04.6) 
1. Python(2.7) 
    * [Biopython](https://biopython.org)   
    * [sklearn](https://scikit-learn.org)  
2. txCdsPredict   

---
### Install
1. Clone this repository to use this tool:  
`git clone https://github.com/LeiCao97/PreLnc.git`  
2. Import python package:  
`pip install biopython`   
`pip install sklearn`    
3. Configure txCdsPredict paths by modifying ~/.bashrc:  
`vi ~/.bashrc`   
Please add:  
 `export PATH=$PATH:'yourpath'/Prelnc/bin/`  
4. Update profile  
`source ~/.bashrc` 

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
`cd PreLnc-master/example`   
`python ../bin/prelnc.py -i example.fa -m ../models/human_model.pkl -r ../models/human_hexamer.tsv -o example.txt`

If you need to model another species, you need to create an hexamer file through a script in CPAT called "make\_hexamer\_tab.py".  
1. prepare the CDS sequence and long non-coding transcripts. All files are in fasta format.   
`python make_hexamer_tab.py -c Human_coding_transcripts_CDS.fa   -n Human_noncoding_transcripts_RNA.fa >Human_Hexamer.tsv`
2. Input high-confident long non-coding transcripts and protein-coding transcripts to build the model (pkl format).
`python createmodel.py -l human_lnc.fa -p human_pc.fa -r human_hexamer.tsv -o human_model.pkl`

***For plants:***
Use the same script as the animals, the script files used are 'prelnc\_plant.py' and 'createmodel\_plant.py'.

--- 
### Output File 
The main output file contains the feature information and prediction results of the transcript.  
***Human and animals (columns describe):***  
1. CDS percentage  
2. Fickett score  
3. Hexamer score.  
4. CDS score of txCdsPredict prediction  
5. CGG   
6. TAG   
7. GC content  
8. Standard deviation of stop codon counts   
9. Isoelectric point   
10. ACG  
11. GGC  
12. Sequence length  
13. CGT  
14. AGC  
15. GAC  
16. GGG  
17. TCA  
18. CAT  
19. Open reading frame integrity   
20. CAG   
21. The possibility of lncRNA:1-lncRNA   

***Plants (columns describe):***  
1. Sequence length  
2. Hexamer score  
3. CDS score of txCdsPredict prediction   
4. Standard deviation of stop codon counts   
5. Fickett score   
6. Isoelectric point          
7. GGG   
8. TAG   
9. CGA   
10. CAA   
11. CDS percentage   
12. GAT   
13. ACT   
14. GC content    
15. ATC   
16. Open reading frame integrity   
17. AAA   
18. AAC   
19. CTA   
20. AGG   
21. GAG   
22. AAT   
23. GAC   
24. CAG   
25. GGA   
26. CAT   
27. AGC   
28. TAT   
29. TAA   
30. ACC   
31. GTA   
32. ATG   
33. TTA       
34. The possibility of lncRNA:1-lncRNA   

---
### Statement: 
The feature extraction method in the scripts mainly comes from from CPAT[1], CPC2[2] and LncRScan-SVM[3].  
The cpmodule package is mainly derived from CPAT and used to calculate Hexamer score.   
The Compiled 'txCdsPredict' is derived from lncRScan-SVM.  

---
### References:  
[1] Wang L , Park H J , Dasari S , et al. CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model[J]. Nucleic Acids Research, 2013, 41(6):e74-e74.  
[2] Kang Y J , Yang D C , Kong L , et al. CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features[J]. Nucleic Acids Research, 2017.  
[3] Sun L , Liu H , Zhang L , et al. lncRScan-SVM: A Tool for Predicting Long Non-Coding RNAs Using Support Vector Machine[J]. PLOS ONE, 2015, 10.





