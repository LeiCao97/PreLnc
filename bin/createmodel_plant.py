#!/usr/bin/python
#coding=utf-8
#Required: Biopython, Statsmodel and txCdsPredict

from __future__ import print_function
import Bio.Alphabet
from Bio import SeqIO
from Bio.SeqUtils import GC,ProtParam
from Bio.Seq import Seq
import math
import numpy
import os
import re
import time
from cpmodule  import orf
from cpmodule  import FrameKmer
import sys, getopt
import statsmodels.api as smf
import pickle

######################### Functions #########################
class Fickett:
    '''
    calculate Fickett TESTCODE for full sequence
    NAR 10(17) 5303-531
    modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar 
    '''
    def __init__(self):
        '''new compiled Fickett look-up table'''
        self.position_parameter  = [1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,0.0]
        self.content_parameter  = [0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.19,0.17,0]
        self.position_probability = {
            "A":[0.51,0.55,0.57,0.52,0.48,0.58,0.57,0.54,0.50,0.36],
            "C":[0.29,0.44,0.55,0.49,0.52,0.60,0.60,0.56,0.51,0.38],
            "G":[0.62,0.67,0.74,0.65,0.61,0.62,0.52,0.41,0.31,0.17],
            "T":[0.51,0.60,0.69,0.64,0.62,0.67,0.58,0.48,0.39,0.24],
            }
        self.position_weight = {"A":0.062,"C":0.093,"G":0.205,"T":0.154}
        self.content_probability = {
            "A":[0.40,0.55,0.58,0.58,0.52,0.48,0.45,0.45,0.38,0.19],
            "C":[0.50,0.63,0.59,0.50,0.46,0.45,0.47,0.56,0.59,0.33],
            "G":[0.21,0.40,0.47,0.50,0.52,0.56,0.57,0.52,0.44,0.23],
            "T":[0.30,0.49,0.56,0.53,0.48,0.48,0.52,0.57,0.60,0.51]
            }
        self.content_weight = {"A":0.084,"C":0.076,"G":0.081,"T":0.055}


    def look_up_position_probability(self,value, base):
        '''
        look up positional probability by base and value
        '''
        if float(value) < 0:
            return None
        for idx,val in enumerate (self.position_parameter):
            if (float(value) >= val):
                return float(self.position_probability[base][idx]) * float(self.position_weight[base])

    def look_up_content_probability(self,value, base):
        '''
        look up content probability by base and value
        '''
        if float(value) < 0:
            return None
        for idx,val in enumerate (self.content_parameter):
            if (float(value) >= val):
                return float(self.content_probability[base][idx]) * float(self.content_weight[base])

    def fickett_value(self,dna):
        '''
        calculate Fickett value from full RNA transcript sequence
        '''
        if len(dna) < 2:
            return 0
        fickett_score=0
        dna=dna
        total_base = len(dna)
        A_content = float(dna.count("A"))/total_base
        C_content = float(dna.count("C"))/total_base
        G_content = float(dna.count("G"))/total_base
        T_content = float(dna.count("T"))/total_base

        phase_0 = dna[::3]
        phase_1 = dna[1::3]
        phase_2 = dna[2::3]
        
        phase_0_A = phase_0.count("A")
        phase_1_A = phase_1.count("A")
        phase_2_A = phase_2.count("A")
        phase_0_C = phase_0.count("C")
        phase_1_C = phase_1.count("C")
        phase_2_C = phase_2.count("C")
        phase_0_G = phase_0.count("G")
        phase_1_G = phase_1.count("G")
        phase_2_G = phase_2.count("G")
        phase_0_T = phase_0.count("T")
        phase_1_T = phase_1.count("T")
        phase_2_T = phase_2.count("T")

        A_content = float(phase_0_A + phase_1_A + phase_2_A)/total_base
        C_content = float(phase_0_C + phase_1_C + phase_2_C)/total_base
        G_content = float(phase_0_G + phase_1_G + phase_2_G)/total_base
        T_content = float(phase_0_T + phase_1_T + phase_2_T)/total_base
        A_position= numpy.max([phase_0_A,phase_1_A,phase_2_A])/(numpy.min([phase_0_A,phase_1_A,phase_2_A]) +1.0)
        C_position= numpy.max([phase_0_C,phase_1_C,phase_2_C])/(numpy.min([phase_0_C,phase_1_C,phase_2_C]) +1.0)
        G_position= numpy.max([phase_0_G,phase_1_G,phase_2_G])/(numpy.min([phase_0_G,phase_1_G,phase_2_G]) +1.0)
        T_position= numpy.max([phase_0_T,phase_1_T,phase_2_T])/(numpy.min([phase_0_T,phase_1_T,phase_2_T]) +1.0)

        fickett_score += self.look_up_content_probability(A_content,"A")
        fickett_score += self.look_up_content_probability(C_content,"C")
        fickett_score += self.look_up_content_probability(G_content,"G")
        fickett_score += self.look_up_content_probability(T_content,"T")
        
        fickett_score += self.look_up_position_probability(A_position,"A")
        fickett_score += self.look_up_position_probability(C_position,"C")
        fickett_score += self.look_up_position_probability(G_position,"G")
        fickett_score += self.look_up_position_probability(T_position,"T")
            
        return fickett_score


class FindCDS:
    '''
    Find the most like CDS in a given sequence 
    The most like CDS is the longest ORF found in the sequence
    When having same length, the upstream ORF is printed
    modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar
    '''
    def __init__(self,seq):
        self.seq = seq
        self.result = (0,0,0,0,0)
        self.longest = 0
        self.basepair = {"A":"T","T":"A","U":"A","C":"G","G":"C","N":"N","X":"X"}

    def _reversecompliment(self):
        return "".join(self.basepair[base] for base in self.seq)[::-1]

    def get_codons(self,frame_number):
        '''
        Record every nucleotide triplet and its coordinate position for input sequence in one frame
        '''
        coordinate = frame_number
        while coordinate + 3 <= len(self.seq):
            yield (self.seq[coordinate:coordinate+3], coordinate)
            coordinate += 3 
    
    def find_longest_in_one(self,myframe,direction,start_codon,stop_codon):
        '''
        find the longest ORF in one reading myframe
        '''
        triplet_got = self.get_codons(myframe)  
        starts = start_codon
        stops = stop_codon
        '''
        Extend sequence by triplet after start codon encountered
        End ORF extension when stop codon encountered
        '''
        while True:
            try: 
                codon,index = triplet_got.next()
            except StopIteration:
                break 
            if codon in starts and codon not in stops:
                '''
                find the ORF start
                '''
                orf_start = index
                end_extension = False
                while True:
                    try: 
                        codon,index = triplet_got.next()
                    except StopIteration:
                        end_extension = True
                        integrity = -1
                    if codon in stops:
                        integrity = 1
                        end_extension = True
                    if end_extension:
                        orf_end = index + 3
                        Length = (orf_end - orf_start)
                        if Length > self.longest:
                            self.longest = Length
                            self.result = [direction,orf_start,orf_end,Length,integrity]
                        if Length == self.longest and orf_start < self.result[1]:
                            '''
                            if ORFs have same length, return the one that if upstream
                            '''
                            self.result = [direction,orf_start,orf_end,Length,integrity]
                        break

    def longest_orf(self,direction,start_codon={"ATG":None}, stop_codon={"TAG":None,"TAA":None,"TGA":None}):
        return_orf = ""
        for frame in range(3):
            self.find_longest_in_one(frame,"+",start_codon,stop_codon)
        return_orf = self.seq[self.result[1]:self.result[2]][:]
        start_coordinate = self.result[1]
        strand_direction = "+"
        orf_integrity = self.result[4]
        '''
        Also check reverse chain if -r is chosen
        '''
        if direction == "-":
            self.seq = self._reversecompliment()
            for frame in range(3):
                self.find_longest_in_one(frame,"-",start_codon,stop_codon)
            if self.result[0] == "-":
                return_orf = self.seq[self.result[1]:self.result[2]][:]
                start_coordinate = self.result[1]
                strand_direction = "-"
                orf_integrity = self.result[4]
        return return_orf,start_coordinate,strand_direction,orf_integrity


## (1) farward frame translations
## result: A translated sequence will be generated
def frame_translation(seq, shift, genetic_code=1):
    from Bio.Seq import translate
    length = len(seq) 
    frame = {} 
    fragment_length = 3 * ((length-shift) // 3) 
    frame = translate(seq[shift:shift+fragment_length], genetic_code) 
    return frame


def extract_feature_from_seq(seq,stt,stp,c_tab,g_tab):
    '''extract features of sequence from fasta entry'''  
    stt_coden = stt.strip().split(',')
    stp_coden = stp.strip().split(',')
    #transtab = maketrans("ACGTNX","TGCANX")
    mRNA_seq = seq.upper()
    mRNA_size = len(seq)
    tmp = orf.ORFFinder(mRNA_seq)
    (CDS_size1, CDS_frame1, CDS_seq1) = tmp.longest_orf(direction="+",start_coden=stt_coden, stop_coden=stp_coden)
    #fickett_score1 = fickett.fickett_value(CDS_seq1)
    hexamer = FrameKmer.kmer_ratio(CDS_seq1,6,3,c_tab,g_tab)
    #return (mRNA_size, CDS_size1, fickett_score1,hexamer)
    return hexamer

def mRNA_translate(mRNA):
    return Seq(mRNA).translate()

def protein_param(putative_seqprot):
    return putative_seqprot.isoelectric_point()

def get_feature(lncFA,pcFA,heDAT,outFILE):
    f_out = open(outFILE+".feature.txt", "w")
    test_count = 0
    transcript_len = dict()
    ### get feature from lnc_fasta
    #### run txCdsPredict on the input FASTA file ####
    cmd = "txCdsPredict " + lncFA + " -anyStart tmp.cds"
    os.system(cmd)
    cds_len = dict()
    cds_score = dict()
    temp=open("tmp.cds")
    for line in temp:
        line_array = line.split()
        id_array = line_array[0].split('|')
        start=0
        end=0
        trans_id = id_array[0]
        pred_start = int(line_array[1])
        pred_end = int(line_array[2])   
        pred_len = pred_end-pred_start
        cds_len[trans_id] = pred_len
        cds_score[trans_id] = float(line_array[5])
    temp.close()
    os.system("rm tmp.cds")
    ############ end of running txCdsPredict ################
    #### extract hexamer (CPAT)####
    #build hexamer table from hexamer frequency file
    coding={}
    noncoding={}    
    start_codons='ATG'
    stop_codons='TAG,TAA,TGA'
    for line in open(heDAT):
        line = line.strip()
        fields = line.split()
        if fields[0] == 'hexamer':continue
        coding[fields[0]] = float(fields[1])
        noncoding[fields[0]] =  float(fields[2])
    ####end extract hexamer ####
    #### extract peptide_length,Fickett_score,ORF_integrity (CPC2) ####
    strand = "+"
    strinfoAmbiguous = re.compile("X|B|Z|J|U",re.I)
    ptU = re.compile("U",re.I)
    fickett_obj = Fickett()
    ####end extract peptide_length,Fickett_score,ORF_integrity (CPC2) ####
    ## read the FASAT file of transcripts
    f_in = open(lncFA)
    for record in SeqIO.parse(f_in, "fasta"):
        ID = record.id
        first = ID.split("|")
        seq = record.seq
        seq_len = len(seq)
        #ACC_mer = seq.count("ACC")*100/float(seq_len-2)
        #ATC_mer = seq.count("ATC")*100/float(seq_len-2)
        #ATG_mer = seq.count("ATG")*100/float(seq_len-2)
        #CAT_mer = seq.count("CAT")*100/float(seq_len-2)
        CTA_mer = seq.count("CTA")*100/float(seq_len-2)
        #GAT_mer = seq.count("GAT")*100/float(seq_len-2)
        #GGT_mer = seq.count("GGT")*100/float(seq_len-2)
        #GTA_mer = seq.count("GTA")*100/float(seq_len-2) 
        #GTC_mer = seq.count("GTC")*100/float(seq_len-2)
        #TAC_mer = seq.count("TAC")*100/float(seq_len-2)
        TGG_mer = seq.count("TGG")*100/float(seq_len-2)
        # Translating nucleotides to peptide sequences according to frame shift
        frame_0 = frame_translation(seq, 0, genetic_code=1)
        stop_count_0 = frame_0.count("*")
        frame_1 = frame_translation(seq, 1, genetic_code=1)
        stop_count_1 = frame_1.count("*")
        frame_2 = frame_translation(seq, 2, genetic_code=1)
        stop_count_2 = frame_2.count("*")
        stop = (stop_count_0, stop_count_1, stop_count_2)
        std_stop = numpy.std(stop)
        seqRNA = ptU.sub("T",str(seq).strip())
        seqRNA = seqRNA.upper()
        seqCDS,start_pos,orf_strand,orf_fullness = FindCDS(seqRNA).longest_orf(strand)
        '''seqCDS:longest ORF'''
        seqprot = mRNA_translate(seqCDS)
        pep_len = len(seqprot) #pep_len = len(seqprot.strip("*"))
        newseqprot = strinfoAmbiguous.sub("",str(seqprot))
        '''exclude ambiguous amio acid X, B, Z, J, Y in peptide sequence'''
        fickett_score = fickett_obj.fickett_value(seqRNA)
        protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot.strip("*")))
        if pep_len > 0:
            isoelectric_point = protein_param(protparam_obj)
        else:
            orf_fullness = -1
            isoelectric_point = 0.0   
        hexamer = extract_feature_from_seq(seq = seq, stt = start_codons,stp = stop_codons,c_tab=coding,g_tab=noncoding)  
        # print features
        #print("%s"%first[0], end='\t', file=f_out)
        print("%d"%seq_len, end='\t', file=f_out)
        print("%0.2f"%GC(seq), end='\t', file=f_out)
        print("%.6f"%std_stop, end='\t', file=f_out)
        cds = cds_len[first[0]]
        #print("%d"%cds, end='\t', file=f_out)
        len_perc = float(cds)/seq_len
        print("%0.2f"%len_perc, end='\t', file=f_out)
        print("%d"%cds_score[first[0]], end='\t', file=f_out)
        #print("%s"%str(pep_len), end='\t', file=f_out)
        print("%s"%str(fickett_score), end='\t', file=f_out)
        #print("%s"%str(isoelectric_point), end='\t', file=f_out)
        print("%s"%str(orf_fullness), end='\t', file=f_out)
        print("%s"%str(hexamer), end='\t', file=f_out)
        #print("%0.4f"%ACC_mer, end = "\t", file=f_out)
        #print("%0.4f"%ATC_mer, end = "\t", file=f_out)
        #print("%0.4f"%ATG_mer, end = "\t", file=f_out)
        #print("%0.4f"%CAT_mer, end = "\t", file=f_out)
        print("%0.4f"%CTA_mer, end = "\t", file=f_out)
        #print("%0.4f"%GAT_mer, end = "\t", file=f_out)
        #print("%0.4f"%GGT_mer, end = "\t", file=f_out)
        #print("%0.4f"%GTA_mer, end = "\t", file=f_out)
        #print("%0.4f"%GTC_mer, end = "\t", file=f_out)
        #print("%0.4f"%TAC_mer, end = "\t", file=f_out)
        print("%0.4f"%TGG_mer, end = "\t", file=f_out)
        print("%d"%1, file=f_out)
    f_in.close()
    ### get feature from pc_fasta
    #### run txCdsPredict on the input FASTA file ####
    cmd = "txCdsPredict " + pcFA + " -anyStart tmp.cds"
    os.system(cmd)
    cds_len = dict()
    cds_score = dict()
    temp=open("tmp.cds")
    for line in temp:
        line_array = line.split()
        id_array = line_array[0].split('|')
        start=0
        end=0
        trans_id = id_array[0]
        pred_start = int(line_array[1])
        pred_end = int(line_array[2])   
        pred_len = pred_end-pred_start
        cds_len[trans_id] = pred_len
        cds_score[trans_id] = float(line_array[5])
    temp.close()
    os.system("rm tmp.cds")
    ############ end of running txCdsPredict ################
    #### extract hexamer (CPAT)####
    #build hexamer table from hexamer frequency file
    coding={}
    noncoding={}    
    start_codons='ATG'
    stop_codons='TAG,TAA,TGA'
    for line in open(heDAT):
        line = line.strip()
        fields = line.split()
        if fields[0] == 'hexamer':continue
        coding[fields[0]] = float(fields[1])
        noncoding[fields[0]] =  float(fields[2])
    ####end extract hexamer ####
    #### extract peptide_length,Fickett_score,ORF_integrity (CPC2) ####
    strand = "+"
    strinfoAmbiguous = re.compile("X|B|Z|J|U",re.I)
    ptU = re.compile("U",re.I)
    fickett_obj = Fickett()
    ####end extract peptide_length,Fickett_score,ORF_integrity (CPC2) ####
    ## read the FASAT file of transcripts
    f_in = open(pcFA)
    for record in SeqIO.parse(f_in, "fasta"):
        ID = record.id
        first = ID.split("|")
        seq = record.seq
        seq_len = len(seq)
        #ACC_mer = seq.count("ACC")*100/float(seq_len-2)
        #ATC_mer = seq.count("ATC")*100/float(seq_len-2)
        #ATG_mer = seq.count("ATG")*100/float(seq_len-2)
        #CAT_mer = seq.count("CAT")*100/float(seq_len-2)
        CTA_mer = seq.count("CTA")*100/float(seq_len-2)
        #GAT_mer = seq.count("GAT")*100/float(seq_len-2)
        #GGT_mer = seq.count("GGT")*100/float(seq_len-2)
        #GTA_mer = seq.count("GTA")*100/float(seq_len-2) 
        #GTC_mer = seq.count("GTC")*100/float(seq_len-2)
        #TAC_mer = seq.count("TAC")*100/float(seq_len-2)
        TGG_mer = seq.count("TGG")*100/float(seq_len-2)
        # Translating nucleotides to peptide sequences according to frame shift
        frame_0 = frame_translation(seq, 0, genetic_code=1)
        stop_count_0 = frame_0.count("*")
        frame_1 = frame_translation(seq, 1, genetic_code=1)
        stop_count_1 = frame_1.count("*")
        frame_2 = frame_translation(seq, 2, genetic_code=1)
        stop_count_2 = frame_2.count("*")
        stop = (stop_count_0, stop_count_1, stop_count_2)
        std_stop = numpy.std(stop)
        seqRNA = ptU.sub("T",str(seq).strip())
        seqRNA = seqRNA.upper()
        seqCDS,start_pos,orf_strand,orf_fullness = FindCDS(seqRNA).longest_orf(strand)
        '''seqCDS:longest ORF'''
        seqprot = mRNA_translate(seqCDS)
        pep_len = len(seqprot) #pep_len = len(seqprot.strip("*"))
        newseqprot = strinfoAmbiguous.sub("",str(seqprot))
        '''exclude ambiguous amio acid X, B, Z, J, Y in peptide sequence'''
        fickett_score = fickett_obj.fickett_value(seqRNA)
        protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot.strip("*")))
        if pep_len > 0:
            isoelectric_point = protein_param(protparam_obj)
        else:
            orf_fullness = -1
            isoelectric_point = 0.0   
        hexamer = extract_feature_from_seq(seq = seq, stt = start_codons,stp = stop_codons,c_tab=coding,g_tab=noncoding)  
        # print features
        #print("%s"%first[0], end='\t', file=f_out)
        print("%d"%seq_len, end='\t', file=f_out)
        print("%0.2f"%GC(seq), end='\t', file=f_out)
        print("%.6f"%std_stop, end='\t', file=f_out)
        cds = cds_len[first[0]]
        #print("%d"%cds, end='\t', file=f_out)
        len_perc = float(cds)/seq_len
        print("%0.2f"%len_perc, end='\t', file=f_out)
        print("%d"%cds_score[first[0]], end='\t', file=f_out)
        #print("%s"%str(pep_len), end='\t', file=f_out)
        print("%s"%str(fickett_score), end='\t', file=f_out)
        #print("%s"%str(isoelectric_point), end='\t', file=f_out)
        print("%s"%str(orf_fullness), end='\t', file=f_out)
        print("%s"%str(hexamer), end='\t', file=f_out)
        #print("%0.4f"%ACC_mer, end = "\t", file=f_out)
        #print("%0.4f"%ATC_mer, end = "\t", file=f_out)
        #print("%0.4f"%ATG_mer, end = "\t", file=f_out)
        #print("%0.4f"%CAT_mer, end = "\t", file=f_out)
        print("%0.4f"%CTA_mer, end = "\t", file=f_out)
        #print("%0.4f"%GAT_mer, end = "\t", file=f_out)
        #print("%0.4f"%GGT_mer, end = "\t", file=f_out)
        #print("%0.4f"%GTA_mer, end = "\t", file=f_out)
        #print("%0.4f"%GTC_mer, end = "\t", file=f_out)
        #print("%0.4f"%TAC_mer, end = "\t", file=f_out)
        print("%0.4f"%TGG_mer, end = "\t", file=f_out)
        print("%d"%0, file=f_out)
    f_in.close()
    f_out.close()

#############################################################
def main():
    lnc_fasta = ''
    pc_fasta = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"h:l:p:r:o:",["Help","lncFASTA=","pcFASTA=","o=outModel"])
    except getopt.GetoptError:
        print("Usage: %s  -l lncTranscripts.fa -p pcTranscripts.fa -r hexamer.csv -o output.pkl" % sys.argv[0])
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("createmodel_plant.py -l/--lncFASTA <FASTA> -p/--pcFASTA <FASTA> -r/--hexamerfile <tsv> -o/--outmodel <pkl>")
            print("Example: createmodel_plant.py -l ath_lnc.fa -p ath_pc.fa -r ath_hexamer.csv -o ath_model.pkl")
            sys.exit()
        elif opt in ("-l", "--lncFASTA"):
            lnc_fasta = arg
        elif opt in ("-p", "--pcFASTA"):
            pc_fasta = arg
        elif opt in ("-r", "--hexamerfile"):
            hexamer_dat = arg
        elif opt in ("-o", "--outmodel"):
            outputfile = arg

################################### main body #####################################
    start_time = time.time()
    # get features
    get_feature(lnc_fasta,pc_fasta,hexamer_dat,outputfile)
    # create model
    fout= open(outputfile, "w")
    data= numpy.loadtxt(outputfile+".feature.txt")
    x=data[:,0:10]
    y=data[:,-1]
    x = smf.add_constant(x,prepend = False,has_constant='add')
    model = smf.Logit(y,x).fit()
    with open(outputfile, 'wb') as f:
        pickle.dump(model, f)
    f.close()
    sys.stderr.write("\n[PreLnc]Model built cost time: %ds\n"%(time.time()-start_time))

if __name__ == '__main__':
    main()
