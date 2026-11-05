from __future__ import absolute_import

import re

import numpy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import ProtParam

from prelnc.config import START_CODONS, STOP_CODONS
from prelnc.legacy.fickett import Fickett
from prelnc.legacy.frame_kmer import kmer_ratio
from prelnc.legacy.orf import FindCDS, ORFFinder


U_PATTERN = re.compile("U", re.I)
AMBIGUOUS_AA_PATTERN = re.compile("X|B|Z|J|U", re.I)
FICKETT = Fickett()


def gc_content(seq):
    sequence = str(seq).upper()
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count * 100.0 / float(len(sequence))


def frame_translation(seq, shift, genetic_code=1):
    from Bio.Seq import translate

    length = len(seq)
    fragment_length = 3 * ((length - shift) // 3)
    return translate(seq[shift:shift + fragment_length], genetic_code)


def transcript_id(record):
    return record.id.split("|")[0]


def tri_nucleotide_percent(seq, mer):
    seq_str = str(seq).upper()
    if len(seq_str) <= 2:
        return 0.0
    return seq_str.count(mer) * 100.0 / float(len(seq_str) - 2)


def build_kmer_features(seq, kmers):
    return dict((kmer, tri_nucleotide_percent(seq, kmer)) for kmer in kmers)


def extract_hexamer_score(seq, coding, noncoding):
    start_coden = START_CODONS.strip().split(",")
    stop_coden = STOP_CODONS.strip().split(",")
    m_rna_seq = str(seq).upper()
    tmp = ORFFinder(m_rna_seq)
    _, _, cds_seq = tmp.longest_orf(direction="+", start_coden=start_coden, stop_coden=stop_coden)
    return kmer_ratio(cds_seq, 6, 3, coding, noncoding)


def sequence_summary(seq):
    frame_0 = frame_translation(seq, 0, genetic_code=1)
    frame_1 = frame_translation(seq, 1, genetic_code=1)
    frame_2 = frame_translation(seq, 2, genetic_code=1)
    stop = (frame_0.count("*"), frame_1.count("*"), frame_2.count("*"))
    std_stop = numpy.std(stop)

    seq_rna = U_PATTERN.sub("T", str(seq).strip()).upper()
    seq_cds, _, _, orf_fullness = FindCDS(seq_rna).longest_orf("+")
    seq_prot = Seq(seq_cds).translate()
    pep_len = len(seq_prot)
    new_seq_prot = AMBIGUOUS_AA_PATTERN.sub("", str(seq_prot))
    fickett_score = FICKETT.value(seq_rna)
    if pep_len > 0:
        protparam_obj = ProtParam.ProteinAnalysis(str(new_seq_prot.strip("*")))
        isoelectric_point = protparam_obj.isoelectric_point()
    else:
        orf_fullness = -1
        isoelectric_point = 0.0

    return {
        "std_stop": std_stop,
        "fickett_score": fickett_score,
        "isoelectric_point": isoelectric_point,
        "orf_integrity": orf_fullness,
        "gc_content": gc_content(seq),
        "sequence_length": len(seq),
    }


def extract_feature_row(record, cds_len, cds_score, hexamer_tables, feature_names):
    seq = record.seq
    transcript = transcript_id(record)
    seq_len = len(seq)
    cds_length = cds_len.get(transcript, 0)
    len_perc = float(cds_length) / seq_len if seq_len else 0.0
    coding, noncoding = hexamer_tables

    summary = sequence_summary(seq)
    kmer_features = build_kmer_features(seq, [name for name in feature_names if len(name) == 3 and name.isupper()])
    row = {
        "transcript_id": transcript,
        "sequence_length": summary["sequence_length"],
        "hexamer_score": extract_hexamer_score(seq, coding, noncoding),
        "cds_score": cds_score.get(transcript, 0.0),
        "std_stop": summary["std_stop"],
        "fickett_score": summary["fickett_score"],
        "isoelectric_point": summary["isoelectric_point"],
        "cds_percentage": len_perc,
        "gc_content": summary["gc_content"],
        "orf_integrity": summary["orf_integrity"],
    }
    row.update(kmer_features)
    return row


def extract_feature_rows(fasta_path, cds_len, cds_score, hexamer_tables, feature_names, progress=None):
    rows = []
    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            rows.append(extract_feature_row(record, cds_len, cds_score, hexamer_tables, feature_names))
            if progress is not None:
                progress.advance()
    return rows

