#!/usr/bin/env python

import sys
import re
import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq


REPAM = re.compile('[A-Z]GG',re.IGNORECASE)
BfuAI5 = Seq('GCGG')
BfuAI3 = Seq('GGAT')
BspQI5 = Seq('AAG')
BspQI3 = Seq('CCG')




def find_pams(seq_rec):
    pams = []
    inseq = str(seq_rec.seq)
    pamiter = REPAM.finditer(inseq)
    for m in pamiter:
        p1 = m.span(0)
        p2 = m.span(1)
        pams.append((p1, p2, inseq[p1:p2]))
    return pams



def pad_to_orf(seq_len, arm_p2, orf_start=0):
    real_len = seq_len + orf_start - arm_p2
    pad_len = 3 - real_len % 3
    if pad_len == 3:
        pad_len = 0
    return pad_len

def find_upstream_arm_sense(genomic_dna, guide_rna, guide_rna_span, armlen = 48, orf_start=0):

    p1, p2 = guide_rna_span
    # Page 20, Figure 20
    cas9_site =  p2 - 3
    arm_p1 = cas9_site - armlen
    arm_p2 = cas9_site
    upstream_arm = genomic_dna.seq[arm_p1:arm_p2]
    print(f"upstream arm: {upstream_arm}")
    # Found upstream arm, now find spacer
    # Spacer is a 3mer whose base is not found in -51 to -49
    # Page 16, item #2
    spacer_base = find_spacer_base(genomic_dna.seq[arm_p1-3:arm_p1])
    print(f"bla bla: {genomic_dna.seq[arm_p1-3:p2]}")
    assert spacer_base, f'No spacer_base found in {genomic_dna.seq[arm_p1-3:arm_p1]}'
    spacer_base *= 3
    upstream_arm = spacer_base + upstream_arm


    # Padding to nearest codon (page 17 iterm #3). This is to ensure
    # integration is in-frame
    # Need to add padding
    # Determine padding length: 0,1, or 2 bases
    pad_len = pad_to_orf(len(genomic_dna.seq), arm_p2, orf_start)
    print(genomic_dna.seq[orf_start:arm_p2])
    print (f"pad_len {pad_len}")
    if pad_len > 0:
        pad_p1, pad_p2 = cas9_site, cas9_site + pad_len
        pad_seq = genomic_dna.seq[pad_p1:pad_p2].complement()
        print (f"padding sequence {pad_seq}")






    return upstream_arm

def find_spacer_base(inseq):
    print(f"inseq {str(inseq).upper()}")
    spacer_base = None
    for i in ['A','T','G','C']:
        if not (i in str(inseq).upper()):
            spacer_base = i
            break
    return spacer_base 

guide_rna = Seq("GCCACAGCATGTTTGGGCAT")
                  
guide_rna_2 = Seq("CCTGTTTGATCACCCGTTAA").reverse_complement()

def get_genomic_dna(sequence_id, source, entrez_email="idoerg@iastate.edu",
        entrez_db="nucleotide", file_format="fasta",gdna_file=None,verbose=0):
    # Source can be: local file, EnsEMBL ID, NCBI ID
    if verbose:
        print("Getting {sequence_id} from {source}")
    if source.upper() == "ENSEMBL":
        ensRest = EnsemblRest()
        genomic_dna_ensembl = ensRest.getSequenceById(id=sequence_id)
        genomic_dna = SeqRecord(id=genomic_dna_ensembl['id'],
                      seq=genomic_dna_ensembl['seq'])
    elif source.upper() == "NCBI" or source.upper() == "ENTREZ":
        assert entrez_email and entrez_db, "Your email and entrez db required"
        Entrez.email = entrez_email
        with (Entrez.efetch(id=sequence_id, db=entrez_db,
            rettype="fasta", retmode="text")) as handle:

            genomic_dna = SeqIO.read(handle,"fasta")
    elif source.upper() == "LOCAL":
        genomic_dna = next(SeqIO.parse(gdna_file, "fasta"))
    else:
        raise ValueError(f"Bad sequence source {source}")
    return genomic_dna


def find_guide_sequence(genomic_dna, guide_rna=guide_rna, strand=1):
    guide_rna = guide_rna.upper()
    genomic_dna.seq = genomic_dna.seq.upper()
    print(type(guide_rna), type(genomic_dna))
    assert len(guide_rna) == 20, f'guide RNA must be 20bp' 
    # Look for the guide RNA on the sense strand
    if strand == 1:
        guide_rna_location = re.search(str(guide_rna), str(genomic_dna.seq))
    else:
        guide_rna_location = re.search(str(guide_rna.reverse_complement()),
            str(genomic_dna))
    if not guide_rna_location:
        raise ValueError(f'Guide RNA {guide_rna} not in genomic DNA')
    guide_rna_span = guide_rna_location.span()

    print(guide_rna, guide_rna_span)
    return (guide_rna, guide_rna_span)
    
def create_parser():
    p = argparse.ArgumentParser(description='Process arguments for pygtag')
    p.add_argument('-q', '--seqid')
    p.add_argument('-u', '--source', default='local')
    p.add_argument('-f', '--file', default=None)
    # p.add_argument('-u', '--source',choices=['ensembl', 'ncbi', 'local'], default='local')
    p.add_argument('-s', '--strand', default=1)
    p.add_argument('--verbose', '-v', action='count', default=0)

    return p

if __name__ == '__main__':
    p = create_parser()
    args = p.parse_args()
    print(args)

    genomic_dna = get_genomic_dna(args.seqid, args.source, gdna_file=args.file,verbose=args.verbose)
    guide_rna, guide_rna_span = find_guide_sequence(genomic_dna, strand=args.strand)
    upstream_arm = find_upstream_arm_sense(genomic_dna, guide_rna,guide_rna_span)
    print(upstream_arm)
    

