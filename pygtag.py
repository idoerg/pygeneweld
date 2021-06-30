
import re
from Bio import SeqIO
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

def find_upstream_arm_sense(seq_rec, pam, armlen = 48):

    # Page 20, Figure 20
    p1, p2, pamstr = pam
    cas9_site = p1 - 3
    assert cas9_site - armlen > 0, f'{p2} is the pam position too close to 5` end'
    arm_p1 = cas9_site - armlen
    upstream_arm = seq_rec[arm_p1:cas9_site]

    # Found upstream arm, now find spacer
    # Spacer is a 3mer whose base is not found in -51 to -49
    # Page 16, item #2
    spacer_base = find_spacer_base(seq_rec[arm_p1-3:arm_p1])
    assert spacer_base, f'No spacer_base found in {seq_rec[arm_p1-3:arm_p1]}'
    spacer_base *= 3
    upstream_arm = spacer_base + upstream_arm


    # Padding to nearest codon (page 17 iterm #3). This is to ensure
    # integration is in-frame


    return upstream_arm

def find_spacer_base(inseq):
    spacer_base = None
    for i in ['A','T','G','C']:
        if not (i in upper(str(inseq)):
            spacer_base = i
            break
    return spacer_base 

guide_rna=Seq("GCCACAGCATGTTTGGGCAT")

def read_sequence(gdna_file, guide_rna):
    genomic_dna = next(SeqIO.parse(gdna_file, "fasta"))
    guide_rna = guide_rna.upper()
    genomic_dna_seq = genomic_dna.seq.upper()
    assert len(guide_rna) == 20, f'guide RNA must be 20bp' 
    guide_rna_location = re.search(str(guide_rna), str(genomic_dna.seq))
    assert guide_rna_location, f'Guide RNA not in genomic DNA'
    guide_rna_span = guide_rna_location.span()

    return guide_rna, genomic_dna_seq, guide_rna_span
    

