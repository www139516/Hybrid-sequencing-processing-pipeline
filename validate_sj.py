"""
Description: This class validates the splicing junctions derived from long reads using SJ generated based on short reads using STAR.
  (FLNC CCSs).
Output: FLNC-validated sequences
Input: FLNC_corrected reads
Other notes: sj fpath generated by STAR
"""
# The codes follow is used for validation of splicing junctions.

# functions.py

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO



def get_sj_from_gtf(gtf_file):
    # Get the splcing junctions using gtf file obtained by long-read sequencing
   
    # 1	gmapZmV4	exon	50866	51217	99.00	-	.	transcript_id "ND-iso-new-7_5p-d_HQ_transcript/11010"; gene_id "ND-iso-new-7_5p-d_HQ_transcript/11010.gene"; gene_name "ND-iso-new-7_5p-d_HQ_transcript/11010";
    # 1	gmapZmV4	exon	51371	51436	100.00	-	.	transcript_id "ND-iso-new-7_5p-d_HQ_transcript/11010"; gene_id "ND-iso-new-7_5p-d_HQ_transcript/11010.gene"; gene_name "ND-iso-new-7_5p-d_HQ_transcript/11010";
    # 1	gmapZmV4	exon	51886	52004	100.00	-	.	transcript_id "ND-iso-new-7_5p-d_HQ_transcript/11010"; gene_id "ND-iso-new-7_5p-d_HQ_transcript/11010.gene"; gene_name "ND-iso-new-7_5p-d_HQ_transcript/11010";
    # 1	gmapZmV4	exon	52137	52294	100.00	-	.	transcript_id "ND-iso-new-7_5p-d_HQ_transcript/11010"; gene_id "ND-iso-new-7_5p-d_HQ_transcript/11010.gene"; gene_name "ND-iso-new-7_5p-d_HQ_transcript/11010";
    # 1	gmapZmV4	exon	52391	52546	100.00	-	.	transcript_id "ND-iso-new-7_5p-d_HQ_transcript/11010"; gene_id "ND-iso-new-7_5p-d_HQ_transcript/11010.gene"; gene_name "ND-iso-new-7_5p-d_HQ_transcript/11010";
    # 1	gmapZmV4	exon	52668	52826	100.00	-	.	transcript_id "ND-iso-new-7_5p-d_HQ_transcript/11010"; gene_id "ND-iso-new-7_5p-d_HQ_transcript/11010.gene"; gene_name "ND-iso-new-7_5p-d_HQ_transcript/11010";

# 1	PacBio	transcript	45948	49822	.	+	.	gene_id "PB.1"; transcript_id "PB.1.1";
# 1	PacBio	exon	45948	46148	.	+	.	gene_id "PB.1"; transcript_id "PB.1.1";
# 1	PacBio	exon	46233	46342	.	+	.	gene_id "PB.1"; transcript_id "PB.1.1";
# 1	PacBio	exon	46451	46633	.	+	.	gene_id "PB.1"; transcript_id "PB.1.1";
# 1	PacBio	exon	47045	47262	.	+	.	gene_id "PB.1"; transcript_id "PB.1.1";
# 1	PacBio	exon	47650	48111	.	+	.	gene_id "PB.1"; transcript_id "PB.1.1";
# 1	PacBio	exon	48200	49247	.	+	.	gene_id "PB.1"; transcript_id "PB.1.1";
# 1	PacBio	exon	49330	49822	.	+	.	gene_id "PB.1"; transcript_id "PB.1.1";
# 1	PacBio	transcript	50935	55715	.	-	.	gene_id "PB.2"; transcript_id "PB.2.1";
# 1	PacBio	exon	50935	51217	.	-	.	gene_id "PB.2"; transcript_id "PB.2.1";
# 1	PacBio	exon	51371	51436	.	-	.	gene_id "PB.2"; transcript_id "PB.2.1";
# 1	PacBio	exon	51886	52004	.	-	.	gene_id "PB.2"; transcript_id "PB.2.1";
# 1	PacBio	exon	52137	52294	.	-	.	gene_id "PB.2"; transcript_id "PB.2.1";
# 1	PacBio	exon	52391	52546	.	-	.	gene_id "PB.2"; transcript_id "PB.2.1";
# 1	PacBio	exon	52668	52826	.	-	.	gene_id "PB.2"; transcript_id "PB.2.1";
# 1	PacBio	exon	52947	53149	.	-	.	gene_id "PB.2"; transcript_id "PB.2.1";
# 1	PacBio	exon	53623	53932	.	-	.	gene_id "PB.2"; transcript_id "PB.2.1";
# 1	PacBio	exon	55222	55715	.	-	.	gene_id "PB.2"; transcript_id "PB.2.1";
# 1	PacBio	transcript	209983	212765	.	-	.	gene_id "PB.3"; transcript_id "PB.3.1";
# 1	PacBio	exon	209983	210382	.	-	.	gene_id "PB.3"; transcript_id "PB.3.1";
# 1	PacBio	exon	210488	210615	.	-	.	gene_id "PB.3"; transcript_id "PB.3.1";
# 1	PacBio	exon	210696	210783	.	-	.	gene_id "PB.3"; transcript_id "PB.3.1";
# 1	PacBio	exon	211005	211143	.	-	.	gene_id "PB.3"; transcript_id "PB.3.1";
# 1	PacBio	exon	211234	211410	.	-	.	gene_id "PB.3"; transcript_id "PB.3.1";
# 1	PacBio	exon	211484	211730	.	-	.	gene_id "PB.3"; transcript_id "PB.3.1";
# 1	PacBio	exon	212055	212164	.	-	.	gene_id "PB.3"; transcript_id "PB.3.1";
# 1	PacBio	exon	212503	212765	.	-	.	gene_id "PB.3"; transcript_id "PB.3.1";


    df_gtf = pd.read_csv(gtf_file, sep='\t', names=['chr', 'genome', 'transcript_class', 'start', 'end', 'score', 'strand', 'dot', 'annotation'])
    # get the id of each isoform
    df_gtf['id'] = df_gtf.annotation.str.split('"').apply(lambda x: x[-2])
    # using a dict to save the start/end points of every SJs, including the start/end points of the chromsome
    sj_long_loc = {} # dictionary, key: isoform id, value: list 
    sj_isoform_loc = [] # temporary vaiable, stare the location of isoform

    # check every id, get the location of SJ if the id of the exon is identical to the id of the current one
    # 
    for i in range(len(df_gtf)-1):
        j = i + 1
        if (df_gtf.id[j] == df_gtf.id[i]) & (str(df_gtf.transcript_class[i]) == 'exon'):
            sj_start = df_gtf.end[i] + 1
            sj_end = df_gtf.start[j] - 1
            sj_isoform_loc.append(str(df_gtf.chr[i]) + '_' +
                                  str(sj_start) + '_' + str(sj_end))
        else:
            sj_long_loc[df_gtf.id[i]] = sj_isoform_loc
            sj_isoform_loc = []

    return sj_long_loc # Returen a dict, store the location of every SJ 


def __sj_validation(short_sj, long_sj):
    
    # check every long_reads' SJ, validate the isoform if only all SJs could be detected in the SJs obtained by short-read sequencing
    sj_loc = list(short_sj.loc[:, 'location'])
    long_sj = long_sj


    # iterate every value of SJ in the long reads
    flags = []
    for value in long_sj.values():
        flag = True # set a flag, its value is true only if all values are true, false otherwise. 
        for i in range(len(value)):
            if value[i] not in sj_loc:
                flag = False
                break

        if flag:
            flags.append(True)
        else:
            flags.append(False)

    return flags

def get_the_seqs(seq_fa, sj_validated, out_dir, prefix):
    # Get the validated seqs id, and extract the sequences from sequence file
    # >PB.1.1|1:45948-49822(+)|m54286_190414_205743/55771290/ccs m54286_190414_205743/55771290/ccs

    sj_validated = sj_validated.loc[sj_validated.validated == True, :]
    id_validated = list(sj_validated.id)
    print(id_validated[:10])
    with open(out_dir + prefix + '.fa', 'w') as fa:
        for seq_record in SeqIO.parse(seq_fa, "fasta"):
            if seq_record.id.split('|')[0] in id_validated:
                fa.write('>{}\n'.format(seq_record.id))
                fa.write(str(seq_record.seq) + '\n')
                
