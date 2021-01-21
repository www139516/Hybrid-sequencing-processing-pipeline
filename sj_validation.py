# The codes follow is used for validation of splicing junctions.

# functions.py

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

# splicing_junctions.py
# import pandas as pd 


class SplicingJunction():

    def __init__(self, sj_star_file):

        # initializing the SplicingJunction class, import the location of SJs calculateds using STAR, and set the name of each column 
        # 1	44948	45665	1	1	0	1	0	66
        # 1	45804	45887	1	1	0	1	0	31
        # 1	46343	46450	1	1	0	3	0	47
        # 1	46634	47044	1	1	0	3	0	67
        # 1	47263	47649	1	1	0	8	0	56

        self.sj_file = pd.read_csv(sj_star_file, sep='\t', names=[
            'chr', 'start', 'end', 'strand', 'in_motif', 'annotated', 'num_uni_cross', 'num_multi_cross', 'max_overhang'])


    def get_sj_loc(self, num_cross=0):
        # add a column named location, store the location of every SJ
        # 
        fil_sj = self.sj_file.loc[self.sj_file.num_uni_cross > num_cross, :]
        fil_sj['location'] = self.sj_file.chr.astype(str) + '_' + self.sj_file.start.astype(str) + '_' + self.sj_file.end.astype(str)
        sj_loc = fil_sj.loc[:, ['location', 'num_uni_cross']]
        return sj_loc
        # print(sj_loc.head())
        # print(sj_loc.dtypes)


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


def sj_validatiing(short_sj, long_sj):
    
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
                

# sj_validation.py

# coding=utf-8

import pandas as pd
from argparse import ArgumentParser

# from splicing_junction import SplicingJunction
# import functions as fun





def main(args):
    # input the splicing_junction file which is produced by star
    with open(args.sj) as star_sj_file:
        short_sj = SplicingJunction (star_sj_file)
        short_sj_loc = short_sj.get_sj_loc()

    with open(args.gtf) as gtf_file:
        sj_long = get_sj_from_gtf(gtf_file)
        sj_validation = sj_validatiing(short_sj_loc, sj_long)

    sj_validation = pd.DataFrame({'id': list(sj_long.keys()), 'validated': sj_validation})

    # print(sj_validation.head())
    # print(sj_validation.validated.describe(include='validated'))

    with open(args.seq) as fa:
        get_the_seqs(fa, sj_validation, args.o, args.p)





if __name__ == "__main__":
    
     
    parser = ArgumentParser(
        description='Validate the splicing junctions identified by long-read sequencing results')
    parser.add_argument('-sj', required=True,
                        help='The splicing junctions file generated by STAR')
    parser.add_argument('-gtf', required=True,
                        help='The gff file generated by collapse_isoforms_by_sam.py')
    parser.add_argument('-seq', required=True,
                        help='The fasta file generated by collapse_isoforms_by_sam.py')
    parser.add_argument(
        '-o', required=True, help='The folder of the output files')
    parser.add_argument('-p', required=True, help='The profix of output files')
    args = parser.parse_args()

    main(args)


