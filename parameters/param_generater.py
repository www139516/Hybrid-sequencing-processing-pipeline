'''
This class is used to generate the parameter file that is used in the program. The file contains all the information that is customed by the user. 
Some of the parameters are set with default values.
Once this file is set, users do not need to provide other infomation, and the program can accomplish all the tasks automatically.
'''

import os



class ParamGenerator:

    def __init__(self):

        # the name of the project, all the files generated from this prgram are stored in the directory under this project.
        self._proj_name = ''
        # The path of the reference genome file
        self._ref_fpath = ''
        # The path of the directory where you want to store the results.
        self._output_dpath = ''
        # The name of the parameter file
        self._param_fname = ''
        # The path of the parameter file
        self._param_fpath = ''

    def fit(self, proj_name, ref_fpath, output_dpath):
        self._proj_name = proj_name
        self._ref_fpath = ref_fpath

        if not output_dpath:
            self._output_dpath = os.path.join(os.getcwd(), 'proj_{}'.format(self._proj_name))
        else:
            self._output_dpath = os.path.join(
                os.path.abspath(output_dpath), 'proj_{}'.format(self._proj_name))
        
        if not os.path.exists(self._output_dpath):
            os.mkdir(self._output_dpath)
        else:
            print('The directory for this project is already exist.')

        self._param_fname = 'param_' + self._proj_name + '.ini'
        self._param_fpath = os.path.join(self._output_dpath, self._param_fname)
        return self

    def generate_param_file(self):
        if os.path.exists(self._param_fpath):
            print('The parameter files is already exist.')
            print('Reading the parameter files....')
        else:
            param_file = open(self._param_fpath, 'a')

            param_file.writelines(
                '; This is the config file containing all the parameters used in this pipeline. \n')
            param_file.writelines('\n\n')

            param_file.writelines('[common_parameters]\n')
            param_file.writelines('; The name of this project. \n proj_name = {} \n'.format(self._proj_name))
            param_file.writelines('; The path of the reference assembly file used in this project. \n ref_assembly = {} \n'.format(self._ref_fpath))
            param_file.writelines('; The path of the directory for the output results. \n output_dir = {} \n'.format(self._output_dpath))
            param_file.writelines('; The threads used in this program, default value is 1. \n threads = 1 \n')
            param_file.writelines('; The directory path of the collapse_isoforms_by_sam.py. \n cibs_dpath = ')
            param_file.writelines('\n\n')

            # ccs movieX.subreads.bam movieX.ccs.bam --noPolish --minPasses 1 
            # lima --isoseq --dump-clips --no-pbi [--peek-guess] -j 24 ccs.bam primers.fasta demux.bam
            # isoseq3 refine --require-polya combined_demux.consensusreadset.xml primers.fasta unpolished.flnc.bam
            param_file.writelines('[step_1_isoseq_parameters] \n')
            param_file.writelines('; input: The path of the input subreads file. \n i_subread_fpath = \n')
            param_file.writelines('; input: The path of the primers.fasta file. \n i_primer_fpath = \n')
            # param_file.writelines(
            #     '; input: the path of the collapse_isoforms_by_sam.py. \n i_collapse_dpath = \n')
            param_file.writelines(
                '; param: Output the initial template derived from the POA (faster but less accurate), default value set in this programe is True \n p_no_polish = True \n')
            param_file.writelines(
                '; param: Minimum number of subreads required to generate CCS, default value set in this program is 1. \n p_min_passes = 1 \n')
            param_file.writelines(
                '; param: set to True if your transcripts have a polyA tail, default value in this program is True. \n p_polya = True \n')
            param_file.writelines('; output: The output of this step, automatically set if you do not set it in this config file. \n o_flnc_out = auto \n')
            param_file.writelines('\n\n')

            param_file.writelines('[step_2_fqtrim_parameters] \n')
            """
                fqtrim -A -l25 -o trimmed.fq.gz exome_reads_1.fastq.gz,exome_reads_2.fastq.gz
            """
            param_file.writelines('; input: The path of the short read file, fastq or fastq.gz format, using comma to separate two paired files if needed. \n i_short_read_file = \n')
            # param_file.writelines('; input: The score system used in short read sequencing. phred33 (default) or phred64. \n phred = \n')
            param_file.writelines('; param: minimal lenth retained after trimming, default value is 40. \n p_lenth = 40 \n')
            param_file.writelines('; param: genomic sequencing (g, means do not trim poly A) or RNA-sequencing (r, default, trim poly A). \n p_sequencing_type = r \n')
            param_file.writelines('; param: quality value for trimming, default is 20. \n p_qvalue = 20 \n')
            param_file.writelines(
                '; onput: The output suffix of the trimmed file. \n o_suffix = trimmed.fq.gz \n')
            param_file.writelines('; output: the path of the dir where you put the output trimmed files. \n o_trimmed_dir = {}/trimmed_files/ \n'.format(self._output_dpath))
            param_file.writelines('\n\n')

            param_file.writelines('[step_3_correction_parameters] \n')
            """
                lordec-correct

                -i|--long_reads <long read FASTA/Q file>
                -2|--short_reads <short read FASTA/Q file(s)>
                -k|--kmer_len <k-mer size>
                -o|--corrected_read_file <output reads file>
                -s|solid_threshold <solid k-mer abundance threshold>
                [-t|--trials <number of paths to try from a k-mer>]
                [-b|--branch <maximum number of branches to explore>]
                [-e|--errorrate <maximum error rate>]
                [-T|--threads <number of threads>]
                [-S|--stat_file <out statistics file>]
                [-c|--complete_search]
                [-a|--abundance-max <abundance max threshold for k-mers>]
                [-O|--out-tmp <GATB graph creation temporary files directory>]
                [-p|--progress]
                [-g|--graph_named_like_output]

            """
            param_file.writelines('; program: using lordec (l, default) or proovread (p). \n program_for_correction = l \n')
            param_file.writelines(
                '; input: The path of the short read file, fastq or fastq.gz format, using comma to separate two paired files if needed. \n i_short_read_fpath = auto \n')
            param_file.writelines('; input: The path of the long read file. \n i_long_read_fpath = auto \n')
            param_file.writelines('; param: kmer length, the recommend values in the lordec website are 17-19 for small genome, 21 for large genome (>2g, default). \n p_kmer_len = 21 \n')
            param_file.writelines('; param: solid threshold, default value is 3. \n p_solid_threshold = 3 \n ')
            param_file.writelines('; param: abundance max, default value is 10000. \n p_abundance_max = 10000 \n')
            param_file.writelines(
                '; param:  the maximum error rate of corrected regions (default is 0.4). \n errorrate = 0.4 \n')
            param_file.writelines('; output: corrected reads file, . \n o_corrected_seq_fpath = auto \n')
            param_file.writelines('; output: stat file. \n o_stat_fpath = auto \n')
            param_file.writelines('; additional prameters for proovread. \n')
            """
                '{path_of_proovread} -l {long_read_seq_fpath} -s {short_read_seq_fpaths} -p {pre} > proovread.log'
            """
            param_file.writelines(
                '; param: Estimated short read coverage, default value is 50X. \n p_prov_coverage = 50 \n')
            # param_file.writelines('; param: use trimmed file (true, default) or untrimmed corrected file (false) for the further process. \n p_trimmed_corrected_file = true')
            param_file.writelines('; output: prefix of output files. \n o_prov_prefix = auto \n\n')

            param_file.writelines('[step_4_alignment_parameters] \n')
            param_file.writelines('; Choose the program for alignment, minimap2 (m, default), or (gmap, g). \n program_for_alignment = m \n')

            """
                minimap2 -d ref.mmi ref.fa # index ref
                minimap2 -ax splice:hq -uf ref.fa iso-seq.fq > aln.sam
            """
            param_file.writelines('; The settings for minimap2 or gmap. \n')
            param_file.writelines('; input: the path of the corrected long read sequencing file. \n i_long_read_fpath = auto \n')
            param_file.writelines('; parameter: DNA (d) or RNA reads (r, default). \n p_dna_or_rna = r \n')
            param_file.writelines('; output: The output sam file generated by minimap2. \n o_out_mm2_file = auto \n')
            param_file.writelines('\n')
            param_file.writelines('; parameters for gmap \n')
            """
                gmap_build -D /home/wangyc/maizeGenome/v4 -d /home/wangyc/gmap_database  /home/wangyc/maizeGenome/v4/v4.fa 
                ./gmap -D /home/wangyc/maizeGenome/v4/gmap_database/gmap_database -d gmap_database -f samse -n 0 -t 8  /home/wangyc/outputs/proovread/DN_pread/DN_pread.untrimmed.fasta > /home/wangyc/outputs/gmap/dn/flnc_ccs_corrected_dn.sam
            """
            param_file.writelines('; parameter: The k-mer used for building gmap index (12-15, default 15). \n gmap_idx_kmer = 15\n')
            param_file.writelines('\n\n')

            param_file.writelines('[step_5_star_parameters] \n')
            """
            # library
            STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /home/wangyc/stargdb --genomeFastaFiles /home/wangyc/maizeGenome/v4/v4.fa > starGD_log

            # splicing junction detection
            STAR --runThreadN 8 --genomeDir /home/wangyc/stargdb --readFilesIn /data3/wyc/raw_data/full_cDNA/RNA-seq/A1_R1.fastq /data3/wyc/raw_data/full_cDNA/RNA-seq/A1_R2.fastq 
            --outFileNamePrefix sj_dn. --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM MD --outFilterMismatchNmax 10 
            --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 200000 200000 --alignMatesGapMax 200000 --chimSegmentMin 15 --chimJunctionOverhangMin 15
            """
            param_file.writelines('; input: the path of the RNA-seq file. \n i_timmed_short_read_fpath = auto \n')
            param_file.writelines('; input: the path of star database. \n i_stargdb = auto \n')
            param_file.writelines(
                '; parameter: out file name prefix. \n p_outfile_prefix = {}_ \n'.format(self._proj_name))
            param_file.writelines(
                '; parameter: out sam type. \n p_out_file_type = BAM \n')
            param_file.writelines(
                '; parameter: outFilterMismatchNmax. \n p_mismatch_max = 10 \n'
            )
            param_file.writelines(
                '; parameter: remove non-cononical (T: true, F: false). \n p_remove_noncononical = T \n')
            param_file.writelines(
                '; parameter: align intron max. \n p_align_intro_max = 200000 \n')
            param_file.writelines(
                '; parameter: align Mates Gap Max. \n p_alignMatesGapMax = 200000 \n')
            param_file.writelines(
                '; parameter: int>=0: minimum length of chimeric segment length, if ==0, no chimeric output. \n chimSegmentMin = 15 \n')
            param_file.writelines(
                '; parameter: int>=0: minimum overhang for a chimeric junction. \n chimJunctionOverhangMin = 15 \n')

            param_file.writelines('\n\n')
            
            param_file.writelines(
                '[step_6_filter_and_sort_sam_file_parameters] \n')
            param_file.writelines(
                '; input: the path of the sam file generated by alignment tool. \n sam_fpath = auto \n')
            param_file.writelines(
                '; parameter: the lowest mapping quality allowed in the output sam file. \n mapping_qual = 30 \n')
            param_file.writelines(
                '; parameter: the highest mismatch number allowed in the output sam file. \n max_mismatch = 0.15 \n\n')

            param_file.writelines(
                '[step_7_validate_SJ_and_output_seqs] \n')
            param_file.writelines(
                '; input: the path of sorted sam file. \n sorted_long_sam_fpath = auto \n'
            )
            param_file.writelines(
                '; input: the path of SJ files generated by STAR. \n SJ_short_fpath = auto \n')
            param_file.writelines(
                '; parameter: Minimum alignment coverage (default: 0.99) \n min_coverage = 0.99 \n')
            param_file.writelines('; parameter: Minimum alignment identity \n min_align_identity = 0.95 \n')
            param_file.writelines('; Max fuzzy junction distance \n fuzzy_junction = 5 \n')
            param_file.writelines('; Maximum allowed 5\' difference if on same exon \n max_allow_5 = 1000 \n')
            param_file.writelines(
                '; Maximum allowed 3\' difference if on same exon \n max_allow_3 = 100 \n')
            param_file.writelines(
                '; Don\'t collapse shorter 5\' transcripts (Default set is false or F, change to true or T if you want to turn on) \n donot_collapse_short5 = F \n')
            
            param_file.close()
        return self._param_fpath
        
