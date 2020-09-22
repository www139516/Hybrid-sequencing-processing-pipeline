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
            self._output_dpath = output_dpath
        
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
        else:
            param_file = open(self._param_fpath, 'a')

            param_file.writelines(
                '; This is the config file containing all the parameters used in this pipeline. \n')
            param_file.writelines('\n\n')

            param_file.writelines('[common_parameters]\n')
            param_file.writelines('; The name of this project. \n proj_name = {} \n'.format(self._proj_name))
            param_file.writelines('; The path of the reference assembly file used in this project. ref_assembly = {} \n'.format(self._ref_fpath))
            param_file.writelines('; The path of the directory for the output results. \n output_dir = {} \n'.format(self._output_dpath))
            param_file.writelines('; The thread used in this program, default value is 1. \n threads = 1 \n')
            param_file.writelines('\n')

            # ccs movieX.subreads.bam movieX.ccs.bam --noPolish --minPasses 1 
            # lima --isoseq --dump-clips --no-pbi [--peek-guess] -j 24 ccs.bam primers.fasta demux.bam
            # isoseq3 refine --require-polya combined_demux.consensusreadset.xml primers.fasta unpolished.flnc.bam
            param_file.writelines('[step_1_isoseq_parameters] \n')
            param_file.writelines('; input: The path of the input subreads file. \n i_subread_fpath = \n')
            param_file.writelines('; input: The path of the primers.fasta file. \n i_primer_fpath = \n')
            param_file.writelines(
                '; param: Output the initial template derived from the POA (faster but less accurate), default value set in this programe is True \n p_no_polish = True \n')
            param_file.writelines(
                '; param: Minimum number of subreads required to generate CCS, default value set in this program is 1. \n p_min_passes = 1 \n')
            param_file.writelines(
                '; param: set to True if your transcripts have a polyA tail, default value in this program is True. \n p_polya = True \n')
            param_file.writelines('; output: The output of this step, automatically set if you do not set it in this config file. \n o_flnc_out = auto \n')
            param_file.writelines('\n')

            param_file.writelines('[step_2_fqtrim_parameters] \n')
            """
                fqtrim -A -l25 -o trimmed.fq.gz exome_reads_1.fastq.gz,exome_reads_2.fastq.gz
            """
            param_file.writelines('; input: The path of the short read file, fastq or fastq.gz format, using comma to separate two paired files if needed. \n i_short_read_file = \n')
            # param_file.writelines('; input: The score system used in short read sequencing. phred33 (default) or phred64. \n phred = \n')
            param_file.writelines('; param: minimal lenth retained after trimming. \n p_lenth = \n')
            param_file.writelines('; param: genomic sequencing (g, means do not trim poly A) or RNA-sequencing (r, default, trim poly A). \n p_sequencing_type = r \n')
            param_file.writelines('; param: quality value for trimming, default is 20. \n p_qvalue = 20 \n')
            param_file.writelines(
                '; onput: The output suffix of the trimmed file. \n o_suffix = trimmed.fq.gz \n')
            param_file.writelines('; output: the path of the dir where you put the output trimmed files. \n o_trimmed_dir = {}/trimmed_files/ \n'.format(self._output_dpath))
            param_file.close()
        return self._param_fpath
        
