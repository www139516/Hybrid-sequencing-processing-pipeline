import os
import subprocess

"""
    # library
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /home/wangyc/stargdb --genomeFastaFiles /home/wangyc/maizeGenome/v4/v4.fa > starGD_log

    # splicing junction detection
    STAR --runThreadN 8 --genomeDir /home/wangyc/stargdb --readFilesIn /data3/wyc/raw_data/full_cDNA/RNA-seq/A1_R1.fastq /data3/wyc/raw_data/full_cDNA/RNA-seq/A1_R2.fastq 
    --outFileNamePrefix sj_dn. --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM MD --outFilterMismatchNmax 10 
    --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 200000 200000 --alignMatesGapMax 200000 --chimSegmentMin 15 --chimJunctionOverhangMin 15
"""


class FindingSJBySTAR():
    def __init__(self):
        self._param_obj = None
        self._lst_i_short_fpath = []
        self._ref_fpath = ''
        self._i_ref_dpath = ''
        self._i_star_ref_dpath = ''
        self._p_thread_num = ''
        self._p_outfile_prefix = ''
        self._p_out_file_type = ''
        self._p_mixmatch_max = ''
        self._p_rmnoncononical = ''
        self._p_intro_max = ''
        # self._mismatches = ''
        self._p_match_gap_max = ''
        self._p_chim_seqment_min = ''
        self._p_chim_junc_overhan_min = ''
        self._o_sj_dpath = ''
        self._o_sj_fpath = ''
        
    def fit(self, param_obj, trimmed_short_read_fpath):
        self._param_obj = param_obj
        self._ref_fpath = self._param_obj.get(
            'common_parameters', 'ref_assembly')
        self._p_thread_num = self._param_obj.get(
            'common_parameters', 'threads')
        self._lst_i_short_fpath = (trimmed_short_read_fpath)

        self._i_star_ref_dpath = self._param_obj.get(
            'step_5_star_parameters', 'i_stargdb'
        )
        if self._i_star_ref_dpath.lower() == 'auto':
            self._i_star_ref_dpath = os.path.join(os.path.dirname(self._ref_fpath), 'stargdb')
        self._p_outfile_prefix = self._param_obj.get(
            'step_5_star_parameters', 'p_outfile_prefix')
        self._p_out_file_type = self._param_obj.get(
            'step_5_star_parameters', 'p_out_file_type')

        self._p_mixmatch_max = self._param_obj.get(
            'step_5_star_parameters', 'p_mismatch_max'
        )

        self._p_rmnoncononical = self._param_obj.get(
            'step_5_star_parameters', 'p_remove_noncononical'
        )
        if self._p_rmnoncononical[0].upper() == 'T':
            self._p_rmnoncononical = 'RemoveNoncanonical'
        else:
            self._p_rmnoncononical = ''

        self._p_intro_max = self._param_obj.get(
            'step_5_star_parameters', 'p_align_intro_max'
        )
        self._p_match_gap_max = self._param_obj.get(
            'step_5_star_parameters', 'p_alignMatesGapMax'
        )
        self._p_chim_seqment_min = self._param_obj.get(
            'step_5_star_parameters', 'chimSegmentMin'
        )
        self._p_chim_junc_overhan_min = self._param_obj.get(
            'step_5_star_parameters', 'chimJunctionOverhangMin'
        )
        self._o_sj_dpath = os.path.join(self._param_obj.get(
            'common_parameters', 'output_dir'), 'sj_sites')
        if not os.path.exists(self._o_sj_dpath):
            os.mkdir(self._o_sj_dpath)
        self._o_sj_fpath = os.path.join(self._o_sj_dpath, '{}SJ.out.tab'.format(self._p_outfile_prefix))

    def __unzip(self):
        lst_unziped_fpaths = []
        for fpath in self._lst_i_short_fpath:
            if fpath.endswith('.gz'):
                print("Unpacking {gz_fpath} ...".format(gz_fpath=fpath))
                cmd_unzip = 'gzip -d {gz_fpath}'.format(gz_fpath=fpath)
                lst_unziped_fpaths.append(fpath[:-3])
                if subprocess.check_call(cmd_unzip, shell=True) != 0:
                    raise SystemCommandError
            else:
                lst_unziped_fpaths.append(fpath)
        return lst_unziped_fpaths
        


    def finding_sj(self):
        print('Finding splicing junctions using short-read sequences...')
        os.chdir(self._o_sj_dpath)
        if not os.path.exists(self._i_star_ref_dpath):
            cmd = 'STAR --runThreadN 8 --runMode genomeGenerate --genomeDir {star_gdb} --genomeFastaFiles {ref_fpath} > starGD_log'.format(
                star_gdb = self._i_star_ref_dpath,
                ref_fpath = self._ref_fpath
            )
            if subprocess.check_call(cmd, shell=True) != 0:
                raise SystemCommandError

        self._lst_i_short_fpath = self.__unzip()
        if len(self._lst_i_short_fpath) > 1:
            short_read_fpaths = ' '.join(self._lst_i_short_fpath)
        else:
            short_read_fpaths = self._lst_i_short_fpath[0]
        cmd_find_sj = 'STAR --runThreadN {thread_num} --genomeDir {star_gdb} --readFilesIn {short_reads_fpath} --outFileNamePrefix {out_prefix} --outSAMtype {out_file_type} SortedByCoordinate --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM MD --outFilterMismatchNmax {mismatchmax}\
        --outFilterIntronMotifs {rm_non_cononical} --alignIntronMax {intro_max} --alignMatesGapMax {mates_gap_max} --chimSegmentMin {chim_segment_min} --chimJunctionOverhangMin {chim_junc_overhan_min}'.format(
            thread_num = self._p_thread_num,
            star_gdb = self._i_star_ref_dpath,
            short_reads_fpath=short_read_fpaths,
            out_prefix = self._p_outfile_prefix,
            out_file_type = self._p_out_file_type,
            mismatchmax = self._p_mixmatch_max,
            rm_non_cononical = self._p_rmnoncononical,
            intro_max = self._p_intro_max,
            mates_gap_max = self._p_match_gap_max,
            chim_segment_min = self._p_chim_seqment_min,
            chim_junc_overhan_min = self._p_chim_junc_overhan_min
        )
        if subprocess.check_call(cmd_find_sj, shell=True) != 0:
            raise SystemCommandError
        print(self._o_sj_dpath)
        return os.path.abspath(self._o_sj_fpath)
