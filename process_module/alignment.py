import os 
import subprocess
import re


class AlignLongReads:
    def __init__(self):
        self._param_obj = None
        self._program_for_alignment = ''
        self._i_corrected_long_fpath = ''
        self._i_ref_fpath = ''
        self._i_ref_dpath = ''
        self._i_ref_idx_fpath = ''
        self._o_mapped_file_dpath = ''
        self._o_mapped_long_fpath = ''
        self._proj_name = ''
        self._p_gmap_idx_kmer = ''
        self._threads = ''

    def fit(self, param_obj, i_correct_long_fpath = '', o_alignment_dpath = 'aligned_seq'):
        self._param_obj = param_obj

        self._proj_name = self._param_obj.get('common_parameters', 'proj_name')
        self._threads = self._param_obj.get('common_parameters', 'threads')
        self._program_for_alignment = self._param_obj.get(
            'step_4_alignment_parameters', 'program_for_alignment')
        self._i_corrected_long_fpath = self._param_obj.get(
            'step_4_alignment_parameters', 'i_long_read_fpath')
        if self._i_corrected_long_fpath.lower() == 'auto':
            self._i_corrected_long_fpath = os.path.abspath(i_correct_long_fpath)
        self._o_mapped_file_dpath = self._param_obj.get(
            'step_4_alignment_parameters', 'o_out_mm2_file')
        if self._o_mapped_file_dpath.lower() == 'auto':
            self._o_mapped_file_dpath = o_alignment_dpath
        if not os.path.exists(self._o_mapped_file_dpath):
            os.mkdir(self._o_mapped_file_dpath)
        self._i_ref_fpath = self._param_obj.get(
            'common_parameters', 'ref_assembly')
        self._i_ref_dpath = os.path.dirname(self._i_ref_fpath)
        self._o_mapped_long_fpath = os.path.join(self._o_mapped_file_dpath, '{}.sam').format(self._proj_name)
        self._p_gmap_idx_kmer = self._param_obj.get(
            'step_4_alignment_parameters', 'gmap_idx_kmer')
        return self

    def __align_use_minimap2(self):
        # check if the index file is exist, if not, creating it.
        mm2_idx_fpath = self._i_ref_fpath + '.mmi'
        if not os.path.exists(mm2_idx_fpath):
            cmd_build_mm2_idx = 'minimap2 -d {ref_mmi} {ref_fa}'.format(ref_mmi = mm2_idx_fpath,
                                                                        ref_fa = self._i_ref_fpath)
            print('The index file of minimap2 is not exist, building it....')
            if subprocess.check_call(cmd_build_mm2_idx, shell=True) != 0:
                raise SystemCommandError

        # using minimap2 to align the corrected reads, generating a sam file
        cmd_align_minimap2 = 'minimap2 -a {ref_fa} {query_fq} > {alignment}'.format(ref_fa = self._i_ref_fpath,
            query_fq = self._i_corrected_long_fpath,
            alignment = self._o_mapped_long_fpath)
        print('Using minimap2 to align the corrected reads....')
        if subprocess.check_call(cmd_align_minimap2, shell=True) != 0:
                raise SystemCommandError
        # mapping cDNA reads

        # mapping DNA reads
        return os.path.abspath(self._o_mapped_long_fpath)

    def __align_use_gmap(self):
        # check if the index file is exist, if not, creating it.
        """
        gmap_build -D /home/wangyc/maizeGenome/v4 -d /home/wangyc/gmap_database  /home/wangyc/maizeGenome/v4/v4.fa
        """
        gmap_idx_dpath = os.path.join(self._i_ref_dpath, 'gmap_database')
        if not os.path.exists(self._i_ref_dpath):
            os.mkdir(self._i_ref_dpath)
            cmd_build_gmap_idx = 'gmap_build -D {ref_dpath} -d {gmap_database_dpath} {ref_fa}'.format(ref_dpath = self._i_ref_dpath)
            if subprocess.check_call(cmd_build_gmap_idx, shell=True) != 0:
                raise SystemCommandError
        # using gmap to align the corrected reads
        cmd_gmap_align = 'gmap -D {gmap_idx_dpath} -d {gmap_idx_dname} -f samse -n 0 -t {threads} {corrected_fasta} > {out_aligned_file}'.format(
            gmap_idx_dpath=gmap_idx_dpath,
            gmap_idx_dname="gmap_database",
            threads = self._threads,
            corrected_fasta=self._i_corrected_long_fpath,
            out_aligned_file=self._o_mapped_long_fpath
        )
        if subprocess.check_call(cmd_gmap_align, shell=True) != 0:
                raise SystemCommandError
        return os.path.abspath(self._o_mapped_long_fpath)

    def align_corrected_reads(self):
        if self._program_for_alignment.lower() == 'm':
            self._o_mapped_long_fpath = self.__align_use_minimap2()
        elif self._program_for_alignment.lower() == 'g':
            self._o_mapped_long_fpath = self.__align_use_gmap()
        return self._o_mapped_long_fpath
