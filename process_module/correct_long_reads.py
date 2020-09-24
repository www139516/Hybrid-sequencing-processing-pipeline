#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/6 10:04
# @Author : Yuancong Wang
# @Site :
# @File : FlncGenerator.py
# @Software: PyCharm
"""
Discription: This program is used to correct the long reading sequences using short read sequences.
software: lordec
"""


import os
import subprocess


class CorrectLongReads:
    def __init__(self):
        self._param_obj = None
        self._program_use = ''
        self._lst_i_short_fpath = []
        self._i_long_fpath = ''
        self._p_kmer_len = ''
        self._p_solid_threshold = ''
        self._p_abundance_max = ''
        self._p_errorrate = ''
        self._p_threads = ''
        self._o_corrected_dpath = ''
        self._o_corrected_fpath = ''
        self._o_stat_fpath = ''
        self._prov_coverage = ''
        self._prov_prefix = ''

    def fit(self, param_obj, flnc_fpath, lst_trimmed_read_fpath):
        self._param_obj = param_obj
        self._program_use = self._param_obj.get(
            'step_3_correction_parameters', 'program_for_correction')
        self._i_long_fpath = self._param_obj.get(
            'step_3_correction_parameters', 'i_long_read_fpath')
        if self._i_long_fpath.lower() == 'auto':
            self._i_long_fpath = flnc_fpath
        self._i_long_fpath = os.path.abspath(self._i_long_fpath)
        tmp_short_fpath = self._param_obj.get(
            'step_3_correction_parameters', 'i_short_read_fpath')
        if tmp_short_fpath.lower() == 'auto':
            self._lst_i_short_fpath = lst_trimmed_read_fpath
        else:
            self._lst_i_short_fpath.append(
                os.path.abspath(tmp_short_fpath.strip().split(',')[0]))
            self._lst_i_short_fpath.append(
                os.path.abspath(tmp_short_fpath.strip().split(',')[1]))
        self._p_kmer_len = self._param_obj.get(
            'step_3_correction_parameters', 'p_kmer_len')
        self._p_solid_threshold = self._param_obj.get(
            'step_3_correction_parameters', 'p_solid_threshold')
        self._p_abundance_max = self._param_obj.get(
            'step_3_correction_parameters', 'p_abundance_max')
        self._p_errorrate = self._param_obj.get(
            'step_3_correction_parameters', '0.4')
        self._p_threads = self._param_obj.get('common_parameters', 'threads')
        self._o_corrected_dpath = os.path.join(self._param_obj.get(
            'common_parameters', 'output_dir'), 'corrected_reads')
        if not os.path.exists(self._o_corrected_dpath):
            os.mkdir(self._o_corrected_dpath)
        self._prov_coverage = self._param_obj.get(
            'step_3_correction_parameters', 'p_prov_coverage')
        self._prov_prefix = self._param_obj.get(
            'step_3_correction_parameters', 'o_prov_prefix')
        if self._prov_prefix.lower() == 'auto':
            self._prov_prefix = 'proovread_corrected'
        
        return self

    def lordec_correct(self):
        print('Using lordec to correct long-read sequences...')
        tmp_fname = os.path.basename(self._i_long_fpath)
        tmp_o_fname = 'lordec_corrected_{}'.format(tmp_fname)
        self._o_corrected_fpath = os.path.join(self._o_corrected_dpath, tmp_o_fname)
        tmp_lor_stat_fname = 'lordec_stats_{}'.format(tmp_fname)
        self._o_stat_fpath = os.path.join(self._o_corrected_dpath, tmp_lor_stat_fname)
        if len(self._lst_i_short_fpath) == 1:
            cmd = 'lordec-correct -i {long_reads} -2 {short_reads} -T {threads} -k {kmer} -o {output} -s {solid_k} -e {errorrate} -a {abundance_max} -S {stat_file}'.format(
                long_reads = self._i_long_fpath,
                short_reads = self._lst_i_short_fpath[0],
                threads = self._p_threads,
                kmer = self._p_kmer_len,
                output = self._o_corrected_dpath,
                solid_k = self._p_solid_threshold,
                errorrate = self._p_errorrate,
                abundance_max = self._p_abundance_max,
                stat_file = self._o_stat_fpath
            )

        else:
            
            cmd = 'lordec-correct -i {long_reads} -2 {short_reads} -T {threads} -k {kmer} -o {output} -s {solid_k} -e {errorrate} -a {abundance_max} -S {stat_file}'.format(
                long_reads=self._i_long_fpath,
                short_reads=','.join(self._lst_i_short_fpath),
                threads=self._p_threads,
                kmer=self._p_kmer_len,
                output=self._o_corrected_dpath,
                solid_k=self._p_solid_threshold,
                errorrate=self._p_errorrate,
                abundance_max=self._p_abundance_max,
                stat_file=self._o_stat_fpath
            )
        print(cmd)
        if subprocess.check_call(cmd, shell=True) != 0:
            raise SystemCommandError

        return self._o_corrected_fpath


        
        

    def proovread_correct(self):
        print('Using proovread to correct long-read sequences...')
        if len(self._lst_i_short_fpath) == 1:
            tmp_short_fpath = self._lst_i_short_fpath[0]
        else:
            tmp_short_fpath = ','.join(self._lst_i_short_fpath)
        cmd = 'proovread -l {long_reads} -s {short_reads} -p {out_prefix} -t {threads} --coverage {coverage} > proovread.log'.format(
            long_reads = self._i_long_fpath,
            short_reads = tmp_short_fpath,
            out_prefix = self._prov_prefix,
            threads = self._p_threads,
            coverage = self._prov_coverage
        )
        print(cmd)
        if subprocess.check_call(cmd, shell=True) != 0:
            raise SystemCommandError
        return self._o_corrected_fpath

    def correct_execute_correct(self):
        pass
