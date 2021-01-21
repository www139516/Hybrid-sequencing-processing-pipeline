#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/6 10:04
# @Author : Yuancong Wang
# @Site :
# @File : FlncGenerator.py
# @Software: PyCharm
"""
Discription: This program is used to trim short read fastq files, return trimmed clean reads files.
"""
import os
import re
import subprocess


class TrimShortReads:
    def __init__(self):
        self._params_obj = None
        self._short_read1_fpath = ''
        self._short_read2_fpath = ''
        self._min_len = ''
        self._is_genom = ''
        self._qvalue = ''
        self._o_suffix = ''
        self._o_dpath = ''
        self._lst_o_fpath = []

    def fit(self, config_File):
        self._params_obj = config_File
        self._min_len = self._params_obj.get(
            'step_2_fqtrim_parameters', 'p_lenth')
        self._is_genom = self._params_obj.get(
            'step_2_fqtrim_parameters', 'p_sequencing_type')
        self._qvalue = self._params_obj.get('step_2_fqtrim_parameters', 'p_qvalue')
        self._o_suffix = self._params_obj.get(
            'step_2_fqtrim_parameters', 'o_suffix')
        self._o_dpath = os.path.abspath(self._params_obj.get(
            'step_2_fqtrim_parameters', 'o_trimmed_dir'))
        if not os.path.exists(self._o_dpath):
            os.mkdir(self._o_dpath)
        lst_input_fq_files = self._params_obj.get(
            'step_2_fqtrim_parameters', 'i_short_read_file').strip().split(',')
        num_input_files = len(lst_input_fq_files)
        if num_input_files == 2:
            self._short_read1_fpath = lst_input_fq_files[0]
            self._short_read2_fpath = lst_input_fq_files[1]
        elif num_input_files == 1:
            self._short_read1_fpath = lst_input_fq_files[0]

        return self

    def cmd_fqtrim(self):
        if not self._short_read2_fpath:
            in_fname = os.path.basename(self._short_read1_fpath)
            out_fname = re.split(r'.f(?:ast)?q', in_fname)[0] + '.{}'.format(self._o_suffix)
            if self._is_genom.lower().strip() == 'g':
                cmd = 'fqtrim -A -l {length} -q {qvalue} --outdir {outdir} -o trimmed.fq.gz {seq_f1}'.format(
                    length = self._min_len,
                    qvalue = self._qvalue,
                    outdir=self._o_dpath,
                    seq_f1=self._short_read1_fpath                 
                )
            elif self._is_genom.lower().strip() == 'r':
                cmd = 'fqtrim -l {length} -q {qvalue} --outdir {outdir} -o trimmed.fq.gz {seq_f1}'.format(
                    length=self._min_len,
                    qvalue=self._qvalue,
                    outdir=self._o_dpath,
                    seq_f1=self._short_read1_fpath
                )
            self._lst_o_fpath.append(os.path.join(self._o_dpath, out_fname))
        else:
            in_fname1 = os.path.basename(self._short_read1_fpath)
            in_fname2 = os.path.basename(self._short_read2_fpath)
            out_fname1 = re.split(r'.f(?:ast)?q', in_fname1)[
                0] + '.{}'.format(self._o_suffix)
            out_fname2 = re.split(r'.f(?:ast)?q', in_fname2)[
                0] + '.{}'.format(self._o_suffix)
            if self._is_genom.lower().strip() == 'g':
                cmd = 'fqtrim -A -l {length} -q {qvalue} --outdir {outdir} -o trimmed.fq.gz {seq_f1},{seq_f2}'.format(
                    length = self._min_len,
                    qvalue = self._qvalue,
                    outdir = self._o_dpath,
                    seq_f1 = self._short_read1_fpath,
                    seq_f2 = self._short_read2_fpath                
                )
            elif self._is_genom.lower().strip() == 'r':
                cmd = 'fqtrim -l {length} -q {qvalue} --outdir {outdir} -o trimmed.fq.gz {seq_f1},{seq_f2}'.format(
                    length = self._min_len,
                    qvalue = self._qvalue,
                    outdir = self._o_dpath,
                    seq_f1 =self._short_read1_fpath,
                    seq_f2 = self._short_read2_fpath
                )
            self._lst_o_fpath.append(os.path.join(self._o_dpath, out_fname1))
            self._lst_o_fpath.append(os.path.join(self._o_dpath, out_fname2))

        print('Trimming short reads...')
        if subprocess.check_call(cmd, shell=True) != 0:
                raise SystemCommandError

        return self._lst_o_fpath
