#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/6 10:04
# @Author : Yuancong Wang
# @Site :
# @File : FlncGenerator.py
# @Software: PyCharm
"""
Description: This class processes the subreads from long-read sequencing to generate full length non-chimeric CCSs
  (FLNC CCSs).
Output: FLNC CCSs file
Input: Subreads.bam file
Other notes: Using IsoSeq3 to process subreads
"""

import os
import subprocess


class FlncGenerator:

    def __init__(self):
        self._params_obj = None
        self._subread_fpath = None
        self._primer_fpath = None
        self._out_dpath = None
        self._movie_name = None
        self._out_ccs_fpath = None
        self._out_flccs_fpath = None
        self._in_flnc_fpath = None
        self._out_flnc_fpath = None
        self._out_flnc_fa_fpath = ''

    def fit(self, param_obj):
        """
        set the values for this class
        :param subread_fpath: the path of the subreads file
        :param primer_fpath: the path of the file containing the primer sequences used in the PacBio-sequencing
        :param out_dpath: the out put directory path
        :return: self
        """

        self._params_obj = param_obj
        self._subread_fpath = os.path.abspath(self._params_obj.get('step_1_isoseq_parameters', 'i_subread_fpath'))
        self._primer_fpath = os.path.abspath(self._params_obj.get('step_1_isoseq_parameters', 'i_primer_fpath'))

        assert os.path.exists(self._subread_fpath), \
            'The subread file is not existed, please make sure that you have the subread file.'
        assert os.path.exists(self._primer_fpath), \
            'The primer file is not existed, please make sure that you have the primer file.'

        out_tmp = self._params_obj.get('step_1_isoseq_parameters', 'o_flnc_out')
        if out_tmp.lower() == 'auto':
            self._out_dpath = os.path.join(os.path.abspath(
                self._params_obj.get('common_parameters', 'output_dir')
                ), 'out_flnc_files')
            if not os.path.exists(self._out_dpath):
                os.mkdir(self._out_dpath)
        else:
            out_tmp = os.path.abspath(out_tmp)
            if not os.path.exists(out_tmp):
                os.mkdir(out_tmp)
            self._out_dpath = out_tmp

        self._movie_name = os.path.basename(self._subread_fpath).split(".")[0]

        return self


    def sub2ccs(self):
        """
        generate ccs based on subreads.bam
        :return: None
        """
        print('Convert subreads to ccs...')
        out_file_name = self._movie_name + '.ccs.bam'
        self._out_ccs_fpath = os.path.join(self._out_dpath, out_file_name)
        if os.path.exists(self._out_ccs_fpath):
            print('{} already exists.'.format(self._out_ccs_fpath))
        else:
            no_polish = self._params_obj.get(
                'step_1_isoseq_parameters', 'p_no_polish')
            min_pass = self._params_obj.get(
                'step_1_isoseq_parameters', 'p_min_passes')
            if no_polish.lower() == 'true':
                cmd_subreads2ccs = 'ccs {input_subreads} {output} '.format(input_subreads=self._subread_fpath,
                                                                        output=self._out_ccs_fpath) + \
                    '--noPolish --minPasses {}'.format(min_pass)
            else:
                cmd_subreads2ccs = 'ccs {input_subreads} {output} '.format(input_subreads=self._subread_fpath,
                                                                           output=self._out_ccs_fpath) + \
                    '--minPasses {}'.format(min_pass)
            print(cmd_subreads2ccs)
            if subprocess.check_call(cmd_subreads2ccs, shell=True) != 0:
                raise SystemCommandError

    def ccs2flccs(self):
        """
        generate full length ccs based on ccs
        :return: None
        """
        print('Convert ccs to flccs...')
        # prefix output file name generated by lima
        out_file_name = self._movie_name + '.flccs.bam'
        self._out_flccs_fpath = os.path.join(self._out_dpath, out_file_name)
        # the actual output file name generated by lima
        in_file_name = self._movie_name + '.flccs.primer_5p--primer_3p.bam'
        self._in_flnc_fpath = os.path.join(self._out_dpath, in_file_name)

        if os.path.exists(self._in_flnc_fpath):
            print('{} already exists.'.format(self._out_flccs_fpath))
        else:
            cmd_ccs2flccs = 'lima {input_ccs} {barcode} {output_flccs} '.format(input_ccs=self._out_ccs_fpath,
                                                                                barcode=self._primer_fpath,
                                                                                output_flccs=self._out_flccs_fpath) + \
                            '--isoseq --peek-guess'
            print(cmd_ccs2flccs)
            if subprocess.check_call(cmd_ccs2flccs, shell=True) != 0:
                raise SystemCommandError

    def flccs2flnc(self):
        """
        generate FLNC CCSs based on FLCCS
        :return: the path of the FLNC CCS file
        """
        print('Convert flccs to flnc ccs...')
        out_file_name = self._movie_name + '.flnc.bam'
        self._out_flnc_fpath = os.path.join(self._out_dpath, out_file_name)
        if os.path.exists(self._out_flnc_fpath):
            print('{} already exists.'.format(self._out_flnc_fpath))
        else:
            # check if need to search polya as a way of identifying full length cDNA.
            polya = self._params_obj.get('step_1_isoseq_parameters', 'p_polya')

            if polya.lower() == 'true':
                cmd_flccs2flnc = "isoseq3 refine {input_flccs} {barcode} {output_flnc} ".format(input_flccs=self._in_flnc_fpath,
                                                                                                barcode=self._primer_fpath,
                                                                                                output_flnc=self._out_flnc_fpath) + \
                                "--require-polya"
            else:
                cmd_flccs2flnc = "isoseq3 refine {input_flccs} {barcode} {output_flnc} ".format(input_flccs=self._in_flnc_fpath,
                                                                                                barcode=self._primer_fpath,
                                                                                                output_flnc=self._out_flnc_fpath) 
            print(cmd_flccs2flnc)                                                                                    
            if subprocess.check_call(cmd_flccs2flnc, shell=True) != 0:
                raise SystemCommandError
            self.__bam2fa()
        return self._out_flnc_fa_fpath

    def __bam2fa(self):
        """
        convert the .bam to .fasta format
        :return: path of the fasta file
        """
        out_file_name = self._movie_name + '.flnc.fasta'
        self._out_flnc_fa_fpath = os.path.join(
            self._out_dpath, out_file_name)
        if not os.path.exists(self._out_flnc_fa_fpath):
            cmd_bam2fa = '''
            samtools view {flnc_bam_path} | awk \
                    '{{OFS="\\t"; print ">"$1"\\n"$10}}' - > {flnc_fa_path}
            '''.format(flnc_bam_path=self._out_flnc_fpath, flnc_fa_path=self._out_flnc_fa_fpath)
            if subprocess.check_call(cmd_bam2fa, shell=True) != 0:
                raise SystemCommandError
        
