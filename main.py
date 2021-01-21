import argparse
import os
from parameters.param_generater import ParamGenerator
from parameters.param_reader import ParamsReader
from process_module.sub2flnc import FlncGenerator
from process_module.trim_fq import TrimShortReads
from process_module.correct_long_reads import CorrectLongReads
from process_module.alignment import AlignLongReads
from process_module.finding_sj_by_star import FindingSJBySTAR
from process_module.sort_sam import SortSamFiles
from process_module.SJ_validation import ValidatingSJ



def main():
    parser = argparse.ArgumentParser(
        description='This is a pipeline that process hybrid sequencing data to gain high confidence full-length reads.')
    parser.add_argument(
        '-p', '--proj', help='Generate the paramater file used in the project.', required=True)
    parser.add_argument('-r', '--ref', help='The path of the refernce assembly used in this project.')
    parser.add_argument(
        '-d', '--dir', help='The directory where the results will be placed. The default directory is current working directory', default='')
    parser.add_argument('-i', '--config', help='The path of the config file used in the corresponding project.')
    parser.add_argument('-s', '--step', default='')
    args = parser.parse_args()
    if not args.step:
        param_gen = ParamGenerator()
        param_gen = param_gen.fit(args.proj, args.ref, args.dir)
        param_fpath = param_gen.generate_param_file()
    else:
        param_fpath = os.path.abspath(args.config)
    param_reader = ParamsReader()
    params_obj = param_reader.fit(param_fpath)
    wkdir = params_obj.get('common_parameters', 'output_dir')
    os.chdir(os.path.abspath(wkdir))

    flnc_fpath = ''
    lst_trimmed_fpath = list()
    
    if '1' in args.step:
        flnc_seq = FlncGenerator()
        flnc_seq = flnc_seq.fit(params_obj)
        flnc_seq.sub2ccs()
        flnc_seq.ccs2flccs()
        flnc_fpath = flnc_seq.flccs2flnc()

    if '2' in args.step:
        trimmed_files = TrimShortReads()
        trimmed_files = trimmed_files.fit(params_obj)
        lst_trimmed_fpath = trimmed_files.cmd_fqtrim()

    if '3' in args.step:
        # if not flnc_fpath:
            # flnc_fpath = os.path.abspath(params_obj.get(
            #     'step_3_correction_parameters', 'i_long_read_fpath'))
        if not lst_trimmed_fpath:
            lst_trimmed_fpath = params_obj.get(
                'step_3_correction_parameters', 'i_short_read_fpath').strip().split(',')
        corrected_file = CorrectLongReads()
        corrected_file = corrected_file.fit(
            params_obj, flnc_fpath, lst_trimmed_fpath)
        corrected_long_read_fpath = corrected_file.execute_correct_program()
    
    if '4' in args.step:
        # input_corrected_long_read_fpath = params_obj.get(
        #     'step_3_correction_parameters', 'o_corrected_seq_fpath')
        input_corrected_long_read_fpath = corrected_long_read_fpath
        if input_corrected_long_read_fpath.lower() == 'auto':
            input_corrected_long_read_fpath = corrected_long_read_fpath
        align_long_reads = AlignLongReads()
        align_long_reads = align_long_reads.fit(params_obj,
            input_corrected_long_read_fpath)
        aligned_reads_fpath = align_long_reads.align_corrected_reads()
        
    if '5' in args.step:
        input_trimmed_short_read_fpath = lst_trimmed_fpath
        
        find_sj = FindingSJBySTAR()
        find_sj.fit(params_obj, input_trimmed_short_read_fpath)
        splicing_junction_fpath = find_sj.finding_sj()
    
    if '6' in args.step:
         
        sort_filter_sam = SortSamFiles()
        sort_filter_sam = sort_filter_sam.fit(params_obj, aligned_reads_fpath)
        sorted_filtered_sam_fpath = sort_filter_sam.sort_filter_sam()
    
    if '7' in args.step:
        validated_seq = ValidatingSJ()
        validated_seq = validated_seq.fit(
            params_obj, corrected_long_read_fpath, sorted_filtered_sam_fpath, splicing_junction_fpath)
        # out_validated_seq_fpath = validated_seq.write_the_validated_seqs()
        print('Completed, the collapsed sequences are stored in {}'.format(out_validated_seq_fpath))

if __name__ == "__main__":
    main()
