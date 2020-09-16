import argparse
from parameters.param_generater import ParamGenerator
from parameters.param_reader import ParamsReader


def main():
    parser = argparse.ArgumentParser(
        description='This is a pipeline that process hybrid sequencing data to gain high confidence full-length reads.')
    parser.add_argument(
        '-p', '--proj', help='Generate the paramater file used in the project.', required=True)
    parser.add_argument('-r', '--ref', help='The path of the refernce assembly used in this project.', required='True')
    parser.add_argument(
        '-d', '--dir', help='The directory where the results will be placed. The default directory is current working directory', default='')
    # parser.add_argument('-s', '--step')
    args = parser.parse_args()
    param_gen = ParamGenerator()
    param_gen = param_gen.fit(args.proj, args.ref, args.dir)
    param_fpath = param_gen.generate_param_file()
    param_reader = ParamsReader()
    params_obj = param_reader.fit(param_fpath)


if __name__ == "__main__":
    main()