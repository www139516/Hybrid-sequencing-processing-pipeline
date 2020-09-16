"""
This file is used to read and store the paramters in a dict, and these paramters will be used in the following steps in the program.
"""

import configparser
import os
import sys


class ParamsReader:
    def __init__(self):
        self._conf_fpath = ''
        self._config = None

    def fit(self, conf_file_path):

        print('config_fpath: {}'.format(conf_file_path))

        if not os.path.exists(conf_file_path):
            print("The config file path is not exist.")
            sys.exit()
        else:
            self._conf_fpath = conf_file_path

        self._config = configparser.ConfigParser()
        self._config.read(self._conf_fpath, encoding='utf-8')
        # print(self._config.sections())
        # print(self._config.items('common_parameters'))

        return self._config
    
