# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 17:31:35 2022

@author: Sangani
"""

#!/bin/bash

import clean_ENCODE_functions as f

file = input('Please enter the file need to be lowercased: ')
new_file = input(' Input the name of file to saved the updated file: ')
f.lower_case(file, new_file)