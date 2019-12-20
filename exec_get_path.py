#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright 2019 Kevin F. Li. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Kevin F. Li
#
"""
Created on Sat Nov  9 15:08:06 2019

@author: kevin

Your sys.path (i.e. the list of directories Python goes through to search for modules and files) is stored in the path attribute of the sys module. Since path is a list, you can use the append method to add new directories to the path.
"""


import os
import matplotlib.pyplot as plt
import sys


def get_path():
    dirpath = os.getcwk()
    print("Current directory is : " + dirpath)
    foldername = os.path.basename(dirpath)
    print("Current folder name is : " + foldername)
    scriptpath = os.path.abspath(os.path.dirname(sys.argv[0]))
    print("Absolute script path is : " + scriptpath)
    lib_path = input ("Enter the directary which store the in-house packages for fast instllation of modules. Otherwise, enter quit: ")
    if lib_path == 'quit':
        return
    else:
        sys.path.append(str(lib_path))
    print("\n An installation-dependent list of packages directories configured at the time virenv Python is used: ")
    return sys.path


if __name__ == '__main__':
    get_path()









