# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 14:34:08 2019

@author: kevin
"""

import numpy as np
import func.thrust_cal as tc
import matplotlib.pyplot as plt


def main(test):    
     num = int(input("""Input number of test you want to run (>0) : """))
     for k in range(num):
         test.append(tc.test_campaign.input_param())
         
     filename = 'test_campaign.txt'
     open(filename, 'w').close()
     for i in range(len(test)):
         try:
             with open(filename, 'a') as f:
                print('\n ==== Test {} ==== {}'.format(i+1,test[i].get_descriptive_name()), file=f)
            
         except:
             with open(filename, 'a') as f:
                 print('\n ==== Test {} ==== {}'.format(i+1,test[i].print_keys_values()), file=f)                 

     f.close()



# ----------------------------------------------------------------------
#   Core
# ----------------------------------------------------------------------
if __name__ == '__main__':
     plt.close('all')
     test = []
     main(test)
     plt.show()

   