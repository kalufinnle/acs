#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Nov 19 19:55:04 2019

@author: Kevin F. Li

Thrust calculation for nozzles where tank pressure and temperature changes as a result of using propellant 
Condition:  <Isentropic nozzle flow>
Assumption: <Upstream pressure decreases as a result of not using regulators>
            <We do not consider too small nozzles, say with chamber size <10 mm and neck size <1 mm, where the effect of boundary layers become predominant>
"""

class test_campaign:
    # class var shared by all instances
    testing_system = 'acs'             

    def __init__(self, test_date, test_title, test_description, a_e, a_c):
        # instance var unique to each instance
        self.date  = test_date           
        self.title = test_title        
        self.doc   = test_description
        self.a_e   = a_e #unit
        self.a_c   = a_c #unit
        self.inputs = {"Test Date": self.date,
                       "Test Title": self.title,
                       "Test Description": self.doc,
                       "Exit Aera of the nozzle": self.a_e,
                       "Choked Aera of the nozzle": self.a_c}
        
        #00000 input here about new param
        
        if self.a_e == self.a_c:
            raise Exception('Exit Area should be larger than Choked area. The value of a_e was: {}'.format(self.a_e))
        if self.a_e > 100*self.a_c:
            raise Exception('Exit Area should NOT exceed 100 times of Choked area. The value of a_e was: {}'.format(self.a_e))

    @classmethod
    def input_param(cls):
        return cls(
                int(input("""=============================
                          \nInput Test Date (e.g.201911) : """)),
                input("""Input Test Title: """),
                input("""Input Test Description: """),
                float(input ("""Input a Exit Area of a Nozzle in numbers : """)),
                float(input ("""Input a Choked Area of a Nozzle in numbers : """))
                #0000000000-add 
               # 1111111111111
        )

        
    def get_descriptive_name(self):
        """Return a neatly formatted descriptive name."""
        long_name = '\r\n |Test Date:  ' + str(self.date) + '\r\n |Test Title: ' + \
                    self.title + '\r\n |Test Description: ' +  \
                    self.doc + '\r\n |Exit Aera of the nozzle: ' + str(self.a_e) + \
                    '\r\n |Choked Aera of the nozzle: ' + str(self.a_c) # 22222222222 #0000000000-add 
        print(long_name)
        return long_name
    
    # redundency for get_descriptive_name()
    def print_keys_values(self):
        for key, value in self.inputs.items(): 
            print ("%s : %s" %(key, value)) 
    
    def algorithm(*array):  #same as *args
        z = 1
        for num in array:
        #Input algorithm from Khalil's code  eg. z*= num  main: multiply = algorithm(3, 5, 10, 6)
            num += 'exp()'
            z = num   
        print(z)
    
    #temperature, pressure, and thrust curve vs time
    def t_p_f_vs_t(self):
        return
    








