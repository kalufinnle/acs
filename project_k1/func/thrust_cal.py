#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Nov 19 19:55:04 2019

@author: Kevin F. Li

1.Thrust calculation for NOZZLES where tank pressure and temperature changes as a result of using propellant 
2.For tank blow-down model (describing how fast a tank depressurize through a small orifice or nozzel attached to a tank)
Assumption: <Isentropic nozzle flow>
            <Absolute pressure is used (psia)>
            <Mass flow rate calculation is suitable for ideal compressible gas>
            <The gas is calorically perfect for Mach less than 5; therefore ratio of specific heat is constant >
            <Isotheraml during on-ground test>
            <Not Isotheral during flight, therefore nitrogen density flutuate due to atmospheric temperature change even nitrogen is not being used>

Condition:  <Upstream pressure decreases as a result of not using regulators>
            <We do not consider too small nozzles, say with chamber size <10 mm and neck size <1 mm, where the effect of boundary layers become predominant>
            <Citical point of nitrogen is at 126.2K/-147.0C/-232.5F && 3.40Mpa/34.0bara/493psia>
Note:  At subsonic and low supersonic Mach numbers (M<5), air is calorically perfect. But under low hypersonic conditions (M>=5), air is calorically imperfect. 
The specific heat capacity changes with the temperature of the flow because of excitation of the vibrational modes of the diatomic nitrogen and oxygen of the atmosphere. 
For simulation illustrating molecular vibration under calorically imperfect condition, use function "calorically_imperf"
"""

class Test_campaign:
    # class var shared by all instances
    testing_system = 'acs'  
    def __init__(self, test_date, test_title, test_description, p_ambient, m_dot, a_c, p_t, t_t, specific_volume, a_e, p_e, t_e, velo_e, m_e):
        # instance var unique to each instance
        self.date  = test_date           
        self.title = test_title        
        self.doc       = test_description
        self.p_ambient = p_ambient #unit
        self.m_dot     = m_dot #
        self.a_c       = a_c #unit
        self.p_t, self.t_t = p_t, t_t #units. Tank params   
        
        self.specific_volume = specific_volume #
        
        self.a_e, self.p_e, self.t_e, self.velo_e, self.m_e  = a_e, p_e, t_e, velo_e, m_e #units. Wrapped all the exit params
        
        self.inputs = {"Test Date": self.date,
                       "Test Title": self.title,
                       "Test Description": self.doc,
                       "Ambient Pressure": self.p_ambient,
                       
                       "Mass Flow Rate": self.m_dot,
                       
                       "Choked Aera of the nozzle": self.a_c,
                       "Tank Starting Pressure":self.p_t,
                       "Tank Starting Temperature":self.t_t,
                       "Specific Volume":self.specific_volume,
                       "Exit Aera of the nozzle": self.a_e,
                       "Exit Pressure":self.p_e,
                       "Exit Temperature":self.t_e,
                       "Exit Velocity":self.velo_e,
                       "Exit Mach number":self.m_e
                       
                      }
        
        #00000 input here about new param
                       
                       "Test Date": self.date,
                       "Test Title": self.title,
                       "Test Description": self.doc,
                       "Ambient Pressure": self.p_ambient,
                       "Mass Flow Rate": self.m_dot,
                       "Choked Aera of the nozzle": self.a_c,
                       "Tank Pressure":self.p_t,
                       "Tank Temperature":self.t_t,
                        "Specific Volume":self.specific_volume,
                        
                       "Exit Aera of the nozzle": self.a_e,
                       "Exit Pressure":self.p_e,
                       "Exit Temperature":self.t_e,
                        "Exit Velocity":self.velo_e,
                       "Exit Mach number":self.m_e
        
        if self.a_e <= self.a_c:
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
                float(input ("""Input a value of Ambient Pressure : """)),
                float(input ("""Input a Exit Area of a Nozzle in numbers : """)),
                
                
                float(input ("""Input a Exit Area of a Nozzle in numbers : """)),
                float(input ("""Input a Choked Area of a Nozzle in numbers : """)),
                float()
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
    
    def calorically_imperf(cv_perf, cp_perf, gam_perf):
        #gam_perf is the ratio of heat capacity for calorically perfect gas condition; theta is the thermal constant equal to 5500 Rankine, T is the static temperature
        #The specific heat capacity at constant volume
        cv = (cv_perf) * (1 + (gam_perf - 1) * [(theta/T)^2 * np.exp(theta/T) /(np.exp(theta/T) -1)^2])
        #The specific heat capacity at constant pressure
        cp = (cp_perf) * (1 + ((gam_perf - 1)/gam_perf) * [(theta/T)^2 * np.exp(theta/T) /(np.exp(theta/T) -1)^2])
        #Ratio of specific heat
        gam = 1 + (gam_perf - 1) / ( 1 + (gam_perf-1) * [(theta/T)^2 * np.exp(theta/T) /(np.exp(theta/T) -1)^2])
        return gam
        
    #temperature, pressure, and thrust curve vs time
    def t_p_f_vs_t(self):
        #m_dot = 
        return
    
    







