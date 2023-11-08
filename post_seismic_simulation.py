#!/usr/bin/env python
# coding: utf-8
# %%

# post earthquake landslide failed. this code is used to simulate the failed hillslope in given time window. 



import numpy as np

import copy
# there has some default input parameter values, such as cov for cohesion and friction angle, thickness, the number of MC simulation times.  
# if we don't assgin a value to run model, the model will use these default value
# we defined "__init__" to recieve input
# we defined "calculate_ky" to create random cohesion and friction angle, and use these values to calculate FS. From all FS values, we select the couples that suit failed condition (FS<1) 
class po_landslide_strength_inversion():
    def __init__(self, data, number_of_iterations=5000,
             cov_c=0.2,
             cov_i_f=0.15,
             tick=2,
             seed=0,
             l=1000,
             **kwds):
        self._seed_generator(seed)
        self.n = int(number_of_iterations)
        self.cov_c= cov_c
        self.cov_i_f= cov_i_f
        self.t=tick
        self.l=l
        self.data=data
    def calculate_ky(self, i):
        self._slope = self.data['slope'][i]
        self._den = self.data['density'][i]
        self._C_mode = self.data['rock_cohesion'][i]
        self._phi_mode = self.data['internal_friction_angle'][i]
        self._m=self.data['m_after'][i]
        self._C = np.random.normal(self._C_mode,self._C_mode*self.cov_c,size=self.n)
        self._phi = np.random.normal(self._phi_mode,self._phi_mode*self.cov_i_f,size=self.n)
        # we use a circulation to ensure the distribution not contain value that <0, because it is meanless
        for i in range(len(self._C)):
            while self._C[i] <= 0:
                self._C[i]=np.random.normal(self._C_mode,self._C_mode*self.cov_c)
        for i in range(len(self._phi)):    
            while self._phi[i] <= 0:
                self._phi[i]=np.random.normal(self._phi_mode,self._phi_mode*self.cov_i_f)
        self._tan_phi=np.tan(self._phi*np.pi/180.)
        self._thi=self.t
        # calculate the FS
        self._FS = self._C/(np.sin(self._slope*np.pi/180.)*self._den*self._thi)+(self._tan_phi/np.tan(self._slope*np.pi/180.))*(1-self._m*9.8/self._den) 
        # filtering cohesion and friction angle according to failed condition (FS<1), and calculate mean
        self._mean_phi=np.mean(self._phi[self._FS<1])
        self._mean_C=np.mean(self._C[self._FS<1])
        # a try to find limit condition, useless in this study 
        if len([x for x in self._FS[self._FS<1]])>=1:
            self._li_phi=np.mean(self._phi[self._FS==np.min([x for x in self._FS[self._FS<1]])])
            self._li_C=np.mean(self._C[self._FS==np.min([x for x in self._FS[self._FS<1]])])
        else:
            self._li_phi=-9999.0
            self._li_C=-9999.0
    # repeat 1000 times (default value) of above steps to output       
    def calculate_inversion_strength(self, **kwds):  
        c=np.random.rand(len(self.data['slope']),self.l)
        fhi=np.random.rand(len(self.data['slope']),self.l)
        c1=np.random.rand(len(self.data['slope']),self.l)
        fhi1=np.random.rand(len(self.data['slope']),self.l)
        for i in range(len(self.data['slope'])):
            for x in range(self.l): 
                self.calculate_ky(i)
                c[i,x]=self._mean_C
                fhi[i,x]=self._mean_phi
                c1[i,x]=self._li_C
                fhi1[i,x]=self._li_phi
        return c, fhi, c1, fhi1
  
    def _seed_generator(self, seed=0):
        """Seed the random-number generator. This method will create the same
        sequence again by re-seeding with the same value (default value is
        zero). To create a sequence other than the default, assign non-zero
        value for seed.
        """
        np.random.seed(seed)    
        


