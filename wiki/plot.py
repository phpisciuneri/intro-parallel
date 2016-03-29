# -*- coding: utf-8 -*-
"""
Created on Tue Mar 04 11:34:20 2014

@author: patrick
"""

import numpy as np
import matplotlib.pyplot as plt

# read points
profile = np.loadtxt( "profile.txt" )
x = profile[:,0]
T = profile[:,1]

plt.plot( x, T, lw=3 )
    
# since we are plotting geometry the aspect ratio needs to be fixed
plt.axes().set_xlim(0,2*np.pi)    
plt.axes().set_ylim(-1,3)    
plt.axes().set_aspect( 'equal' )

# label figure
plt.xlabel( 'x', fontsize=16, fontweight='bold' )
plt.ylabel( 'T', fontsize=16, fontweight='bold' )

# plot
plt.box('on')
plt.grid('on')
plt.show()