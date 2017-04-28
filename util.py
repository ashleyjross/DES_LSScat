import os
import numpy as np
def write_catalog(data,filename):
    if os.path.isfile(filename):
        print 'Removing previous',filename
        os.remove(filename)    
    data.write(filename,format="fits")
    print 'Wrote',filename
def radec2thphi(ra,dec):
    theta = (90.0 - dec)*np.pi/180.
    phi = ra*np.pi/180.
    return theta,phi
