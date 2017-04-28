import numpy as np
import healpy as hp
from math import *
from astropy.table import Table
from optparse import OptionParser
import sys
import os
import sample_selection
from util import *
'''
This script will make a sub-selection from an input fits or ascii file, according to one of the selections in sample_selection.
Optionally, it will perform a masking step according to a bad regions mask and depth map.
Run this for every z bin you want to make (include all bins in sample_selection.py).
You will want to modify the path to inputs and outputs, as well as mask/map filenames in 
the LSS sample class definition.
Usage: mk_LSS_sample.py [options]
Options:
  -h, --help            show this help message and exit
  -f FILENAME, --file=FILENAME
                        Read data from FILENAME
  --mask                Mask catalog
  --depth=DEPTH_LIMIT   Depth limit for masking
  --sample=SAMPLE_NAME  Sample name
  --quiet               Remove verbosity
'''

class lss_sample:
    def __init__(self,filename,dl,sample):
        inputdir = '/data1/des/y1a1/gold/y1a1_gold_1.0_masks/' #directory for maps
        outdir = '/scratch/sevilla/' #directory for output data
        self.inputfile = filename #change to whatever input file is
        self.outputfile = outdir+'lss_sample.fits' #change to whatever input file is
        self.badmask = inputdir+'y1a1_gold_1.0_wide_badmask_4096.fit'
        self.depthmap = inputdir+'y1a1_gold_1.0.2_wide_auto_nside4096_i_10sigma.fits.gz'
        self.depth_limit = dl
        self.sample_name = sample
        self.data = Table.read(self.inputfile)
        for oldname in self.data.colnames:
            if oldname != oldname.upper():
                self.data.rename_column(oldname,oldname.upper())

    def dump_info(self):
        print "Input file name",self.inputfile
        print "Output file name",self.outputfile
        print "Bad regions file name",self.badmask
        print "Depth map file name",self.depthmap
        print "Depth limit",self.depth_limit
        print "Sample name",self.sample_name
        print "Input catalog size",len(self.data)

    def mask(self):
        if not os.path.isfile(self.badmask):
            print 'Bad mask file',self.badmask,'not found.'
            return 2
        if not os.path.isfile(self.depthmap):
            print 'Depth map file',self.badmask,'not found.'            
            return 2
        bad_mask = hp.read_map(self.badmask,verbose=False)
        depth_map = hp.read_map(self.depthmap,verbose=False)
        ra = self.data['RA']
        dec = self.data['DEC']
        theta,phi = radec2thphi(ra,dec)
        if len(bad_mask) != len(depth_map):
            print 'Incompatible mask and depth maps being used, skipping masking.'
            return 3
        nside = hp.npix2nside(len(bad_mask))
        pix = hp.ang2pix(nside,theta,phi,nest=False)
        good = np.empty(len(self.data),dtype=bool)
        for counter,i in enumerate(pix):
            if counter%10000000 == 0:
                print 'Masked',counter,'objects'
            if (bad_mask[i] == 0) & (depth_map[i] > self.depth_limit):
                good[counter] = True 
        print 'Masked from',len(self.data),'to',len(self.data[good])
        self.data = self.data[good]
        masked_filename = os.path.splitext(self.inputfile)[0]+'_masked.fits'
        write_catalog(self.data,masked_filename)
        return 0

    def cut(self):
        mask = sample_selection.sample_cuts(self.data,self.sample_name)
        print 'Cut from',len(self.data),'to',len(self.data[mask])
        self.data = self.data[mask]
        return 0

    def write(self):
        write_catalog(self.data,self.outputfile)
        return 0

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]. Use %prog --help to see available options."
    parser = OptionParser(usage=usage)
    parser.set_defaults(filename='/data1/des/y1a1/gold/y1a1_gold_lss_sample_wadatpc_masked.fits',
                        toggle_mask=False,depth_limit=22.0,sample_name="OFFICIAL_RED",toggle_quiet=False)
    parser.add_option("-f","--file", action="store", type="string", dest="filename", help="Read data from FILENAME")
    parser.add_option("--mask", action="store_true", dest="toggle_mask", help="Mask catalog")
    parser.add_option("--depth", action="store", type="float", dest="depth_limit", help="Depth limit for masking")
    parser.add_option("--sample", action="store", type="string", dest="sample_name", help="Sample name")
    parser.add_option("--quiet", action="store_true", dest="toggle_quiet", help="Remove verbosity")
    ### Parse command line
    (options, args) = parser.parse_args()
    ### Arguments
    if len(args) !=0 :
        errormsg = "Incorrect number of arguments!"
        parser.error(errormsg)
        return 1
    if not options.toggle_quiet:
        print 'Reading catalog from file:',options.filename
    if os.path.isfile(options.filename):
        lss = lss_sample(options.filename,options.depth_limit,options.sample_name)
        if not options.toggle_quiet:
            lss.dump_info()
    else:
        print "Bad input catalog filename",options.filename
        return 2
    if options.toggle_mask:
        if not options.toggle_quiet:
            print 'Masking...'
        err = lss.mask() #create masked sample cutting at depth limit
        if err != 0:
            return err
    if not options.toggle_quiet:
        print 'Selecting...'
    err = lss.cut() #cut sample
    err = lss.write() #write new catalog
    print 'DONE'
    return 0

if __name__ == '__main__':
    sys.exit(main())
