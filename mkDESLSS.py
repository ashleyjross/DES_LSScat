import fitsio
from healpix import thphi2radec,radec2thphi#,healpix,ang2pix_ring,pix2ang_ring# #comes from AJR's healpix routines
try:
    import healpy as hp
    hpm = True
except:
    print 'no healpy, this will cause problems '
    hpm = False
#from lambert import lb2rd
from math import *
import numpy
from numpy import loadtxt as load
from numpy import array,zeros
from cubspline import *
import pylab as plb
from random import random

inputdir = '/Users/ashleyross/DESY1/' #directory for input data
maskdir = '/Users/ashleyross/DESY1/Y1masks/'
#inputfile = 'BAO_red_v1622.ssv' #this one was used for VF
#inputfile = 'BAO_red_v1622.ssv' #this one was used for VF
#inputfile = 'Y1_LSS_AUTO_MOF_BPZv1.BPZ_mof_orig_sva1prior_priormag_mof_i.fits'
inputfile = 'LSS_Y1_PZ_merge_v1.1.csv'
outf = inputfile.split('.')
#define how to read data from input file
exten = inputfile.split('.')[-1]
print exten
if exten == 'csv':
    spls = ','
if exten == 'ssv':
    spls = ' '
if exten == 'dat':
    spls = ' '
if exten == 'fits':
    spls = ' '
    exten = 'dat'
footmask = maskdir+'y1a1_gold_1.0.2_wide_footprint_4096.fit'
#badmask = maskdir+'y1a1_gold_1.0.2_wide_badmask_4096.fit'
badmask = inputdir+'y1a1_gold_1.0.3_wide_badmask_4096.fits.gz'
depthmap = maskdir+'y1a1_gold_1.0.2_wide_auto_nside4096_i_10sigma.fits'
badammask = inputdir+'Y1A1NEW_COADD_SPT_AIRMASSge1_band_z_nside4096_oversamp4_AIRMASS__mean.fits.gz'
#fracdetmap = inputdir+'Y1A1_SPT_and_S82_IMAGE_SRC_band_i_nside4096_oversamp4_count__fracdet.fits.gz'
fracdetmap = maskdir +'Y1A1_WIDE_frac_combined_griz_o.4096_t.32768_EQU.fits.gz'
maskredcut = maskdir + 'Y1LSSmask_2rz23.7.dat' #depth mask for other bands involved in color cuts
starmap = 'y1a1_gold_1.0.2_stars_nside0512.fits'
outdir = '/Users/ashleyross/DESY1/test/'
maskout = outdir+'Y1LSSmask_v2_redlimcut_il22_seeil4.0_4096ring.dat' #mask created by Jack/Martin including redlimcut

photoz = 'MEAN_Z_DNF_MOF' #column name in data base for what is used to split data by redshift, etc.
photozmc = 'Z_MC_DNF_MOF' #column name for mc draw for photoz from its pdf
#photoz = 'MEAN_Z_BPZv1_MOF_ORIG' #column name in data base for what is used to split data by redshift, etc.
#photozmc = 'Z_MC_BPZv1_MOF_ORIG' #column name for mc draw for photoz from its pdf
#photoz = 'mean_z_bpz' #column name in data base for what is used to split data by redshift, etc.
#photozmc = 'z_mc_bpz' #column name for mc draw for photoz from its pdf
dl = 22
fracd = .8
mf = 'LSSmask_mof_redlimcut_il'
#seeic = 4. #was used for i band cut
seeic = 3.7 #testing for z-band cute
syscut = 'seezl'+str(seeic)
#syscut = 'seeil'+str(seeic) #use if constructing full mask from scratch
#syscut = '_seeil4.0_am' #use if reconstructing from previous full mask file
#syscut = '' #use if not cutting on seeing
#mf='goldFoot_fracdet4TEST'+str(fracd)+'wide_iautodepth' #used a few times downstream
header = open(inputdir+inputfile).readline().split(spls)
test = ''#'test' #toggle this to change whether it is the test file being worked on

def maskY1_tot(dl=22,fracd=.8):
    #combine masks
    np = 4096*4096*12
    amm = fitsio.read(badammask)
    fam = []
    frc = []
    for i in range(0,np):
        fam.append(0)
        frc.append(0)
    for i in range(0,len(amm)):
        fam[amm[i]['PIXEL']] = -1
    mrc = open(maskredcut).readlines()
    for i in range(0,len(mrc)):
        p = int(mrc[i])
        frc[p] = 1
    if hpm == True:
        mg = hp.read_map(badmask)
        dpm = hp.read_map(depthmap)
        fm = hp.read_map(footmask)
        fdm = hp.read_map(fracdetmap)
    
    #fdf = fitsio.read(fracdetmap)
    #for i in range(0,len(fdf)):
    #    fdm[fdf[i]['PIXEL']] = fdf[i]['SIGNAL']

    fo = open(outdir+'Y1'+mf+str(dl)+'4096ring.dat','w')
    print outdir+'Y1'+mf+str(dl)+'4096ring.dat'
    n = 0
    for i in range(0,np):
        if dpm[i] > dl and mg[i] <= 3 and fm[i] >= 1 and fdm[i] > .8 and fam[i] != -1 and frc[i] == 1 and int(mg[i]) & 512 <=0:
            th,phi = hp.pix2ang(4096,i)
            ra,dec = thphi2radec(th,phi)
            if ra > 50 or dec < -30:
                n += fdm[i]
                fo.write(str(i)+' '+str(fdm[i])+'\n')
    fo.close()
    ndeg = 360*360./pi*n/(12.*4096*4096)
    print n,ndeg
    return True

def maskY1_totfilepam():
    np = 4096*4096*12
    amm = fitsio.read(badammask)
    fam = []
    for i in range(0,np):
        fam.append(0)
    for i in range(0,len(amm)):
        fam[amm[i]['PIXEL']] = -1
    mout = open(maskout).readlines()
    for i in range(0,len(mout)):
        p = int(mout[i].split()[0])
        if fam[p] != -1:
            fam[p] = 1
    fo = open(outdir+'Y1'+mf+str(dl)+'_seeil'+str(seeic)+'_am4096ring.dat','w')
    for i in range(0,len(fam)):
        if fam[i] == 1:
            fo.write(str(i)+'\n')
    fo.close()

def maskY1see(seeic=seeic):
    d = numpy.loadtxt(outdir+'Y1'+mf+str(dl)+'4096ring.dat').transpose()
    f = fitsio.read(inputdir+'Y1A1NEW_COADD_SPT_band_i_FWHM/Y1A1NEW_COADD_SPT_band_i_nside4096_oversamp4_FWHM_MEAN_coaddweights3_mean.fits.gz')
    seemap = mkseemap()
#    np = 12*4096**2
#    for i in range(0,np):
#        seemap.append(0)
#    for i in range(0,len(f)):
#        p = f[i]['PIXEL']
#        seemap[p] = f[i]['SIGNAL']
    fo = open(outdir+'Y1'+mf+str(dl)+'seeil'+str(seeic)+'4096ring.dat','w')
    print outdir+'Y1'+mf+str(dl)+'seeil'+str(seeic)+'4096ring.dat'
    for i in range(0,len(d[0])):
        if seemap[int(d[0][i])] < seeic:
            fo.write(str(d[0][i])+' '+str(d[1][i])+'\n')
    fo.close()
    return True

def maskY1seez(seeic=seeic):
    d = numpy.loadtxt(outdir+'Y1'+mf+str(dl)+'4096ring.dat').transpose()
    #f = fitsio.read(inputdir+'Y1A1NEW_COADD_SPT_band_z/Y1A1NEW_COADD_SPT_band_z_nside4096_oversamp4_FWHM_MEAN_coaddweights3_mean.fits.gz')
    seemap = mkseemap(bnd='z')
    #    np = 12*4096**2
    #    for i in range(0,np):
    #        seemap.append(0)
    #    for i in range(0,len(f)):
    #        p = f[i]['PIXEL']
    #        seemap[p] = f[i]['SIGNAL']
    fo = open(outdir+'Y1'+mf+str(dl)+'seezl'+str(seeic)+'4096ring.dat','w')
    print outdir+'Y1'+mf+str(dl)+'seezl'+str(seeic)+'4096ring.dat'
    for i in range(0,len(d[0])):
        if seemap[int(d[0][i])] < seeic:
            fo.write(str(d[0][i])+' '+str(d[1][i])+'\n')
    fo.close()
    return True


def maskd(res,dl=22):
    #degrade mask
    f = open(outdir+'Y1'+mf+str(dl)+syscut+'4096ring.dat').readlines()
    mo = []
    for i in range(0,12*res*res):
        mo.append(0)
    frac = (res/4096.)**2.
    for i in range(0,len(f)):
        ln = f[i].split()
        p = int(float(ln[0]))
        mv = float(ln[1])
        if mv > 0:
            th,phi = hp.pix2ang(4096,p)
            a = frac*mv
            po = hp.ang2pix(res,th,phi)
            mo[po] += a
    fo = open(outdir+'Y1'+mf+str(dl)+syscut+str(res)+'ring.dat','w')
    for i in range(0,len(mo)):
        if mo[i] > 0:
            fo.write(str(i)+' '+str(mo[i])+'\n')
    fo.close()
    return True

def mksampfitsBPZ():
    import fitsio
    f = fitsio.read(inputdir+inputfile)
    maskf = open(outdir+'Y1'+mf+str(dl)+syscut+'4096ring.dat')
    print outdir+'Y1'+mf+str(dl)+syscut+'4096ring.dat'
    np = 12*4096*4096
    mask = []
    for i in range(0,np):
        mask.append(0)
    for line in maskf:
        pix = int(float(line.split()[0]))
        mask[pix] = 1
    outf = inputfile.split('.')
    fo = open(outdir+'Y1red'+outf[0]+photoz+test+'.dat','w')
    for i in range(0,len(f)):
        z = f[i]['MEAN_Z']
        magi = f[i]['mag_auto_i']
        if magi < 19.0 + 3.0*z and z > 0.6 and z < 1.1:
            ra,dec = f[i]['ra'],f[i]['dec']
            th,phi = radec2thphi(ra,dec)
            p = hp.ang2pix(4096,th,phi)
            if mask[p] == 1:
                fo.write(str(ra)+' '+str(dec)+' '+str(z)+' '+str(f[i]['Z_MC'])+' '+str(f[i]['coadd_objects_id'])+'\n')
    fo.close()
    return True

class mksample_merge:
    def __init__(self):
        self.findcolumns()
    # find column numbers for important quantities EVERYONE DOUBLE CHECK THESE!
    def findcolumns(self):
        self.col_ra = findcol('RA')
        self.col_dec = findcol('DEC')
        self.col_gmag_auto = findcol('MAG_AUTO_G')
        self.col_rmag_auto = findcol('MAG_AUTO_R')
        self.col_imag_auto = findcol('MAG_AUTO_I')
        self.col_zmag_auto = findcol('MAG_AUTO_Z')
        self.col_pz = findcol(photoz)
        self.col_pzmc = findcol(photozmc)
        self.col_id = findcol('COADD_OBJECTS_ID')
    
    def createmaskedY1redsimp(self):
        #mask file and apply imag cut
        maskf = open(outdir+'Y1'+mf+str(dl)+syscut+'4096ring.dat')
        print outdir+'Y1'+mf+str(dl)+syscut+'4096ring.dat'
        np = 12*4096*4096
        mask = []
        for i in range(0,np):
            mask.append(0)
        for line in maskf:
            pix = int(float(line.split()[0]))
            mask[pix] = 1
        outf = inputfile.split('.')
        fo = open(outdir+'Y1red'+outf[0]+photoz+syscut+'.dat','w')
        fbo = open(outdir+'bad'+photoz+'.dat','w')
        f = open(inputdir+inputfile)
        f.readline() #assumes input file has one line header
        nb = 0
        ng = 0
        for line in f:
            ln = line.split(spls)
            ra,dec = float(ln[self.col_ra]),float(ln[self.col_dec])
            th,phi = radec2thphi(ra,dec)
            p = hp.ang2pix(4096,th,phi)
            if mask[p] == 1:
                #magi = float(ln[self.col_imag_auto])
                try:
                    pz = float(ln[self.col_pz])
                    #ng += 1
                    #fo.write(ln[self.col_ra]+' '+ln[self.col_dec]+' '+str(pz)+' '+ln[self.col_pzmc].strip('\n')+' '+ln[self.col_id]+'\n')
                except:
                    #print ln[self.col_pz]
                    nb += 1
                    print nb
                    fbo.write(str(ra)+' '+str(dec)+'\n')
                    pz = 0
                    pass
                #if magi < 19.0 + 3.0*pz and pz > .6 and pz < 1.:
                if pz > .6 and pz < 1.:
                    ng += 1
                    fo.write(ln[self.col_ra]+' '+ln[self.col_dec]+' '+str(pz)+' '+ln[self.col_pzmc].strip('\n')+' '+ln[self.col_id]+'\n')

        print ng
        fo.close()
        return True


class mksample:
    def __init__(self):
        self.findcolumns()
    # find column numbers for important quantities EVERYONE DOUBLE CHECK THESE!
    def findcolumns(self):
        self.col_ra = findcol('ra')
        self.col_dec = findcol('dec')
        self.col_gmag_auto = findcol('mag_auto_g')
        self.col_rmag_auto = findcol('mag_auto_r')
        self.col_imag_auto = findcol('mag_auto_i')
        self.col_zmag_auto = findcol('mag_auto_z')
        self.col_spreadmod = findcol('spread_model_i')
        self.col_spreaderr = findcol('spreaderr_model_i')
        self.col_pz = findcol(photoz)
        self.col_pzmc = findcol(photozmc)
        self.col_type = findcol('t_b')
        self.col_rfpz = findcol('photo_z_rf')
        self.col_id = findcol('coadd_objects_id')



    def creatmaskedY1sampgoldc3pz(self,dl):
        #mg = hp.read_map(badmask)
        #dpm = hp.read_map(depthmap)
        #fo = open(outdir+'Y1gold3pz_m'+str(dl)+'.dat','w')
        maskf = open(outdir+'Y1'+mf+str(dl)+'4096ring.dat')
        np = 12*4096*4096
        mask = []
        for i in range(0,np):
            mask.append(0)
        for line in maskf:
            pix = int(line.split()[0])
            mask[pix] = 1
        outf = inputfile.split('.')
        fo = open(outdir+outf[0]+'_m22.'+outf[1],'w')
        f = open(inputdir+inputfile)
        f.readline() #assumes input file has one line header
        for line in f:
            ln = line.split(spls)
            im = float(ln[self.col_imag_auto])
            if im < dl:
                ra,dec = float(ln[self.col_ra]),float(ln[self.col_dec])
                th,phi = radec2thphi(ra,dec)
                p = hp.ang2pix(4096,th,phi)
                if mask[p] == 1:
                    fo.write(line)
        fo.close()
        return True

    def mkRedsampradecz(self,dl=22,ccut=0.9,sccut=.007):
        f = open(outdir+outf[0]+'_m22.'+outf[1])
        fo = open(outdir+'Y1red'+outf[0]+photoz+'.dat','w')
        n = 0
        for line in f:
            ln = line.split(spls)
            g = float(ln[self.col_gmag_auto])
            r = float(ln[self.col_rmag_auto])
            i = float(ln[self.col_imag_auto])
            z = float(ln[self.col_zmag_auto])
            gmr = g-r
            rmi = r-i
            imz = i-z
            if gmr > -1 and gmr < 3 and rmi > -1 and rmi < 2.5 and imz > -1 and imz < 2: #crazy colors
                sc = float(ln[self.col_spreadmod])+5/3.*float(ln[self.col_spreaderr]) #star/galaxy separation criteria
            #sc = float(ln[self.col_spreadmod])+1.6*float(ln[self.col_spreaderr])
            #sc = 1.
                if sc > sccut and i > 17.5: #s/g cut
                    cred = imz+.6*rmi

                    if cred > ccut: #red color cut
                        zp = float(ln[self.col_pz])
                    #if zp >= 0.6 and zp <= 1.0: #BAO range
                        fo.write(ln[self.col_ra]+' '+ln[self.col_dec]+' '+str(zp)+' '+ln[self.col_pzmc]+'\n')
                        n += 1
#        r =  read_in_file( outdir+outf[0]+'_m22.'+outf[1], 4, [])
#        i =  read_in_file( outdir+outf[0]+'_m22.'+outf[1], 5, [])
#        z =  read_in_file( outdir+outf[0]+'_m22.'+outf[1], 6, [])
#        red_relation = 0.6*(r-i) + (i-z)
#        sindex = numpy.where( red_relation > 0.9 )[0]
#        ra = read_in_file( outdir+outf[0]+'_m22.'+outf[1], 1, sindex)
#        
#        print n,len(ra)
        return True

    def mkRedsampradecztest(self,dlr=22,ccut=0.9,acut=0.6,sccut=.007,islp=0,rmimin=-1,zsld=100,rmimax=2.5,rmzmax=10,tpmin=-1,tpmax=10,gmin=0,zmmax=30,drfmax=100):
        f = open(outdir+outf[0]+'_m22.'+outf[1])
        fo = open(outdir+'Y1red'+outf[0]+photoz+'test.dat','w')
        n = 0
        for line in f:
            ln = line.split(spls)
            g = float(ln[self.col_gmag_auto])
            r = float(ln[self.col_rmag_auto])
            i = float(ln[self.col_imag_auto])
            z = float(ln[self.col_zmag_auto])
            gmr = g-r
            rmi = r-i
            imz = i-z
            rmz = r-z
            tp = float(ln[self.col_type])
            if gmr > -1 and gmr < 3 and rmi > rmimin and rmi < rmimax and imz > -1 and imz < 2 and rmz < rmzmax and tp>tpmin and tp < tpmax and z < zmmax: #crazy colors
                #sc = float(ln[self.col_spreadmod])+5/3.*float(ln[self.col_spreaderr]) #star/galaxy separation criteria
                sc = float(ln[self.col_spreadmod])+1.6*float(ln[self.col_spreaderr])
                #sc = 1.
                if sc > sccut and i > 17.5 and i < dlr and g > gmin: #s/g cut
                    cred = imz+acut*rmi
                    
                    if cred > ccut: #red color cut
                        zp = float(ln[self.col_pz])
                        drf = 0
                        if drfmax != 100:
                            drf = abs(zp-float(ln[self.col_rfpz]))
                        
                        #if zp >= 0.6 and zp <= 1.0: #BAO range
                        if drf < drfmax and i < zsld+zp*islp:
                            fo.write(ln[self.col_ra]+' '+ln[self.col_dec]+' '+str(zp)+' '+ln[self.col_pzmc]+'\n')
                        n += 1
        #        r =  read_in_file( outdir+outf[0]+'_m22.'+outf[1], 4, [])
        #        i =  read_in_file( outdir+outf[0]+'_m22.'+outf[1], 5, [])
        #        z =  read_in_file( outdir+outf[0]+'_m22.'+outf[1], 6, [])
        #        red_relation = 0.6*(r-i) + (i-z)
        #        sindex = numpy.where( red_relation > 0.9 )[0]
        #        ra = read_in_file( outdir+outf[0]+'_m22.'+outf[1], 1, sindex)
        #        
        #        print n,len(ra)
        return True

    def createmaskedY1simp(self):
        maskf = open(outdir+'Y1'+mf+str(dl)+syscut+'4096ring.dat')
        print outdir+'Y1'+mf+str(dl)+syscut+'4096ring.dat'
        np = 12*4096*4096
        mask = []
        for i in range(0,np):
            mask.append(0)
        for line in maskf:
            pix = int(float(line.split()[0]))
            mask[pix] = 1
        outf = inputfile.split('.')
        fo = open(outdir+outf[0]+'_m22.'+outf[1],'w')
        f = open(inputdir+inputfile)
        f.readline() #assumes input file has one line header
        for line in f:
            ln = line.split(spls)
            ra,dec = float(ln[self.col_ra]),float(ln[self.col_dec])
            th,phi = radec2thphi(ra,dec)
            p = hp.ang2pix(4096,th,phi)
            if mask[p] == 1:
                fo.write(line)
        fo.close()
        return True

    def mkRedsampradeczsimp(self):
        f = open(outdir+outf[0]+'_m22.'+outf[1])
        fo = open(outdir+'Y1red'+outf[0]+photoz+'.dat','w')
        n = 0
        for line in f:
            ln = line.split(spls)
            zp = float(ln[self.col_pz])
            #if zp >= 0.6 and zp <= 1.0: #BAO range
            fo.write(ln[self.col_ra]+' '+ln[self.col_dec]+' '+str(zp)+' '+ln[self.col_pzmc]+' '+ln[self.col_id]+'\n')
            n += 1
        return True


def mknz():
#calculate dN/dz for each based on the MC photoz in the red file
    d = numpy.loadtxt(outdir+'Y1red'+outf[0]+photoz+'.dat').transpose()
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        zw = str(zl)+str(zh)
        nz = numpy.zeros((200)) #200 bins for range 0 < z < 2
        ntb = 0
        std = 0
        stdno = 0
        nno = 0
        for j in range(0,len(d[2])):
            z = d[2][j]
            if z >= zl and z < zh:
                zmc = d[3][j]+.01*random()-.005
                if zmc > 0 and zmc < 2:
                    zind = int(100*zmc)
                    nz[zind] += 1.
                    ntb += 1.
                    std += (zmc-z)**2.
                    if abs(zmc-z) < (1.+(zh+zl)/2.)*.1:
                        nno += 1.
                        stdno += (zmc-z)**2.
                else:
                    print z,zmc
        print ntb,sqrt(std/ntb)/(1.+(zh+zl)/2.),nno/ntb,sqrt(stdno/nno)/(1.+(zh+zl)/2.)
        norm = ntb/100.
        fo = open(outdir+'dNdzY1red'+outf[0]+photoz+zw+'.dat','w')
        for j in range(0,len(nz)):
            fo.write(str(.005+j/100.)+' '+str(nz[j]/norm)+'\n')
        fo.close()
    for i in range(0,4):
        zl = 0.6+0.1*i
        zh = zl+.1
        zw = str(zl)+str(zh)
        nz = numpy.zeros((200)) #200 bins for range 0 < z < 2
        ntb = 0
        std = 0
        stdno = 0
        nno = 0
        for j in range(0,len(d[2])):
            z = d[2][j]
            if z >= zl and z < zh:
                zmc = d[3][j]+.01*random()-.005
                if zmc > 0 and zmc < 2:
                    zind = int(100*zmc)
                    nz[zind] += 1.
                    ntb += 1.
                    std += (zmc-z)**2.
                    if abs(zmc-z) < (1.+(zh+zl)/2.)*.1:
                        nno += 1.
                        stdno += (zmc-z)**2.
                else:
                    print z,zmc
        print ntb,sqrt(std/ntb)/(1.+(zh+zl)/2.),nno/ntb,sqrt(stdno/nno)/(1.+(zh+zl)/2.)
        norm = ntb/100.
        fo = open(outdir+'dNdzY1red'+outf[0]+photoz+zw+'.dat','w')
        for j in range(0,len(nz)):
            fo.write(str(.005+j/100.)+' '+str(nz[j]/norm)+'\n')
        fo.close()
    return True

def mknztest():
    #calculate dN/dz for each based on the MC photoz in the red test file
    d = numpy.loadtxt(outdir+'Y1red'+outf[0]+photoz+'test.dat').transpose()
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        zw = str(zl)+str(zh)
        nz = numpy.zeros((200)) #200 bins for range 0 < z < 2
        ntb = 0
        std = 0
        stdno = 0
        nno = 0
        for j in range(0,len(d[2])):
            z = d[2][j]
            if z >= zl and z < zh:
                zmc = d[3][j]+.01*random()-.005
                if zmc > 0 and zmc < 2:
                    zind = int(100*zmc)
                    nz[zind] += 1.
                    ntb += 1.
                    std += (zmc-z)**2.
                    if abs(zmc-z) < (1.+(zh+zl)/2.)*.1:
                        nno += 1.
                        stdno += (zmc-z)**2.
                else:
                    print z,zmc
        print ntb,sqrt(std/ntb)/(1.+(zh+zl)/2.),nno/ntb,sqrt(stdno/nno)/(1.+(zh+zl)/2.)
        norm = ntb/100.
        fo = open(outdir+'dNdzY1red'+outf[0]+photoz+'test'+zw+'.dat','w')
        for j in range(0,len(nz)):
            fo.write(str(.005+j/100.)+' '+str(nz[j]/norm)+'\n')
        fo.close()
    return True


def plotstfits():
#plots the stellar contamination best-fits
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.clf()
    btl = [0.967,0.964,0.931,0.916,0.941,0.978,0.963,0.961]
    mtl = [0.0672,0.0746,0.144,0.173,0.123,0.0476,0.0782,0.0798]
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        zw = str(zl)+str(zh)
        pp = PdfPages(outdir+'Y1red'+outf[0]+zw+photoz+'vnst.pdf')
        b,m = dolinbestfit('Y1red'+outf[0]+zw+photoz+'vstarsjackerr')
        bt,mt = btl[i],mtl[i]
        d = numpy.loadtxt(outdir+'Y1red'+outf[0]+zw+photoz+'vstarsjackerr.dat').transpose()
        dh = numpy.loadtxt(outdir+'Y1red'+outf[0]+zw+photoz+'vstarshist.dat').transpose()
        ds = numpy.loadtxt(outdir+'Y1red'+outf[0]+zw+photoz+'vstarssimp.dat').transpose()
        dt = numpy.loadtxt(outdir+'plotfile512'+str(i)+'.txt').transpose()
        plt.errorbar(d[0],d[1],d[2],fmt='ko',markersize=10,elinewidth=2)
        plt.errorbar(dt[0],dt[1],dt[2],fmt='rd',markersize=10,elinewidth=2)
        plt.plot(dh[0],dh[1],'bs')
        plt.plot(ds[0],ds[1],'g^')
        plt.plot(d[0],m*d[0]+b,'k-',linewidth=3)
        # plt.plot(d[0],mt*d[0]+bt,'r--',linewidth=3)
        plt.xlabel(r'stellar density (arcmin$^{-2}$)',size=16)
        plt.ylabel('galaxy density/mean galaxy density',size=16)
        if b < 1:
            per = str(100.*(1.-b))[:3]
        else:
            per = str(100.*(1.-b))[:4]
        plt.text(0.2,1.2,per+'% stellar contamination',size=16)
        plt.text(0.25,1.25,str(zl)+' < z <'+str(zh),size=16)
        #plt.title('stellar contamination for red sample using BPZ and 0.007 cut')
        plt.ylim(0.9,1.3)
        pp.savefig()
        pp.close()
        plt.clf()
    return True

def plotseefits(xlim=4.):
    #plots the stellar contamination best-fits
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.clf()
    #btl = [0.967,0.964,0.931,0.916,0.941,0.978,0.963,0.961]
    #mtl = [0.0672,0.0746,0.144,0.173,0.123,0.0476,0.0782,0.0798]
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        zw = str(zl)+str(zh)
        pp = PdfPages(outdir+'Y1red'+outf[0]+zw+photoz+'wstvseei.pdf')
        b,m = dolinbestfit('Y1red'+outf[0]+zw+photoz+'wstvseeijackerr',xlim=xlim)
        # bt,mt = btl[i],mtl[i]
        d = numpy.loadtxt(outdir+'Y1red'+outf[0]+zw+photoz+'wstvseeijackerr.dat').transpose()
        #dh = numpy.loadtxt(outdir+'Y1red'+outf[0]+zw+photoz+'vstarshist.dat').transpose()
        #ds = numpy.loadtxt(outdir+'Y1red'+outf[0]+zw+photoz+'vstarssimp.dat').transpose()
        #dt = numpy.loadtxt(outdir+'plotfile512'+str(i)+'.txt').transpose()
        plt.errorbar(d[0],d[1],d[2],fmt='ko',markersize=10,elinewidth=2)
        # plt.errorbar(dt[0],dt[1],dt[2],fmt='rd',markersize=10,elinewidth=2)
        #plt.plot(dh[0],dh[1],'bs')
        #plt.plot(ds[0],ds[1],'g^')
        plt.plot(d[0],m*d[0]+b,'k-',linewidth=3)
        # plt.plot(d[0],mt*d[0]+bt,'r--',linewidth=3)
        plt.xlabel(r'i-band seeing (pixels)',size=16)
        plt.ylabel('galaxy density/mean galaxy density',size=16)
            #if b < 1:
            #per = str(100.*(1.-b))[:3]
            #else:
            #per = str(100.*(1.-b))[:4]
            #plt.text(0.2,1.2,per+'% stellar contamination',size=16)
        plt.text(3.,1.25,str(zl)+' < z <'+str(zh),size=16)
        #plt.title('stellar contamination for red sample using BPZ and 0.007 cut')
        plt.ylim(0.8,1.3)
        pp.savefig()
        pp.close()
        plt.clf()
    return True


def plotvsdepth():
    #plots the stellar contamination best-fits
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.clf()
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        zw = str(zl)+str(zh)
        pp = PdfPages(outdir+'Y1red'+outf[0]+zw+photoz+test+'vdepthi.pdf')
        #b,m = dolinbestfit('Y1red'+outf[0]+zw+photoz+'vstarsjackerr')
        d = numpy.loadtxt(outdir+'Y1red'+outf[0]+zw+photoz+test+'vdepthjackerr.dat').transpose()
        plt.errorbar(d[0],d[1],d[2],fmt='ko',markersize=10,elinewidth=2)
        #plt.plot(d[0],m*d[0]+b,'k-',linewidth=3)
        plt.xlabel(r'i band depth',size=16)
        plt.ylabel('galaxy density/mean galaxy density',size=16)
        #if b < 1:
        #    per = str(100.*(1.-b))[:3]
        #else:
        #    per = str(100.*(1.-b))[:4]
        #plt.text(0.2,1.2,per+'% stellar contamination',size=16)
        plt.text(22.15,1.15,str(zl)+' < z <'+str(zh),size=16)
        if test == '':
            plt.title('depth relationship for red sample using BPZ and 0.007 cut')
        else:
            plt.title('depth relationship for red test sample')
        plt.ylim(0.8,1.2)
        pp.savefig()
        pp.close()
        plt.clf()
    return True

def plotvssys(sys):
    #plots the stellar contamination best-fits
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.clf()
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        zw = str(zl)+str(zh)
        pp = PdfPages(outdir+'Y1red'+outf[0]+zw+photoz+test+'v'+sys+'.pdf')
        #b,m = dolinbestfit('Y1red'+outf[0]+zw+photoz+'vstarsjackerr')
        d = numpy.loadtxt(outdir+'Y1red'+outf[0]+zw+photoz+test+'v'+sys+'jackerr.dat').transpose()
        plt.errorbar(d[0],d[1],d[2],fmt='ko',markersize=10,elinewidth=2)
        #plt.plot(d[0],m*d[0]+b,'k-',linewidth=3)
        plt.xlabel(sys,size=16)
        plt.ylabel('galaxy density/mean galaxy density',size=16)
        #if b < 1:
        #    per = str(100.*(1.-b))[:3]
        #else:
        #    per = str(100.*(1.-b))[:4]
        #plt.text(0.2,1.2,per+'% stellar contamination',size=16)
        #plt.text(22.15,1.15,str(zl)+' < z <'+str(zh),size=16)
        if test == '':
            plt.title(sys+' relationship for red sample using BPZ and 0.007 cut, '+str(zl)+' < z <'+str(zh))
        else:
            plt.title(sys+' relationship for red test sample, '+str(zl)+' < z <'+str(zh))
        plt.ylim(0.8,1.2)
        pp.savefig()
        pp.close()
        plt.clf()
    return True


def writestfit():
    fo = open(outdir+'Y1redstfits'+outf[0]+photoz+'.dat','w')
    fo.write('#zcen intercept slope\n')
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        zw = str(zl)+str(zh)
        b,m = dolinbestfit('Y1red'+outf[0]+zw+photoz+syscut+test+'vstarsjackerr')
        fo.write(str(zl+.025)+' '+str(b)+' '+str(m)+'\n')
    fo.close()
    return True

def writeseefit(xlim=4.):
    fo = open(outdir+'Y1redseeifits'+outf[0]+photoz+'.dat','w')
    fo.write('#zcen intercept slope\n')
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        zw = str(zl)+str(zh)
        b,m = dolinbestfit('Y1red'+outf[0]+zw+photoz+test+'wstvseeijackerr',xlim=xlim)
        fo.write(str(zl+.025)+' '+str(b)+' '+str(m)+'\n')
    fo.close()
    return True


def addseeweight():
    #adds a weight for seeing as a function of redshift
    #1st step is to fit the slope at each redshift
    meansee = 3.45
    seemap = mkseemap()
    cl = []
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        zw = str(zl)+str(zh)
        b,m = dolinbestfit('Y1red'+outf[0]+zw+photoz+test+'wstvseeijackerr',xlim=4.0)
        print b
        cl.append(m)
    zl = [0.59,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975,1.01] #low and high values needed to allow spline to work over full range
    stl = [cl[0],cl[0],cl[1],cl[2],cl[3],cl[4],cl[5],cl[6],cl[7],cl[7]]
    stsp = spline(zl,stl)
    f = open(outdir+'Y1red'+outf[0]+photoz+test+'wst.dat')
    fo = open(outdir+'Y1red'+outf[0]+photoz+test+'wstsee.dat','w')
    nlsee = 0
    wseemax = 1
    seemin = 100
    for line in f:
        ln = line.split()
        z = float(ln[2])
        zmc = float(ln[3])
        w = float(ln[4])
        if z >= 0.6 and z <= 1.0:
            stc = Splev(stsp,z)
            bc = 1.-meansee*stc
            
            ra,dec = float(ln[0]),float(ln[1])
            th,phi = radec2thphi(ra,dec)
            p = hp.ang2pix(4096,th,phi)
            see = seemap[p]
            wsee = 1.
            if see > 0:
                wsee = 1./(bc+stc*see)
            else:
                nlsee += 1.
                if see < seemin:
                    seemin = see
            if wsee > wseemax:
                wseemax = wsee
                seewmak = see
            fo.write(str(ra)+' '+str(dec)+' '+str(z)+' '+str(zmc)+' '+str(w)+' '+str(wsee)+' '+str(w*wsee)+'\n')
    fo.close()
    print nlsee,wseemax,seewmak,seemin
#    putngalvnstar3jackN(zmin=0.75,zmax=0.8,smin=2,smax=4.,res=4096,t=.2,wm='wstsee',pz=photoz,smd='seei',bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    from matplotlib import pyplot as plt
#    from matplotlib.backends.backend_pdf import PdfPages
#    plt.clf()
#    zw = '0.750.8'
#    d = numpy.loadtxt(outdir+'Y1red'+outf[0]+zw+photoz+test+'wstvseeijackerr.dat').transpose()
#    dw = numpy.loadtxt(outdir+'Y1red'+outf[0]+zw+photoz+test+'wstseevseeijackerr.dat').transpose()
#    plt.errorbar(d[0],d[1],d[2],fmt='ko',markersize=10,elinewidth=2)
#    plt.errorbar(dw[0],dw[1],dw[2],fmt='rd',markersize=10,elinewidth=2)
#    plt.ylabel('galaxy density/mean galaxy density',size=16)
#    plt.show()
    return True

def addseeweightone(xlim=4.5):
    #adds a weight for seeing as a function of redshift
    #1st step is to fit the slope at each redshift
    meansee = 3.44
    seemap = mkseemap()
    zl = 0.6
    zh = 1.
    zw = str(zl)+str(zh)
    b,m = dolinbestfit('Y1red'+outf[0]+zw+photoz+'wst'+syscut+test+'vseeijackerr',xlim=xlim)
    print b
    cl = m
    f = open(outdir+'Y1red'+outf[0]+photoz+syscut+test+'wst.dat')
    fo = open(outdir+'Y1red'+outf[0]+photoz+syscut+test+'wstsee.dat','w')
    fo.write('#ra dec zmean zmc wstar*wsee wsee wstar coadd_objects_id\n')
    nlsee = 0
    wseemax = 1
    seemin = 100
    for line in f:
        ln = line.split()
        z = float(ln[2])
        zmc = float(ln[3])
        w = float(ln[4])
        if z >= 0.6 and z <= 1.0:
            bc = 1.-meansee*m
            
            ra,dec = float(ln[0]),float(ln[1])
            th,phi = radec2thphi(ra,dec)
            p = hp.ang2pix(4096,th,phi)
            see = seemap[p]
            wsee = 1.
            if see > 0:
                wsee = 1./(bc+m*see)
            else:
                nlsee += 1.
                if see < seemin:
                    seemin = see
            if wsee > wseemax:
                wseemax = wsee
                seewmak = see
            fo.write(str(ra)+' '+str(dec)+' '+str(z)+' '+str(zmc)+' '+str(w*wsee)+' '+str(wsee)+' '+str(w)+' '+str(int(ln[-1]))+'\n')
    fo.close()
    print nlsee,wseemax,seewmak,seemin
    return True



def addstweight():
#adds a weight for stellar contamination as a function of redshift
#1st step is to get contamination in each z bin
    cl = []
    ml = []
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        zw = str(zl)+str(zh)
        b,m = dolinbestfit('Y1red'+outf[0]+zw+photoz+syscut+test+'vstarsjackerr',xlim=1.8)
        print b
        cl.append(1.-b)
        ml.append(m)
#    cl[0] += .01
#    cl[1] += .01
#    cl[2] += .03
#    cl[3] += .03
#    cl[4] += .035
#    cl[5] += .035
#    cl[6] += .035
#    cl[7] += .035
    zl = [0.59,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975,1.01] #low and high values needed to allow spline to work over full range
    stl = [cl[0],cl[0],cl[1],cl[2],cl[3],cl[4],cl[5],cl[6],cl[7],cl[7]]
    stsp = spline(zl,stl)
    mtl = [ml[0],ml[0],ml[1],ml[2],ml[3],ml[4],ml[5],ml[6],ml[7],ml[7]]
    msp = spline(zl,mtl)
#test spline
#    from matplotlib import pyplot as plt
#    zt = numpy.arange(.6,1.,.001)
#    stt = []
#    for i in range(0,len(zt)):
#        #stt.append(Splev(zt[i],stsp))
#        stt.append(Splev(stsp,zt[i]))
#    plt.plot(zt,stt,'k-',zl,stl,'ko')
#    plt.show()

#now calculate average number of stars
    mask = open(outdir+'Y1'+mf+str(dl)+syscut+'512ring.dat').readlines()
    ml = []
    for i in range(0,12*512*512):
        ml.append(0.0)
        np = 0
    for i in range(0,len(mask)):
        ln = mask[i].split()
        p = int(ln[0])
        fr = float(ln[1])
        ml[p] = fr
    starl = hp.read_map(inputdir+starmap)
    nst = 0
    summl = 0
    for i in range(0,len(starl)):
        if starl[i] >= 0:
            nst += starl[i]*ml[i]
            summl += ml[i]
    stave = nst/summl
    print summl,sum(ml)
    print stave
    f = open(outdir+'Y1red'+outf[0]+photoz+syscut+test+'.dat')
    fo = open(outdir+'Y1red'+outf[0]+photoz+syscut+test+'wst.dat','w')
    for line in f:
        ln = line.split()
        z = float(ln[2])
        zmc = float(ln[3])
        if z >= 0.6 and z <= 1.0:
            stc = Splev(stsp,z)
            slp = stc/stave
           #slp = Splev(msp,z)
            ra,dec = float(ln[0]),float(ln[1])
            th,phi = radec2thphi(ra,dec)
            p = hp.ang2pix(512,th,phi)
            if ml[p] == 0:
                print 0
            nst = starl[p]
            w = 1.
            if nst >= 0:
                w = 1./((1.-stc)+slp*nst)
            fo.write(str(ra)+' '+str(dec)+' '+str(z)+' '+str(zmc)+' '+str(w)+' '+str(int(float(ln[-1])))+'\n')
    fo.close()
    return True

def addgdepthweight(gmap):
    #adds a weight for g band depth as a function of redshift
    #1st step is to get coefficients in each z bin
    cl = []
    ml = []
    for i in range(0,4):
        zl = 0.6+0.1*i
        zh = zl+.1
        zw = str(zl)+str(zh)
        b,m = dolinbestfit('Y1red'+outf[0]+zw+photoz+'wstsee'+syscut+test+'vdepthgjackerr',xlim=30)
        print b
        cl.append(1.-b)
        ml.append(m)
    zl = [0.58,0.65,0.75,0.85,0.95,1.02] #low and high values needed to allow spline to work over full range
    stl = [cl[0],cl[0],cl[1],cl[2],cl[3],cl[3]]
    stsp = spline(zl,stl)
    mtl = [ml[0],ml[0],ml[1],ml[2],ml[3],ml[3]]
    msp = spline(zl,mtl)
    #test spline
    from matplotlib import pyplot as plt
    zt = numpy.arange(.6,1.,.001)
    stt = []
    for i in range(0,len(zt)):
        #stt.append(Splev(zt[i],stsp))
        stt.append(Splev(stsp,zt[i]))
    plt.plot(zt,stt,'k-',zl,stl,'ko')
    plt.show()

    meang = 23.67
    f = open(outdir+'Y1red'+outf[0]+photoz+syscut+test+'wstsee.dat')
    f.readline()
    fo = open(outdir+'Y1red'+outf[0]+photoz+syscut+test+'wstseegdepth.dat','w')
    fo.write('#ra dec z z_mc wtot wdepth wsee wst id\n')
    for line in f:
        ln = line.split()
        z = float(ln[2])
        zmc = float(ln[3])
        if z >= 0.6 and z <= 1.0:
            stc = Splev(stsp,z)
            slp = stc/meang
            #slp = Splev(msp,z)
            ra,dec = float(ln[0]),float(ln[1])
            th,phi = radec2thphi(ra,dec)
            p = hp.ang2pix(4096,th,phi)
            nst = gmap[p]
            w = 1.
            if nst >= 0:
                w = 1./((1.-stc)+slp*nst)
            fo.write(str(ra)+' '+str(dec)+' '+str(z)+' '+str(zmc)+' '+str(w*float(ln[4]))+' '+ str(w)+' '+ln[5]+' '+ln[6]+' '+str(int(float(ln[-1])))+'\n')
    fo.close()
    return True


def addblindweight(res=512,blindf=.4):
    #adds a weight for stellar contamination as a function of redshift
    #1st step is to get contamination in each z bin
    from random import random
    cl = []
    #odll = []
    odl = mkodensred(0.6,1.,res)
    fo = open(outdir+'blindfac.dat','w')
    for i in range(0,8):
        zl = 0.6+0.05*i
        zh = zl+.05
        #odll.append(mkodensred(zl,zh,res))
        fac = 1.-blindf/2.+blindf#*random()
        fo.write(str(fac)+'\n')
        cl.append(fac)
    zl = [0.59,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975,1.01] #low and high values needed to allow spline to work over full range
    stl = [cl[0],cl[0],cl[1],cl[2],cl[3],cl[4],cl[5],cl[6],cl[7],cl[7]]
    stsp = spline(zl,stl)
    #print len(odll)
    f = open(outdir+'Y1red'+outf[0]+photoz+test+'wst.dat')
    fo = open(outdir+'Y1red'+outf[0]+photoz+test+'wstBLIND.dat','w')
    wbave = 0
    n = 0
    for line in f:
        ln = line.split()
        z = float(ln[2])
        wb = 1.
        if z >= 0.6 and z < 1.:
            zb = int((z-.6)/.05)
            if zb < 0 or zb > 7:
                print z,zb
            stc = Splev(stsp,z)
            ra,dec = float(ln[0]),float(ln[1])
            th,phi = radec2thphi(ra,dec)
            p = hp.ang2pix(512,th,phi)
            #od = odll[zb][p]
            od = odl[p]
            if od == -1:
                print z,zb,p
            #wb = stc-stc/(od+1.)+1./(od+1.)
            #if zb == 0:
            wb = (stc*od+1.)/(od+1.)
            wbave += wb
            n += 1.
            wt = wb*float(ln[-1])
            fo.write(str(ra)+' '+str(dec)+' '+str(z)+' '+str(wt)+'\n')
    fo.close()
    print wbave/n
    return True



def mkgalmapY1red(res,zmin,zmax,wm='',wo=False):
    gl = []
    for i in range(0,12*res*res):
        gl.append(0)
    #if inputfile =='Y1_LSS_AUTO_MOF_BPZv1.BPZ_mof_orig_sva1prior_priormag_mof_i.fits':
    #    f = open(outdir+'Y1_LSS_AUTO_MOF_BPZv1_m22.'+exten)
    #else:
    f = open(outdir+'Y1red'+outf[0]+photoz+syscut+test+wm+'.dat')
    ngt = 0
    w = 1.
    for line in f:
        if line[0] != '#':
            ln = line.split()
            ra,dec,z = float(ln[0]),float(ln[1]),float(ln[2])
            if wm != '':
                w = float(ln[4])
            if z >= zmin and z < zmax:
                th,phi = radec2thphi(ra,dec)
                p = hp.ang2pix(res,th,phi)
                gl[p] += w
                ngt += w
    print len(gl),ngt
    if wo == True:
        fo = open(outdir+'galy1redmap'+str(zmin)+str(zmax)+str(res)+'.dat','w')
        for i in range(0,len(gl)):
            fo.write(str(gl[i])+'\n')
        fo.close()
    return gl

def mkodensred(zmin=0,zmax=2,res=512,wm=''):
    gall = mkgalmapY1red(res,zmin,zmax,wm=wm)
    mask = open(outdir+'Y1'+mf+str(dl)+str(res)+'ring.dat').readlines()
    ml = []
    odl = []
    for i in range(0,12*res*res):
        ml.append(0.0)
        odl.append(0)
        np = 0
    ng = 0
    for i in range(0,len(mask)):
        ln = mask[i].split()
        p = int(ln[0])
        if res == 4096:
            fr = 1.
        else:
            fr = float(ln[1])
            ml[p] = fr
        np += fr
        ng += gall[p]
    print np,ng
    ave = ng/np
    print ave
    for i in range(0,len(ml)):
        if ml[i] > 0:
            odl[i] = gall[i]/ml[i]/ave-1.
    return odl

def testreweight(fac=1.2,zmin=.6,zmax=1.,res=512):
    od = mkodensred(zmin,zmax)
    gall = mkgalmapY1red(res,zmin,zmax)
    wave = 0
    nt = 0
    for i in range(0,len(gall)):
        if gall[i] > 0:
            w = (fac*od[i]+1.)/(od[i]+1.)
            wave += gall[i]*w
            nt += gall[i]
    print wave/nt,nt,wave
    return True


def ngalvnstarsysl(gl,sl,dl,res=512,t=.2,nbin=10,smin=50,smax=600,smd='stars',decmin=-70,decmax=5):
    from healpix import thphi2radec,pix2ang_ring
    mask = open(outdir+'Y1'+mf+str(dl)+syscut+str(res)+'ring.dat').readlines()
    ml = []
    for i in range(0,12*res*res):
        ml.append(0.0)
        np = 0
    for i in range(0,len(mask)):
        ln = mask[i].split()
        p = int(float(ln[0])+.0001)
            #if res == 4096:
            #fr = 1.
            #else:
        fr = float(ln[1])
        ml[p] = fr
        np += 1
    print np
    sysl = []
    for i in range(0,nbin*2):
        sysl.append([])
    min = 1000
    max = -100
    nb = 0
    for i in range(0,len(sl)):
        m = ml[i]
        if m > 0:
            th,phi = pix2ang_ring(res,i)
            dec = -180./pi*th+90.
            if dec < decmin or dec > decmax:
                m = 0
                nb += 1
            if m > t:
                if sl[i] > max and sl[i] < smax:
                    max = sl[i]
                if sl[i] < min and sl[i] > smin:
                    min = sl[i]
    print 'min,max,nb'
    print min,max,nb
    gtl = []
    pl = []
    for i in range(0,nbin):
        gtl.append(0)
        pl.append(0)
    gt = 0
    pt = 0
    slt = 0
    print len(gl)
    for i in range(0,len(gl)):
        m = ml[i]
        if m > 0:
            th,phi = pix2ang_ring(res,i)
            dec = -180./pi*th+90.
            if dec < decmin or dec > decmax:
                m = 0
            if m > t:
                slv = sl[i]
                if slv >= min and slv <= max:
                    pt += ml[i]
                    gt += gl[i]

                    slt += sl[i]
                    bin = int((slv-min)/(max*1.00001-min)*nbin)
                    if bin < nbin:
                        sysl[bin*2].append(gl[i])
                        sysl[bin*2+1].append(ml[i])
                    gtl[bin] += slv*ml[i]
                    pl[bin] += ml[i]
    bcl = []
    for i in range(0,nbin):
        bcl.append(gtl[i]/pl[i])
    print gt,pt,t,slt
    print slt/pt
    return sysl,gt/pt,min,max,bcl

def testseemask():
    f = fitsio.read(inputdir+'Y1A1NEW_COADD_SPT_band_i_FWHM/Y1A1NEW_COADD_SPT_band_i_nside4096_oversamp4_FWHM_MEAN_coaddweights3_mean.fits.gz')
    starl = []
    np = 12*4096**2
    for i in range(0,np):
        starl.append(0)
    for i in range(0,len(f)):
        p = f[i]['PIXEL']
        starl[p] = f[i]['SIGNAL']
    n = 0
    mask = open(outdir+'Y1'+mf+str(dl)+syscut+'4096ring.dat').readlines()
    for i in range(0,len(mask)):
        ln = mask[i].split()
        p = int(float(ln[0])+.0001)
        if starl[p] > 4.0:
            n += 1.
    print n
    return True



def mkseemap(bnd='i'):
    f = fitsio.read(inputdir+'Y1A1NEW_COADD_SPT_band_'+bnd+'/Y1A1NEW_COADD_SPT_band_'+bnd+'_nside4096_oversamp4_FWHM_MEAN_coaddweights3_mean.fits.gz')
    starl = []
    np = 12*4096**2
    for i in range(0,np):
        starl.append(0)
    for i in range(0,len(f)):
        p = f[i]['PIXEL']
        starl[p] = f[i]['SIGNAL']
    f = fitsio.read(inputdir+'Y1A1NEW_COADD_STRIPE82/nside4096_oversamp4/Y1A1NEW_COADD_STRIPE82_band_'+bnd+'_nside4096_oversamp4_FWHM_MEAN_coaddweights3_mean.fits.gz')
    for i in range(0,len(f)):
        p = f[i]['PIXEL']
        starl[p] = f[i]['SIGNAL']

    return starl

def putngalvnstar3jackN(zmin=0,zmax=2,smin=0.0,smax=1.8,njack=20,res=512,t=.8,nbin=10,dl=22,stellcon=0,pz=photoz,smd='stars',bnd='i',decmin=-70,decmax=5,wm='',mapin=False):
    gall = mkgalmapY1red(res,zmin,zmax,wm=wm)
    if smd == 'stars':
        starl = hp.read_map(inputdir+'y1a1_gold_1.0.2_stars_nside0512.fits')
    try:
        print mapin[0]
        starl = mapin
    except:
    #if mapin == False:
        if smd == 'depth':
            starl = hp.read_map(depthmap)
        if smd == 'seei':
            f = fitsio.read(inputdir+'Y1A1NEW_COADD_SPT_band_i_FWHM/Y1A1NEW_COADD_SPT_band_i_nside4096_oversamp4_FWHM_MEAN_coaddweights3_mean.fits.gz')
            starl = []
            np = 12*4096**2
            for i in range(0,np):
                starl.append(0)
            for i in range(0,len(f)):
                p = f[i]['PIXEL']
                starl[p] = f[i]['SIGNAL']
#   else:
#    starl = mapin
    sysl,mt,smin,smax,bcl = ngalvnstarsysl(gall,starl,dl,res=res,t=t,nbin=nbin,smin=smin,smax=smax,smd=smd,decmin=decmin,decmax=decmax)
    print mt
    sl = []
    ml = []
    print len(sysl)
    for i in range(0,2*nbin,2):
        print len(sysl[i])
    for i in range(0,2*nbin,2):
        std = 0
        ng = 0
        np = 0
        print len(sysl[i]),len(sysl[i+1])
        for j in range(0,len(sysl[i])):
            ng += sysl[i][j]
            np += sysl[i+1][j]
        mean = ng/np/mt
        print i,mean
        ml.append(mean)
        if len(sysl[i]) < njack:
            for k in range(0,len(sysl[i])):
                std += (sysl[i][k]/sysl[i+1][k]/mt-mean)**2.
                std = sqrt(std/(len(sysl[i])-1.))/sqrt(len(sysl[i])-1.)
        else:
            jkf = len(sysl[i])/njack
            for k in range(0,njack):
                ng = 0
                np = 0
                min = jkf*k
                max = jkf*(k+1)
                for j in range(0,len(sysl[i])):
                    if j < min or j >= max:
                        ng += sysl[i][j]
                        np += sysl[i+1][j]
                mj = ng/np/mt
                std += (mj-mean)**2.
            std = sqrt((njack-1.)/float(njack)*std)
        print i,mean,std
        sl.append(std)
    rw = ''
    rw += pz
    rw += wm
    fo = open(outdir+'Y1red'+outf[0]+str(zmin)+str(zmax)+rw+syscut+test+'v'+smd+'jackerr.dat','w')
    for i in range(0,nbin):
        slb = sl[i]
        st = i*(smax-smin)/float(nbin)+(smax-smin)/nbin/2.+smin
        #fo.write(str(st)+' '+str(ml[i])+' '+str(slb)+'\n')
        fo.write(str(bcl[i])+' '+str(ml[i])+' '+str(slb)+'\n')
    fo.close()
    return True

def ngalvnstarsimp(zmin=0,zmax=2,smin=0,smax=2.,njack=20,res=512,t=.8,nbin=10,dl=22,stellcon=0,pz=photoz,smd='stars',bnd='i',decmin=-70,decmax=5,wm=''):
    gall = mkgalmapY1red(res,zmin,zmax,wm=wm)
    if smd == 'stars':
        starl = hp.read_map(inputdir+'y1a1_gold_1.0.2_stars_nside0512.fits')
    
    if smd == 'depth':
        starl = hp.read_map(depthmap)
    if smd == 'seei':
        f = fitsio.read(inputdir+'Y1A1NEW_COADD_SPT_band_i_FWHM/Y1A1NEW_COADD_SPT_band_i_nside4096_oversamp4_FWHM_MEAN_coaddweights3_mean.fits.gz')
        starl = []
        np = 12*4096**2
        for i in range(0,np):
            starl.append(0)
        for i in range(0,len(f)):
            p = f[i]['PIXEL']
            starl[p] = f[i]['SIGNAL']

    mask = open(outdir+'Y1'+mf+str(dl)+str(res)+'ring.dat').readlines()
    ml = []
    for i in range(0,12*res*res):
        ml.append(0.0)
        np = 0
    for i in range(0,len(mask)):
        ln = mask[i].split()
        p = int(ln[0])
        if res == 4096:
            fr = 1.
        else:
            fr = float(ln[1])
        ml[p] = fr
        np += 1
    print np
    gtl = []
    pixl = []
    for i in range(0,nbin):
        gtl.append(0)
        pixl.append(0)
    min = 1000
    max = -100
    nb = 0
    for i in range(0,len(starl)):
        m = ml[i]
        if m > t:
            if starl[i] > max and starl[i] < smax:
                max = starl[i]
            if starl[i] < min and starl[i] > smin:
                min = starl[i]
    print 'min,max,nb'
    print min,max,nb
    gt = 0
    pt = 0
    slt = 0
    print len(gall)
    for i in range(0,len(gall)):
        m = ml[i]
        if m > t:
            slv = starl[i]
            if slv >= min and slv <= max:
                pt += ml[i]
                gt += gall[i]

                slt += starl[i]
                bin = int((slv-min)/(max*1.00001-min)*nbin)
                if bin < nbin:
                    gtl[bin] += gall[i]
                    pixl[bin] += ml[i]
    rw = ''
    rw += pz
    print gt,pt
    print gtl,pixl
    fo = open(outdir+'Y1red'+outf[0]+str(zmin)+str(zmax)+rw+test+'v'+smd+'simp.dat','w')
    for i in range(0,nbin):
        
        st = i*(max-min)/float(nbin)+(max-min)/nbin/2.+min
        fo.write(str(st)+' '+str(gtl[i]/pixl[i]*pt/gt)+'\n')
    fo.close()
    return True

def ngalvnstardivsee(zmin=0,zmax=2,smin=0,smax=2.,split=3.5,njack=20,res=512,t=.8,nbin=10,dl=22,stellcon=0,pz=photoz,smd='stars',bnd='i',decmin=-70,decmax=5,wm=''):
    gall = mkgalmapY1red(res,zmin,zmax,wm=wm)
    if smd == 'stars':
        starl = hp.read_map(inputdir+'y1a1_gold_1.0.2_stars_nside0512.fits')
    seel = []
    f = fitsio.read(inputdir+'Y1A1NEW_COADD_SPT_band_i_FWHM/Y1A1NEW_COADD_SPT_band_i_nside4096_oversamp4_FWHM_MEAN_coaddweights3_mean.fits.gz')
    np = 12*4096**2
    for i in range(0,np):
        seel.append(0)
    for i in range(0,len(f)):
        p = f[i]['PIXEL']
        seel[p] = f[i]['SIGNAL']

    mask = open(outdir+'Y1'+mf+str(dl)+str(res)+'ring.dat').readlines()
    ml = []
    for i in range(0,12*res*res):
        ml.append(0.0)
        np = 0
    for i in range(0,len(mask)):
        ln = mask[i].split()
        p = int(ln[0])
        if res == 4096:
            fr = 1.
        else:
            fr = float(ln[1])
            ml[p] = fr
            np += 1
    print np
    gtl = []
    pixl = []
    gth = []
    pixh = []
    for i in range(0,nbin):
        gtl.append(0)
        pixl.append(0)
        gth.append(0)
        pixh.append(0)
    min = 1000
    max = -100
    nb = 0
    for i in range(0,len(starl)):
        m = ml[i]
        if m > t:
            if starl[i] > max and starl[i] < smax:
                max = starl[i]
            if starl[i] < min and starl[i] > smin:
                min = starl[i]
    print 'min,max,nb'
    print min,max,nb
    gt = 0
    pt = 0
    slt = 0
    print len(gall)
    for i in range(0,len(gall)):
        m = ml[i]
        if m > t:
            slv = starl[i]
            if slv >= min and slv <= max:
                pt += ml[i]
                gt += gall[i]
                
                slt += starl[i]
                bin = int((slv-min)/(max*1.00001-min)*nbin)
                if bin < nbin:
                    th,phi = hp.pix2ang(res,i)
                    p4 = hp.ang2pix(4096,th,phi)
                    seev = seel[p4]
                    if seev < split:
                        gtl[bin] += gall[i]
                        pixl[bin] += ml[i]
                    else:
                        gth[bin] += gall[i]
                        pixh[bin] += ml[i]

    rw = ''
    rw += pz
    print gt,pt
    print gtl,pixl
    fo = open(outdir+'Y1red'+outf[0]+str(zmin)+str(zmax)+rw+test+'v'+smd+'splseei'+str(split)+'.dat','w')
    for i in range(0,nbin):
        st = i*(max-min)/float(nbin)+(max-min)/nbin/2.+min
        fo.write(str(st)+' '+str(gtl[i]/pixl[i]*pt/gt)+' '+str(pixl[i])+' '+str(gth[i]/pixh[i]*pt/gt)+' '+str(pixh[i])+'\n')
    fo.close()
    return True


def ngalvnstarhist(zmin=0,zmax=2,njack=20,res=512,t=.8,nbin=10,dl=22,stellcon=0,pz=photoz,smd='stars',bnd='i',decmin=-70,decmax=5,wm=''):
    gall = mkgalmapY1red(res,zmin,zmax,wm=wm)
    if smd == 'stars':
        starl = hp.read_map(inputdir+'y1a1_gold_1.0.2_stars_nside0512.fits')
    
    if smd == 'depth':
        starl = hp.read_map(depthmap)
    if smd == 'seei':
        f = fitsio.read(inputdir+'Y1A1NEW_COADD_SPT_band_i_FWHM/Y1A1NEW_COADD_SPT_band_i_nside4096_oversamp4_FWHM_MEAN_coaddweights3_mean.fits.gz')
        starl = []
        np = 12*4096**2
        for i in range(0,np):
            starl.append(0)
        for i in range(0,len(f)):
            p = f[i]['PIXEL']
            starl[p] = f[i]['SIGNAL']

    mask = open(outdir+'Y1'+mf+str(dl)+str(res)+'ring.dat').readlines()
    ml = []
    for i in range(0,12*res*res):
        ml.append(0.0)
        np = 0
    for i in range(0,len(mask)):
        ln = mask[i].split()
        p = int(ln[0])
        if res == 4096:
            fr = 1.
        else:
            fr = float(ln[1])
            ml[p] = fr
            np += 1
    print np
    starlh = []
    mlh = []
    gallh = []
    for i in range(0,len(ml)):
        if ml[i] > t:
            if starl[i] > 0:
                starlh.append(starl[i])
                mlh.append(ml[i])
                gallh.append(gall[i])
            #else:
            #    print i,ml[i],starl[i]
    histg = numpy.histogram(starlh,weights=gallh)
    histm = numpy.histogram(starlh,weights=mlh)
    print histg[0],histm[0]
    rw = ''
    rw += pz
    pt = sum(mlh)
    gt = sum(gallh)
    print gt,pt,pt/gt
    fo = open(outdir+'Y1red'+outf[0]+str(zmin)+str(zmax)+rw+test+'v'+smd+'hist.dat','w')
    smin = min(starlh)
    smax = max(starlh)
    print smin,smax
    for i in range(0,nbin):
        
        st = i*(smax-smin)/float(nbin)+(smax-smin)/nbin/2.+smin
        fo.write(str(st)+' '+str(histg[0][i]/histm[0][i]*pt/gt)+'\n')
    fo.close()
    return True



def findcol(name):
    for i in range(0,len(header)):
        if header[i].strip('\n').strip('\r') == name:
            return i
    print name+' not found'
    print name,header[-1].strip('\r').strip('\n')
    a = header[-1].strip('\n').strip('\r')
    for i in range(0,len(a)):
        print a[i]
    print len(header[-1].strip('\r').strip('\n'))
    if name == str(header[-1].strip('\r').strip('\n')):
        print 'yup'
    return 'Error, column not found!!!'

def Splev(spl,x):
    return splint(spl[1],spl[2],spl[0],x)

def dolinbestfit(file,xlim=2.2,xmin=0):
    from optimize import fmin
    f = load(outdir+file+'.dat').transpose()
    lf = linfit(f[0],f[1],f[2],xlim=xlim,xmin=xmin)
    inl = numpy.array([1.,0])
    b0,m0 = fmin(lf.chilin,inl)
    return b0,m0

class linfit:
    def __init__(self,xl,yl,el,xlim=1000.,xmin=-1000):
        self.xl = xl
        self.yl = yl
        self.el = el
        self.xlim = xlim
        self.xmin = xmin
    def chilin(self,bml):
        chi = 0
        b = bml[0]
        m = bml[1]
        for i in range(0,len(self.xl)):
            if self.xl[i] < self.xlim and self.xl[i] > self.xmin:
                y = b+m*self.xl[i]
                chi += (self.yl[i]-y)**2./self.el[i]**2.
        return chi

def read_in_file(fname, icol, index ):
    data  = plb.loadtxt( fname, usecols=(icol,), unpack=True, skiprows=1)
    
    if len(index) < 1:
        return  data
    else:
        return data[index]

def masksyscut(file,syscut):
    maskf = open(outdir+'Y1goldFoot_fracdet'+str(fracd)+'wide_iautodepth'+str(dl)+syscut+'4096ring.dat')
    np = 12*4096*4096
    mask = []
    for i in range(0,np):
        mask.append(0)
    for line in maskf:
        pix = int(float(line))
        mask[pix] = 1
    f = open(outdir+file+'.dat')
    fo = open(outdir+file+syscut+'.dat','w')
    for line in f:
        ln = line.split()
        ra,dec = float(ln[0]),float(ln[1])
        th,phi = radec2thphi(ra,dec)
        p = hp.ang2pix(4096,th,phi)
        if mask[p] == 1:
            fo.write(line)
    fo.close()
    return True



if __name__ == '__main__':
#do everything
    #maskY1_tot()
    #maskY1see()
    #maskd(512)
    #m = mksample_merge() #for new, merged file, initialize class, gets column names
    #m.createmaskedY1redsimp()
#    m = mksample() #for old file, initialize class, gets column names
#    m.createmaskedY1simp() #only masks sample
#    m.mkRedsampradeczsimp() #only write out relevant columns from simple masking
#    mksampfitsBPZ() #was used with fits file, should be depreciated


#stellar density 1D
#    putngalvnstar3jackN(.6,.65)
#    putngalvnstar3jackN(.65,.7)
#    putngalvnstar3jackN(.7,.75)
#    putngalvnstar3jackN(.75,.8)
#    putngalvnstar3jackN(.8,.85)
#    putngalvnstar3jackN(.85,.9)
#    putngalvnstar3jackN(.9,.95)
#    putngalvnstar3jackN(.95,1.)
#    writestfit()
#    addstweight()
#add seeing weight, fit over full range
#    wst = 'wst'
#    sys = 'seei'
#    smin = 2.5
#    smax = 4.7
#    seemap = mkseemap()
#    putngalvnstar3jackN(zmin=0.6,zmax=1.,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    addseeweightone(xlim=smax)
#testing g-band depth weights
    wst = 'wstsee'
    sys = 'depthg'
    smin = 21
    smax = 25
    gmap = hp.read_map(inputdir+'y1a1_gold_1.0.2_wide_auto_nside4096_g_10sigma.fits.gz')
#    putngalvnstar3jackN(zmin=0.6,zmax=.7,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='g',decmin=-70,decmax=5,mapin=gmap)
#    putngalvnstar3jackN(zmin=0.7,zmax=.8,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='g',decmin=-70,decmax=5,mapin=gmap)
#    putngalvnstar3jackN(zmin=0.8,zmax=.9,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='g',decmin=-70,decmax=5,mapin=gmap)
#    putngalvnstar3jackN(zmin=0.9,zmax=1.,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='g',decmin=-70,decmax=5,mapin=gmap)
    addgdepthweight(gmap)

#below are tests that are sometimes nice to do
#test seeing weights
#    wst = 'wstsee'
#    sys = 'seei'
#    smin = 2.5
#    smax = 4.7
#    seemap = mkseemap()
#    putngalvnstar3jackN(zmin=0.6,zmax=1.,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='i',decmin=-70,decmax=5,mapin=seemap)

#test seeing
#    sys = 'seei'
#    smin = 1.5
#    smax = 6.0
#    seemap = mkseemap()
#    putngalvnstar3jackN(zmin=0.6,zmax=1.,smin=smin,smax=smax,res=4096,t=.2,wm='',pz=photoz,smd=sys,bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    sys = 'seez'
#    smin = 1.5
#    smax = 6.0
#    seemap = mkseemap(bnd='z')
#    putngalvnstar3jackN(zmin=0.6,zmax=1.,smin=smin,smax=smax,res=4096,t=.2,wm='',pz=photoz,smd=sys,bnd='z',decmin=-70,decmax=5,mapin=seemap)


#test after weighting
#    putngalvnstar3jackN(.75,.8,wm='wst')
#    putngalvnstar3jackN(.75,.8,wm='wstsee')
#stellar density 1D, histogram method
#    ngalvnstarhist(.6,.65)
#    ngalvnstarhist(.65,.7)
#    ngalvnstarhist(.7,.75)
#    ngalvnstarhist(.75,.8)
#    ngalvnstarhist(.8,.85)
#    ngalvnstarhist(.85,.9)
#    ngalvnstarhist(.9,.95)
#    ngalvnstarhist(.95,1.)

# simple stellar density 1D
#    ngalvnstarsimp(.6,.65)
#    ngalvnstarsimp(.65,.7)
#    ngalvnstarsimp(.7,.75)
#    ngalvnstarsimp(.75,.8)
#    ngalvnstarsimp(.8,.85)
#    ngalvnstarsimp(.85,.9)
#    ngalvnstarsimp(.9,.95)
#    ngalvnstarsimp(.95,1.)

# simple stellar density 1D, split by seeing
#    ngalvnstardivsee(.6,.65)
#    ngalvnstardivsee(.65,.7)
#    ngalvnstardivsee(.7,.75)
#    ngalvnstardivsee(.75,.8)
#    ngalvnstardivsee(.8,.85)
#    ngalvnstardivsee(.85,.9)
#    ngalvnstardivsee(.9,.95)
#    ngalvnstardivsee(.95,1.)


#add blinding weight
#    addblindweight()
#testreweight()
#plot stellar density fits
#    plotstfits()

#find nz from MC redshifts, print dispersions based on MC
#    mknz()
#test depth
#    wst = 'wst'
#    putngalvnstar3jackN(zmin=0.6,zmax=0.65,smin=22,smax=23.,res=4096,t=.2,wm=wst,pz=photoz,smd='depth',bnd='i',decmin=-70,decmax=5)
#    putngalvnstar3jackN(zmin=0.65,zmax=0.7,smin=22,smax=23.,res=4096,t=.2,wm=wst,pz=photoz,smd='depth',bnd='i',decmin=-70,decmax=5)
#    putngalvnstar3jackN(zmin=0.7,zmax=0.75,smin=22,smax=23.,res=4096,t=.2,wm=wst,pz=photoz,smd='depth',bnd='i',decmin=-70,decmax=5)
#    putngalvnstar3jackN(zmin=0.75,zmax=0.8,smin=22,smax=23.,res=4096,t=.2,wm=wst,pz=photoz,smd='depth',bnd='i',decmin=-70,decmax=5)
#    putngalvnstar3jackN(zmin=0.8,zmax=0.85,smin=22,smax=23.,res=4096,t=.2,wm=wst,pz=photoz,smd='depth',bnd='i',decmin=-70,decmax=5)
#    putngalvnstar3jackN(zmin=0.85,zmax=0.9,smin=22,smax=23.,res=4096,t=.2,wm=wst,pz=photoz,smd='depth',bnd='i',decmin=-70,decmax=5)
#    putngalvnstar3jackN(zmin=0.9,zmax=0.95,smin=22,smax=23.,res=4096,t=.2,wm=wst,pz=photoz,smd='depth',bnd='i',decmin=-70,decmax=5)
#    putngalvnstar3jackN(zmin=0.95,zmax=1.,smin=22,smax=23.,res=4096,t=.2,wm=wst,pz=photoz,smd='depth',bnd='i',decmin=-70,decmax=5)
#plot depth
#    plotvsdepth()
#test any potential systematics
#   residual systematics
#    putngalvnstar3jackN(.6,.7,wm='wstsee')
#    putngalvnstar3jackN(.7,.8,wm='wstsee')
#    putngalvnstar3jackN(.8,.9,wm='wstsee')
#    putngalvnstar3jackN(.9,1.,wm='wstsee')
#    testseemask()
#    seemap = mkseemap()
#    putngalvnstar3jackN(zmin=0.6,zmax=.7,smin=2.5,smax=4.5,res=4096,t=.2,wm='wstsee',pz=photoz,smd='seei',bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    putngalvnstar3jackN(zmin=0.7,zmax=.8,smin=2.5,smax=4.,res=4096,t=.2,wm='wstsee',pz=photoz,smd='seei',bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    putngalvnstar3jackN(zmin=0.8,zmax=.9,smin=2.5,smax=4.,res=4096,t=.2,wm='wstsee',pz=photoz,smd='seei',bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    putngalvnstar3jackN(zmin=0.9,zmax=1.,smin=2.5,smax=4.,res=4096,t=.2,wm='wstsee',pz=photoz,smd='seei',bnd='i',decmin=-70,decmax=5,mapin=seemap)

#    putngalvnstar3jackN(zmin=0.65,zmax=0.7,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    putngalvnstar3jackN(zmin=0.7,zmax=0.75,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    putngalvnstar3jackN(zmin=0.75,zmax=0.8,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    putngalvnstar3jackN(zmin=0.8,zmax=0.85,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    putngalvnstar3jackN(zmin=0.85,zmax=0.9,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    putngalvnstar3jackN(zmin=0.9,zmax=0.95,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    putngalvnstar3jackN(zmin=0.95,zmax=1.,smin=smin,smax=smax,res=4096,t=.2,wm=wst,pz=photoz,smd=sys,bnd='i',decmin=-70,decmax=5,mapin=seemap)
#    writeseefit()
#    plotseefits()
#    addseeweight()
#plot sys
#    plotvssys(sys)
#test different red selection
#    cc = 1.7 #fiducial is 0.9
#    acut = 2. #fiducial is 0.6
#    dlr = 22. #fiducial is 22
#    zmmax = 30. #fiducial is 30
#    rmimin = -1 #fiducial is -1
#    rmimax = 2.5 #fiducial is 2.5
#    rmzmax = 10 #fiducial is 10
#    tpmin = -1 #fidcucial is -1
#    tpmax = 10. #fiducial is 10
#    gmin = 0 #fiducial is 0
#    drfmax = 100 #fiducial is 100
#    zsld = 19. #fiducial is 100
#    islp = 3. #fiducial is 0
#    m = mksample()
#    m.mkRedsampradecztest(ccut=cc,dlr=dlr,rmimin=rmimin,rmimax=rmimax,tpmin=tpmin,tpmax=tpmax,gmin=gmin,rmzmax=rmzmax,zmmax=zmmax,drfmax=drfmax,zsld=zsld,acut=acut,islp=islp)
#    mknz()
