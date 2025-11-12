#!/usr/bin/env python
#---------------------------------------------------------------------
#  plot_des_results_fast.py
#  
#  Copyright 2013 Matias Carrasco Kind <mcarras2@illinois.edu>
#  

from pylab import *
from numpy import *
import sys,getopt
import matplotlib.cm as cm
from scipy import optimize
from utils import *
import tawala as tw
a=tw.Tawala()
a.setVersion('scipy', '0.11.0', '0')
a.setVersion('matplotlib', '1.0.0', '1')




def isNaN(num): return num != num


###################################
#########    MAIN     #############
###################################

def usage():
    print
    print "Usage:: "
    print "-----------------------------------------------------------------------"
    print "plot_des_results_fast --t1 <test_file>  --err1 <err_cut> "
    print "-----------------------------------------------------------------------"
    print "* error cuts are optional"
    print "* order doesn't matter"
    print "------------------------------"
    print "-> MORE OPTIONAL PARAMETERS <-"
    print "------------------------------"
    print "* --output <output_root> "
    print "* the output root creates a file with summary of results for each test"
    print "* --name1 <name test 1>"
    print "* --w1 <weights file test 1>"
    print "* --frac1 <fraction between 0 and 1 of best selected galaxies (smallest errors) test 1>"
    print "* --only_all yes Compute only whole range metrics (for .pz files)"
    print
    print "* See README file for more information"
    sys.exit(0)

if len(sys.argv) < 2 or sys.argv[1]=='-h' or sys.argv[1]=='--help':
    usage()

t1file = ''
v1file = ''
w1file=''
err1=1.e10
outfile1=''
outroot=''
name1_input=''
frac1=-1.

oall='no'

try:
    opts, args = getopt.getopt(sys.argv[1:],"h",["t1=","err1=","output=","name1=","w1=","frac1=","only_all="])
except getopt.GetoptError:
    usage()
for opt, arg in opts:
    if opt in ("--t1"):
        t1file= arg
    elif opt in ("--err1"):
        err1=float(arg)
    elif opt in ("--output"):
        outroot=arg
    elif opt in ("--name1"):
        name1_input=arg
    elif opt in ("--w1"):
        w1file=arg
    elif opt in ("--frac1"):
        frac1=float(arg)
    elif opt in ("--only_all"):
        oall='yes'
if outroot!='':
    outfile1=outroot+'_test1.stats'


if t1file[-3:]=='.pz':
    zmid_t1,zs_t1,fid_t1,number_t1,BIG_t1=read_pdf_file(t1file)
    zp_t1,ep_t1=get_vals_pdf(zmid_t1,BIG_t1)
else:
    zp_t1,ep_t1,zs_t1,fid_t1,number_t1=loadtxt(t1file,unpack=True)
    
if w1file!='': 
    fid_w1,number_w1,vals_w1=loadtxt(w1file,unpack=True)
    NW1=len(fid_w1)
    mask=check_fid(fid_w1,number_w1,fid_t1,number_t1)
    fid_w1=fid_w1[mask]
    number_w1=number_w1[mask]
    vals_w1=vals_w1[mask]

if frac1>0:
    Dz=ep_t1
    sDz=argsort(Dz)
    keep_n=int(frac1*len(Dz))-1
    mask=sDz[0:keep_n]
    zp_t1=zp_t1[mask]
    ep_t1=ep_t1[mask]
    zs_t1=zs_t1[mask]
    fid_t1=fid_t1[mask]
    number_t1=number_t1[mask]
    if t1file[-3:]=='.pz':BIG_t1=BIG_t1[mask]
    err1=1.e10
    if w1file!='':
        fid_w1=fid_w1[mask]
        number_w1=number_w1[mask]
        vals_w1=vals_w1[mask]

werror=where(ep_t1 < err1)[0]
ep_t1=ep_t1[werror]
zp_t1=zp_t1[werror]
zs_t1=zs_t1[werror]
fid_t1=fid_t1[werror]
number_t1=number_t1[werror]


if w1file!='': 
    fid_w1=fid_w1[werror]
    number_w1=number_w1[werror]
    vals_w1=vals_w1[werror]
    NTOT1=1.*NW1
else:
    vals_w1=ones(len(zp_t1))
    vals_w1/=sum(vals_w1)
    NTOT1=1.*len(vals_w1)
    

minz=0
maxz=1.5

zx1=zs_t1
zy1=zp_t1
e1=ep_t1
fid1=fid_t1
nu1=number_t1
if t1file[-3:]=='.pz':
    name1='Test 1 (pdf)'
    zM1=zmid_t1
    B1=BIG_t1
else: name1='Test 1'
if name1_input!='' : name1=name1_input


err1_all=(zy1-zx1)/e1
wzero=where((e1>0.))[0]
bias_norm1=sum(err1_all[wzero]*vals_w1[wzero])/sum(vals_w1[wzero])
sigma_norm1=sqrt(sum(vals_w1[wzero]*(err1_all[wzero]-bias_norm1)**2)/sum(vals_w1[wzero]))
if t1file[-3:]=='.pz':
    D_KS1=KS_pdf(zmid_t1,BIG_t1,zx1,vals_w1)
else:
    D_KS1=KS(zy1,zx1,vals_w1)
    
dz_bin=0.1

if t1file[-3:]=='.pz':
    B1p=metrics_pdf_w(zx1,zM1,B1,vals_w1,name1,0,1.5,dz_bin,NTOT1,mode=1,oall=oall)
else:
    B1p=metrics_w(zx1,zy1,vals_w1,name1,0,1.5,dz_bin,NTOT1,mode=1)
    

print
print '----------------------'
print '       TEST 1         '
print '----------------------'
print 'N_total         : %d  ' % B1p.ntotal2
print 'Bias            : %.5f' % B1p.mean_all
print 'Median          : %.5f' % B1p.median_all
print 'Sigma           : %.5f' % B1p.sigma_all
print 'Sigma68         : %.5f' % B1p.sigma68_all
print 'Frac2           : %.5f' % B1p.frac2_all
print 'Frac3           : %.5f' % B1p.frac3_all
print 'Bias_Norm       : %.5f' % bias_norm1
print 'Sigma_Norm      : %.5f' % sigma_norm1
print 'D_N_poisson     : %.5f' % B1p.n_poisson_rms
print 'D_KS_test       : %.5f' % D_KS1
print 

if outroot!='':
    F1=open(outfile1,'w')
    F1.write('#bin  zmid    Nin   bias    median   sigma  sigma68  frac2   frac3  bias_norm    sigma_norm  n_poisson   D_KS_Nz\n')
    F1.write('#Whole sample\n')
    line='%2d  %.4f  %5d  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n' % (\
    0,0.75,B1p.ntotal2,B1p.mean_all,B1p.median_all,B1p.sigma_all,\
    B1p.sigma68_all,B1p.frac2_all,B1p.frac3_all,bias_norm1,sigma_norm1,B1p.n_poisson_rms,D_KS1)
    F1.write(line)
    F1.write('#Binned sample\n')
    for j in xrange(len(B1p.bins)):
        line='%2d  %.4f  %.2f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n' % (\
            j+1,B1p.bins[j],B1p.n[j],B1p.mean[j],B1p.median[j],B1p.sigma[j],\
            B1p.sigma68[j],B1p.frac2[j],B1p.frac3[j],0.,0.,B1p.n_poisson[j],0.)
        F1.write(line)
    F1.close()

