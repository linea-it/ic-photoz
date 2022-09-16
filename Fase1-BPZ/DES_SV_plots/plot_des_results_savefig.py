#!/usr/bin/env python
#---------------------------------------------------------------------
#  plot_des_results.py
#  
#  Copyright 2013 Matias Carrasco Kind <mcarras2@illinois.edu>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import tawala as tw
a=tw.Tawala()
a.setVersion('scipy', '0.11.0', '0')
a.setVersion('matplotlib', '1.0.0', '1')



from pylab import *
from numpy import *
import sys,getopt
import matplotlib.cm as cm
from scipy import optimize
from utils import *
try:
    from Tkinter import *
    import tkFont
    notk=False
except ImportError:
    notk=True
    

class Parameter:
    def __init__(self, value): self.value = value
    def set(self, value): self.value = value 
    def __call__(self): return self.value

def fit(function, parameters, x, y):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(x)
    if x is None: x = arange(y.shape[0])
    p = [param() for param in parameters]
    optimize.leastsq(f, p)
    
def isNaN(num): return num != num

def onpick2(event):
    if event.mouseevent.inaxes:
        s3B.cla()
        s3B.set_xlim(x1,x2)
        x=event.mouseevent.xdata
        y=event.mouseevent.ydata
        cur_ax=event.mouseevent.inaxes
        if cur_ax.is_first_col() :inx=1
        else : inx=2
        if inx==1: zx=zx1;zy=zy1
        if inx==2: zx=zx2;zy=zy2
        dx=array(x-zx[event.ind],dtype=float)
        dy=array(y-zy[event.ind],dtype=float)
        dd=hypot(dx,dy)
        indmin=dd.argmin()
        dataind=event.ind[indmin]
        print dataind
        current1.set_visible(True)
        current1.set_data(zx1[dataind],zy1[dataind])
        current2.set_visible(True)
        current2.set_data(zx2[dataind],zy2[dataind])
        if t1file[-3:]=='.pz':
            t1_pdf=BIG_t1[dataind]
            s3B.plot(zmid_t1,t1_pdf,'b-',label=name1)
            maxy1=max(t1_pdf)*1.05
        else:
            s3B.plot([zy1[dataind],zy1[dataind]],[0,1.],'b-',label=name1)
            maxy1=1
        if t2file[-3:]=='.pz':
            t2_pdf=BIG_t2[dataind]
            s3B.plot(zmid_t2,t2_pdf,'g-',label=name2)
            maxy2=max(t2_pdf)*1.05
        else:
            s3B.plot([zy2[dataind],zy2[dataind]],[0,1],'g-',label=name2)
            maxy2=1
        s3B.plot([zx1[dataind],zx1[dataind]],[0,1],'r--',label='z spec')
        maxy3=min((maxy1,maxy2))
        s3B.set_ylim(0,maxy3)
        if inx==1:
            s3B.plot([0,0],[0.,0],label='field_id ='+str(int(fid_t1[dataind])),visible=False)
            s3B.plot([0,0],[0.,0],label='number ='+str(int(number_t1[dataind])),visible=False)
            sw='weight = %.6f' % vals_w1[dataind]
            s3B.plot([0,0],[0.,0],label=sw,visible=False)
        if inx==2:
            s3B.plot([0,0],[0.,0],label='field_id ='+str(int(fid_t2[dataind])),visible=False)
            s3B.plot([0,0],[0.,0],label='number ='+str(int(number_t2[dataind])),visible=False)
            sw='weight = %.6f' % vals_w2[dataind]
            s3B.plot([0,0],[0.,0],label=sw,visible=False)
        s3B.set_xlim(minz,maxz)
        s3B.set_xlabel('Redshift', fontsize=15)
        s3B.set_title(str(dataind))
        s3B.legend(loc=0)
        cf=gcf()
        cf.canvas.draw()


def on_key(event):
    if event.key=='Q':
        close('all')
        root.destroy()
    if event.key=='q':
        cf=gcf()
        close(cf)
        
def on_key3(event):
    if event.key=='h':
        s3B.cla()
        s3B.plot(xx1,pmf1,'ko',label=name1)
        fit1=r'$\mu=$%5.3f, $\sigma=$%5.3f' % (mu1(),sigma1()) 
        fit2=r'$\mu=$%5.3f, $\sigma=$%5.3f' % (mu2(),sigma2()) 
        s3B.plot(xg,gauss1(xg),'b-',label=fit1,lw=1.5)#/max(y)*max(G[0]))
        s3B.set_xlabel('$\Delta z/\sigma_{\Delta z}$',fontsize=15)
        s3B.set_ylabel('Number density',fontsize=15)
        if t2file !='':
            s3B.plot(xx2,pmf2,'rs',label=name2)
            s3B.plot(xg,gauss2(xg),'g-',label=fit2,lw=1.5)
        s3B.legend(loc=2,prop={'size':12})
        s3B.set_xlim(-6,6)
        cf=gcf()
        cf.canvas.draw()


def on_key2(event):
    global count_o,mode_p
    if event.key=='p':
        mode_p='zp'
        cspec1a.set_visible(False)
        cspec2a.set_visible(False) 
        cphot1a.set_visible(True)  
        cphot2a.set_visible(True)  

        s1M.set_xticklabels('')
        s1M.legend((cphot1a,cphot2a),(name1+ext1p,name2+ext2p),loc=0)
        cspec1b.set_visible(False)
        if t2file!='':cspec2b.set_visible(False) 
        cphot1b.set_visible(True)  
        if t2file!='':cphot2b.set_visible(True)

        s2M.set_xticklabels('')
        
        cspec1c.set_visible(False) 
        cspec2c.set_visible(False) 

        cphot1c.set_visible(True)  
        cphot2c.set_visible(True)

        s3M.set_xticklabels('')
        s3M.legend((cphot1c,cphot2c),(ext1p_c,ext2p_c),loc=0)
        cspec1d.set_visible(False) 
        cspec2d.set_visible(False) 
        cphot1d.set_visible(True)  
        cphot2d.set_visible(True)
        s2M.legend((cphot1d,cphot2d),(ext1p_d,ext2p_d),loc=0)

        cspec1e.set_visible(False) 
        cspec2e.set_visible(False) 
        cphot1e.set_visible(True)  
        cphot2e.set_visible(True)
        cspec1f.set_visible(False) 
        cspec2f.set_visible(False) 
        cphot1f.set_visible(True)  
        cphot2f.set_visible(True)
        
        cspec1g.set_visible(False) 
        if t2file!='':cspec2g.set_visible(False) 
        cphot1g.set_visible(True)  
        if t2file!='':cphot2g.set_visible(True)
        cspec1h.set_visible(False) 
        if t2file!='':cspec2h.set_visible(False) 
        cphot1h.set_visible(True)  
        if t2file!='':cphot2h.set_visible(True)
        
        s8M.set_xticklabels('')
        s6M.set_xlabel(r'$z_{\rm phot}$',fontsize=20)
        s9M.set_xlabel(r'$z_{\rm phot}$',fontsize=20)
        s6M.legend((cphot1e,cphot2e),(ext1p_e,ext2p_e),loc=0)
        s5M.legend((cphot1f,cphot2f),(ext1p_f,ext2p_f),loc=0)
        if t2file !='':
            s7M.legend((cphot1g,cphot2g),(ext1p_g,ext2p_g),loc=0)
            s9M.legend((cphot1b,cphot2b),(ext1p_b,ext2p_b),loc=0)
            s8M.legend((cphot1h,cphot2h),(ext1p_h,ext2p_h),loc=0)
        else:
            s7M.legend((cphot1g,),(ext1p_g,),loc=0)
            s9M.legend((cphot1b,),(ext1p_b,),loc=0)
            s8M.legend((cphot1h,),(ext1p_h,),loc=0)

        figM.canvas.draw()
    if event.key=='t':
        mode_p='zs'
        cspec1a.set_visible(True) 
        cspec2a.set_visible(True) 
        cphot1a.set_visible(False)
        cphot2a.set_visible(False)

        s1M.set_xticklabels('')
        s1M.legend((cspec1a,cspec2a),(name1+ext1s,name2+ext2s),loc=0)
        
        cspec1b.set_visible(True) 
        if t2file!='':cspec2b.set_visible(True) 
        cphot1b.set_visible(False)
        if t2file!='':cphot2b.set_visible(False)

        s2M.set_xticklabels('')
        
        cspec1c.set_visible(True) 
        cspec2c.set_visible(True) 
        cphot1c.set_visible(False)
        cphot2c.set_visible(False)

        s3M.set_xticklabels('')
        s3M.legend((cspec1c,cspec2c),(ext1s_c,ext2s_c),loc=0)
        cspec1d.set_visible(True) 
        cspec2d.set_visible(True) 
        cphot1d.set_visible(False)
        cphot2d.set_visible(False)
        s2M.legend((cspec1d,cspec2d),(ext1s_d,ext2s_d),loc=0)

        s8M.set_xticklabels('')
        cspec1e.set_visible(True) 
        cspec2e.set_visible(True) 
        cphot1e.set_visible(False)
        cphot2e.set_visible(False)
        cspec1f.set_visible(True) 
        cspec2f.set_visible(True) 
        cphot1f.set_visible(False)
        cphot2f.set_visible(False)
        
        cspec1g.set_visible(True) 
        if t2file!='':cspec2g.set_visible(True) 
        cphot1g.set_visible(False)
        if t2file!='':cphot2g.set_visible(False)
        cspec1h.set_visible(True) 
        if t2file!='':cspec2h.set_visible(True) 
        cphot1h.set_visible(False)
        if t2file!='':cphot2h.set_visible(False)
        s9M.set_xlabel(r'$z_{\rm spec}$',fontsize=20)
        s6M.set_xlabel(r'$z_{\rm spec}$',fontsize=20)
        s6M.legend((cspec1e,cspec2e),(ext1s_e,ext2s_e),loc=0)
        s5M.legend((cspec1f,cspec2f),(ext1s_f,ext2s_f),loc=0)
        if t2file !='':
            s7M.legend((cspec1g,cspec2g),(ext1s_g,ext2s_g),loc=0)
            s9M.legend((cspec1b,cspec2b),(ext1s_b,ext2s_b),loc=0)
            s8M.legend((cspec1h,cspec2h),(ext1s_h,ext2s_h),loc=0)
        else:
            s7M.legend((cspec1g,),(ext1s_g,),loc=0)
            s9M.legend((cspec1b,),(ext1s_b,),loc=0)
            s8M.legend((cspec1h,),(ext1s_h,),loc=0)
        figM.canvas.draw()


###################################
#########    MAIN     #############
###################################

def usage():
    print
    print "Usage:: "
    print "----------------------------------------------------------------------------------------------------------------------------"
    print "plot_des_results --t1 <test_file> --v1 <val_file> --t2 <test_file 2> --v2 <val_file 2> --err1 <err_cut> --err2 <err_cut 2>"
    print "----------------------------------------------------------------------------------------------------------------------------"
    print "* second set of files is optional "
    print "* error cuts are optional"
    print "* if validation file is not passed, test file is repeated"
    print "* order doesn't matter"
    print "------------------------------"
    print "-> MORE OPTIONAL PARAMETERS <-"
    print "------------------------------"
    print "* --output <output_root> "
    print "* the output root creates a file with summary of results for each test"
    print "* --name1 <name test 1>"
    print "* --name2 <name test 2>"
    print "* --w1 <weights file test 1>"
    print "* --w2 <weights file test 2>"
    print "* --frac1 <fraction between 0 and 1 of best selected galaxies (smallest errors) test 1>"
    print "* --frac2 <fraction between 0 and 1 of best selected galaxies (smallest errors) test 2>"
    print
    print "* See README file for more information"
    sys.exit(0)

if len(sys.argv) < 2 or sys.argv[1]=='-h' or sys.argv[1]=='--help':
    usage()

t1file = ''
v1file = ''
t2file= ''
v2file=''
w1file=''
w2file=''
err1=1.e10
err2=1.e10
outfile1=''
outfile2=''
outroot=''
name1_input=''
name2_input=''
frac1=-1.
frac2=-1.

try:
    opts, args = getopt.getopt(sys.argv[1:],"h",["t1=","v1=","t2=","v2=","err1=","err2=","output=","name1=","name2=","w1=","w2=","frac1=","frac2="])
except getopt.GetoptError:
    usage()
for opt, arg in opts:
    if opt in ("--t1"):
        t1file= arg
        v1file=arg
    elif opt in ("--v1"):
        v1file=arg
    elif opt in ("--t2"):
        t2file=arg
        v2file=arg
    elif opt in ("--v2"):
        v2file=arg
    elif opt in ("--err1"):
        err1=float(arg)
    elif opt in ("--err2"):
        err2=float(arg)
    elif opt in ("--output"):
        outroot=arg
    elif opt in ("--name1"):
        name1_input=arg
    elif opt in ("--name2"):
        name2_input=arg
    elif opt in ("--w1"):
        w1file=arg
    elif opt in ("--w2"):
        w2file=arg
    elif opt in ("--frac1"):
        frac1=float(arg)
    elif opt in ("--frac2"):
        frac2=float(arg)
if outroot!='':
    outfile1=outroot+'_test1.stats'
    if t2file!='': outfile2=outroot+'_test2.stats'


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
    
if t1file==v1file:
    zp_v1=zp_t1
    zs_v1=zs_t1
    ep_v1=ep_t1
    fid_v1=fid_t1
    number_v1=number_t1
    vals_w2=ones(len(zs_v1))
    vals_w2/=sum(vals_w2)
    if t1file[-3:]=='.pz':
        zmid_v1=zmid_t1
        zs_v1=zs_t1
        fid_v1=fid_t1
        number_v1=number_t1
        BIG_v1=BIG_t1
else:
    if v1file[-3:]=='.pz':
        zmid_v1,zs_v1,fid_v1,number_v1,BIG_v1=read_pdf_file(v1file)
        zp_v1,ep_v1=get_vals_pdf(zmid_v1,BIG_v1)
    else:
        zp_v1,ep_v1,zs_v1,fid_v1,number_v1=loadtxt(v1file,unpack=True)
    werror=where(ep_v1 < err1)[0]
    zp_v1=zp_v1[werror]
    zs_v1=zs_v1[werror]
    ep_v1=ep_v1[werror]
    fid_v1=fid_v1[werror]
    number_v1=number_v1[werror]
    vals_w2=ones(len(zs_v1))
    vals_w2/=sum(vals_w2)
    NTOT2=1.*len(vals_w2)

    


if t2file!='':
    if t2file[-3:]=='.pz':
        zmid_t2,zs_t2,fid_t2,number_t2,BIG_t2=read_pdf_file(t2file)
        zp_t2,ep_t2=get_vals_pdf(zmid_t2,BIG_t2)
    else:
        zp_t2,ep_t2,zs_t2,fid_t2,number_t2=loadtxt(t2file,unpack=True)
        
    if w2file!='': 
        fid_w2,number_w2,vals_w2=loadtxt(w2file,unpack=True)
        NW2=len(vals_w2)
        mask=check_fid(fid_w2,number_w2,fid_t2,number_t2)
        fid_w2=fid_w2[mask]
        number_w2=number_w2[mask]
        vals_w2=vals_w2[mask]
        
    if frac2>0:
        Dz=ep_t2
        sDz=argsort(Dz)
        keep_n=int(frac2*len(Dz))-1
        mask=sDz[0:keep_n]
        zp_t2=zp_t2[mask]
        ep_t2=ep_t2[mask]
        zs_t2=zs_t2[mask]
        fid_t2=fid_t2[mask]
        number_t2=number_t2[mask]
        if t2file[-3:]=='.pz':BIG_t2=BIG_t2[mask]
        err2=1.e10
        if w2file!='':
            fid_w2=fid_w2[mask]
            number_w2=number_w2[mask]
            vals_w2=vals_w2[mask]
    

    werror=where(ep_t2 < err2)[0]
    zp_t2=zp_t2[werror]
    zs_t2=zs_t2[werror]
    ep_t2=ep_t2[werror]
    fid_t2=fid_t2[werror]
    number_t2=number_t2[werror]


    if w2file!='':
        fid_w2=fid_w2[werror]
        number_w2=number_w2[werror]
        vals_w2=vals_w2[werror]
        NTOT2=1.*NW2
    else:
        vals_w2=ones(len(zp_t2))
        vals_w2/=sum(vals_w2)
        NTOT2=1.*len(vals_w2)
        
    if t2file==v2file:
        zp_v2=zp_t2
        zs_v2=zs_t2
        ep_v2=ep_t2
        fid_v2=fid_t2
        number_v2=number_t2
        if t2file[-3:]=='.pz':
            zmid_v2=zmid_t2
            zs_v2=zs_t2
            fid_v2=fid_t2
            number_v2=number_t2
            BIG_v2=BIG_t2
    else:
        if v2file[-3:]=='.pz':
            zmid_v2,zs_v2,fid_v2,number_v2,BIG_v2=read_pdf_file(v2file)
            zp_v2,ep_v2=get_vals_pdf(zmid_v2,BIG_v2)
        else:
            zp_v2,ep_v2,zs_v2,fid_v2,number_v2=loadtxt(v2file,unpack=True)
        werror=where(ep_v2 < err2)[0]
        zp_v2=zp_v2[werror]
        zs_v2=zs_v2[werror]
        ep_v2=ep_v2[werror]
        fid_v1=fid_v1[werror]
        number_v1=number_v1[werror]


#Window with commands help
if not notk:
    root = Tk()
    root.title("Commands help")
    myfont=tkFont.Font(size=12)
    text = Text(bg='black',fg='white',font=myfont)
    helpfile=file("help.txt") 
    helptxt= helpfile.read() 
    helpfile.close() 
    text.insert(0.0,helptxt)
    text.pack(expand=1, fill=BOTH) 
    text.config(state=DISABLED)

minz=0
maxz=1.5


#PLOT ERRORS!!
xg=linspace(-10,10,10000)
yg=1./sqrt(2.*pi)*exp(-0.5*xg*xg)
def gauss1(x): return 1./(sqrt(2.*pi)*sigma1())*exp(-0.5*((x-mu1())/sigma1())**2)
def gauss2(x): return 1./(sqrt(2.*pi)*sigma2())*exp(-0.5*((x-mu2())/sigma2())**2)
mu1=Parameter(0.)
sigma1=Parameter(1.)
mu2=Parameter(0.)
sigma2=Parameter(1.)


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

if t2file!='':
    zx2=zs_t2
    zy2=zp_t2
    e2=ep_t2
    fid2=fid_t2
    nu2=number_t2
    if t2file[-3:]=='.pz':
        name2='Test 2 (pdf)'
        zM2=zmid_t2
        B2=BIG_t2
    else:
        name2='Test 2'

    if name2_input!='' : name2=name2_input
else:
    zx2=zs_v1
    zy2=zp_v1
    e2=ep_v1
    fid2=fid_v1
    nu2=number_v1
    if v1file[-3:]=='.pz':
        name2='Val 1 (pdf)'
        zM2=zmid_v1
        B2=BIG_v1
    else:
        name2='Val 1'
    if v1file==t1file : name2=name1

err1=(zy1-zx1)/e1
we1=where(abs(err1)<=5.)[0] #ONLY USED FOR FITTING NORMAL DISTRIBUTION
err1=err1[we1]
ww1=vals_w1[we1]
ww1=ww1/sum(ww1)
err1_all=(zy1-zx1)/e1
wzero=where((e1>0.))[0]
bias_norm1=sum(err1_all[wzero]*vals_w1[wzero])/sum(vals_w1[wzero])
sigma_norm1=sqrt(sum(vals_w1[wzero]*(err1_all[wzero]-bias_norm1)**2)/sum(vals_w1[wzero]))
if t1file[-3:]=='.pz':
    D_KS1=KS_pdf(zmid_t1,BIG_t1,zx1,vals_w1)
else:
    D_KS1=KS(zy1,zx1,vals_w1)


err2=(zy2-zx2)/e2
we2=where(abs(err2)<=5.)[0] #ONLY USED FOR FITTING NORMAL DISTRIBUTION
err2=err2[we2]
ww2=vals_w2[we2]
ww2=ww2/sum(ww2)
err2_all=(zy2-zx2)/e2
wzero=where((e2>0.))[0]
bias_norm2=sum(err2_all[wzero]*vals_w2[wzero])/sum(vals_w2[wzero])
sigma_norm2=sqrt(sum(vals_w2[wzero]*(err2_all[wzero]-bias_norm2)**2)/sum(vals_w2[wzero]))
if t2file[-3:]=='.pz':
    D_KS2=KS_pdf(zmid_t2,BIG_t2,zx2,vals_w2)
else:
    D_KS2=KS(zy2,zx2,vals_w2)





G1=histogram(err1,bins=35,normed=True,weights=ww1)
G2=histogram(err2,bins=35,normed=True,weights=ww2)
pmf1=G1[0]*1.
bins1=G1[1]
pmf2=G2[0]*1.
bins2=G2[1]
xx1=0.5*(bins1[1:]+bins1[:-1])
xx2=0.5*(bins2[1:]+bins2[:-1])

fit(gauss1,[mu1,sigma1],xx1,pmf1)
fit(gauss2,[mu2,sigma2],xx2,pmf2)

#Figure 1, true plots and PDFs

rms=0.5

figB=figure(1,figsize=(16,12))
s1B=figB.add_subplot(2,2,1)
s2B=figB.add_subplot(2,2,2)
s3B=figB.add_subplot(2,2,3)
s4B=figB.add_subplot(2,2,4)

x1=y1=min(zx1)
x2=y2=max(zx1)


if t2file !='':s1B.plot(zx1,zy1,'w.',picker=5)
s1B.plot(zx1,zy1,'k.')
if t2file !='':current1,=s1B.plot(zx1[0],zy1[0],'yo',ms=10,alpha=0.8,visible=False)
s1B.set_xlabel(r'$z_{\rm spec}$',fontsize=18)
s1B.set_title(name1)
s1B.set_ylabel(r'$z_{\rm phot}$',fontsize=18)
s1B.plot([x1,x2],[x1,x2],'r-',lw=1.5)
s1B.set_xlim(minz,maxz)
s1B.set_ylim(minz,maxz)


if t2file !='':s2B.plot(zx2,zy2,'w.',picker=5)
s2B.plot(zx2,zy2,'k.')
if t2file !='':current2,=s2B.plot(zx2[0],zy2[0],'yo',ms=10,alpha=0.8,visible=False)
s2B.set_xlabel(r'$z_{\rm spec}$',fontsize=18)
s2B.set_title(name2)
s2B.set_ylabel(r'$z_{\rm phot}$',fontsize=18)
s2B.plot([x1,x2],[y1,y2],'r-',lw=1.5)
s2B.set_xlim(minz,maxz)
s2B.set_ylim(minz,maxz)


s3B.plot(xx1,pmf1,'ko',label=name1)
fit1=r'$\mu=$%5.3f, $\sigma=$%5.3f' % (mu1(),sigma1()) 
fit2=r'$\mu=$%5.3f, $\sigma=$%5.3f' % (mu2(),sigma2()) 
s3B.plot(xg,gauss1(xg),'b-',label=fit1,lw=1.5)#/max(y)*max(G[0]))
s3B.set_xlabel('$\Delta z/\sigma_{\Delta z}$',fontsize=15)
s3B.set_ylabel('Number density',fontsize=15)
if t2file !='':
    s3B.plot(xx2,pmf2,'rs',label=name2)
    s3B.plot(xg,gauss2(xg),'g-',label=fit2,lw=1.5)
s3B.legend(loc=2,prop={'size':12})
s3B.set_xlim(-6,6)



hbins=30

hist_bins=linspace(minz,maxz,hbins+1)
dB=hist_bins[1]-hist_bins[0]

s4B.hist(zs_t1,bins=hist_bins,normed=True,weights=vals_w1,color='gray',histtype='stepfilled',label=r'$z_{spec}$',lw=1.5,alpha=0.5)
if t1file[-3:]=='.pz':
    hist_pdf1=zeros(hbins)
    dz=zM1[1]-zM1[0]
    for i in xrange(len(B1)):
        for j in xrange(hbins):
            wb=where((zM1>hist_bins[j])&(zM1<=hist_bins[j+1]))[0]
            if shape(wb)[0] > 0 :
                hist_pdf1[j]+=sum(B1[i][wb])*vals_w1[i]
    hist_pdf1/=sum(hist_pdf1)
    hist_pdf1/=dB
    s4B.step(hist_bins[:-1],hist_pdf1,color='blue',label=name1,where='post',lw=2)
else:
    s4B.hist(zp_t1,bins=hist_bins,normed=True,weights=vals_w1,color='blue',histtype='step',label=name1,lw=2)

if t2file[-3:]=='.pz':
    hist_pdf2=zeros(hbins)
    dz=zM2[1]-zM2[0]
    for i in xrange(len(B2)):
        for j in xrange(hbins):
            wb=where((zM2>hist_bins[j])&(zM2<=hist_bins[j+1]))[0]
            if shape(wb)[0]>0: hist_pdf2[j]+=sum(B2[i][wb])*vals_w2[i]
    hist_pdf2/=sum(hist_pdf2)
    hist_pdf2/=dB
    s4B.step(hist_bins[:-1],hist_pdf2,color='green',label=name2,where='post',lw=2)
else:
    if t2file !='':
        s4B.hist(zp_t2,bins=hist_bins,normed=True,weights=vals_w2,color='green',histtype='step',label=name2,lw=2)
            
s4B.set_xlabel('Redshift')
s4B.set_ylabel('N(z)')
s4B.set_xlim(minz,maxz)
s4B.legend(loc=0)


figB.canvas.mpl_connect('pick_event',onpick2)
figB.canvas.mpl_connect('key_press_event',on_key)
figB.canvas.mpl_connect('key_press_event',on_key3)

#Figure 2, metrics

modephot=1
if modephot==0 : xtt=r'$z_{spec}$'
if modephot==1 : xtt=r'$z_{phot}$'


dz_bin=0.1


if t1file[-3:]=='.pz':
    B1p=metrics_pdf_w(zx1,zM1,B1,vals_w1,name1,0,1.5,dz_bin,NTOT1,mode=1)
    B1s=metrics_pdf_w(zx1,zM1,B1,vals_w1,name1,0,1.5,dz_bin,NTOT1,mode=0)
else:
    B1p=metrics_w(zx1,zy1,vals_w1,name1,0,1.5,dz_bin,NTOT1,mode=1)
    B1s=metrics_w(zx1,zy1,vals_w1,name1,0,1.5,dz_bin,NTOT1,mode=0)
    

if t1file==v1file:
    B2p=B1p
    B2s=B1s
    V1p=B1p
    V1s=B1s
else:
    if v1file[-3:]=='.pz':
        V1p=metrics_pdf(zs_v1,zmid_v1,BIG_v1,'vali 1',0,1.5,dz_bin,mode=1)
        B2p=V1p
        V1s=metrics_pdf(zs_v1,zmid_v1,BIG_v1,'vali 1',0,1.5,dz_bin,mode=0)
        B2s=V1s
    else:
        V1s=metrics(zs_v1,zp_v1,'vali 1',0,1.5,dz_bin,mode=0)
        V1p=metrics(zs_v1,zp_v1,'vali 1',0,1.5,dz_bin,mode=1)
        B2s=V1s
        B2p=V1p

D1s=diff_metrics(B1s,V1s)
D1p=diff_metrics(B1p,V1p)
if t2file!='':
    if t2file[-3:]=='.pz':
        B2p=metrics_pdf_w(zx2,zM2,B2,vals_w2,name2,0,1.5,dz_bin,NTOT2,mode=1)
        B2s=metrics_pdf_w(zx2,zM2,B2,vals_w2,name2,0,1.5,dz_bin,NTOT2,mode=0)
    else:
        B2p=metrics_w(zx2,zy2,vals_w2,name2,0,1.5,dz_bin,NTOT2,mode=1)
        B2s=metrics_w(zx2,zy2,vals_w2,name2,0,1.5,dz_bin,NTOT2,mode=0)
    if t2file==v2file:
        V2p=B2p
        V2s=B2s
    else:     
        if v2file[-3:]=='.pz':
            V2p=metrics_pdf(zs_v2,zmid_v2,BIG_v2,'vali 2',0,1.5,dz_bin,mode=1)
            V2s=metrics_pdf(zs_v2,zmid_v2,BIG_v2,'vali 2',0,1.5,dz_bin,mode=0)
        else:
            V2p=metrics(zs_v2,zp_v2,'vali 2',0,1.5,dz_bin,mode=1)
            V2s=metrics(zs_v2,zp_v2,'vali 2',0,1.5,dz_bin,mode=0)    
    
    D2s=diff_metrics(B2s,V2s)
    D2p=diff_metrics(B2p,V2p)



figM=figure(2,figsize=(16,12))
figM.subplots_adjust(hspace=0.15)
s1M=figM.add_subplot(4,2,1)
s2M=figM.add_subplot(4,2,5)
s3M=figM.add_subplot(4,2,3)
#s4M=figM.add_subplot(2,2,4)
s5M=figM.add_subplot(4,2,2)
s6M=figM.add_subplot(4,2,7)
s7M=figM.add_subplot(4,2,4)

s8M=figM.add_subplot(4,2,6)
s9M=figM.add_subplot(4,2,8)



cspec1a,=s1M.plot(B1s.bins[B1s.nzero],B1s.mean[B1s.nzero],'b->',lw=1.5,visible=False)
cspec2a,=s1M.plot(B2s.bins[B2s.nzero],B2s.mean[B2s.nzero],'g-o',lw=1.5,visible=False)
cphot1a,=s1M.plot(B1p.bins[B1p.nzero],B1p.mean[B1p.nzero],'b->',lw=1.5)
cphot2a,=s1M.plot(B2p.bins[B2p.nzero],B2p.mean[B2p.nzero],'g-o',lw=1.5)

ext1p=', rms = %.4f' % B1p.mean_rms
ext2p=', rms = %.4f' % B2p.mean_rms
ext1s=', rms = %.4f' % B1s.mean_rms
ext2s=', rms = %.4f' % B2s.mean_rms

s1M.plot([minz,maxz],[0,0],'k--')
s1M.legend((cphot1a,cphot2a),(name1+ext1p,name2+ext2p),loc=0)
s1M.set_ylabel('Bias')
s1M.set_xlim(minz,maxz)
s1M.set_xticklabels('')



cspec1c,=s3M.plot(B1s.bins[B1s.nzero],B1s.sigma[B1s.nzero],'b->',lw=1.5,visible=False)
cspec2c,=s3M.plot(B2s.bins[B2s.nzero],B2s.sigma[B2s.nzero],'g-o',lw=1.5,visible=False)
cphot1c,=s3M.plot(B1p.bins[B1p.nzero],B1p.sigma[B1p.nzero],'b->',lw=1.5)
cphot2c,=s3M.plot(B2p.bins[B2p.nzero],B2p.sigma[B2p.nzero],'g-o',lw=1.5)

s3M.plot([minz,maxz],[.12,0.12],'r--',lw=1.3)
s3M.set_xlim(minz,maxz)

s3M.set_ylabel(r'$\sigma$',fontsize=18)
s3M.set_xlim(minz,maxz)
s3M.set_xticklabels('')

ext1p_c='rms = %.4f' % B1p.sigma_rms
ext2p_c='rms = %.4f' % B2p.sigma_rms
ext1s_c='rms = %.4f' % B1s.sigma_rms
ext2s_c='rms = %.4f' % B2s.sigma_rms

s3M.legend((cphot1c,cphot2c),(ext1p_c,ext2p_c),loc=0)



cspec1d,=s2M.plot(B1s.bins[B1s.nzero],B1s.sigma68[B1s.nzero],'b->',lw=1.5,visible=False)
cspec2d,=s2M.plot(B2s.bins[B2s.nzero],B2s.sigma68[B2s.nzero],'g-o',lw=1.5,visible=False)
cphot1d,=s2M.plot(B1p.bins[B1p.nzero],B1p.sigma68[B1p.nzero],'b->',lw=1.5)
cphot2d,=s2M.plot(B2p.bins[B2p.nzero],B2p.sigma68[B2p.nzero],'g-o',lw=1.5)


s2M.plot([minz,maxz],[.12,0.12],'r--',lw=1.3)
s2M.set_ylabel(r'$\sigma_{68}$',fontsize=18)
s2M.set_xlim(minz,maxz)
s2M.set_xticklabels('')

ext1p_d='rms = %.4f' % B1p.sigma68_rms
ext2p_d='rms = %.4f' % B2p.sigma68_rms
ext1s_d='rms = %.4f' % B1s.sigma68_rms
ext2s_d='rms = %.4f' % B2s.sigma68_rms

s2M.legend((cphot1d,cphot2d),(ext1p_d,ext2p_d),loc=0)


cspec1e,=s6M.plot(B1s.bins[B1s.nzero],B1s.frac2[B1s.nzero],'b->',lw=1.5,visible=False)
cspec2e,=s6M.plot(B2s.bins[B2s.nzero],B2s.frac2[B2s.nzero],'g-o',lw=1.5,visible=False)
cphot1e,=s6M.plot(B1p.bins[B1p.nzero],B1p.frac2[B1p.nzero],'b->',lw=1.5)
cphot2e,=s6M.plot(B2p.bins[B2p.nzero],B2p.frac2[B2p.nzero],'g-o',lw=1.5)

s6M.plot([minz,maxz],[.1,0.1],'r--',lw=1.3)
s6M.set_ylabel(r'$frac (>2 \sigma)$',fontsize=18)
s6M.set_xlabel(r'$z_{\rm phot}$',fontsize=20)
s6M.set_xlim(minz,maxz)

ext1p_e='rms = %.4f' % B1p.frac2_rms
ext2p_e='rms = %.4f' % B2p.frac2_rms
ext1s_e='rms = %.4f' % B1s.frac2_rms
ext2s_e='rms = %.4f' % B2s.frac2_rms

s6M.legend((cphot1e,cphot2e),(ext1p_e,ext2p_e),loc=0)


cspec1f,=s5M.plot(B1s.bins[B1s.nzero],B1s.frac3[B1s.nzero],'b->',lw=1.5,visible=False)
cspec2f,=s5M.plot(B2s.bins[B2s.nzero],B2s.frac3[B2s.nzero],'g-o',lw=1.5,visible=False)
cphot1f,=s5M.plot(B1p.bins[B1p.nzero],B1p.frac3[B1p.nzero],'b->',lw=1.5)
cphot2f,=s5M.plot(B2p.bins[B2p.nzero],B2p.frac3[B2p.nzero],'g-o',lw=1.5)

s5M.plot([minz,maxz],[.015,0.015],'r--',lw=1.3)
s5M.set_ylabel(r'$frac (>3 \sigma)$',fontsize=18)
s5M.set_xlim(minz,maxz)
s5M.set_xticklabels('')

ext1p_f='rms = %.4f' % B1p.frac3_rms
ext2p_f='rms = %.4f' % B2p.frac3_rms
ext1s_f='rms = %.4f' % B1s.frac3_rms
ext2s_f='rms = %.4f' % B2s.frac3_rms

s5M.legend((cphot1f,cphot2f),(ext1p_f,ext2p_f),loc=0)



cspec1g,=s7M.plot(D1s.bins,D1s.biasdiff,'b->',lw=1.5,visible=False)
if t2file !='': cspec2g,=s7M.plot(D2s.bins,D2s.biasdiff,'g-o',lw=1.5,visible=False)
cphot1g,=s7M.plot(D1p.bins,D1p.biasdiff,'b->',lw=1.5)
if t2file!='': cphot2g,=s7M.plot(D2p.bins,D2p.biasdiff,'g-o',lw=1.5)

s7M.plot([minz,maxz],[0.001,0.001],'r--')
s7M.set_ylabel('$\Delta$ bias',fontsize=18)
s7M.set_xlim(minz,maxz)
s7M.set_xticklabels('')

ext1p_g='rms = %.4f' % D1p.biasrms
ext1s_g='rms = %.4f' % D1s.biasrms
if t2file !='':
    ext2p_g='rms = %.4f' % D2p.biasrms
    ext2s_g='rms = %.4f' % D2s.biasrms


if t2file !='':
    s7M.legend((cphot1g,cphot2g),(ext1p_g,ext2p_g),loc=0)
else:
    s7M.legend((cphot1g,),(ext1p_g,),loc=0)



cspec1h,=s8M.plot(D1s.bins,D1s.sigmadiff,'b->',lw=1.5,visible=False)
if t2file !='': cspec2h,=s8M.plot(D2s.bins,D2s.sigmadiff,'g-o',lw=1.5,visible=False)
cphot1h,=s8M.plot(D1p.bins,D1p.sigmadiff,'b->',lw=1.5)
if t2file!='': cphot2h,=s8M.plot(D2p.bins,D2p.sigmadiff,'g-o',lw=1.5)

s8M.plot([minz,maxz],[0.003,0.003],'r--')
s8M.set_ylabel('$\Delta \sigma$',fontsize=18)
s8M.set_xlim(minz,maxz)
s8M.set_xticklabels('')

ext1p_h='rms = %.4f' % D1p.sigmarms
ext1s_h='rms = %.4f' % D1s.sigmarms
if t2file !='':
    ext2p_h='rms = %.4f' % D2p.sigmarms
    ext2s_h='rms = %.4f' % D2s.sigmarms


if t2file !='':
    s8M.legend((cphot1h,cphot2h),(ext1p_h,ext2p_h),loc=0)
else:
    s8M.legend((cphot1h,),(ext1p_h,),loc=0)


cspec1b,=s9M.plot(D1s.bins,D1s.outfrac2diff,'b->',lw=1.5,visible=False)
if t2file !='': cspec2b,=s9M.plot(D2s.bins,D2s.outfrac2diff,'g-o',lw=1.5,visible=False)
cphot1b,=s9M.plot(D1p.bins,D1p.outfrac2diff,'b->',lw=1.5)
if t2file!='': cphot2b,=s9M.plot(D2p.bins,D2p.outfrac2diff,'g-o',lw=1.5)

s9M.plot([minz,maxz],[0.01,0.01],'r--')
s9M.set_ylabel('$\Delta(> 2 \sigma)$',fontsize=18)
s9M.set_xlim(minz,maxz)
s9M.set_xlabel(r'$z_{\rm phot}$',fontsize=20)

ext1p_b='rms = %.4f' % D1p.outfrac2rms
ext1s_b='rms = %.4f' % D1s.outfrac2rms
if t2file !='':
    ext2p_b='rms = %.4f' % D2p.outfrac2rms
    ext2s_b='rms = %.4f' % D2s.outfrac2rms


if t2file !='':
    s9M.legend((cphot1b,cphot2b),(ext1p_b,ext2p_b),loc=0)
else:
    s9M.legend((cphot1b,),(ext1p_b,),loc=0)




figM.canvas.mpl_connect('key_press_event', on_key)
figM.canvas.mpl_connect('key_press_event', on_key2)


print

print '----------------------'
print '       TEST 1         '
print '----------------------'
print 'N_total_1         : %d  ' % B1p.ntotal2
print 'Bias_1            : %.5f' % B1p.mean_all
print 'Median_1          : %.5f' % B1p.median_all
print 'Sigma_1           : %.5f' % B1p.sigma_all
print 'Sigma68_1         : %.5f' % B1p.sigma68_all
print 'Frac2_1           : %.5f' % B1p.frac2_all
print 'Frac3_1           : %.5f' % B1p.frac3_all
print 'Bias_Norm_1       : %.5f' % bias_norm1
print 'Sigma_Norm_1      : %.5f' % sigma_norm1
print 'D_N_poisson_1     : %.5f' % B1p.n_poisson_rms
print 'D_KS_test_1       : %.5f' % D_KS1
if v1file==t1file: print '** Cross metrics ** [NO VALID] with itself'
else: print '** Cross metrics **'
print 'D_bias_1          : %.5f'% D1p.biasrms
print 'D_sigma_1         : %.5f'% D1p.sigmarms
print 'D_sigma68_1       : %.5f'% D1p.sigma68rms
print 'D_frac2_1         : %.5f'% D1p.outfrac2rms
print 'D_frac3_1         : %.5f'% D1p.outfrac3rms
print 'covar_1           : %.5f'% D1p.covrms
print 
print
print



if t2file !='':
    print '----------------------'
    print '       TEST 2 '
    print '----------------------'
    print 'N_total_2         : %d  ' % B2p.ntotal2
    print 'Bias_2            : %.5f' % B2p.mean_all
    print 'Median_2          : %.5f' % B2p.median_all
    print 'Sigma_2           : %.5f' % B2p.sigma_all
    print 'Sigma68_2         : %.5f' % B2p.sigma68_all
    print 'Frac2_2           : %.5f' % B2p.frac2_all
    print 'Frac3_2           : %.5f' % B2p.frac3_all
    print 'Bias_Norm_2       : %.5f' % bias_norm2
    print 'Sigma_Norm_2      : %.5f' % sigma_norm2
    print 'D_N_poisson_2     : %.5f' % B2p.n_poisson_rms
    print 'D_KS_test_2       : %.5f' % D_KS2
    if v2file==t2file: print '** Cross metrics ** [NO VALID] with itself'
    else: print '** Cross metrics **'
    print 'D_bias_2          : %.5f'% D2p.biasrms
    print 'D_sigma_2         : %.5f'% D2p.sigmarms
    print 'D_sigma68_2       : %.5f'% D2p.sigma68rms
    print 'D_frac2_2         : %.5f'% D2p.outfrac2rms
    print 'D_frac3_2         : %.5f'% D2p.outfrac3rms
    print 'covar_2           : %.5f'% D2p.covrms




if outroot!='':
    F1=open(outfile1,'w')
    F1.write('#bin  zmid    Nin   bias    median   sigma  sigma68  frac2   frac3  bias_norm    sigma_norm  n_poisson   D_KS_Nz  D_bias  D_sigma  D_sigma68  D_frac2  D_frac3\n')
    F1.write('#Whole sample\n')
    line='%2d  %.4f  %5d  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n' % (\
    0,0.75,B1p.ntotal2,B1p.mean_all,B1p.median_all,B1p.sigma_all,\
    B1p.sigma68_all,B1p.frac2_all,B1p.frac3_all,bias_norm1,sigma_norm1,B1p.n_poisson_rms,D_KS1,D1p.biasrms,D1p.sigmarms,D1p.sigma68rms,D1p.outfrac2rms,D1p.outfrac3rms)
    F1.write(line)
    F1.write('#Binned sample\n')
    for j in xrange(len(B1p.bins)):
        line='%2d  %.4f  %.2f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n' % (\
            j+1,B1p.bins[j],B1p.n[j],B1p.mean[j],B1p.median[j],B1p.sigma[j],\
            B1p.sigma68[j],B1p.frac2[j],B1p.frac3[j],0.,0.,B1p.n_poisson[j],0.,D1p.biasdiff2[j],D1p.sigmadiff2[j],D1p.sigma68diff2[j],D1p.outfrac2diff2[j],D1p.outfrac3diff2[j])
        F1.write(line)
    F1.close()
    
if outroot!='' and t2file!='':
    F1=open(outfile2,'w')
    F1.write('#bin  zmid    Nin   bias    median   sigma  sigma68  frac2   frac3  bias_norm    sigma_norm  n_poisson  D_KS_Nz  D_bias  D_sigma  D_sigma68  D_frac2  D_frac3\n')
    F1.write('#Whole sample\n')
    line='%2d  %.4f  %5d  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n' % (\
    0,0.75,B2p.ntotal2,B2p.mean_all,B2p.median_all,B2p.sigma_all,\
    B2p.sigma68_all,B2p.frac2_all,B2p.frac3_all,bias_norm2,sigma_norm2,B2p.n_poisson_rms,D_KS2,D2p.biasrms,D2p.sigmarms,D2p.sigma68rms,D2p.outfrac2rms,D2p.outfrac3rms)
    F1.write(line)
    F1.write('#Binned sample\n')
    for j in xrange(len(B2p.bins)):
        line='%2d  %.4f  %.2f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n' % (\
            j+1,B2p.bins[j],B2p.n[j],B2p.mean[j],B2p.median[j],B2p.sigma[j],\
            B2p.sigma68[j],B2p.frac2[j],B2p.frac3[j],0.,0.,B1p.n_poisson[j],0.,D2p.biasdiff2[j],D2p.sigmadiff2[j],D2p.sigma68diff2[j],D2p.outfrac2diff2[j],D2p.outfrac3diff2[j])
        F1.write(line)
    F1.close()

#show()
figB.savefig(t1file+'_figB.png')
figM.savefig(t1file+'_figM.png')
if not notk:
    root.mainloop() 












