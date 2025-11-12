#  utils.py
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
from numpy import *
import matplotlib.pyplot as plt



def percentile(N, percent, key=lambda x:x): #copied from TPZ
    """
    N - is a list of values.
    percent - a float value from 0.0 to 1.0.
    """
    if (len(N) <= 0):
        print '\n getperc warning: n <= 0, val : ', n, 0.
        return 0.
    N=sort(N)
    k = (len(N)-1) * percent
    f = floor(k)
    c = ceil(k)
    if f == c: return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1
    
    
def perc_pdf(vals,ws,percent):
    svals=argsort(vals)
    vals=vals[svals]
    ws=ws[svals]
    stemp=cumsum(ws)
    stemp/=stemp.max()
    return interp(percent,stemp,vals)



def getperc(arr, dd): #copied from the idl code
    n = len(arr)
    if (n <= 0):
        val = 0.
        print
        print 'getperc warning: n <= 0, val : ', n, val
        return val
    sn=argsort(arr)
    arr=arr[sn]
    if (n == 1):
        val = arr[0]
        return val
    x = 1. + float(n-1)*dd
    i = floor(x)
    i1 = max([1,i])
    i2 = min([n,i+1])

    if (i1 == i2):
        if (i2 == n):
            i1 = i2-1
        else:
            i2 = i1 + 1

    x1 = float(i1)
    val1 = arr[i1-1]
    x2 = float(i2)
    val2 = arr[i2-1]
    val = (val2-val1)/(x2-x1) * (x-x1) + val1
    return val

class metrics: 
    def __init__(self,zs,zb,name,zmin,zmax,dz,mode=1): #mode 0: zspec mode 1: zphot bins
        zbins=arange(zmin,zmax+dz,dz)
        zmid=0.5*(zbins[1:]+zbins[:-1])
        self.bins=zmid
        self.zs=zs
        self.zb=zb
        self.n=zeros(len(zmid))
        self.name=name
        self.mean=zeros(len(zmid))
        self.median=zeros(len(zmid))
        self.sigma=zeros(len(zmid))
        self.sigma68=zeros(len(zmid))
        self.frac2=zeros(len(zmid))
        self.frac3=zeros(len(zmid))
        self.n_poisson=zeros(len(zmid))
        top= self.name+' : %d galaxies, mode %d:' %(len(zs),mode)
        print top
        self.ntotal=len(zs)
        for i in xrange(len(zmid)):
            wts=where((zs >= zbins[i])&(zs < zbins[i+1]))
            wtp=where((zb >= zbins[i])&(zb < zbins[i+1]))
            if mode ==0 :wt=wts
            if mode ==1 :wt=wtp
            deltaz=zb-zs # zphot-zspec
            if shape(wt)[1]==0: continue #empty bins
            tempz=sort(deltaz[wt])
            self.n[i]=shape(wt)[1]
            self.n_poisson[i]=(1.*shape(wtp)[1]-1.*shape(wts)[1])/sqrt(1.*shape(wts)[1])
            self.sigma68[i]=0.5*(percentile(tempz,0.84)-percentile(tempz,0.16))
            self.mean[i]=mean(deltaz[wt])
            self.median[i]=median(deltaz[wt])
            self.sigma[i]=std(deltaz[wt],ddof=1)
            w2=where(abs(deltaz[wt]-self.mean[i]) > 2.*self.sigma[i])
            w3=where(abs(deltaz[wt]-self.mean[i]) > 3.*self.sigma[i])
            self.frac2[i]=1.*shape(w2)[1]/(1.*shape(wt)[1])
            self.frac3[i]=1.*shape(w3)[1]/(1.*shape(wt)[1])
        nonzero=where(self.n>0)[0]

        self.nzero=nonzero
        self.mean_rms=sqrt(mean(self.mean[nonzero]**2))
        self.median_rms=sqrt(mean(self.median[nonzero]**2))
        self.sigma_rms=sqrt(mean(self.sigma[nonzero]**2))
        self.sigma68_rms=sqrt(mean(self.sigma68[nonzero]**2))
        self.frac2_rms=sqrt(mean(self.frac2[nonzero]**2))
        self.frac3_rms=sqrt(mean(self.frac3[nonzero]**2))
        self.n_poisson_rms=sqrt(mean(self.n_poisson[nonzero]**2))
        if mode ==0 :wt=where((zs >= zmin)&(zs <= zmax))
        if mode ==1 :wt=where((zb >= zmin)&(zb <= zmax))
        dz_all=zb[wt]-zs[wt]
        if mode == 0: self.mean_all2=mean(dz_all/(1.+zs[wt]))
        if mode == 1: self.mean_all2=mean(dz_all/(1.+zb[wt]))
        self.mean_all=mean(dz_all)
        self.median_all=median(dz_all)
        self.sigma_all=std(dz_all,ddof=1)
        self.sigma68_all=0.5*(percentile(dz_all,0.84)-percentile(dz_all,0.16))
        w2_all=where(abs(dz_all-self.mean_all) > 2.*self.sigma_all)
        w3_all=where(abs(dz_all-self.mean_all) > 3.*self.sigma_all)
        self.frac2_all=1.*shape(w2_all)[1]/(1.*len(dz_all))
        self.frac3_all=1.*shape(w3_all)[1]/(1.*len(dz_all))


class diff_metrics:
    def __init__(self,test,val):
        nze=intersect1d(test.nzero,val.nzero)
        self.nze=nze
        self.test=test
        self.vali=val
        self.bins=test.bins[nze]
        self.biasdiff=test.mean[nze]/(1.+test.bins[nze])-val.mean[nze]/(1.+val.bins[nze])
        self.sigmadiff = test.sigma[nze]-val.sigma[nze]
        self.sigma68diff = test.sigma68[nze]-val.sigma68[nze]
        self.outfrac2diff = test.frac2[nze]-val.frac2[nze]
        self.outfrac3diff = test.frac3[nze]-val.frac3[nze]
        
        self.biasdiff2=zeros(len(test.bins))
        self.sigmadiff2=zeros(len(test.bins))
        self.sigma68diff2=zeros(len(test.bins))
        self.outfrac2diff2=zeros(len(test.bins))
        self.outfrac3diff2=zeros(len(test.bins))
        
        self.biasdiff2[nze] =self.biasdiff
        self.sigmadiff2[nze] = self.sigmadiff
        self.sigma68diff2[nze] =self.sigma68diff
        self.outfrac2diff2[nze] = self.outfrac2diff
        self.outfrac3diff2[nze] = self.outfrac3diff
        
        self.biasrms = sqrt(mean(self.biasdiff**2))
        self.sigmarms = sqrt(mean(self.sigmadiff**2))
        self.sigma68rms = sqrt(mean(self.sigma68diff**2))
        self.outfrac2rms = sqrt(mean(self.outfrac2diff**2))
        self.outfrac3rms = sqrt(mean(self.outfrac3diff**2))
        nbin = len(self.biasdiff)
        ncov = nbin*(nbin-1)/2
        covdiff2 = zeros(ncov)
        icov = -1
        for i in  xrange(nbin-2+1):
            for j in xrange(i+1, nbin-1+1):
                icov = icov + 1
                covdiff2[icov] = abs(self.biasdiff[i]*self.biasdiff[j])
        
        self.covrms = sqrt(mean(covdiff2))

def get_vals_pdf(zmid,BPDF):
    zphot=zeros(len(BPDF))
    err=zeros(len(BPDF))
    for i in xrange(len(BPDF)):
        pdf=BPDF[i]
        if sum(pdf)>0: pdf/=sum(pdf)
        zphot[i]=sum(pdf*zmid)
        err[i]=sqrt(sum(pdf*(zmid-zphot[i])**2))
        if isnan(err[i]) : err[i]=10.
    return zphot,err

def read_pdf_file(filename):
    F=open(filename,'r')
    ALL=F.readlines()
    F.close()
    L0=ALL[0].strip()
    L0=L0.split(' ')
    if L0[0]!='': Ng=int(L0[0])
    zmid=[]
    Ng=len(ALL)-1
    for i in xrange(1,len(L0)):
        if L0[i]!='': zmid.append(float(L0[i]))
    zmid=array(zmid)
    zs=zeros(Ng)
    field_id=zeros(Ng,dtype='int')
    number=zeros(Ng,dtype='int')
    win=where((zmid>=0.)&(zmid<=1.5))[0]
    zmid=zmid[win]
    BIG=zeros((Ng,len(zmid)))
    for i in xrange(Ng):
        L=ALL[i+1].strip()
        L=L.split(' ')
        temp=[]
        for j in xrange(len(L)):
            if L[j]!='': temp.append(float(L[j]))
        zs[i]=temp[0]
        field_id[i]=int(temp[1])
        number[i]=int(temp[2])
        tpdf=array(temp[3:])
        BIG[i]=tpdf[win]
        if isnan(BIG[i][0]): BIG[i]=0.
        if sum(BIG[i]) >0 : BIG[i]/=sum(BIG[i])
    return zmid,zs,field_id,number,BIG


class metrics_pdf:
    def __init__(self,zs,z_pdf,B_pdf,name,zmin,zmax,dz,mode=1): #mode 0: zspec mode 1: zphot bins
        zbins=arange(zmin,zmax+dz,dz)
        zmid=0.5*(zbins[1:]+zbins[:-1])
        self.bins=zmid
        self.zs=zs
        self.zb=sum(z_pdf*B_pdf[:],axis=1)
        self.n=zeros(len(zmid))
        self.name=name
        self.mean=zeros(len(zmid))
        self.median=zeros(len(zmid))
        self.sigma=zeros(len(zmid))
        self.sigma68=zeros(len(zmid))
        self.frac2=zeros(len(zmid))
        self.frac3=zeros(len(zmid))
        self.n_poisson=zeros(len(zmid))
        top= self.name+' : %d galaxies (PDF),mode %d' %(len(zs),mode)
        print top
        self.ntotal=len(zs)
        
        for i in xrange(len(zmid)):
            if mode==1:
                wt=where((z_pdf >= zbins[i])&(z_pdf < zbins[i+1]))[0]
                wts=where((zs >= zbins[i])&(zs < zbins[i+1]))
                for j in xrange(self.ntotal):
                    pdf=B_pdf[j]
                    self.n[i]+=sum(pdf[wt])
                    self.mean[i]+=sum(pdf[wt]*(z_pdf[wt]-zs[j]))
                self.mean[i]/=self.n[i]
                self.n_poisson[i]=(1.*self.n[i]-1.*shape(wts)[1])/sqrt(1.*shape(wts)[1])
                list_dz=array([])
                list_pdf=array([])
                for j in xrange(self.ntotal):
                    pdf=B_pdf[j]
                    self.sigma[i]+=sum(pdf[wt]*((z_pdf[wt]-zs[j])-self.mean[i])**2)
                    list_dz=concatenate((list_dz,(z_pdf[wt]-zs[j])))
                    list_pdf=concatenate((list_pdf,pdf[wt]))
                self.sigma[i]/=self.n[i]
                self.sigma[i]=sqrt(self.sigma[i])
                self.sigma68[i]=0.5*(perc_pdf(list_dz,list_pdf,0.84)-perc_pdf(list_dz,list_pdf,0.16))
                self.median[i]=perc_pdf(list_dz,list_pdf,0.5)
                w2=0.
                w3=0.
                for j in xrange(self.ntotal):
                    pdf=B_pdf[j]
                    ddz=abs((z_pdf[wt]-zs[j])-self.mean[i])
                    w2t=where(ddz>2.*self.sigma[i])[0]
                    if len(w2t)>0 : w2+=sum(pdf[wt[w2t]])
                    w3t=where(ddz>3.*self.sigma[i])[0]
                    if len(w3t)>0 : w3+=sum(pdf[wt[w3t]])
                self.frac2[i]=w2/self.n[i]
                self.frac3[i]=w3/self.n[i]
            if mode==0:
                wts=where((zs >= zbins[i])&(zs < zbins[i+1]))
                self.n[i]=1.*shape(wts)[1]
                Nin=int(self.n[i])
                if Nin >0: wt=wts[0]
                else: continue
                for ij in xrange(Nin):
                    j=wt[ij]
                    pdf=B_pdf[j]
                    self.mean[i]+=sum(pdf*(z_pdf-zs[j]))
                self.mean[i]/=self.n[i]
                list_dz=array([])
                list_pdf=array([])
                for ij in xrange(Nin):
                    j=wt[ij]
                    pdf=B_pdf[j]
                    self.sigma[i]+=sum(pdf*((z_pdf-zs[j])-self.mean[i])**2)
                    list_dz=concatenate((list_dz,(z_pdf-zs[j])))
                    list_pdf=concatenate((list_pdf,pdf))
                self.sigma[i]/=self.n[i]
                self.sigma[i]=sqrt(self.sigma[i])
                self.sigma68[i]=0.5*(perc_pdf(list_dz,list_pdf,0.84)-perc_pdf(list_dz,list_pdf,0.16))
                self.median[i]=perc_pdf(list_dz,list_pdf,0.5)
                w2=0.
                w3=0.
                for ij in xrange(Nin):
                    j=wt[ij]
                    pdf=B_pdf[j]
                    ddz=abs((z_pdf-zs[j])-self.mean[i])
                    w2t=where(ddz>2.*self.sigma[i])[0]
                    if len(w2t)>0 : w2+=sum(pdf[w2t])
                    w3t=where(ddz>3.*self.sigma[i])[0]
                    if len(w3t)>0 : w3+=sum(pdf[w3t])
                self.frac2[i]=w2/self.n[i]
                self.frac3[i]=w3/self.n[i]
            
        
        nonzero=where(self.n>0)[0]
        self.nzero=nonzero
        self.mean_rms=sqrt(mean(self.mean[nonzero]**2))
        self.median_rms=sqrt(mean(self.median[nonzero]**2))
        self.sigma_rms=sqrt(mean(self.sigma[nonzero]**2))
        self.sigma68_rms=sqrt(mean(self.sigma68[nonzero]**2))
        self.frac2_rms=sqrt(mean(self.frac2[nonzero]**2))
        self.frac3_rms=sqrt(mean(self.frac3[nonzero]**2))
        self.n_poisson_rms=sqrt(mean(self.n_poisson[nonzero]**2))
        
        self.n_all=0.
        self.mean_all=0.
        self.sigma_all=0.
        
        for j in xrange(self.ntotal):
            pdf=B_pdf[j]
            self.n_all+=sum(pdf)
            self.mean_all+=sum(pdf*(z_pdf-zs[j]))
        self.mean_all/=self.n_all
        list_dz=array([])
        list_pdf=array([])
        for j in xrange(self.ntotal):
            pdf=B_pdf[j]
            self.sigma_all+=sum(pdf*((z_pdf-zs[j])-self.mean_all)**2)
            list_dz=concatenate((list_dz,(z_pdf-zs[j])))
            list_pdf=concatenate((list_pdf,pdf))
        self.sigma_all/=self.n_all
        self.sigma_all=sqrt(self.sigma_all)
        self.sigma68_all=0.5*(perc_pdf(list_dz,list_pdf,0.84)-perc_pdf(list_dz,list_pdf,0.16))
        self.median_all=perc_pdf(list_dz,list_pdf,0.5)
        w2=0.
        w3=0.
        for j in xrange(self.ntotal):
            pdf=B_pdf[j]
            ddz=abs((z_pdf-zs[j])-self.mean_all)
            w2t=where(ddz>2.*self.sigma_all)[0]
            if len(w2t)>0 : w2+=sum(pdf[w2t])
            w3t=where(ddz>3.*self.sigma_all)[0]
            if len(w3t)>0 : w3+=sum(pdf[w3t])
        self.frac2_all=w2/self.n_all
        self.frac3_all=w3/self.n_all
 
class metrics_w: 
    def __init__(self,zs,zb,w_s,name,zmin,zmax,dz,ntt,mode=1): #mode 0: zspec mode 1: zphot bins
        zbins=arange(zmin,zmax+dz,dz)
        zmid=0.5*(zbins[1:]+zbins[:-1])
        self.bins=zmid
        self.zs=zs
        self.zb=zb
        self.n=zeros(len(zmid))
        self.name=name
        self.mean=zeros(len(zmid))
        self.median=zeros(len(zmid))
        self.sigma=zeros(len(zmid))
        self.sigma68=zeros(len(zmid))
        self.frac2=zeros(len(zmid))
        self.frac3=zeros(len(zmid))
        self.n_poisson=zeros(len(zmid))
        top= self.name+' : %d galaxies, mode %d:' %(len(zs),mode)
        print top
        self.ntotal=1.*ntt
        self.ntotal2=len(zs)
        for i in xrange(len(zmid)):
            wts=where((zs >= zbins[i])&(zs < zbins[i+1]))
            wtp=where((zb >= zbins[i])&(zb < zbins[i+1]))
            if mode ==0 :wt=wts
            if mode ==1 :wt=wtp
            deltaz=zb-zs # zphot-zspec
            if shape(wt)[1]==0: continue #empty bins
            tempz=sort(deltaz[wt])
            self.n[i]=sum(w_s[wt])*1.
            self.n_poisson[i]=(sum(w_s[wtp]*self.ntotal)-1.*sum(w_s[wts]*self.ntotal))/sqrt(1.*sum(w_s[wts]*self.ntotal))
            list_dz=deltaz[wt]
            list_pdf=w_s[wt]
            self.sigma68[i]=0.5*(perc_pdf(list_dz,list_pdf,0.84)-perc_pdf(list_dz,list_pdf,0.16))
            self.median[i]=perc_pdf(list_dz,list_pdf,0.5)
            self.mean[i]=sum(deltaz[wt]*w_s[wt])/sum(w_s[wt])
            self.sigma[i]=sqrt(sum(w_s[wt]*(deltaz[wt]-self.mean[i])**2)/self.n[i])
        
            w2=where(abs(deltaz[wt]-self.mean[i]) > 2.*self.sigma[i])[0]
            w3=where(abs(deltaz[wt]-self.mean[i]) > 3.*self.sigma[i])[0]
            self.frac2[i]=1.*sum(w_s[wt][w2])/self.n[i]
            self.frac3[i]=1.*sum(w_s[wt][w3])/self.n[i]

        nonzero=where(self.n>0)[0]

        self.nzero=nonzero
        self.mean_rms=sqrt(mean(self.mean[nonzero]**2))
        self.median_rms=sqrt(mean(self.median[nonzero]**2))
        self.sigma_rms=sqrt(mean(self.sigma[nonzero]**2))
        self.sigma68_rms=sqrt(mean(self.sigma68[nonzero]**2))
        self.frac2_rms=sqrt(mean(self.frac2[nonzero]**2))
        self.frac3_rms=sqrt(mean(self.frac3[nonzero]**2))
        self.n_poisson_rms=sqrt(mean(self.n_poisson[nonzero]**2))
        if mode ==0 :wt=where((zs >= zmin)&(zs <= zmax))[0]
        if mode ==1 :wt=where((zb >= zmin)&(zb <= zmax))[0]
        dz_all=zb[wt]-zs[wt]
        wsa=w_s[wt]
        if mode == 0:
            temp=dz_all/(1.+zs[wt])
            self.mean_all2=sum(wsa*temp)/sum(wsa)
        if mode == 1:
            temp=dz_all/(1.+zb[wt])
            self.mean_all2=sum(wsa*temp)/sum(wsa)
            
        self.mean_all=sum(dz_all*wsa)/sum(wsa)
        self.median_all=perc_pdf(dz_all,wsa,0.5)
        self.sigma68_all=0.5*(perc_pdf(dz_all,wsa,0.84)-perc_pdf(dz_all,wsa,0.16))
        self.sigma_all=sqrt(sum(wsa*(dz_all-self.mean_all)**2)/sum(wsa))
 
        w2_all=where(abs(dz_all-self.mean_all) > 2.*self.sigma_all)
        w3_all=where(abs(dz_all-self.mean_all) > 3.*self.sigma_all)
        self.frac2_all=1.*sum(wsa[w2_all])/sum(wsa)
        self.frac3_all=1.*sum(wsa[w3_all])/sum(wsa)
        self.n=self.n*self.ntotal

class metrics_pdf_w:
    def __init__(self,zs,z_pdf,B_pdf,w_s,name,zmin,zmax,dz,ntt,mode=1,oall='no'): #mode 0: zspec mode 1: zphot bins
        zbins=arange(zmin,zmax+dz,dz)
        zmid=0.5*(zbins[1:]+zbins[:-1])
        self.bins=zmid
        self.zs=zs
        self.zb=sum(z_pdf*B_pdf[:],axis=1)
        self.n=zeros(len(zmid))
        self.n2=zeros(len(zmid))
        self.name=name
        self.mean=zeros(len(zmid))
        self.median=zeros(len(zmid))
        self.sigma=zeros(len(zmid))
        self.sigma68=zeros(len(zmid))
        self.frac2=zeros(len(zmid))
        self.frac3=zeros(len(zmid))
        self.n_poisson=zeros(len(zmid))
        top= self.name+' : %d galaxies (PDF),mode %d' %(len(zs),mode)
        print top
        self.ntotal=1.*ntt
        self.ntotal2=len(zs)
        
        if oall=='no':
            for i in xrange(len(zmid)):
                if mode==1:
                    wt=where((z_pdf >= zbins[i])&(z_pdf < zbins[i+1]))[0]
                    wts=where((zs >= zbins[i])&(zs < zbins[i+1]))[0]
                    for j in xrange(self.ntotal2):
                        pdf=B_pdf[j]
                        D_Z=z_pdf[wt]-zs[j]
                        self.n[i]+=w_s[j]
                        self.n2[i]+=sum(pdf[wt]*w_s[j])
                        if sum(pdf[wt]>0): self.mean[i]+=w_s[j]*sum(pdf[wt]*D_Z)/sum(pdf[wt])
                    self.mean[i]/=self.n[i]
                    self.n_poisson[i]=(1.*self.n2[i]*self.ntotal-sum(w_s[wts]*self.ntotal))/sqrt(sum(w_s[wts]*self.ntotal))
                    list_dz=[]
                    list_pdf=[]
                    for j in xrange(self.ntotal2):
                        pdf=B_pdf[j]
                        D_Z=z_pdf[wt]-zs[j]
                        if sum(pdf[wt])>0.: self.sigma[i]+=w_s[j]*(sum(pdf[wt]*(D_Z-self.mean[i]))/sum(pdf[wt]))**2
                        tsum=sum(pdf[wt]*D_Z)
                        if sum(pdf[wt])>0. : tsum=tsum/sum(pdf[wt])
                        list_dz.append(tsum)
                        list_pdf.append(w_s[j])
                    list_dz=array(list_dz)
                    list_pdf=array(list_pdf)
                    self.sigma[i]/=self.n[i]
                    self.sigma[i]=sqrt(self.sigma[i])
                    self.sigma68[i]=0.5*(perc_pdf(list_dz,list_pdf,0.84)-perc_pdf(list_dz,list_pdf,0.16))
                    self.median[i]=perc_pdf(list_dz,list_pdf,0.5)
                    w2=0.
                    w3=0.
                    for j in xrange(self.ntotal2):
                        pdf=B_pdf[j]
                        D_Z=z_pdf[wt]-zs[j]
                        ddz=sum(pdf[wt]*abs(D_Z-self.mean[i]))
                        if sum(pdf[wt]) >0.: ddz=ddz/sum(pdf[wt])                        
                        if ddz > 2.*self.sigma[i]: w2+=w_s[j]
                        if ddz > 3.*self.sigma[i]: w3+=w_s[j]

                        #w2t=where(ddz>2.*self.sigma[i])[0]
                        #if len(w2t)>0 : w2+=sum(pdf[wt[w2t]])
                        #w3t=where(ddz>3.*self.sigma[i])[0]
                        #if len(w3t)>0 : w3+=sum(pdf[wt[w3t]])
                    self.frac2[i]=w2/self.n[i]
                    self.frac3[i]=w3/self.n[i]
                if mode==0:
                    wts=where((zs >= zbins[i])&(zs < zbins[i+1]))[0]
                    self.n[i]=1.*sum(w_s[wts])
                    Nin=int(self.n[i])
                    Nin2=shape(wts)[0]
                    if Nin2 >0: wt=wts
                    else: continue
                    for ij in xrange(Nin2):
                        j=wt[ij]
                        pdf=B_pdf[j]
                        D_Z=z_pdf-zs[j]
                        self.mean[i]+=w_s[j]*sum(pdf*D_Z)
                    self.mean[i]/=self.n[i]
                    list_dz=[]
                    list_pdf=[]
                    for ij in xrange(Nin2):
                        j=wt[ij]
                        pdf=B_pdf[j]
                        D_Z=z_pdf-zs[j]
                        self.sigma[i]+=w_s[j]*sum(pdf*(D_Z-self.mean[i]))**2
                        tsum=sum(pdf*D_Z)
                        list_dz.append(tsum)
                        list_pdf.append(w_s[j])
                    list_dz=array(list_dz)
                    list_pdf=array(list_pdf)
                    self.sigma[i]/=self.n[i]
                    self.sigma[i]=sqrt(self.sigma[i])
                    self.sigma68[i]=0.5*(perc_pdf(list_dz,list_pdf,0.84)-perc_pdf(list_dz,list_pdf,0.16))
                    self.median[i]=perc_pdf(list_dz,list_pdf,0.5)
                    w2=0.
                    w3=0.
                    for ij in xrange(Nin2):
                        j=wt[ij]
                        pdf=B_pdf[j]
                        D_Z=z_pdf-zs[j]
                        ddz=sum(pdf*abs(D_Z-self.mean[i]))
                        if ddz > 2.*self.sigma[i]: w2+=w_s[j]
                        if ddz > 3.*self.sigma[i]: w3+=w_s[j]
                        #w2t=where(ddz>2.*self.sigma[i])[0]
                        #if len(w2t)>0 : w2+=sum(pdf[w2t])
                        #w3t=where(ddz>3.*self.sigma[i])[0]
                        #if len(w3t)>0 : w3+=sum(pdf[w3t])
                    self.frac2[i]=w2/self.n[i]
                    self.frac3[i]=w3/self.n[i]
            
        if oall=='yes':
            for i in xrange(len(zmid)):
                if mode==1:
                    wt=where((z_pdf >= zbins[i])&(z_pdf < zbins[i+1]))[0]
                    wts=where((zs >= zbins[i])&(zs < zbins[i+1]))[0]
                    for j in xrange(self.ntotal2):
                        pdf=B_pdf[j]*w_s[j]
                        self.n2[i]+=sum(pdf[wt])
                    self.n_poisson[i]=(1.*self.n2[i]*self.ntotal-sum(w_s[wts]*self.ntotal))/sqrt(sum(w_s[wts]*self.ntotal))

        nonzero=where(self.n>0)[0]
        self.nzero=nonzero
        self.mean_rms=sqrt(mean(self.mean[nonzero]**2))
        self.median_rms=sqrt(mean(self.median[nonzero]**2))
        self.sigma_rms=sqrt(mean(self.sigma[nonzero]**2))
        self.sigma68_rms=sqrt(mean(self.sigma68[nonzero]**2))
        self.frac2_rms=sqrt(mean(self.frac2[nonzero]**2))
        self.frac3_rms=sqrt(mean(self.frac3[nonzero]**2))
        self.n_poisson_rms=sqrt(mean(self.n_poisson[nonzero]**2))
        
        self.n_all=0.
        self.mean_all=0.
        self.sigma_all=0.
        
        for j in xrange(self.ntotal2):
            pdf=B_pdf[j]
            D_Z=z_pdf-zs[j]
            self.n_all+=w_s[j]
            self.mean_all+=sum(pdf*D_Z)*w_s[j]
        self.mean_all/=self.n_all
        list_dz=[]
        list_pdf=[]
        for j in xrange(self.ntotal2):
            pdf=B_pdf[j]
            D_Z=z_pdf-zs[j]
            self.sigma_all+=w_s[j]*sum(pdf*(D_Z-self.mean_all))**2
            tsum=sum(pdf*D_Z)
            list_dz.append(tsum)
            list_pdf.append(w_s[j])
        list_dz=array(list_dz)
        list_pdf=array(list_pdf)
        self.sigma_all/=self.n_all
        self.sigma_all=sqrt(self.sigma_all)
        self.sigma68_all=0.5*(perc_pdf(list_dz,list_pdf,0.84)-perc_pdf(list_dz,list_pdf,0.16))
        self.median_all=perc_pdf(list_dz,list_pdf,0.5)
        w2=0.
        w3=0.
        for j in xrange(self.ntotal2):
            pdf=B_pdf[j]
            D_Z=z_pdf-zs[j]
            ddz=sum(pdf*abs(D_Z-self.mean_all))
            if ddz>2.*self.sigma_all : w2+=w_s[j]
            if ddz>3.*self.sigma_all : w3+=w_s[j]
            #w2t=where(ddz>2.*self.sigma_all)[0]
            #if len(w2t)>0 : w2+=sum(w2t)
            #w3t=where(ddz>3.*self.sigma_all)[0]
            #if len(w3t)>0 : w3+=sum(w3t)
        self.frac2_all=w2/self.n_all
        self.frac3_all=w3/self.n_all
        self.n=self.n*self.ntotal
        
def check_fid(fw,nw,ft,nt):
    mask=[]
    for i in xrange(len(ft)):
        wt=where((fw==ft[i])&(nw==nt[i]))[0]
        if shape(wt)[0]==1: mask.append(wt[0])
    return array(mask)

def edfw(data,ws,x):
    n=1.*len(data)
    ws=ws/sum(ws)
    w=where((data<=x),ws,0.)
    return sum(w*n)/n
    
    
def edfw_pdf(zm,BIG,ws,x):
    n=1.*len(BIG)
    ws=ws/sum(ws)
    zin=where(zm<=x)[0]
    w=0.
    for j in xrange(len(BIG)):
        pdf=BIG[j]
        if sum(BIG[j]) > 0: pdf=pdf/sum(BIG[j])
        w+=sum(pdf[zin])*ws[j]*n
    return w/n


def KS(zphot,zspec,weights):
    
    zz=linspace(0,1.5,1000)
    Nd=len(zz)
    ecf_zs=zeros(Nd)
    ecf_zp=zeros(Nd)
    for iz in xrange(len(zz)):
        ecf_zs[iz]=edfw(zspec,weights,zz[iz])
        ecf_zp[iz]=edfw(zphot,weights,zz[iz])
    Dk=max(abs(ecf_zs-ecf_zp))
    return Dk
        
def KS_pdf(zm,BP,zspec,weights):
    
    zz=linspace(0,1.5,1000)
    Nd=len(zz)
    ecf_zs=zeros(Nd)
    ecf_zp=zeros(Nd)
    for iz in xrange(len(zz)):
        ecf_zs[iz]=edfw(zspec,weights,zz[iz])
        ecf_zp[iz]=edfw_pdf(zm,BP,weights,zz[iz])
    Dk=max(abs(ecf_zs-ecf_zp))
    return Dk


