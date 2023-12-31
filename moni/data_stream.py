# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 22:06:28 2018

@author: Dell
"""
from obspy import Stream,Trace,read,UTCDateTime
import h5py
import obspy
def stream_to_list(datain=Stream(),nsmp=2048,starttime=-9999.0):
    data1=[];
    
    for i in range(0,len(datain)):
        st=datain[i];
        dt=0.0
        if type(starttime)!=float:
            dt=st.stats.starttime-starttime
        idt=round(dt/st.stats.delta)
        if idt>=0:
            data=[0,]*idt+st.data.tolist()
        else:
            data=st.data.tolist()[-idt:]
        dn=nsmp-len(data);
        if(dn<=0):
            data1.append(data[0:nsmp])
        else:
            #print('nsmp is small, zero filled,nsmp,npts',nsmp,st.stats.npts,i)
            data1.append(data[0:]+[0,]*dn)
        #print(len(data1[0]),idt,dn,st.stats.npts)
    data=[[[data1[ir][j],data1[ir+1][j],data1[ir+2][j]] for j in range(nsmp)] for ir in range(0,len(datain),3)]
    return data;
def load_sac(saclist=['D:/Italy/201610.20to24continue_sac/2016-08-20T00.00.00.000000ZIV.ASSB.HHE.SAC','HHN.sac','HHZ.sac'],
             win=[],
             freq=[2.0,8.0],
             nsmpout=1024):
    for i in range(0,len(saclist)):
        try:
            st_tmp=read(saclist[i]);
            st_tmp.interpolate(sampling_rate=20)
            st_tmp=st_tmp.slice(starttime=win[0],endtime=win[0]+win[1],keep_empty_traces=True);
            st_tmp.detrend("demean")
            st_tmp.detrend("linear")
            st_tmp.filter('bandpass',freqmin=freq[0],freqmax=freq[1])
            st_tmp.taper(max_percentage=0.00001)
            #print(st_tmp[0].stats)
        except:
            st_tmp=Stream(Trace(np.array([0,]*100)))
            print('no data in',saclist[i]);
        if i==0:
            st=st_tmp;
        else:
            st=st+st_tmp;
    data1event=stream_to_list(datain=st,nsmp=nsmpout)
    data1event=np.clip(np.nan_to_num(np.array(data1event)),-1.0e10,1.0e10).tolist()
    print('st',len(st))
    return data1event


def h5py_getObspyTrace(trcn,hf):
    dataset=hf.get(trcn)
    #print(dataset.attrs['delta'],dataset.attrs['sampling_rate'])
    if dataset==None:
        return dataset
    tr = obspy.Trace(data=np.array(dataset))
    tr.stats.starttime = UTCDateTime(dataset.attrs['starttime'])
    tr.stats.delta = dataset.attrs['delta']
    tr.stats.channel = dataset.attrs['channel']
    tr.stats.station = dataset.attrs['station']
    tr.stats.network = dataset.attrs['network']
    return tr

def h5py_getObspyStream(trcns,hf):
    trcs=[]
    for i,trcn in enumerate(trcns):
        tr=h5py_getObspyTrace(trcn,hf)
        if tr==None:
            tr=Trace(np.array([0,]*30))
        trcs.append(tr)
    return obspy.Stream(trcs)

def h5py_getObspyStream_data(trcns,hf,prep={'freq':[2.0,8.0],'nsmpout':1024},win=[]):
    freq=prep['freq']
    nsmpout=prep['nsmpout']
    nbrok=0
    for i,trcn in enumerate(trcns):
        if type(trcn)==str:
            st_tmp=h5py_getObspyStream(trcns=[trcn],hf=hf)
        else:
            st_tmp=Stream(trcn.copy())
        try:
            #print(win,st_tmp[0].stats.starttime,st_tmp[0].stats.delta,st_tmp[0].stats.npts )
            st_tmp.filter('bandpass',freqmin=freq[0],freqmax=freq[1])
            st_tmp.interpolate(sampling_rate=20)
            #print(win,st_tmp[0].stats.starttime,st_tmp[0].stats.delta,st_tmp[0].stats.npts)
            st_tmp=st_tmp.slice(starttime=win[0],endtime=win[0]+win[1],keep_empty_traces=True);
            st_tmp.detrend("demean")
            st_tmp.detrend("linear")
            st_tmp.taper(max_percentage=0.00001)
            #print(st_tmp[0].stats)
        except:
            st_tmp=Stream(Trace(np.array([0,]*100)))
            print('532no data in',trcn,st_tmp[0].data);
            nbrok=nbrok+1
        if i==0:
            st=st_tmp;
        else:
            st=st+st_tmp;
    data1event=stream_to_list(datain=st,nsmp=nsmpout,starttime=win[0])
    data1event=np.clip(np.nan_to_num(np.array(data1event)),-1.0e10,1.0e10).tolist()
    #print('st',len(st))
    return {'data':data1event,'nbrok':nbrok}


def h5py_getObspyStream_data1(trcns,hf,freq=[2.0,8.0],nsmpout=1024,win=[]):
    paz_wa = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
                'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
    wa_amp=[]
    nbrok=0
    for i,trcn in enumerate(trcns):
        st_tmp=h5py_getObspyStream(trcns=[trcn],hf=hf)
        try:
            st_tmp.filter('bandpass',freqmin=freq[0],freqmax=freq[1])
            st_tmp_mag=st_tmp.copy().slice(starttime=win[0]+30.0,endtime=win[0]+30.0+35.0,keep_empty_traces=True);
            st_tmp_mag.simulate(paz_remove = None, paz_simulate = paz_wa)
            st_tmp_mag_amp=np.max(np.abs(st_tmp_mag[0].data))
            wa_amp.append(st_tmp_mag_amp)
            #print('554',st_tmp[0].data.shape)
            
            st_tmp.interpolate(sampling_rate=20)
            st_tmp=st_tmp.slice(starttime=win[0],endtime=win[0]+win[1],keep_empty_traces=True);
            st_tmp.detrend("demean")
            st_tmp.detrend("linear")
            st_tmp.taper(max_percentage=0.00001)

        except Exception as e:
            print('no data in',e,trcn,st_tmp[0].data.shape,st_tmp_mag[0].data.shape);
            st_tmp=Stream(Trace(np.array([0,]*100)))
            wa_amp.append(0.0)
            nbrok=nbrok+1
        if i==0:
            st=st_tmp;
        else:
            st=st+st_tmp;
    data1event=stream_to_list(datain=st,nsmp=nsmpout,starttime=win[0])
    data1event=np.clip(np.nan_to_num(np.array(data1event)),-1.0e10,1.0e10).tolist()
    
    return {'wa_amp':np.array(wa_amp),'data':data1event,'nbrok':nbrok}
    

def h5py_get(trcn,hf):
    chn=[]
    data=[]
    tmpdata=hf.get(trcn+'.HHE')
    if tmpdata!=None:
        chn.append('.HHE')
        data.append(tmpdata)
    else:
        tmpdata=hf.get(trcn+'.HH2')
        if tmpdata!=None:
            chn.append('.HH2')
            data.append(tmpdata)
        else:
            chn.append('None')
    
    tmpdata=hf.get(trcn+'.HHN')
    if tmpdata!=None:
        chn.append('.HHN')
        data.append(tmpdata)
    else:
        tmpdata=hf.get(trcn+'.HH1')
        if tmpdata!=None:
            chn.append('.HH1')
            data.append(tmpdata)
        else:
            chn.append('None')
    
    tmpdata=hf.get(trcn+'.HHZ')
    if tmpdata!=None:
        chn.append('.HHZ')
        data.append(tmpdata)
    else:
        chn.append('None')
    return {'channel':chn,'data':data}


def testdata():
    
    with open('D:/loca_aug/plot/tmp.txt') as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        stxy=np.array([[float(i[0]),float(i[1])] for i in txtlines])
        stnam=[i[3]+'.'+i[4] for i in txtlines]
        tmp1,tmp2=coordinate2model_stn(stnx=stxy[:,0],stny=stxy[:,1])
        stxy[:,0]=tmp1
        stxy[:,1]=tmp2
        stns=[{'X':stxy[i,0],'Y':stxy[i,1],'stnam':stnam[i]} for i in range(len(stnam))] #'E'->x
        stns=sort_stn(stns=stns)
        print(stns)
        stxy=[[i['X'],i['Y']] for i in stns]
        stnam=[i['stnam'] for i in stns]
        
    with open('D:/loca_aug/plot/tmp1.txt') as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        #eve_samp=[[i[0],float(i[4]),float(i[5]),float(i[6]),float(i[7])] for i in txtlines]
        ev_t0=[i[0] for i in txtlines]
        ev_xyz=np.array([[float(i[1]),float(i[2]),float(i[3])*0.001] for i in txtlines])
        tmp1,tmp2=coordinate2model_stn(stnx=ev_xyz[:,0],stny=ev_xyz[:,1])
        ev_xyz[:,0]=tmp1
        ev_xyz[:,1]=tmp2
        #print(tmp1,tmp2)
        #ev_mag=[float(i[7]) for i in txtlines]
        #eve_samp=[[ev_t0[i],ev_xyz[i,0],ev_xyz[i,1],ev_xyz[i,2],ev_mag[i]] for i in range(len(ev_t0))]
    h5obj=h5py.File('D:/loca_aug/data/Ridgecrest/Ridgecrest_abvm2.5.h5','r')
    
    xdata=[]
    srcxyz=[]
    stns1=[]
    for i in range(len(txtlines[0:500])):
        prelist=txtlines[i][0].replace(':','.');
        saclist=[]
        for j in range(len(stns)):
            #saclist.append(prelist+stnam[j])
            saclist=saclist+[prelist+stnam[j]+'.HHE',
                             prelist+stnam[j]+'.HHN',
                             prelist+stnam[j]+'.HHZ']
        tmp=h5py_getObspyStream_data(trcns=saclist,hf=h5obj,win=[UTCDateTime(ev_t0[i])-15,30.0])
        if tmp['nbrok']>12:
            continue
        print(ev_t0[i],txtlines[i][0],ev_xyz[i])
        #tmp1,tmp2=rotate(data_E=xdata_tmp[:,:,0],data_N=xdata_tmp[:,:,1])
        #xdata_tmp[:,:,0]=tmp1
        #xdata_tmp[:,:,1]=tmp2
        xdata_tmp=check_and_norm_data(tmp['data'],stn=stxy)
        xdata.append(xdata_tmp['xdata'])
        stns1.append(xdata_tmp['stns'])
        srcxyz.append(ev_xyz[i])
        
    sio.savemat('samples_data_ridgecrest.mat',{'xdata':np.array(xdata),
                                          'srcxyz':np.array(srcxyz),
                                          'stns':np.array(stns1)})
    

def testdata_conti():
    
    with open('D:/loca_aug/plot/stations_ridge.txt') as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        stxy=np.array([[float(i[0]),float(i[1])] for i in txtlines])
        stnam=[i[3]+'.'+i[4] for i in txtlines]
        tmp1,tmp2=coordinate2model_stn(stnx=stxy[:,0],stny=stxy[:,1])
        stxy[:,0]=tmp1
        stxy[:,1]=tmp2
        stns=[{'X':stxy[i,0],'Y':stxy[i,1],'stnam':stnam[i]} for i in range(len(stnam))] #'E'->x
        stns=sort_stn(stns=stns)
        print(stns)
        stxy=[[i['X'],i['Y']] for i in stns]
        stnam=[i['stnam'] for i in stns]
        
    h5obj=h5py.File('D:/loca_aug/data/Ridgecrest/20190704170000.h5','r')
    
    xdata=[]
    srcxyz=[]
    stns1=[]
    prelist='2019-07-04T17:00:00.000000Z'.replace(':','.');
    saclist=[]
    for j in range(len(stns)):
        #saclist.append(prelist+stnam[j])
        saclist=saclist+[prelist+stnam[j]+'.HHE',
                         prelist+stnam[j]+'.HHN',
                         prelist+stnam[j]+'.HHZ']
    for i in range(500):
        ev_t0=UTCDateTime('2019-07-04T17:00:00.000000Z')+i*0.5
        tmp=h5py_getObspyStream_data(trcns=saclist,hf=h5obj,win=[ev_t0,30.0])
        if tmp['nbrok']>12:
            continue
        print(ev_t0)
        #tmp1,tmp2=rotate(data_E=xdata_tmp[:,:,0],data_N=xdata_tmp[:,:,1])
        #xdata_tmp[:,:,0]=tmp1
        #xdata_tmp[:,:,1]=tmp2
        xdata_tmp=check_and_norm_data(tmp['data'],stn=stxy)
        xdata.append(xdata_tmp['xdata'])
        stns1.append(xdata_tmp['stns'])
        
    sio.savemat('test_data_ridgecrest20190704170000.mat',{'xdata':np.array(xdata),
                                          'srcxyz':np.array(srcxyz),
                                          'stns':np.array(stns1)})
    
def check_and_norm_data(data,stns,stnr=[[0.0,82.0],[0.0,100.0]],nsmpout=1024,max_num_st=12):
    nr=max_num_st#len(data)
    nsmp=len(data[0])
    maxval_st=[max([max([abs(i[0]),abs(i[1]),abs(i[2])]) for i in j]) for j in data]
    maxval=max(maxval_st)
    dx=stnr[0][1]-stnr[0][0]
    dy=stnr[1][1]-stnr[1][0]
    data1=[]
    stns1=[]
    for idx, j in enumerate(data):
        if maxval_st[idx]<0.1E-23:
            continue
        stn=stns[idx]
        stnx=(stn['X']-stnr[0][0])/dx
        stny=(stn['Y']-stnr[1][0])/dy
        if stnx<0.0 or stnx>1.0 or stny<0.0 or stny>1.0:
            stnx=0.0;stny=0.0
            continue
        stns1.append(stn)
        data1.append([[float(i[0])/(maxval+0.01e-100),
                       float(i[1])/(maxval+0.01e-100),
                       float(i[2])/(maxval+0.01e-100), 
                       stnx,stny] for i in j]+[[0.0,0.0,0.0,0,0]]*(nsmpout-nsmp))
    data1=data1+[[[0.0,0.0,0.0,0.0,0.0]]*nsmpout]*(nr-len(data1))
    stns1=stns1+[{}]*(nr-len(stns1))
    return {'xdata':np.array(data1),'stns':stns1,'maxval':maxval,'maxval_st':maxval_st}
from math import cos,sin
import math
from distaz import DistAz
def rotate(data_E,data_N,rt_clockwise=0.0): #rt_clockwise is "data rotate clokwisly" equal to "coordinate rotate anticlockwisly" 
    rt_clockwise1=rt_clockwise*math.pi/180.0
    tmp1=data_E*cos(rt_clockwise1)+data_N*sin(rt_clockwise1)
    tmp2=-data_E*sin(rt_clockwise1)+data_N*cos(rt_clockwise1)
    return tmp1,tmp2

def coordinate2model_stn(stnx,stny,
                         data_org={'xy':[-10631,3980],'scale':[111.19*0.813,111.19],'rt_clockwise':0.0},
                         model_org=[82.0/2,100.0/2],
                         ):
    rt_clockwise=data_org['rt_clockwise']
    stnx1=stnx*data_org['scale'][0]-data_org['xy'][0]
    stny1=stny*data_org['scale'][1]-data_org['xy'][1]
    tmp1,tmp2=rotate(stnx1,stny1,rt_clockwise=rt_clockwise)
    tmp1=tmp1+model_org[0]
    tmp2=tmp2+model_org[1]
    return tmp1,tmp2
    

def sort_stn(stns,org=[0.0,0.0]): #X-->E
    for i in range(len(stns)):
        stn=stns[i]
        tmp=stn['X']-org[0]
        tmp1=stn['Y']-org[1]
        dist=(tmp*tmp+tmp1*tmp1)#**0.5+1.0/stn[0]+1.0/stn[1]
        for j in range(i,len(stns)):
            stn1=stns[j]
            tmp=stn1['X']-org[0]
            tmp1=stn1['Y']-org[1]
            dist_tmp=(tmp*tmp+tmp1*tmp1)#**0.5+1.0/stn1[0]+1.0/stn1[1]
            if dist>dist_tmp:
                tmp=stns[i]
                stns[i]=stns[j]
                stns[j]=tmp
                dist=dist_tmp
    return stns

def set_saclist(prelist,stnam,chn):
     saclist=[]
     for j in range(len(stnam)):
         stnam1=stnam[j].split('.')
         if len(stnam1)==2:
             saclist=saclist+[prelist+stnam[j]+'.'+chn[0],
                              prelist+stnam[j]+'.'+chn[1],
                              prelist+stnam[j]+'.'+chn[2]]
         else:
             saclist=saclist+[prelist+'.'.join(stnam1[0:2]+[stnam1[2]]),
                              prelist+'.'.join(stnam1[0:2]+[stnam1[3]]),
                              prelist+'.'.join(stnam1[0:2]+[stnam1[4]])]
         #saclist.append(prelist+stnam[j])
         
     return saclist

import json
class DataStream(): 
    def __init__(self,datajs_file='data_info.json',modeljs_file='model_info.json'):
        with open(datajs_file,'r') as f:
            datajs=json.load(f)
        with open(modeljs_file,'r') as f:
            modeljs=json.load(f)
        self.h5obj=h5py.File(datajs['h5file'],'r')
        with open(datajs['stations_file']) as f:
            txtlines=[i.strip().split() for i in f.readlines()];
            st_geo=np.array([[float(i[0]),float(i[1])] for i in txtlines])
            stnam=[i[3]+'.'+i[4] for i in txtlines]
            tmp1,tmp2=coordinate2model_stn(stnx=st_geo[:,0],stny=st_geo[:,1],
                                           data_org=datajs['data_org'],model_org=modeljs['model_org'])
            stns=[{'X':tmp1[i],'Y':tmp2[i],'geo':st_geo[i],'stnam':stnam[i]} for i in range(len(stnam))] #'E'->x
            stns=sort_stn(stns=stns)
            print(stns)
            stxy=[[i['X'],i['Y']] for i in stns]
            st_geo=[[i['geo'][0],i['geo'][1]] for i in stns]
            stnam=[i['stnam'] for i in stns]
        self.stnam=stnam
        self.stxy=stxy
        self.st_geo=st_geo
        self.stns=stns
        self.win_len=modeljs['twin_len']
        self.modeljs=modeljs
        self.datajs=datajs
        self.set_datastream()
        
    def set_datastream(self):
        prelist=self.datajs['prelist'];
        if type(prelist)!=str:
            return None
        datajs=self.datajs
        sacnams=set_saclist(prelist=prelist, stnam=self.stnam, chn=[datajs['E'],datajs['N'],datajs['Z']])
        self.saclist=h5py_getObspyStream(trcns=sacnams,hf=self.h5obj)           
        
    def get_waveform_events(self,datatype='DetecLoca'):
        stns1=[]
        xdata=[]
        maxvals=[]
        if datatype=='DetecLoca':
            prep={'freq':self.modeljs['filter'],'nsmpout':self.modeljs['model_input_size'][1]}
        if datatype=='Mag':
            prep={'freq':self.modeljs['filter_mag'],'nsmpout':self.modeljs['model_input_size'][1]}
            
        datajs=self.datajs
        rt_clockwise=self.datajs['data_org']['rt_clockwise']
        prelist=datajs['prelist']
        tbegin=datajs['tbegin']
        chn=[datajs['E'],datajs['N'],datajs['Z']]
        stnam=self.stnam
        tbegins=[]
        for i in range(len(datajs['prelist'])):
            saclist=set_saclist(prelist[i], stnam, chn)
            ev_t0=UTCDateTime(tbegin[i]);
            tmp=h5py_getObspyStream_data(trcns=saclist,hf=self.h5obj,win=[ev_t0,self.win_len],prep=prep)
            if tmp['nbrok']>12:
                continue
            print(ev_t0)
            tmp['data']=np.array(tmp['data'])
            tmp1,tmp2=rotate(data_E=tmp['data'][:,:,0],data_N=tmp['data'][:,:,1],rt_clockwise=rt_clockwise)
            tmp['data'][:,:,0]=tmp1
            tmp['data'][:,:,1]=tmp2
            xdata_tmp=check_and_norm_data(tmp['data'],stns=self.stns)
            if datatype=='DetecLoca':
                xdata.append(xdata_tmp['xdata'])
            elif datatype=='Mag':
                xdata.append(xdata_tmp['xdata'][0:])
            stns1.append(xdata_tmp['stns'])
            maxvals.append(xdata_tmp['maxval']*self.datajs['ampscale'])
            tbegins.append(ev_t0)
        self.windata={'xdata':np.array(xdata),'stns':stns1,'maxvals':maxvals,'tbegins':tbegins}
        return self.windata
        
    def get_data(self):
        datajs=self.datajs
        if type(datajs['prelist'])==str:
            windata=self.get_waveform_conti(tbegin=UTCDateTime(datajs['tbegin']),
                                           nwin=datajs['nwin'],
                                           dt=datajs['dtwin'],datatype='DetecLoca')
            windata_mag=self.get_waveform_conti(tbegin=UTCDateTime(datajs['tbegin']),
                                           nwin=datajs['nwin'],
                                           dt=datajs['dtwin'],datatype='Mag')
            return {'windata':windata,'windata_mag':windata_mag}
        else:
            windata=self.get_waveform_events(datatype='DetecLoca')
            windata_mag=self.get_waveform_events(datatype='Mag')
            return {'windata':windata,'windata_mag':windata_mag}
   
        
    def get_waveform_conti(self,tbegin,nwin,dt=1.0, datatype='DetecLoca'):
        stns1=[]
        xdata=[]
        maxvals=[]
        if datatype=='DetecLoca':
            prep={'freq':self.modeljs['filter'],'nsmpout':self.modeljs['model_input_size'][1]}
        if datatype=='Mag':
            prep={'freq':self.modeljs['filter_mag'],'nsmpout':self.modeljs['model_input_size'][1]}
        if type(tbegin)==str:
            tbegin=UTCDateTime(tbegin)
        tbegins=[]
        rt_clockwise=self.datajs['data_org']['rt_clockwise']
        for i in range(nwin):
            ev_t0=tbegin+i*dt;
            tmp=h5py_getObspyStream_data(trcns=self.saclist,hf=self.h5obj,win=[ev_t0,self.win_len],prep=prep)
            if tmp['nbrok']>12:
                continue
            #print(ev_t0)
            tmp['data']=np.array(tmp['data'])
            tmp1,tmp2=rotate(data_E=tmp['data'][:,:,0],data_N=tmp['data'][:,:,1],rt_clockwise=rt_clockwise)
            tmp['data'][:,:,0]=tmp1
            tmp['data'][:,:,1]=tmp2
            xdata_tmp=check_and_norm_data(tmp['data'],stns=self.stns)
            if datatype=='DetecLoca':
                xdata.append(xdata_tmp['xdata'])
            elif datatype=='Mag':
                xdata.append(xdata_tmp['xdata'][0:])
            stns1.append(xdata_tmp['stns'])
            maxvals.append(xdata_tmp['maxval']*self.datajs['ampscale'])
            tbegins.append(ev_t0)
        self.windata={'xdata':np.array(xdata),'stns':stns1,'maxvals':maxvals,'tbegins':tbegins}
        return self.windata
    
    def save_waveform_conti(self,file_name):
        sio.savemat(file_name, self.windata)
    
#####################################################################
import numpy as np
def img2xyz(xr=[0.25,0.01,24],yr=[-0.2,0.013,32],zr=[3.07,0.01,18],imgs=[]):
     xyz=[]
     for i in range(0,len(imgs)):
        mvalue=np.amax(imgs[i])
        idx=np.where(imgs[i]==mvalue)
        xyz=xyz+[[xr[0]+xr[1]*idx[0][0],yr[0]+yr[1]*idx[1][0],zr[0]+zr[1]*idx[2][0],mvalue,i,idx[0][0],idx[1][0],idx[2][0]]]
     return xyz;
      
def output_result1(r=[],imgs=[],namout='test_xyz.txt'):
     xyz=img2xyz(xr=r[0],yr=r[1],zr=r[2],imgs=imgs)
     with open(namout,'w') as f:
         for i in xyz:
             f.write("{:.6f} {:.6f} {:.6f} {:.6f}\n".format(i[0],i[1],i[2],i[3]))
     return xyz

from keras.models import *
import scipy.io as sio
#from dataio import load_xdata_norm,lonlat2xy
def lonlat2xy_events(lonlats):
    xys=[]
    for tmp in lonlats:
        xy=lonlat2xy(lonlat=tmp)
        xys.append(xy)
    return np.array(xys)

def predict():
    print('load testing samples.')
    r=[[-2667.0,0.5,96],[-4145.0,0.5,192],[0.0,0.5,64],25,(10.0/10.0)**2.0]
    r1=[[r[0][0],r[0][0]+r[0][1]*r[0][2]],[r[1][0],r[1][0]+r[1][1]*r[1][2]],[0.0,15.0]]
    print('r,r1:',r,r1)
    #wave_test,loca_true=sgydata.load_sgylist_xyz1(shuffle='false',sgylist=['./waveform_data/','testing_samples.txt'],
    #                            sgyr=[0,-1,1],xr=r[0],yr=r[1],zr=r[2],r=0.05,shiftdata=[list(range(20,50))+list(range(-200,-20)),0])
    #loca_true=np.reshape(loca_true,(len(loca_true),80,128,30)) 
    
    print('save true labels (true location images).')
    #sio.savemat('test_truelabel.mat', {'loca_true':loca_true})
    #xdata,srcxyz=load_xdata_norm(pklfile='../samples_data_truewave.pkl'); srcxyz=lonlat2xy_events(srcxyz)
    xdata=sio.loadmat('tmp.mat')['xdata'];print(xdata);srcxyz=[0]
    wave_test=xdata[0:40]
    
    print('load trained network model.')
    model=load_model('./FCNloca.hdf5')
#    model=load_model('D:/cnnloca/shift_filter/model20_50_-200_-20a1/unet.hdf5')

    
    print('location prediction.')
    loca_predict = model.predict(wave_test, batch_size=1, verbose=1)
    print('output location results.')
    xyz_predict=output_result1(r=r,imgs=loca_predict,namout='test_xyz.txt')
    
    print('save predicted location images.')
    sio.savemat('test_predictedlabel.mat', {'loca_predict':loca_predict,'xyzrange':r[0:3],'wave_test':wave_test,
                                            'srcxyz':srcxyz,'xyz_predict':xyz_predict})

    print('end predict')

if __name__ == '__main__':
    #predict()
    #testdata()
    #testdata_conti()
    datastream=DataStream()
    datastream.get_waveform_conti('2019-07-04T17:00:00.000000Z', 10)
    datastream.save_waveform_conti('test_data_ridgecrest20190704170000.mat')

