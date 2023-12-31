# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 13:01:27 2023

@author: huawei
"""

from distaz import DistAz
import scipy.io as sio
import numpy as np
import random
import pickle
from math import sin, cos,pi
from obspy import read,UTCDateTime
from obspy import Stream,Trace
import obspy
from math import log,atan2
import h5py
from obspy.taup import TauPyModel

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
    paz_wa = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
                'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
    wa_amp=[]
    for i in range(0,len(saclist)):
        try:
            st_tmp=read(saclist[i]);
            st_tmp.filter('bandpass',freqmin=freq[0],freqmax=freq[1])
            st_tmp_mag=st_tmp.copy().slice(starttime=win[0]+30.0,endtime=win[0]+30.0+35.0,keep_empty_traces=True);
            st_tmp_mag.simulate(paz_remove = None, paz_simulate = paz_wa)
            st_tmp_mag_amp=np.max(np.abs(st_tmp_mag[0].data))
            wa_amp.append(st_tmp_mag_amp)
            st_tmp.interpolate(sampling_rate=20)
            st_tmp=st_tmp.slice(starttime=win[0],endtime=win[0]+win[1],keep_empty_traces=True);
            st_tmp.detrend("demean")
            st_tmp.detrend("linear")
            st_tmp.taper(max_percentage=0.00001)
            #print(st_tmp[0].stats)
        except:
            st_tmp=Stream(Trace(np.array([0,]*100)))
            wa_amp.append(0.0)
            print('no data in',saclist[i]);
        if i==0:
            st=st_tmp;
        else:
            st=st+st_tmp;
    data1event=stream_to_list(datain=st,nsmp=nsmpout)
    data1event=np.clip(np.nan_to_num(np.array(data1event)),-1.0e10,1.0e10).tolist()
    print('st',len(st))
    return {'data':data1event,'wa_amp':[[wa_amp[i],wa_amp[i+1],wa_amp[i+2]] for i in range(0,len(wa_amp),3)]}

def stream_to_list_(datain=Stream(),nsmp=2048,starttime=-9999.0):
    data1=[];
    info=[]
    for i in range(0,len(datain)):
        st=datain[i];
        info.append({'cmpaz':st.stats['sac']['cmpaz'],'cmpinc':st.stats['sac']['cmpinc']})
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
    info=[[info[ir],info[ir+1],info[ir+2]] for ir in range(0,len(info),3)]
    return {'data':data,'info':info};
def load_sac_rt(saclist=['D:/Italy/201610.20to24continue_sac/2016-08-20T00.00.00.000000ZIV.ASSB.HHE.SAC','HHN.sac','HHZ.sac'],
             win=[],
             freq=[2.0,8.0],
             nsmpout=1024):
    paz_wa = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
                'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
    wa_amp=[]
    for i in range(0,len(saclist)):
        try:
            st_tmp=read(saclist[i]);
            st_tmp.differentiate()
            st_tmp.filter('bandpass',freqmin=freq[0],freqmax=freq[1])
            st_tmp1=st_tmp.copy()#.filter('bandpass',freqmin=freq[0],freqmax=freq[1])
            st_tmp_mag=st_tmp1.slice(starttime=win[0]+30.0,endtime=win[0]+30.0+35.0,keep_empty_traces=True);
            st_tmp_mag.simulate(paz_remove = None, paz_simulate = paz_wa)
            st_tmp_mag_amp=np.max(np.abs(st_tmp_mag[0].data))
            wa_amp.append(st_tmp_mag_amp)
            st_tmp.interpolate(sampling_rate=20)
            st_tmp=st_tmp.slice(starttime=win[0],endtime=win[0]+win[1],keep_empty_traces=True);
            st_tmp.detrend("demean")
            st_tmp.detrend("linear")
            st_tmp.taper(max_percentage=0.00001)
            #print(st_tmp[0].stats)
        except:
            st_tmp=Stream(Trace(np.array([0,]*100)))
            wa_amp.append(0.0)
            print('no data in',saclist[i]);
        if i==0:
            st=st_tmp;
        else:
            st=st+st_tmp;
    tmp=stream_to_list_(datain=st,nsmp=nsmpout)
    data1event=tmp['data']
    info=tmp['info']
    data1event=np.clip(np.nan_to_num(np.array(data1event)),-1.0e10,1.0e10)
    for i in range(0,len(data1event)):
        tmpE=data1event[i,:,0]*sin(info[i][0]['cmpaz'])+data1event[i,:,1]*sin(info[i][1]['cmpaz'])
        tmpN=data1event[i,:,0]*cos(info[i][0]['cmpaz'])+data1event[i,:,1]*cos(info[i][1]['cmpaz'])
        data1event[i,:,0]=tmpE
        data1event[i,:,1]=tmpN
    print('st',len(st))
    return {'data':data1event.tolist(),'wa_amp':[[wa_amp[i],wa_amp[i+1],wa_amp[i+2]] for i in range(0,len(wa_amp),3)]}

def value2key(src,srcr=[[0.0,0.5],[0.0,0.5]]):
    szr=srcr[1]
    srr=srcr[0]
    iz=int((src[1]-szr[0])/szr[1]+0.5)
    ir=int((src[0]-srr[0])/srr[1]+0.5)
    key='r{:04d}z{:04d}'.format(ir,iz)
    return key

def check_h5_data(h5obj,nam):
    dataset=h5py_get(nam,h5obj)
    #print(dataset,nam)
    if not ('None' in dataset['channel']):
        
        data=dataset['data']
        dmax0=np.array(data[0]).max()
        dmax1=np.array(data[1]).max()
        dmax2=np.array(data[2]).max()
        dlen0=len(data[0])
        dlen1=len(data[1])
        dlen2=len(data[2])
        if dmax0>0.1e-23 and dmax1>0.1e-23 and dmax2>0.1e-23 and dlen0>200 and dlen1>200 and dlen2>200:
            return {'flg':0,'dataset':dataset}
    return {'flg':-9999,'dataset':dataset}

def h5py_get(trcn,hf):
    chn=[]
    data=[]
    chns1=['.HHE','.BHE','.HNE','.BNE']
    chns2=['.HHN','.BHN','.HNN','.BNN']
    chns3=['.HHZ','.BHZ','.HNZ','.BNZ']
    for i,chn1 in enumerate(chns1):
        tmpdata=hf.get(trcn+chn1)
        #print(trcn+chn1)
        if tmpdata!=None:
            chn.append(chn1)
            data.append(tmpdata)
            break
    if chn==[]:
        chn.append('None')
    for i,chn2 in enumerate(chns2):
        tmpdata=hf.get(trcn+chn2)
        if tmpdata!=None:
            chn.append(chn2)
            data.append(tmpdata)
            break
    if len(chn)<2:
        chn.append('None')
    for i,chn3 in enumerate(chns3):
        tmpdata=hf.get(trcn+chn3)
        if tmpdata!=None:
            chn.append(chn3)
            data.append(tmpdata)
            break
    if len(chn)<3:
        chn.append('None')
    return {'channel':chn,'data':data}


def dist_az(staeast_km,stanorth_km,eveast_km,evnorth_km):
    east_km=staeast_km-eveast_km
    north_km=stanorth_km-evnorth_km
    dist=east_km*east_km+north_km*north_km;dist=dist**0.5
    az=atan2(east_km,north_km)*180/pi
    if az<0:
        az=360+az
    return {'dist':dist,'az':az}
import copy
def gen_1sample(st_samps,stns,src):
        sample=[]
        num=0
        for i,stn in enumerate(stns):
            #print(stn,i)
            #tmp=DistAz(stalon=stn[0],stalat=stn[1],evtlon=src[0],evtlat=src[1])
            #dist=tmp.getDistanceKm()
            #az=tmp.getAz()
            tmp=dist_az(staeast_km=stn[0],stanorth_km=stn[1],eveast_km=src[0],evnorth_km=src[1])
            dist=tmp['dist']
            az=tmp['az']
            srcz=src[2]
            sk=value2key(src=[dist,srcz])
            #print(sk,dist,az)
            try:
                #print(st_samps[sk])
                #st_sam=random.sample(st_samps[sk],1)[0]#.deepcopy()
                st_sam=copy.deepcopy(random.sample(st_samps[sk],1)[0])
                #print(az,st_sam['azimuth'])
                st_sam['daz']=az-st_sam['azimuth']
                #print(st_sam['daz'])
                num=num+1
            except:
                st_sam={}
            sample.append(st_sam)
       # print(src,sample,num)
        return {'srcxyz':src,'samp':sample,'numstn':num}
    
def sort_stn(stns,org=[0.0,0.0]):
    for i in range(len(stns)):
        stn=stns[i]
        tmp=stn[0]-org[0]
        tmp1=stn[1]-org[1]
        dist=(tmp*tmp+tmp1*tmp1)#**0.5+1.0/stn[0]+1.0/stn[1]
        for j in range(i,len(stns)):
            stn1=stns[j]
            tmp=stn1[0]-org[0]
            tmp1=stn1[1]-org[1]
            dist_tmp=(tmp*tmp+tmp1*tmp1)#**0.5+1.0/stn1[0]+1.0/stn1[1]
            if dist>dist_tmp:
                tmp=stns[i]
                stns[i]=stns[j]
                stns[j]=tmp
                dist=dist_tmp
    return stns

def merge_trace_info(nam0,nam1,outnam):
    with open(nam0,'rb') as f:
        st_samps1=pickle.load(f)
    with open(nam1,'rb') as f:
        st_samps0=pickle.load(f)
    for key in st_samps1:
        if key in st_samps0:
            st_samps0[key]=st_samps0[key]+copy.deepcopy(st_samps1[key])#+copy.deepcopy(st_samps1[key]) \
        else:
            st_samps0[key]=copy.deepcopy(st_samps1[key])#+copy.deepcopy(st_samps1[key]) \
    with open(outnam,'wb') as f:
        pickle.dump(st_samps0,f)

def cata2stsamp(cata,stns):
    from pytool import cal_taupPtraveltime
    ttmodel=TauPyModel('prem')
    st_samps={}
    info=[]
    for i,ev in enumerate(cata):
        for j,st in enumerate(stns):
            srcxy=[ev['geo'][0],ev['geo'][1],ev['dep']]
            rxy=[st['geo'][0],st['geo'][1]]
            ctime=str(UTCDateTime(ev['t0']))
            #print(srcxy,rxy)
            #tmp=client.distaz(stalon=rxy[0],stalat=rxy[1],evtlon=srcxy[0],evtlat=srcxy[1])
            tmp=DistAz(stalon=rxy[0],stalat=rxy[1],evtlon=srcxy[0],evtlat=srcxy[1])
            dist=tmp.getDistanceKm()
            srcz=srcxy[2]
            pstime=cal_taupPtraveltime(tauPmodel=ttmodel,source_depth_in_km=srcz,distance_in_degree=tmp.getDelta())
            info.append([dist,srcz])
            skey=value2key(src=[dist,srcz],srcr=[[0.0,0.5],[0.0,0.5]])
            if skey in st_samps:
                st_samps[skey].append({'time':ctime,'distance':dist,'ptime':pstime['ptime'],'stime':pstime['stime'],
                                      'azimuth':tmp.getAz(),'stnam':st['stnam']+'.'+st['channel'][0]+'.'+
                                      st['channel'][1]+'.'+st['channel'][2],
                                      'srcz':srcz,'stz':st['dep'],'mag':ev['mag']})
            else:
                st_samps[skey]=[{'time':ctime,'distance':dist,'ptime':pstime['ptime'],'stime':pstime['stime'],
                                      'azimuth':tmp.getAz(),'stnam':st['stnam']+'.'+st['channel'][0]+'.'+
                                      st['channel'][1]+'.'+st['channel'][2],
                                      'srcz':srcz,'stz':st['dep'],'mag':ev['mag']}]
    return st_samps

def get_waamp(stream,win=['',30],prep={'freq':[1,9]}):
    paz_wa = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
                'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
    wa_amp=[]
    freq=prep['freq']
    for i,st_tmp in enumerate(stream):
        try:
            st_tmp=Stream(st_tmp.copy())
            st_tmp.filter('bandpass',freqmin=freq[0],freqmax=freq[1])
            #print(type(st_tmp),st_tmp,win)
            st_tmp_mag=st_tmp.slice(starttime=win[0],endtime=win[0]+win[1],keep_empty_traces=True);
            st_tmp_mag.simulate(paz_remove = None, paz_simulate = paz_wa)
            st_tmp_mag_amp=np.max(np.abs(st_tmp_mag[0].data))
            wa_amp.append(st_tmp_mag_amp)
        except:
            st_tmp=Stream(Trace(np.array([0,]*100)))
            wa_amp.append(0.0)
            #print('no data in',st_tmp);
    return {'wa_amp':[[wa_amp[i],wa_amp[i+1],wa_amp[i+2]] for i in range(0,len(wa_amp),3)]}

def get_data_from_stream(stream,win=['',30],nsmpout=1024,prep={'freq':[],'samp_rate':20}):
    freq=prep['freq']
    samp_rate=prep['samp_rate']
    nbrok=0
    for i,st_tmp in enumerate(stream):
        try:
            st_tmp=Stream(st_tmp.copy())
            if freq!=[]:
                st_tmp.filter('bandpass',freqmin=freq[0],freqmax=freq[1])
            st_tmp.interpolate(sampling_rate=samp_rate)
            #print(type(st_tmp),st_tmp,win)
            st_tmp=st_tmp.slice(starttime=win[0],endtime=win[0]+win[1],keep_empty_traces=True);
            st_tmp.detrend("demean")
            st_tmp.detrend("linear")
            st_tmp.taper(max_percentage=0.00001)
        except Exception as e:
            print('no data in',e,i);
            st_tmp=Stream(Trace(np.array([0,]*100)))
            nbrok=nbrok+1
        if i==0:
            st=st_tmp;
        else:
            st=st+st_tmp;
    data1event=stream_to_list(datain=st,nsmp=nsmpout)
    data1event=np.clip(np.nan_to_num(np.array(data1event)),-1.0e10,1.0e10).tolist()
    return {'data':data1event,'nbrok':nbrok}

def st_samples_info(nam):
    import pickle
    with open(nam,'rb') as f:
        st_samps=pickle.load(f)
    num_st=0
    keys=st_samps.keys()
    for i,key in enumerate(keys):
        num_st=num_st+len(st_samps[key])
    print('number of seismograms:',num_st)
 
def st_samples_extr_trace(nam,dist,z):
    import pickle
    with open(nam,'rb') as f:
        st_samps=pickle.load(f)
    key=value2key(src=[dist,z])
    return st_samps[key]
    
def h5py_load_dataset(h5file,rg=[0,-1,1]):
    h5obj=h5py.File(h5file,'r')
    data=[]
    for i,key in enumerate(list(h5obj.keys())[rg[0]:rg[1]:rg[2]]):
        dataset=h5obj.get(key)
        data.append(np.array(dataset))
    return {'data':np.array(data)}

        
    
















