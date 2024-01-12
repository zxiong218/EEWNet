# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 22:06:28 2018

@author: Dell
"""
from obspy import Stream,Trace,read,UTCDateTime
import h5py
import obspy
import scipy.io as sio
from obspy.taup import TauPyModel
from data_stream import DataStream,rotate
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
        st_tmp=h5py_getObspyStream(trcns=[trcn],hf=hf)
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
        
    
def check_and_norm_data(data,stn,stnr=[[0.0,82.0],[0.0,100.0]],nsmpout=1024,max_num_st=12):
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
        stnx=(stn[idx][0]-stnr[0][0])/dx
        stny=(stn[idx][1]-stnr[1][0])/dy
        if stnx<0.0 or stnx>1.0 or stny<0.0 or stny>1.0:
            stnx=0.0;stny=0.0
            continue
        stns1.append([stn[idx][0],stn[idx][1]])
        data1.append([[float(i[0])/(maxval+0.01e-100),
                       float(i[1])/(maxval+0.01e-100),
                       float(i[2])/(maxval+0.01e-100), 
                       stnx,stny] for i in j]+[[0.0,0.0,0.0,stnx,stny]]*(nsmpout-nsmp))
    data1=data1+[[[0.0,0.0,0.0,0.0,0.0]]*nsmpout]*(nr-len(data1))
    return {'xdata':np.array(data1),'stnxy':np.array(stns1)}
from math import cos,sin
import math
from distaz import DistAz
def coordinate2model_stn(stnx,stny,data_org={'org':[-10680,4005],'x':111.19*0.816,'y':111.19},
                         model_org=[82.0/2,100.0/2],
                         rt_clockwise=0.0):
    stnx1=stnx*data_org['x']-data_org['org'][0]
    stny1=stny*data_org['y']-data_org['org'][1]
    tmp1,tmp2=rotate(stnx1,stny1,rt_clockwise=rt_clockwise)
    tmp1=tmp1+model_org[0]
    tmp2=tmp2+model_org[1]
    return tmp1,tmp2

def coordinate2model(lon,lat,center={'lon':0,'lat':0},rt_clockwise=0.0):
    tmp=DistAz(stalon=center['lon'],stalat=center['lat'],evtlon=center['lon']+1,evtlat=center['lat'])
    scale=tmp.getDistanceKm()
    data_org={'org':[center['lon']*scale,center['lat']*111.19],'x':scale,'y':111.19}
    tmp1,tmp2=coordinate2model_stn(lon,lat,data_org=data_org,rt_clockwise=rt_clockwise)
    return tmp1,tmp2

def coordinate2geo(ex,ey,center={'lon':0,'lat':0},rt_clockwise_in2model=0.0,model_org=[82.0/2,100.0/2]):
    tmp=DistAz(stalon=center['lon'],stalat=center['lat'],evtlon=center['lon']+1,evtlat=center['lat'])
    scale=tmp.getDistanceKm()
    data_org={'org':[center['lon']*scale,center['lat']*111.19],'x':scale,'y':111.19}
    sx=ex-model_org[0]
    sy=ey-model_org[1]
    sx,sy=rotate(sx,sy,rt_clockwise=-rt_clockwise_in2model)
    sx=sx+data_org['org'][0]
    sy=sy+data_org['org'][1]
    sx=sx/data_org['x']
    sy=sy/data_org['y']
    return sx,sy

def coordinate2geo1(ex,ey,coordinate=
                    {'data_org':{'lonlat':[0,0],'xy':[0*111.19*0.813,0*111.19],
                                                  'scale':[111.19*0.813,111.19],'rt_clockwise':0.0},
                                      'model_org':[82.0/2,100.0/2]}):
    data_org=coordinate['data_org']
    model_org=coordinate['model_org']
    sx=ex-model_org[0]
    sy=ey-model_org[1]
    sx,sy=rotate(sx,sy,rt_clockwise=-data_org['rt_clockwise'])
    sx=sx+data_org['xy'][0]
    sy=sy+data_org['xy'][1]
    sx=sx/data_org['scale'][0]
    sy=sy/data_org['scale'][1]
    return sx,sy

def coordinate2geo_test():
    data=sio.loadmat('D:/loca_aug/loca_net/test/test_predictedlabel.mat');
    print(data['srcxyz'][0,0],data['xyz_predict'])
    #sx,sy=coordinate2geo(data['srcxyz'][:,0],data['srcxyz'][:,1],center={'lon':13.25,'lat':42.75},rt_clockwise_in2model=160)
    #sx_p,sy_p=coordinate2geo(data['xyz_predict'][:,0],data['xyz_predict'][:,1],center={'lon':13.25,'lat':42.75},rt_clockwise_in2model=160)
    sx,sy=coordinate2geo(data['srcxyz'][:,0],data['srcxyz'][:,1],center={'lon':-117.7105,'lat':36.0194},rt_clockwise_in2model=0)
    sx_p,sy_p=coordinate2geo(data['xyz_predict'][:,0],data['xyz_predict'][:,1],center={'lon':-117.7105,'lat':36.0194},rt_clockwise_in2model=0)

    print(sx,sy)
    import matplotlib.pyplot as plt
    plt.plot(sx,sy,'*')
    plt.plot(sx_p,sy_p,'*r')
    with open('test.txt','w') as f:
        for i in range(len(sx_p)):
            tmp1=data['xyz_predict'][i,0]-data['srcxyz'][i,0];tmp2=data['xyz_predict'][i,1]-data['srcxyz'][i,1];disterr=(tmp1*tmp1+tmp2*tmp2)**0.5;
            f.write("{} {} {} {} {} {} {} {}\n".format(sx_p[i],sy_p[i],data['xyz_predict'][i,2],data['xyz_predict'][i,3],
                    sx[i],sy[i],data['srcxyz'][i,2],disterr))
    
    

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
    
#####################################################################
#####################################################################
import numpy as np
def img2xyz(xr=[0.25,0.01,24],yr=[-0.2,0.013,32],zr=[3.07,0.01,18],imgs=[]):
     xyz=[]
     for i in range(0,len(imgs)):
        mvalue=np.amax(imgs[i])
        idx=np.where(imgs[i]==mvalue)
        #xyz=xyz+[[xr[0]+xr[1]*idx[0][0],yr[0]+yr[1]*idx[1][0],zr[0]+zr[1]*idx[2][0],mvalue,i,idx[0][0],idx[1][0],idx[2][0]]]
        xyz.append({'xyz':[xr[0]+xr[1]*idx[0][0],yr[0]+yr[1]*idx[1][0],zr[0]+zr[1]*idx[2][0],mvalue],
                    'others':[i,idx[0][0],idx[1][0],idx[2][0]]})
     return xyz;
      
def output_result1(r=[],imgs=[],namout='test_xyz.txt'):
     xyz=img2xyz(xr=r[0],yr=r[1],zr=r[2],imgs=imgs)
    # with open(namout,'w') as f:
    #     for i in xyz:
    #         xyz1=i['xyz']
    #         f.write("{:.6f} {:.6f} {:.6f} {:.6f}\n".format(xyz1[0],xyz1[1],xyz1[2],xyz1[3]))
     return xyz
 
def result_loca(r,imgs,coordinate):
    re=output_result1(r=r,imgs=imgs)
    for tmp in re:
        xyz=tmp['xyz']
        lon,lat=coordinate2geo1(xyz[0],xyz[1],coordinate=coordinate)
        tmp['geo']=[lon,lat]
    return re
 
def img2detec(r=[],imgs=[]):
    detec=[]
    for i in range(len(imgs)):
        mvalue=np.amax(imgs[i])
        idx=np.where(imgs[i]==mvalue)[1][0]
        detec=detec+[[idx,mvalue,i]]
    return detec

def result_detect(r=[],imgs=[],namout='test_detec.txt'):
    detec=img2detec(r=r,imgs=imgs)
    #with open(namout,'w') as f:
    #    for i in detec:
    #        f.write('{} {:.6f}\n'.format(i[0],i[1]))
    return detec

 
def img2mag(img=[],mag_range=[],amp=-1.0):
    img1=img[:,0]
    mrt=np.max(img1)
    idx=np.where(mrt==img1)[0][0]
    #print(mrt,idx)
    mag=idx*mag_range[1]+mag_range[0]+amp*1.0
    return [mag,mrt]
def result_mag(imgs=[],mag_range=[],amps=[]):
    num=len(imgs)
    mag=[]
    for i in range(num):
        tmp=img2mag(img=imgs[i],mag_range=mag_range,amp=np.math.log10(amps[i]+0.1E-30))
        #print(tmp)
        mag.append(tmp)
   # with open('test_mag.txt','w') as f:
   #     for i in range(num):
   #         f.write('{} {} {} \n'.format(mag[i][0],mag[i][1],np.math.log10(amps[i]+0.1E-30)))
    return mag

def mag_trcs2event(mags):
    mag_valid=[]
    for i in range(len(mags)):
        if mags[i][1]>0.6:
            mag_valid.append(mags[i][0])
    if mag_valid==[]:
        mag_valid=-9999
    return {'mag_valid':mag_valid,'mag_mean':np.mean(mag_valid),'mag_median':np.median(mag_valid),'mags':mags}
        
def mag_data1event(event_data,dists,dist_range=[0.0,110.0]):
    data=[]
    amps=[]
    for i in range(len(event_data)):
        tmp=np.max(np.abs(event_data[i,:,0:3]))
        tmp1=[np.max(np.abs(event_data[i,:,0])),np.max(np.abs(event_data[i,:,1]))]
        if tmp<0.1E-30:
            continue
        amps.append((tmp1[0]+tmp1[1])/2.0)
        data.append([[j[0]/tmp,j[1]/tmp,j[2]/tmp,(dists[i]-dist_range[0])/dist_range[1]] for j in event_data[i]])
    return {'data':np.array(data),'amps':np.array(amps)}

from keras.models import load_model
#from dataio import load_xdata_norm,lonlat2xy
def lonlat2xy_events(lonlats):
    xys=[]
    for tmp in lonlats:
        xy=lonlat2xy(lonlat=tmp)
        xys.append(xy)
    return np.array(xys)
from operator import itemgetter
def sort_stn_xy(stn=[[0,1],[8,9]]):
    stnx=[(i,stn[i][0]) for i in range(len(stn))]
    stny=[(i,stn[i][1]) for i in range(len(stn))]
    stnx=sorted(stnx,key=itemgetter(1))
    stny=sorted(stny,key=itemgetter(1))
    return {'stnx':stnx,'stny':stny}
def sort_data(data,stn):
    #stnxy=[[i['X'],i['Y']] for i in stn if i != {}]+[[0.0,0.0]]*len(stn)
    tmp=sort_stn_xy(stn=stn)
    stnx=tmp['stnx']
    stny=tmp['stny']
    data1=[]
    for i in range(len(stnx)):
        idx=stnx[i][0]
        idy=stny[i][0]
        #data1.append(data[idx]+data[idy])
        data1.append([data[idx][j]+data[idy][j] for j in range(len(data[idx]))])
    return data1

def geo2dist(src_geo,stns_geo):
    dists=[]
    deltas=[]
    for stn in stns_geo:
        tmp=DistAz(stalat=stn[1], stalon=stn[0], evtlat=src_geo[1], evtlon=src_geo[0])
        dists.append(tmp.getDistanceKm())
        deltas.append(tmp.getDelta())
    return {'dists':dists,'deltas':deltas}
    


def predict():
    print('load testing samples.')
    #r=[[16.0,0.5208,96],[0.0,0.5208,192],[0.0,0.35,64],(5.5)**2.0,(5.0/5*10.0/10.0)**2.0]
    r=[[16.0,0.5208,96],[0.0,0.5208,192],[-6.0,0.45,64],(5.5)**2.0,(5.0/5*10.0/10.0)**2.0]
    r1=[[r[0][0],r[0][0]+r[0][1]*r[0][2]],[r[1][0],r[1][0]+r[1][1]*r[1][2]],[0.0,15.0]]
    print('r,r1:',r,r1)
    #wave_test,loca_true=sgydata.load_sgylist_xyz1(shuffle='false',sgylist=['./waveform_data/','testing_samples.txt'],
    #                            sgyr=[0,-1,1],xr=r[0],yr=r[1],zr=r[2],r=0.05,shiftdata=[list(range(20,50))+list(range(-200,-20)),0])
    #loca_true=np.reshape(loca_true,(len(loca_true),80,128,30)) 
    
    print('save true labels (true location images).')
    #sio.savemat('test_truelabel.mat', {'loca_true':loca_true})
    #xdata,srcxyz=load_xdata_norm(pklfile='../samples_data_truewave.pkl'); srcxyz=lonlat2xy_events(srcxyz)
    data=sio.loadmat('samples_data_japan.mat');print('ss',data['xdata'][0].shape,data['stns'][0][0].shape);xdata=[sort_data(data['xdata'][i].tolist(),data['stns'][0][i].tolist()+[[0,0]]*(12-len(data['stns'][0][i]))) for i in range(len(data['xdata']))];srcxyz=data['srcxyz']
    #data=sio.loadmat('samples_data_ridgecrest.mat');xdata=data['xdata'];srcxyz=data['srcxyz']
    wave_test=np.array(xdata[0:500])
    
    print('load trained network model.',wave_test)
    #model=load_model('../tmp2/checkpoint/bak02-0.04.hdf5')
    #model=load_model('../loca/train_dep-6_good/09-0.04.hdf5')
    model=load_model('../loca/tmp/09-0.04.hdf5')
    #model=load_model('../tmp2/FCNloca_epo7.hdf5')
#    model=load_model('D:/cnnloca/shift_filter/model20_50_-200_-20a1/unet.hdf5')

    
    print('location prediction.')
    loca_predict = model.predict(wave_test, batch_size=1, verbose=1)
    print('output location results.')
    xyz_predict=output_result1(r=r,imgs=loca_predict,namout='test_xyz.txt')
    
    print('save predicted location images.')
    sio.savemat('test_predictedlabel.mat', {'loca_predict':loca_predict[0:10],'xyzrange':r[0:3],'wave_test':wave_test[0:10],
                                            'srcxyz':srcxyz,'xyz_predict':xyz_predict,'stns':data['stns']})

    print('end predict')
    
def cal_taupPtraveltime(tauPmodel,source_depth_in_km,distance_in_degree):
    arrivals=tauPmodel.get_travel_times(source_depth_in_km=source_depth_in_km, 
                                   distance_in_degree=distance_in_degree, phase_list=["P","p","S","s"]);
    i = 0
    pi = 0
    si = 0
    while(i<len(arrivals)):
        arr = arrivals[i]
        i = i + 1
        if((arr.name == 'P' or arr.name == 'p') and pi == 0):
           # pname = arr.name
            p_time = arr.time
            #p_ray_param = arr.ray_param*2*ny.pi/360
            #p_hslowness = -1*(p_ray_param/111.19)/math.tan(arr.takeoff_angle*math.pi/180)
            pi = 1

        if((arr.name == 'S' or arr.name == 's') and si == 0):
            #sname = arr.name
            s_time = arr.time
            #s_ray_param = arr.ray_param*2*ny.pi/360
            #s_hslowness = -1*(s_ray_param/111.19)/math.tan(arr.takeoff_angle*math.pi/180)
            si = 1
            if(pi == 1 and si == 1):
                break
    return {'ptime':p_time,'stime':s_time}

def cal_taupPtraveltimes(tauPmodel,source_depth_in_km,distance_in_degree):
    ptimes=[]
    stimes=[]
    for i in range(len(source_depth_in_km)):
        tmp=cal_taupPtraveltime(tauPmodel, source_depth_in_km=source_depth_in_km[i], 
                                distance_in_degree=distance_in_degree[i])
        ptimes.append(tmp['ptime'])
        stimes.append(tmp['stime'])    
    return {'ptimes':ptimes,'stimes':stimes}

    
def time_convert(tbegin,tdetec,ptimes,winlen=30.0):
        t0win=tbegin
        tpa=tdetec
        torg=tpa-min(ptimes)
        dtp=(t0win+winlen)-tpa
        #print(dtp,ptimes,tpa,t0win,winlen)
        dtps=dtp-(ptimes-min(ptimes))
        return {'torg':torg,'dtps':dtps}
    
def mag_trcs2event_prep(dtps,mags,maglen=2.0):
    mags1=[]
    mags1_idx=[]
    for i in range(len(mags)):
        if dtps[i]>maglen:
            mags1.append(mags[i][0])
            mags1_idx.append(i)
    if mags1==[]:
        return {'mags':[],'mags_idx':[],'mean':-9999.0,'median':-9999.0}
    mags1_mean=np.mean(mags1)
    mags1_median=np.median(mags1)
    return {'mags':mags1,'mags_idx':mags1_idx,'mean':mags1_mean,'median':mags1_median}

import json
class NetModel():
    def __init__(self,modeljs_file='model_info.json',datajs_file='data_info.json'):
        with open(modeljs_file,'r') as f:
            modeljs=json.load(f)
        print('loading stream...')
        self.datastream=DataStream(modeljs_file=modeljs_file,datajs_file=datajs_file)
        datajs=self.datastream.datajs
        print('loading model...')
        self.model_dl=load_model(modeljs['model_dl_file'])
        self.model_m=load_model(modeljs['model_m_file'])
        self.xyz_range=modeljs['xyz_range']
        self.mag_range=modeljs['mag_range']
        
        self.coordinate={'data_org':datajs['data_org'],'model_org':modeljs['model_org']}
        self.sample_dt=modeljs['sample_dt']
        self.modeljs=modeljs
        self.datajs=datajs
        self.dl_re={}
        self.dlm_re={}
        print('load data from stream...')
        self.windata=self.datastream.get_data()
        print('load data end.')
 
    def predict_dl(self):
        xyz_predict=[]
        detec_predict=[]
        #tbegin=self.datastream.datajs['tbegin']
        windata=self.windata['windata']
        xdata=[]
        for i in range(len(windata['xdata'])):
            stn=windata['stns'][i]
            stn=[[i['X'],i['Y']] for i in stn if i != {}]
            stn=stn+[[0.0,0.0]]*(12-len(stn))
            tmp=windata['xdata'][i].tolist()
            tmp=sort_data(tmp, stn)
            xdata.append(tmp)
        xdata=np.array(xdata)
        #xdata=[sort_data(windata['xdata'][i].tolist(), windata['stns'][i].tolist())
        #                     for i in range(len(windata['xdata']))];xdata=np.array(xdata)
        loca_data=xdata
        detec_data=xdata[:,:,:,0:3]
        dl_predict=self.model_dl.predict([loca_data,detec_data],batch_size=1, verbose=1)
        xyz_predict=result_loca(r=self.xyz_range,imgs=dl_predict[0],coordinate=self.coordinate)
        detec_predict=result_detect(imgs=dl_predict[1])
        self.dl_re={'xyz_predict':xyz_predict,'detec_predict':detec_predict,'windata':windata,
                    'dl_predict':dl_predict}
        return self.dl_re
    def dl_re_check(self,minpdf=0.0):
        detec_predict=self.dl_re['detec_predict']
        xyz_predict=self.dl_re['xyz_predict']
        windata=self.dl_re['windata']
        tbegin=windata['tbegins']
        dl_re_c=[]
        for dr in detec_predict:
            if dr[1]>minpdf:
                idx=dr[2]
                xyz=xyz_predict[idx]['xyz']
                geo=xyz_predict[idx]['geo']
                tb=UTCDateTime(tbegin[idx])
                t0=[tb+dr[0]*self.sample_dt,dr[1],tb,idx]
                d_img=self.dl_re['dl_predict'][1][idx]
                l_img=self.dl_re['dl_predict'][0][idx]
                dl_re_c.append({'xyz':xyz,'geo':geo,'t0':t0,'xdata':windata['xdata'][idx],
                                'stns':windata['stns'][idx],'maxval':windata['maxvals'][idx],
                                'd_img':d_img,'l_img':l_img})
        self.dl_re_c=dl_re_c
        return self.dl_re_c  
    
    def predict_mag1event(self,xdata,amps):
        mag_img=self.model_m.predict(xdata, batch_size=1, verbose=0)
        mag_predict=result_mag(mag_range=self.mag_range,imgs=mag_img,amps=amps)
        mag_predict=mag_trcs2event(mags=mag_predict)
        mag_predict['mag_img']=mag_img
        return mag_predict

    
    def predict_m_dlrec(self):
        self.dlm_re=self.dl_re_c
        ttmodel=TauPyModel('prem')
        windata=self.windata['windata_mag']
        for dl in self.dlm_re:
            geo=dl['geo']
            idx=dl['t0'][3]
            stns=windata['stns'][idx]
            xdata=windata['xdata'][idx]
            maxval=windata['maxvals'][idx]
            stns_geo=[[i['geo'][0],i['geo'][1]] for i in stns if i != {}]
            if stns_geo==[] or dl['t0'][1]<0.1:
               stns_geo_num=len(stns_geo)
               dl['mag_prep']={'mags':[[-999.0,0.0]],'mags_idx':[0],'mean':-9999.0,'median':-9999.0}
               dl['mag']={'mag_valid':-9999,'mag_mean':-9999,'mag_median':-9999,'mags':[[-999.0,0.0]]*stns_geo_num,'mag_img':np.array([[[0]]*self.mag_range[2]]*stns_geo_num)}
               dl['torg']=UTCDateTime('1980-01-01T00:00:00.00Z')
               dl['dtps']=np.array([-999]*stns_geo_num)
               #print('no data for magnitude predict in predict_m_dlrec()')
               continue
            tmp=geo2dist(src_geo=geo, stns_geo=stns_geo);
            dists=tmp['dists']
            pstimes=cal_taupPtraveltimes(tauPmodel=ttmodel,source_depth_in_km=[dl['xyz'][2]]*len(tmp['deltas']),
                                        distance_in_degree=tmp['deltas'])
            tmp_time=time_convert(tbegin=dl['t0'][2], tdetec=dl['t0'][0], ptimes=pstimes['ptimes'])
            dists=dists+[0]*(12-len(dists))
            tmp=mag_data1event(event_data=xdata, dists=dists, dist_range=self.modeljs['dist_range'])
            mag=self.predict_mag1event(xdata=tmp['data'], amps=(tmp['amps']*maxval))
            tmp2=mag_trcs2event_prep(dtps=tmp_time['dtps'], mags=mag['mags'])
            dl['mag_prep']=tmp2
            dl['mag']=mag
            dl['torg']=tmp_time['torg']
            dl['dtps']=tmp_time['dtps']
        return self.dlm_re
    
    def predict_all(self):
        self.predict_dl()
        self.dl_re_check()
        self.predict_m_dlrec()
    
    def save_dlm(self,filename='test_dlm.csv'):
        with open(filename,'w') as f:
            for dl in self.dlm_re:
                mags_s=['[{:.3f},{:.3f}]'.format(i[0],i[1]) for i in dl['mag']['mags']]
                dtps_s=['{:.3f}'.format(i) for i in dl['dtps']]
                f.write('{}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {}, {:.6f}, {:.6f}, {}, {}\n'.format(str(dl['t0'][0]),dl['t0'][1],
                        dl['xyz'][0],dl['xyz'][1],dl['xyz'][2],dl['xyz'][3],
                        dl['geo'][0],dl['geo'][1],dl['t0'][2],
                        dl['mag_prep']['median'],dl['mag_prep']['mean'],mags_s,dtps_s))
    
    def save_data(self,filename='test_dlm_imgs.mat'):
        l_imgs=[]
        wave_test=[]
        xyzs=[]
        geos=[]
        stns=[]
        d_imgs=[]
        d_t0=[]
        d_t0win=[]
        mag_imgs=[]
        mag_mags=[]
        mag_dtps=[]
        if 'output_3dimags_range' in self.datajs:
            outrg=self.datajs['output_3dimags_range']
        else:
            outrg=[0,100,1]
        for dl in self.dlm_re[outrg[0]:outrg[1]:outrg[2]]:
            l_imgs.append(dl['l_img'])
            wave_test.append(dl['xdata'][:,:,0:3])
            xyzs.append(dl['xyz'])
            geos.append(dl['geo'])
            stns.append(dl['stns'])
            d_imgs.append(dl['d_img'])
            d_t0.append(dl['t0'][0])
            d_t0win.append(dl['t0'][1])
            mag_imgs.append(dl['mag']['mag_img'])
            mag_mags.append(dl['mag']['mags'])
            mag_dtps.append(dl['dtps'])
        sio.savemat(filename, {'l_imgs':l_imgs,'xyzrange':self.modeljs['xyz_range'][0:3],'wave_test':wave_test,
                                  'xyzs':xyzs,'geos':geos,'stns':stns,'d_imgs':d_imgs,'d_t0':d_t0,
                                  'd_t0win':d_t0win,'mag_imgs':mag_imgs,'mag_mags':mag_mags,'mag_dtps':mag_dtps})
        
    def save_re(self,filename='test_dlm.pkl'):
        import pickle
        z_range=self.modeljs['xyz_range'][2]
        twin_nsmp=self.modeljs['twin_nsmp']
        dlm_re=[]
        for dl in self.dlm_re:
            iz=int((dl['xyz'][2]-z_range[0])/z_range[1]+0.5)
            dlm_re.append({'d_t0':str(dl['t0'][0]),'d_pdf':dl['t0'][1],'xyz':dl['xyz'],'geo':dl['geo'],'d_t0win':str(dl['t0'][2]),
                           'mag_median':dl['mag_prep']['median'],'mag_mean':dl['mag_prep']['mean'],'mag_mags':dl['mag']['mags'],
                           'mag_dtps':dl['dtps'],'d_img':dl['d_img'],'mag_img':dl['mag']['mag_img'],'img2d':dl['l_img'][:,:,iz],
                           'wave_test':dl['xdata'][:,0:twin_nsmp,0:3],'torg':str(dl['torg'])})
        with open(filename,'wb') as f:
            pickle.dump(dlm_re, f)
    
class MagModel():
    def __init__(self,modeljs_file='model_info.json'):
        with open(modeljs_file,'r') as f:
            modeljs=json.load(f)
        self.datastream=DataStream()
        self.model_m=load_model(modeljs['model_m_file'])
        self.mag_range=modeljs['mag_range']
        self.sample_dt=modeljs['sample_dt']
        self.modeljs=modeljs
    def predict_mags(self,xdata,amps):
        mag_img=self.model_m.predict(xdata, batch_size=1, verbose=1)
        mag_predict=result_mag(mag_range=self.mag_range,imgs=mag_img,amps=amps)
        return {'mag_and_pdf':mag_predict,'mag_img':mag_img}
    def predict_mags_from_waves(self,waves,dists):
        tmp=mag_data1event(event_data=waves, dists=dists, dist_range=self.modeljs['dist_range'])
        xdata=tmp['xdata'];amps=tmp['amps']
        mag_re=self.predict_mags(xdata, amps)
        return mag_re
    def relative_mags_from_amps(self,amps,mags):
        return np.array(mags)-np.log10(np.array(amps))
    
        
    
        
        
    
    


if __name__ == '__main__':
    import sys
    #import getopt
    arg_len=len(sys.argv)
    if arg_len<2:
        print('usage: netmodel function_name ...')
        print('python netmodel.py moni model_info.json data_info.json 1')
        sys.exit()
    else:
        arg_fun=sys.argv
    print('ss',arg_len,arg_fun)
    if arg_len==5 and arg_fun[1]=='moni':
        import os
        os.environ['CUDA_VISIBLE_DEVICES']=arg_fun[4]
        netmodel=NetModel(modeljs_file=arg_fun[2],datajs_file=arg_fun[3])
        datajs=netmodel.datajs
        netmodel.predict_all()
        netmodel.save_dlm(filename=datajs['out_nam']+'_dlm.csv')
        netmodel.save_data(filename=datajs['out_nam']+'_dlm.mat')
        netmodel.save_re(filename=datajs['out_nam']+'_dlm.pkl')
    


