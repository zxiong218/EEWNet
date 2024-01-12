# -*- coding: utf-8 -*-
"""
Created on Tue May 16 11:19:57 2023

@author: huawei
"""
from distaz import DistAz
import numpy as np
from obspy import Stream,Trace,UTCDateTime
import obspy
from math import cos,sin
import math
from obspy.signal.trigger import recursive_sta_lta, trigger_onset, classic_sta_lta, plot_trigger
def nam2utctime(nam):
    nam1=nam[0:13]+':'+nam[14:16]+':'+nam[17:27]
    utctime=UTCDateTime(nam1)
    return {'str':nam1,'time':utctime}
    
def merge_rdlm_pkl(nams,outnam):
    pkls=[]
    for i in range(len(nams)):
        tmp_pkl=read_pkl(nams[i])
        pkls=pkls+tmp_pkl
    write_pkl(pkls,outnam)
def reshape_rdlms(rdims):
    rdims0=[]
    for i0 in range(len(rdims[0])):
        rdims_tmp=[]
        for i1 in range(len(rdims)):
            rdims_tmp.append(rdims[i1][i0])
        rdims0.append(rdims_tmp)
    return {'rdims':rdims0}
def read_pkls(nams):
    pkls=[]
    for i in range(len(nams)):
        tmp_pkl=read_pkl(nams[i])
        pkls.append(tmp_pkl)
    return {'pkls':pkls}
def read_h5_info(nam):
    import h5py
    h5obj=h5py.File(nam, 'r')
    keylist=list(h5obj.keys())
    data1=np.array(h5obj.get(keylist[0]))
    print('data[0],shape,min,max',data1[0],data1.shape,np.max(data1),np.min(data1))
    print('key:',keylist[0:10],'...',keylist[-10:])
    print('length of keys',len(keylist))
    print('data attrs',dict(h5obj.get(keylist[0]).attrs))
    h5obj.close()
def show_h5_data(nam,key):
    import h5py
    import scipy.io as sio
    h5obj=h5py.File(nam, 'r')
    data=h5obj.get(key)
    print(np.array(data),data)
    sio.savemat('tmp.mat', {'data':np.array(data)})
    write_sacs([np.array(data)],dt=data.attrs['delta'])
    h5obj.close()
def read_h5_datasets(h5obj,keys):
    data=[]
    for j in range(len(keys)):
        tmpdata=h5obj.get(keys[j])
        data.append(list(tmpdata))
    return {'data':np.array(data)}
def write_pkl(dic,outnam):
    import pickle
    with open(outnam,'wb') as f:
        pickle.dump(dic,f)
def read_pkl(nam):
    import pickle
    with open(nam,'rb') as f:
        st_samps=pickle.load(f)  
    return st_samps
def distaz(evlon,evlat,stlon,stlat):
    tmp=DistAz(stalat=stlat, stalon=stlon, evtlat=evlat, evtlon=evlon)
    return {'dist':tmp.getDistanceKm(),'az':tmp.getAz()}
    
def read_catalog(nam):
    with open(nam) as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        evs=[[i[0],float(i[1]),float(i[2]),float(i[3]),float(i[4])]+i[5:] for i in txtlines]
    return {'evs':evs}

def read_stns(nam):
    with open(nam) as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        stns=[[float(i[0]),float(i[1]),float(i[2])]+i[3:] for i in txtlines]
    return {'stns':stns}

def read_catalog_dic(nam):
    with open(nam) as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        evs=[{'t0':i[0],'geo':[float(i[1]),float(i[2])],'dep':float(i[3]),'mag':float(i[4])} for i in txtlines]
    return {'evs':evs}

def read_stns_dic(nam, channel=['E','N','Z']):
    with open(nam) as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        stns=[{'geo':[float(i[0]),float(i[1])],'dep':float(i[2]),'stnam':i[3]+'.'+i[4],'channel':channel} for i in txtlines]
    return {'stns':stns}

def read_json(nam):
    import json
    with open(nam,'r') as f:
        js=json.load(f)
    return js

def write_sacs(data,dt=0.05,starttime='2023-01-01T00:00:00.00',nam='test'):
    for i in range(len(data)):
        trace = Trace(data=data[i])  # 创建Trace对象
        trace.stats.station = f"st{i+1}"  # 设置台站名称
        trace.stats.sampling_rate = dt  # 设置采样率
        trace.stats.starttime = UTCDateTime(starttime)
        trace.write(f"test_st{i+1}.SAC", format="SAC")  # 输出为SAC文件

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

def geo2dist(src_geo,stns_geo):
    dists=[]
    deltas=[]
    for stn in stns_geo:
        tmp=DistAz(stalat=stn[1], stalon=stn[0], evtlat=src_geo[1], evtlon=src_geo[0])
        dists.append(tmp.getDistanceKm())
        deltas.append(tmp.getDelta())
    return {'dists':dists,'deltas':deltas}

def cal_pstimes(evs,stns,tauPmodel=None):
    if tauPmodel==None:
       from obspy.taup import TauPyModel
       ttmodel=TauPyModel('prem')
    else:
       ttmodel=tauPmodel
    pstimes=[]
    for j,ev in enumerate(evs):
        #print(ev)
        tmp=geo2dist(src_geo=ev[1:3], stns_geo=stns)
        pstimes1=cal_taupPtraveltimes(tauPmodel=ttmodel,source_depth_in_km=[ev[2]]*len(stns),
                                distance_in_degree=tmp['deltas'])
        pstimes1['ptimes_min']=np.min(pstimes1['ptimes'])
        pstimes.append(pstimes1)
    return pstimes

def cal_magnitude_1st(data,dist,scale=1.0,sampling_rate=20):
    paz_wa = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
                'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
    data1=data.T
    st=data2stream(data1, delta=1.0/sampling_rate)
    st.simulate(paz_remove = None, paz_simulate = paz_wa)
    waamp=[np.max(np.abs(st[i].data))*scale for i in range(3)]
    r=dist
    mag_r=3.0+1.110*math.log10(r/100.0)+0.00189*(r-100.0)
    mag_amp=math.log10((waamp[0]+waamp[1])/2)
    return {'ml':mag_r+mag_amp}

def cal_magnitude_1event(data,dist,scale=1.0,sampling_rate=20):
    mags=[]
    for i, trc in enumerate(data):
        tmp=cal_magnitude_1st(data=trc, dist=dist[i],scale=scale,sampling_rate=sampling_rate)
        mags.append(tmp['ml'])
    return {'mag_mean':np.mean(mags),'mag_median':np.median(mags),'mags':mags}
    
def cal_mag_ml(wa_amps,dists):
    mags=[]
    mags_valid=[]
    for i in range(len(dists)):
        wa_amp=wa_amps[i]
        if abs(wa_amp[0]-0.0)<0.01E-35 or abs(wa_amp[1]-0.0)<0.01E-35 :
           mags.append(-9999)
           continue
        amp=(wa_amp[0]+wa_amp[1])/2+1.0E-100
        r=dists[i]
        mag_r=3.0+1.110*np.math.log10(r/100.0)+0.00189*(r-100.0)
        mag=np.math.log10(amp)+mag_r
        mags.append(mag)
        mags_valid.append(mag)
    return {'mags':mags,'mag_mean':np.mean(mags_valid),'mag_median':np.median(mags_valid)}

def data2stream(data,delta,starttime='2023-01-01T00:00:00.000'):
    st=[]
    for i,trc in enumerate(data):
        trc1=obspy.Trace(trc)
        trc1.stats.starttime=UTCDateTime(starttime)
        trc1.stats.delta=delta
        #tr.stats.channel
        st.append(trc1)
    return obspy.Stream(st)
        

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

def h5py_read(hf,trcns,freq=[2.0,8.0],nsmpout=1024,win=[]):
    paz_wa = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
                'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
    wa_amp=[]
    nbrok=0
    for i,trcn in enumerate(trcns):
        st_tmp=h5py_getObspyStream(trcns=[trcn],hf=hf)
        try:
            st_tmp.filter('bandpass',freqmin=freq[0],freqmax=freq[1])
            st_tmp_mag=st_tmp.copy().slice(starttime=win[0],endtime=win[0]+win[1],keep_empty_traces=True);
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
            print('no data in',e,trcn);
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

def h5py_trcns(prelist,stnams,chns):
    trcns=[]
    for i,stnam in enumerate(stnams):
        trcns=trcns+[prelist+stnam+chns[0],prelist+stnam+chns[1],prelist+stnam+chns[2]]
    return trcns

def sac2h5(sacnams,h5file,mode='w'):
    import h5py
    from obspy import read
    import numpy as np
    h5obj=h5py.File(h5file,mode)
    #print(sacnams)
    wave=read(sacnams[0])
    for i in range(1,len(sacnams)):
        wave=wave+read(sacnams[i])
        
    for j in range(0,len(wave)):
        #print(wave[j])
        trcn=sacnams[j]
        tmp=h5obj.create_dataset(trcn,data=np.array(wave[j].data))
        tmp.attrs['network']=wave[j].stats.network
        tmp.attrs['station']=wave[j].stats.station
        tmp.attrs['channel']=wave[j].stats.channel
        tmp.attrs['starttime']=str(wave[j].stats.starttime)
        tmp.attrs['endtime']=str(wave[j].stats.endtime)
        tmp.attrs['sampling_rate']=wave[j].stats.sampling_rate
        tmp.attrs['delta']=wave[j].stats.delta
        tmp.attrs['calib']=wave[j].stats.calib
        tmp.attrs['npts']=wave[j].stats.npts
        tmp.attrs['location']=wave[j].stats.location
    h5obj.close()

def ltasta(data,smp_rate,thred=[15,10]):
    df=smp_rate
    cft = recursive_sta_lta(data, int(0.1 * df), int(2.5 * df))
    on_of = trigger_onset(cft, thred[0], thred[1])
    #print('onof',on_of,data)
    if on_of==[]:
        on_of=[[99990,99990]]
    return {'on_of':on_of,'tp':on_of[0][0]/df,'cft':cft}

def pick_1event(data,smp_rate):
    pt=[]
    for i, data0 in enumerate(data):
      tmp=ltasta(data=data0, smp_rate=smp_rate)['tp']
      pt.append(tmp)
    #print(pt)
    return {'tp':pt,'tpmin':np.min(pt)}

def ltasta_zx(data,sn=10,ln=20,rg=[0,-1]):
    data1=data*data
    rt=[]
    for i in range(len(data1)-sn):
        sta=np.mean(data1[i:i+sn])
        if i>=ln:
            lta=np.mean(data1[i-ln:i])
            rt.append(sta/(lta+0.1E-20))
        else:
            rt.append(1.0)
    rt=rt+[0.0]*(len(data1)-len(rt))
    rt=np.array(rt)
    rtmax=max(rt[rg[0]:rg[1]])
    imax=np.where(rt[rg[0]:rg[1]]==rtmax)[0][0]
    return {'imax':imax+rg[0],'rtmax':rtmax,'rt':rt}

def pick_refine_1event(data,smp_rate,rgs=[[0,10.0]]):
    pt=[]
    rgs1=[[int(i[0]*smp_rate+0.5),int(i[1]*smp_rate+0.5)] for i in rgs]
    #rgs1=[[i[0],i[1]] if i[0]>0 and i[0]<len(data[0]) else [0, i[1]] for i in rgs1]
    pt_rt=[]
    for i in range(len(rgs1)):
        dt=rgs1[i][1]-rgs1[i][0]
        if rgs1[i][0]<0:
            rgs1[i][0]=0
            rgs1[i][1]=dt
            print('warning: pick out of range.')
        elif rgs1[i][0]>=len(data[0]):
            rgs1[i][0]=-2
            rgs1[i][1]=-1
            print('warning: pick out of range. 1')
    #print('in pytool',rgs1,'dl',len(data),rgs)
    for i, data0 in enumerate(data):
      tmp=ltasta_zx(data=data0,rg=rgs1[i])
      tp=tmp['imax']/smp_rate
      pt.append(tp)
      pt_rt.append(tmp['rtmax'])
    #print(pt)
    return {'tp':pt,'tp_rt':pt_rt,'tpmin':np.min(pt)}

def ttgrid(tps,stns_geo,rg):
    xr=rg['xr']
    yr=rg['yr']
    zr=rg['zr']
    from obspy.taup import TauPyModel
    ttmodel=TauPyModel('prem')
    obj=10E20
    tps_min=np.min(tps)
    for ix in range(xr[2]):
        x=xr[0]+xr[1]*ix
        for iy in range(yr[2]):
            y=yr[0]+yr[1]*iy
            for iz in range(zr[2]):
                z=zr[0]+zr[1]*iz
                tmp=geo2dist(src_geo=[x,y], stns_geo=stns_geo)
                ptimes=cal_taupPtraveltimes(tauPmodel=ttmodel, source_depth_in_km=[z]*len(stns_geo), 
                                            distance_in_degree=tmp['deltas'])
                ptimes['ptimes_min']=np.min(ptimes['ptimes'])
                #print(tps,ptimes,x,y,z,obj)
                tmp=(tps-tps_min)-(ptimes['ptimes']-ptimes['ptimes_min'])
                tmp=np.abs(tmp)
                obj_tmp=np.sum(tmp[np.where(tmp<6.0)])
                #print('tmp',tmp[np.where(tmp<6.0)],obj_tmp,z)
                if obj>obj_tmp:
                    obj=obj_tmp
                    x0=x
                    y0=y
                    z0=z
                    ptimes0=ptimes
    return {'loca':[x0,y0,z0],'syn_ptimes':ptimes0['ptimes']-ptimes0['ptimes_min'],
            'real_ptimes':(tps-tps_min),'obj':obj}

def stda_degree2km(center=[135.622-0.05, 34.844+0.22],stnr=[[0,82],[0,100]]):
    tmp=DistAz(stalon=center[0],stalat=center[1],evtlon=center[0]+1,evtlat=center[1])
    dist=tmp.getDistanceKm()
    tmp=stnr[0][1]/dist
    tmp1=stnr[1][1]/111.19
    dgr=[[center[0]-tmp/2.0,center[0]+tmp/2.0],[center[1]-tmp1/2.0,center[1]+tmp1/2.0]]
    return {'lonmin':dgr[0][0],'lonmax':dgr[0][1],'latmin':dgr[1][0],'latmax':dgr[1][1],'scale':[dist,111.19],
            'center':center,'center_km':[center[0]*dist,center[1]*111.19]}
def rotate(data_E,data_N,rt_clockwise=0.0):#rt_clockwise is "data rotate clokwisly" equal to "coordinate rotate anticlockwisly" 
    rt_clockwise1=rt_clockwise*math.pi/180.0
    tmp1=data_E*cos(rt_clockwise1)+data_N*sin(rt_clockwise1)
    tmp2=-data_E*sin(rt_clockwise1)+data_N*cos(rt_clockwise1)
    return tmp1,tmp2
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
def coordinate2model_stn(stnx,stny,data_org={'org':[-10680,4005],'x':111.19*0.816,'y':111.19},
                         model_org=[82.0/2,100.0/2],
                         rt_clockwise=0.0):
    stnx1=stnx*data_org['x']-data_org['org'][0]
    stny1=stny*data_org['y']-data_org['org'][1]
    tmp1,tmp2=rotate(stnx1,stny1,rt_clockwise=rt_clockwise)
    tmp1=tmp1+model_org[0]
    tmp2=tmp2+model_org[1]
    return tmp1,tmp2

def coordinate2model(lon,lat,coordinate=
                    {'data_org':{'lonlat':[0,0],'xy':[0*111.19*0.813,0*111.19],
                                                  'scale':[111.19*0.813,111.19],'rt_clockwise':0.0},
                                      'model_org':[82.0/2,100.0/2]}):
    data_org=coordinate['data_org']
    model_org=coordinate['model_org']
    x1=lon*data_org['scale'][0]-data_org['xy'][0]
    y1=lat*data_org['scale'][1]-data_org['xy'][1]
    tmp1,tmp2=rotate(x1,y1,rt_clockwise=data_org['rt_clockwise'])
    tmp1=tmp1+model_org[0]
    tmp2=tmp2+model_org[1]
    return tmp1,tmp2

def coordinate_info(model_json,data_json):
    import json
    with open(model_json,'r') as f:
        modeljs=json.load(f)
    with open(data_json,'r') as f:
        datajs=json.load(f)
    coordinate={'data_org':datajs['data_org'],'model_org':modeljs['model_org']}
    return coordinate
                    

def inbox(box,loca):
    if loca[0]>=box[0][0] and loca[0]<=box[0][1] and loca[1]>=box[1][0] and loca[1]<=box[1][1]:
        return True
    else:
        return False

def monitor_data_org(stns,bias=10,boxsize=[82.0,100.0]):
    stns_center=[np.mean([i[0] for i in stns]),np.mean([i[1] for i in stns])]
    xscale=DistAz(stalon=stns_center[0],stalat=stns_center[1],evtlon=stns_center[0]+1,evtlat=stns_center[1]).getDistanceKm()
    scale=[xscale,111.19]
    xy=[stns_center[0]*scale[0],stns_center[1]*scale[1]]
    box=[[xy[0]-boxsize[0]/2.0,xy[0]+boxsize[0]/2.0],[xy[1]-boxsize[1]/2.0,xy[1]+boxsize[1]/2.0]]
    stns1=[]
    stns1_out=[]
    stns1_xy=[]
    stns1_out_xy=[]
    for i, stn in enumerate(stns):
        xy1=[stn[0]*scale[0],stn[1]*scale[1]]
        if inbox(box,xy1):
            stns1.append(stn)
            stns1_xy.append(xy1)
        else:
            stns1_out.append(stn)
            stns1_out_xy.append(xy1)
    new_center_xy=[np.mean([i[0] for i in stns1_xy]),np.mean([i[1] for i in stns1_xy])]
    tmp=(new_center_xy[0]-xy[0]);tmp1=(new_center_xy[1]-xy[1]);dist=(tmp*tmp+tmp1*tmp1)**0.5
    if dist>bias:
        print('warning: the station deleted make bias for the monitoring center in monitor_data_org(), bias',dist) 
    stn_box_km=[[xy[0]-boxsize[0]/2,xy[0]+boxsize[0]/2],[xy[1]-boxsize[1]/2,xy[1]+boxsize[1]/2]]
    stn_box_geo=[[stn_box_km[0][0]/scale[0],stn_box_km[0][1]/scale[0]],
                 [stn_box_km[1][0]/scale[1],stn_box_km[1][1]/scale[1]]]
    return {'stns':stns1,'stns1_out':stns1_out,'bias':dist,'stn_box_geo':stn_box_geo,'stn_box_km':stn_box_km,
            'data_org':{'lonlat':stns_center,'xy':xy,'scale':scale,'rt_clockwise':0.0}}   
    
def monitor_select_stations(stns,expect_center,station_range=[[0,82],[0,100]],bias=15):
    dists=geo2dist(src_geo=expect_center, stns_geo=stns)['dists']
    maxdist=((station_range[0][1]**2.0+station_range[1][1]**2.0)**0.5)*0.8
    stns1=[]
    for i in range(len(dists)):
        for j in range(i,len(dists)):
            if dists[i]>dists[j]:
                tmp=dists[i]
                dists[i]=dists[j]
                dists[j]=tmp
                tmp=stns[i]
                stns[i]=stns[j]
                stns[j]=tmp
    stns1=[]
    for i in range(len(dists)):
        if dists[i]<maxdist and len(stns1)<12:
            stns1.append(stns[i])
    tmp=monitor_data_org(stns=stns1)
    if inbox(box=tmp['stn_box_geo'],loca=expect_center):
        return tmp
    else:
        return None

def align_mag_img(mag_img, imag):
    mag_img1=list(mag_img)
    maxv = max(mag_img1)  # 获取mag_img中的最大值
    maxv_index = mag_img1.index(maxv)  # 获取最大值在mag_img中的位置
    # 计算需要平移的位移量
    shift = -int((maxv_index - imag) % len(mag_img1))
    # 对mag_img进行平移
    aligned_vector = mag_img1[-shift:] + mag_img1[:-shift]
    return {'mag_img':aligned_vector}

def align_mag_imgs(mag_imgs, imags):
    mag_imgs1=[]
    for i in range(len(imags)):
        tmp=align_mag_img(mag_imgs[i], imags[i])['mag_img']
        mag_imgs1.append(tmp)
    return {'mag_imgs':np.array(mag_imgs1)}
                  
if __name__ == '__main__':
    import sys
    #import getopt
    arg_len=len(sys.argv)
    if arg_len<2:
        print('usage: pytool function_name ...')
        print('pytool show_h5_info test.h5')
        print('pytool h5data test.h5 idx1')
        print('pytool samp-st_samples_info st_samples.pkl')
        print('pytool samp-st_samples_extr_trace st_samples.pkl 40 8')
        sys.exit()
    else:
        arg_fun=sys.argv
    print('ss',arg_len,arg_fun)
    if arg_len==3 and arg_fun[1]=='show_h5_info':
        read_h5_info(arg_fun[2])
    elif arg_len==4 and arg_fun[1]=='h5data':
        show_h5_data(arg_fun[2], arg_fun[3])
    elif arg_len==3 and arg_fun[1]=='samp-st_samples_info':
        from samp import st_samples_info
        st_samples_info(arg_fun[2])               
    elif arg_len==5 and arg_fun[1]=='samp-st_samples_extr_trace':
        from samp import st_samples_extr_trace
        stsam=st_samples_extr_trace(arg_fun[2],float(arg_fun[3]),float(arg_fun[4]))               
        print(stsam)
    
    


    
        
    
    
    
    
    
    
    
    
        
