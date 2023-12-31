# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 20:23:52 2023

@author: huawei
"""

import numpy as np
def check_moni_setting(stns_geo,data_org,bias=10,boxsize=[82.0,100.0]):
    from pytool import inbox
    from distaz import DistAz
    stns_center=[np.mean([i[0] for i in stns_geo]),np.mean([i[1] for i in stns_geo])]
    xscale=DistAz(stalon=stns_center[0],stalat=stns_center[1],evtlon=stns_center[0]+1,evtlat=stns_center[1]).getDistanceKm()
    scale=[xscale,111.19]
    xy=[stns_center[0]*scale[0],stns_center[1]*scale[1]]
    scale1=data_org['scale']
    xy1=data_org['xy']
    box=[[xy1[0]-boxsize[0]/2.0,xy1[0]+boxsize[0]/2.0],[xy1[1]-boxsize[1]/2.0,xy1[1]+boxsize[1]/2.0]]
    stns_out=[]
    stns_out_id=[]
    stn_E=1000.0
    stn_W=1000.0
    stn_N=1000.0
    stn_S=1000.0
    for i,stn in enumerate(stns_geo):
        stn_xy=[stn[0]*scale1[0],stn[1]*scale1[1]]
        tmp_W=stn_xy[0]-box[0][0]
        tmp_E=stn_xy[0]-box[0][1]
        tmp_S=stn_xy[1]-box[1][0]
        tmp_N=stn_xy[1]-box[1][1]
        if tmp_W<stn_W:
            stn_W=tmp_W
        if tmp_E<stn_E:
            stn_E=tmp_E
        if tmp_N<stn_N:
            stn_N=tmp_N
        if tmp_S<stn_S:
            stn_S=tmp_S
        if not inbox(box, loca=stn_xy):
            stns_out.append(stn)
            stns_out_id.append(i)
    data_org0={"lonlat":stns_center,"xy":xy,"scale":scale,"rt_clockwise":0.0}
    #tmp_bias=((xy[0]-xy1[0])**2+(xy[1]-xy1[1])**2)**0.5
    tmp_bias=DistAz(stalon=stns_center[0],stalat=stns_center[1],evtlon=data_org['lonlat'][0],evtlat=data_org['lonlat'][1]).getDistanceKm()
    if len(stns_out)!=0:
        return {'flg':1,'flg_str':'stns_out','stns_out':stns_out,'stns_out_id':stns_out_id,'data_org':data_org0,
                'stn_W':stn_W,'stn_E':stn_E,'stn_N':stn_N,'stn_S':stn_S}
    elif tmp_bias>bias:
        return {'flg':2,'flg_str':'center_bias','data_org':data_org0,
                'stn_W':stn_W,'stn_E':stn_E,'stn_N':stn_N,'stn_S':stn_S,'bias':tmp_bias}
    else:
        return {'flg':0,'flg_str':'pass','data_org':data_org0,
                'stn_W':stn_W,'stn_E':stn_E,'stn_N':stn_N,'stn_S':stn_S}


def check_js_files(datajs_file,modeljs_file,bias=10):
    import json
    with open(datajs_file,'r') as f:
        datajs=json.load(f)
    with open(modeljs_file,'r') as f:
        modeljs=json.load(f)
    with open(datajs['stations_file']) as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        stns_geo=np.array([[float(i[0]),float(i[1])] for i in txtlines])
    station_range=modeljs['station_range']
    boxsize=[station_range[0][1]-station_range[0][0],station_range[1][1]-station_range[1][0]]
    tmp=check_moni_setting(stns_geo, datajs['data_org'],bias=bias,boxsize=boxsize)
    if tmp['flg']==0:
        print('Checking pass..')
    elif tmp['flg']==1:
        print('some stations are out of monitoring range, listed as follows:',
              tmp['stns_out'])
    elif tmp['flg']==2:
        print('warning: the station coverage may be not good; you should set data_org again')
        
def check_renew_datajsfiles(datajs_file,modeljs_file,orgsetting_file,out_station_file,bias=10,center_type='expected'):
    import json
    with open(datajs_file,'r') as f:
        datajs=json.load(f)
    with open(modeljs_file,'r') as f:
        modeljs=json.load(f)
    with open(orgsetting_file) as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        stns_geo=np.array([[float(i[1]),float(i[2])] for i in txtlines[6:]])
        data_org_expe={'lonlat':[float(txtlines[0][1]),float(txtlines[0][2])],
                     'scale':[float(txtlines[0][3]),float(txtlines[0][4])],
                     'xy':[float(txtlines[0][1])*float(txtlines[0][3]),float(txtlines[0][2])*float(txtlines[0][4])],
                     'rt_clockwise':0.0}
        data_org_mean={'lonlat':[float(txtlines[3][1]),float(txtlines[3][2])],
                     'scale':[float(txtlines[3][3]),float(txtlines[3][4])],
                     'xy':[float(txtlines[3][1])*float(txtlines[3][3]),float(txtlines[3][2])*float(txtlines[3][4])],
                     'rt_clockwise':0.0}
        if center_type=='expected':
            data_org=data_org_expe
        elif center_type=='mean':
            data_org=data_org_mean
        else:
            print('Not recognized center_type in check_renew_jsfile()')
    station_range=modeljs['station_range']
    boxsize=[station_range[0][1]-station_range[0][0],station_range[1][1]-station_range[1][0]]
    tmp=check_moni_setting(stns_geo, data_org,bias=bias,boxsize=boxsize)
    if tmp['flg']==0:
        datajs['data_org']=data_org
        datajs['stations_file']=out_station_file
        with open(out_station_file,'w') as f:
            for i, txt in enumerate(txtlines[6:]):
                f.write('{} {} {} {} {}\n'.format(txt[1],txt[2],txt[3],txt[4].split('.')[0],txt[4].split('.')[1]))
        with open(datajs_file,'w') as f:
            json.dump(datajs,f,indent=4)
        print('Checking pass..\n','{} file renewed...'.format(datajs_file))
    elif tmp['flg']==1:
        print('some stations are out of monitoring range, listed as follows:',
              tmp['stns_out'])
    elif tmp['flg']==2:
        print('warning: the station coverage may be not good; you should set data_org again',tmp['bias'])
        
def remove_stn_bymedian(stns):
    from distaz import DistAz
    stns_center=[np.median([i[0] for i in stns]),np.median([i[1] for i in stns])]
    dist=0
    for i,stn in enumerate(stns):
        tmp=DistAz(stalon=stns_center[0],stalat=stns_center[1],evtlon=stn[0],evtlat=stn[1]).getDistanceKm()
        if tmp>dist:
            dist=tmp
            dist_id=i
    stns1=[stns[i] for i in range(len(stns)) if i !=dist_id]
    return stns1

def adjust_data_org(stbias,data_org,bias=10,bias_err=0.01):
    xy=[data_org['xy'][0],data_org['xy'][1]]
    if stbias['stn_W']<0 and stbias['stn_E']>0:
        dx0=abs(stbias['stn_W'])+bias_err;dx1=abs(stbias['stn_E'])
        if dx0<dx1:
            xy[0]=xy[0]-dx0
        else:
            return {'flg':1,'flg_str':'err_W'}
    if stbias['stn_W']>0 and stbias['stn_E']<0:
        dx0=abs(stbias['stn_E'])+bias_err;dx1=abs(stbias['stn_W'])
        if dx0<dx1:
            xy[0]=xy[0]+dx0
        else:
            return {'flg':1,'flg_str':'err_E'}
    if stbias['stn_S']<0 and stbias['stn_N']>0:
        dx0=abs(stbias['stn_S'])+bias_err;dx1=abs(stbias['stn_N'])
        if dx0<dx1:
            xy[1]=xy[1]-dx0
        else:
            return {'flg':1,'flg_str':'err_S'}
    if stbias['stn_S']>0 and stbias['stn_N']<0:
        dx0=abs(stbias['stn_N'])+bias_err;dx1=abs(stbias['stn_S'])
        if dx0<dx1:
            xy[1]=xy[1]+dx0
        else:
            return {'flg':1,'flg_str':'err_N'}
    if stbias['stn_W']<0 and stbias['stn_E']<0:
        return {'flg':1,'flg_str':'err_WE'}
    if stbias['stn_S']<0 and stbias['stn_N']<0:
        return {'flg':1,'flg_str':'err_SN'}
    scale=data_org['scale']
    lonlat=[xy[0]/scale[0],xy[1]/scale[1]]
    data_org1={"lonlat":lonlat,"xy":xy,"scale":scale,"rt_clockwise":0.0}
    return {'flg':0,'flg_str':'done','data_org':data_org1}
        
        
def recommend_moni_range(stns,bias=10,boxsize=[82.0,100.0],min_num_stns=3):
    from distaz import DistAz
    stns_center=[np.mean([i[0] for i in stns]),np.mean([i[1] for i in stns])]
    xscale=DistAz(stalon=stns_center[0],stalat=stns_center[1],evtlon=stns_center[0]+1,evtlat=stns_center[1]).getDistanceKm()
    scale=[xscale,111.19]
    xy=[stns_center[0]*scale[0],stns_center[1]*scale[1]]
    data_org0={"lonlat":stns_center,"xy":xy,"scale":scale,"rt_clockwise":0.0}
    lenstns=len(stns)
    stns1=stns
    for i in range(lenstns-min_num_stns+1):
        tmp=check_moni_setting(stns_geo=stns1, data_org=data_org0, bias=bias, boxsize=boxsize)
        if tmp['flg']==0:
            return {'flg':0,'flg_str':'all_stns_used','data_org':data_org0,'stns':stns1}
        if tmp['flg']==1:
            tmp1=adjust_data_org(stbias=tmp, data_org=data_org0)
            if tmp1['flg']==0:
                data_org0=tmp1['data_org']
                return {'flg':0,'flg_str':'adjust_org','data_org':data_org0,'stns':stns1}
            else:
                if len(stns1)>min_num_stns:
                   stns1=remove_stn_bymedian(stns1)
                   stns_center=[np.mean([i[0] for i in stns1]),np.mean([i[1] for i in stns1])]
                   xscale=DistAz(stalon=stns_center[0],stalat=stns_center[1],evtlon=stns_center[0]+1,evtlat=stns_center[1]).getDistanceKm()
                   scale=[xscale,111.19]
                   xy=[stns_center[0]*scale[0],stns_center[1]*scale[1]]
                   data_org0={"lonlat":stns_center,"xy":xy,"scale":scale,"rt_clockwise":0.0}
                else:
                   return {'flg':1,'flg_str':'failure'}
    return {'flg':1,'flg_str':'failure'}

def recommend_moni_range_json(datajs_file,modeljs_file,bias=10):
    import json
    with open(datajs_file,'r') as f:
        datajs=json.load(f)
    with open(modeljs_file,'r') as f:
        modeljs=json.load(f)
    with open(datajs['stations_file']) as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        stns=[[float(i[0]),float(i[1])]+i[2:] for i in txtlines]
    station_range=modeljs['station_range']
    boxsize=[station_range[0][1]-station_range[0][0],station_range[1][1]-station_range[1][0]]
    tmp=recommend_moni_range(stns=stns,bias=bias,boxsize=boxsize)
    if tmp['flg']==1:
        print('Error:recommend monitoring range setting failure. You can set it by your own')
    else:
        boxsize_stn=[modeljs['station_range'][0][1],modeljs['station_range'][1][1]]
        boxsize_ev=[modeljs['xyz_range'][0][1]*modeljs['xyz_range'][0][2],
                    modeljs['xyz_range'][1][1]*modeljs['xyz_range'][1][2]]
        image_moni_range(tmp['stns'], tmp['data_org'], boxsize_stn, boxsize_ev)
        print("the suggestion settings and stations could be: ", tmp)
        
def image_moni_range(stns,data_org,boxsize_stn,boxsize_ev):
    import matplotlib.pyplot as plt
    xy=data_org['xy']
    scale=data_org['scale']
    box_stn=[[xy[0]-boxsize_stn[0]/2.0,xy[0]+boxsize_stn[0]/2.0],[xy[1]-boxsize_stn[1]/2.0,xy[1]+boxsize_stn[1]/2.0]]
    box_ev=[[xy[0]-boxsize_ev[0]/2.0,xy[0]+boxsize_ev[0]/2.0],[xy[1]-boxsize_ev[1]/2.0,xy[1]+boxsize_ev[1]/2.0]]
    box_stn_geo=np.array([[box_stn[0][0]/scale[0],box_stn[0][1]/scale[0]],[box_stn[1][0]/scale[1],box_stn[1][1]/scale[1]]])
    box_ev_geo=np.array([[box_ev[0][0]/scale[0],box_ev[0][1]/scale[0]],[box_ev[1][0]/scale[1],box_ev[1][1]/scale[1]]])
    stns_geo=np.array([[float(i[0]),float(i[1])] for i in stns])
    plt.plot([box_stn_geo[0,0],box_stn_geo[0,1]],[box_stn_geo[1,0],box_stn_geo[1,0]],'red')
    plt.plot([box_stn_geo[0,0],box_stn_geo[0,1]],[box_stn_geo[1,1],box_stn_geo[1,1]],'red')
    plt.plot([box_stn_geo[0,1],box_stn_geo[0,1]],[box_stn_geo[1,0],box_stn_geo[1,1]],'red')
    plt.plot([box_stn_geo[0,0],box_stn_geo[0,0]],[box_stn_geo[1,0],box_stn_geo[1,1]],'red')
    
    plt.plot([box_ev_geo[0,0],box_ev_geo[0,1]],[box_ev_geo[1,0],box_ev_geo[1,0]],'red')
    plt.plot([box_ev_geo[0,0],box_ev_geo[0,1]],[box_ev_geo[1,1],box_ev_geo[1,1]],'red')
    plt.plot([box_ev_geo[0,1],box_ev_geo[0,1]],[box_ev_geo[1,0],box_ev_geo[1,1]],'red')
    plt.plot([box_ev_geo[0,0],box_ev_geo[0,0]],[box_ev_geo[1,0],box_ev_geo[1,1]],'red')
    
    plt.plot(stns_geo[:,0],stns_geo[:,1],'.r')
    for i in range(len(stns)):
        plt.text(stns[i][0], stns[i][1], '.'.join(stns[i][3:5]))
    plt.savefig('tmp.png')

def image_moni_range_json(datajs_file,modeljs_file):
    import json
    with open(datajs_file,'r') as f:
        datajs=json.load(f)
    with open(modeljs_file,'r') as f:
        modeljs=json.load(f)
    with open(datajs['stations_file']) as f:
        txtlines=[i.strip().split() for i in f.readlines()];
        stns=[[float(i[0]),float(i[1])]+i[2:] for i in txtlines]
    boxsize_stn=[modeljs['station_range'][0][1],modeljs['station_range'][1][1]]
    boxsize_ev=[modeljs['xyz_range'][0][1]*modeljs['xyz_range'][0][2],
                modeljs['xyz_range'][1][1]*modeljs['xyz_range'][1][2]]
    image_moni_range(stns=stns, data_org=datajs['data_org'], boxsize_stn=boxsize_stn, boxsize_ev=boxsize_ev)
    

if __name__=='__main__':
    import sys
    #import getopt
    arg_len=len(sys.argv)
    if arg_len<2:
        print('usage: moni_check function_name ...')
        print('moni_check check_json datajs_file.json modeljs_file.json 12')
        print('moni_check recommend_setting datajs_file.json modeljs_file.json 12')
        print('moni_check image_moni datajs_file.json modeljs_file.json')
        print('moni_check check_renew_datajson datajs_file.json modeljs_file.json orgsetting.txt 12')
        sys.exit()
    else:
        arg_fun=sys.argv
    print('ss',arg_len,arg_fun)
    if arg_len==5 and arg_fun[1]=='check_json':
        check_js_files(datajs_file=arg_fun[2], modeljs_file=arg_fun[3],bias=float(arg_fun[4]))
    if arg_len==5 and arg_fun[1]=='recommend_setting':
        recommend_moni_range_json(datajs_file=arg_fun[2], modeljs_file=arg_fun[3],bias=float(arg_fun[4]))
    if arg_len==4 and arg_fun[1]=='image_moni':
        image_moni_range_json(datajs_file=arg_fun[2], modeljs_file=arg_fun[3])
    if arg_len==6 and arg_fun[1]=='check_renew_datajson':
        check_renew_datajsfiles(datajs_file=arg_fun[2], modeljs_file=arg_fun[3], orgsetting_file=arg_fun[4], 
                                out_station_file='tmp_station.txt',bias=float(arg_fun[5]))
        
    
    
    
    
    
