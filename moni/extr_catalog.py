# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 10:26:43 2023

@author: huawei
"""
from obspy import UTCDateTime
import numpy as np
def mean_t0(t0):
    if len(t0)==0:
        return 0
    dt0=[t0[i]-t0[0] for i in range(len(t0))]
    return t0[0]+np.mean(dt0)
    
def mag_prep(mags,dtps):
    mags0=[]
    mags1=[]
    mags2=[]
    mags3=[]
    mags0_idx=[]
    mags1_idx=[]
    mags2_idx=[]
    mags3_idx=[]
    mags_idx=[]
    for i in range(len(mags)):
        if mags[i][1]>0.6 and dtps[i]>2.0:
            mags1.append(mags[i][0])
            mags1_idx.append(i)
        if mags[i][1]>0.6 and dtps[i]>6.0:
            mags2.append(mags[i][0])
            mags2_idx.append(i)
        if mags[i][1]>0.6 and dtps[i]>10.0:
            mags3.append(mags[i][0])
            mags3_idx.append(i)
        if mags[i][1]>0.6 and dtps[i]>0.0:
            mags0.append(mags[i][0])
            mags0_idx.append(i)
    if len(mags3)>4:
        mag_mean=np.mean(mags3)
        mag_median=np.median(mags3)
        mags_idx=mags3_idx
    elif len(mags2)>4:
        mag_mean=np.mean(mags2)
        mag_median=np.median(mags2)
        mags_idx=mags2_idx
    elif len(mags1)>0:
        mag_mean=np.mean(mags1)
        mag_median=np.median(mags1)
        mags_idx=mags1_idx
    elif len(mags0)>0:
        mag_mean=np.mean(mags0)
        mag_median=np.median(mags0)
        mags_idx=mags0_idx
    else:
        return {'mag_mean':-9999,'mag_median':-9999,'mags_idx':mags_idx,
                'mags1_idx':mags1_idx,'mags2_idx':mags2_idx,'mags3_idx':mags3_idx}
    return {'mag_mean':mag_mean,'mag_median':mag_median,'mags_idx':mags_idx,
            'mags1_idx':mags1_idx,'mags2_idx':mags2_idx,'mags3_idx':mags3_idx}

def cata_criterion(dlm,criterion={'d_pdf_min':0.6,'l_pdf_min':0.4,'l_pdf_good':0.65,'dtps_good':[10,23],
                                  'same_event_dxyz_dt0':[10,3]}):
    d_pdf_min=criterion['d_pdf_min']
    l_pdf_min=criterion['l_pdf_min']
    l_pdf_good=criterion['l_pdf_good']
    same_event_dxyz_dt0=criterion['same_event_dxyz_dt0']
    dtps_good=criterion['dtps_good']
    dls=[]
    for i,dl in enumerate(dlm):
        mag_tmp=mag_prep(mags=dl['mag_mags'], dtps=dl['mag_dtps'])
        dl['mag_mean']=mag_tmp['mag_mean'];dl['mag_median']=mag_tmp['mag_median']
        if dl['d_pdf']>=d_pdf_min and dl['xyz'][3]>=l_pdf_min:
            if len(dls)==0:
                dls.append(dl)
                continue
            dl1=dls[-1]
            dt0=UTCDateTime(dl['d_t0'])-UTCDateTime(dl1['d_t0']);dt0=abs(dt0)
            tmp=(dl['xyz'][0]-dl1['xyz'][0])**2+(dl['xyz'][1]-dl1['xyz'][1])**2+\
                (dl['xyz'][2]-dl1['xyz'][2])**2
            dxyz=tmp**0.5
            if dxyz>same_event_dxyz_dt0[0] or dt0>same_event_dxyz_dt0[1]:
                dls.append(dl)
            else:
                dwin=max(dl['mag_dtps'])
                dwin_=max(dl1['mag_dtps'])
                if dwin>dtps_good[0] and dwin<dtps_good[1] and dl['xyz'][3]>l_pdf_good:
                #if dl['xyz'][3]>dl1['xyz'][3]:
                    #if dl1['xyz'][3]<0.6:
                     #   dls[-1]=dl
                     #   continue
                    if dwin_<=dtps_good[0] or dwin_>=dtps_good[1]:
                        dls[-1]=dl
                        continue
                    #if dl1['mag_median']<dl['mag_median']:
                        #dls[-1]=dl
                    #    continue
                if dl1['xyz'][3]<l_pdf_good and dl['xyz'][3]>dl1['xyz'][3]:
                    dls[-1]=dl
                    continue
#def dlm1event_trigeron(dlm1event_pre,dl,criterion={'d_pdf_min':0.6,'l_pdf_min':0.4}):
#    d_pdf_min=criterion['d_pdf_min']
#    l_pdf_min=criterion['l_pdf_min']

def slt_pdf_criterion(dlm,criterion={'d_pdf_min':0.7,'l_pdf_min':0.6}):
    d_pdf_min=criterion['d_pdf_min']
    l_pdf_min=criterion['l_pdf_min']
    dl_pdf_good=[]
    for i, dl in enumerate(dlm):
        if dl['d_pdf']>d_pdf_min and dl['xyz'][3]>l_pdf_min:
            dl_pdf_good.append(dl)
    return {'dl_pdf_good':dl_pdf_good}

def stati_from_dlm(dlm,criterion={'d_pdf_min':0.7,'l_pdf_min':0.6,'l_pdf_good':0.67,'dtps_good':[10,23],
                                  'same_event_dxyz_dt0':[10,3],'min_num_dlm1event':2}):
    same_event_dxyz_dt0=criterion['same_event_dxyz_dt0']
    min_num_dlm1event=criterion['min_num_dlm1event']
    dlm=slt_pdf_criterion(dlm,criterion=criterion)['dl_pdf_good']
    cata=[]; dlm1event=[dlm[0]]; j=0
    for i,dl in enumerate(dlm[1:]):
        dl1=dlm1event[-1]
        dt0=UTCDateTime(dl['d_t0'])-UTCDateTime(dl1['d_t0']);dt0=abs(dt0)
        tmp=(dl['xyz'][0]-dl1['xyz'][0])**2+(dl['xyz'][1]-dl1['xyz'][1])**2+\
            (dl['xyz'][2]-dl1['xyz'][2])**2
        dxyz=tmp**0.5
        if dxyz<same_event_dxyz_dt0[0] and dt0<same_event_dxyz_dt0[1]:
            dlm1event.append(dl)
            j=j+1
        else:
            if len(dlm1event) > min_num_dlm1event:
                dlm_best=slt_1event_dlm(dlm1event, criterion)['dlm_best']
                cata.append(dlm_best)
            dlm1event=[dl]
            j=0
    return {'catalog':cata}

def stati_catalog():
    import pickle
    with open('test_dlm_all.pkl','rb') as f:
        dlm=pickle.load(f)
    dls=stati_from_dlm(dlm)['catalog']
    print(dls)
    with open('test_dlm_cat.txt','w') as f:
        for i,dl in enumerate(dls):
            mags_s=['[{:.3f},{:.3f}]'.format(i[0],i[1]) for i in dl['mag_mags']]
            dtps_s=['{:.3f}'.format(i) for i in dl['mag_dtps']]
            mag_tmp=mag_prep(mags=dl['mag_mags'], dtps=dl['mag_dtps'])
            dl['mag_mean']=mag_tmp['mag_mean'];dl['mag_median']=mag_tmp['mag_median']
            f.write('{} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {} {:.6f} {:.6f} {} {} {}\n'.format(
                dl['d_t0'],dl['d_pdf'],dl['xyz'][0],dl['xyz'][1],dl['xyz'][2],dl['xyz'][3],dl['geo'][0],dl['geo'][1],
                dl['d_t0win'],dl['mag_median'],dl['mag_mean'],mags_s,dtps_s,dl['reflag']))

def slt_1event_dlm(dlm1event,criterion={'d_pdf_min':0.7,'l_pdf_min':0.6,'l_pdf_good':0.67,'dtps_good':[10,23],
                                  'same_event_dxyz_dt0':[10,3],'min_num_dlm1event':2}):
    l_pdf_good=criterion['l_pdf_good']
    dtps_good=criterion['dtps_good']
    dlm_good=[]
    dlm_pdf=dlm1event[0]
    for i,dl in enumerate(dlm1event):
        dwin=max(dl['mag_dtps'])
        if dwin>dtps_good[0] and dwin<dtps_good[1] and dl['xyz'][3]>l_pdf_good:
            dlm_good.append(dl)
        if dlm_pdf['xyz'][3]<dl['xyz'][3]:
            dlm_pdf=dl
    if len(dlm_good)==0:
        dlm_best=dlm_pdf
        dlm_best['reflag']='pdf_max'
    else:
        dlm_tmp=dlm_good[0]
        for i,dl in enumerate(dlm_good):
            if dlm_tmp['xyz'][3]<dl['xyz'][3]:
                dlm_tmp=dl
        dlm_best=dlm_tmp
        dlm_best['reflag']='overall'
    return {'dlm_best':dlm_best,'dlm_good':dlm_good,'dlm_pdf':dlm_pdf}
        

def dlm1event_criterion(dlm1event,criterion={'d_pdf_min':0.6,'l_pdf_min':0.4,'l_pdf_good':0.65,'dtps_good':[10,23]}):
    d_pdf_min=criterion['d_pdf_min']
    l_pdf_min=criterion['l_pdf_min']
    l_pdf_good=criterion['l_pdf_good']
    dtps_good=criterion['dtps_good']
   # dt0=criterion['dt0']
    dls=[]
    dl_good=[]
    dl_best=[]
    dl_pdf_best=[]
    for i,dl in enumerate(dlm1event):
        mag_tmp=mag_prep(mags=dl['mag_mags'], dtps=dl['mag_dtps'])
        dl['mag_mean']=mag_tmp['mag_mean'];dl['mag_median']=mag_tmp['mag_median']
        if dl['d_pdf']>=d_pdf_min and dl['xyz'][3]>=l_pdf_min:
            dls.append(dl)
            dwin=max(dl['mag_dtps'])
            if dwin>dtps_good[0] and dwin<dtps_good[1] and dl['xyz'][3]>l_pdf_good:
                dl_good.append(dl)
                if len(dl_best)>0:
                    if dl['xyz'][3]>dl_best['xyz'][3]:
                        dl_best=dl
                else:
                    dl_best=dl
        if dl_pdf_best!=[]:
            if dl_pdf_best['xyz'][3]<dl['xyz'][3]:
                dl_pdf_best=dl
        else:
            dl_pdf_best=dl
    if dl_best==[]:
        dl_best=dl_pdf_best        
    return {'dl_pdf_good':dls,'dl_best':dl_best,'dl_good':dl_good}

def extr_catalog():
    import pickle
    from obspy import UTCDateTime
    with open('test_dlm.pkl','rb') as f:
        dlm=pickle.load(f)
    #print(dlm[0])
    j=0
    j1=0
    trg_on=0
    trg_off=1
    minpdf=0.6
    minpdf1=0.4
    dt0=999.0
    xyz=[]
    mag=[]
    dls=[]
    for i,dl in enumerate(dlm):
        mag_tmp=mag_prep(mags=dl['mag_mags'], dtps=dl['mag_dtps'])
        dl['mag_mean']=mag_tmp['mag_mean'];dl['mag_median']=mag_tmp['mag_median']
        if dl['d_pdf']>=minpdf and dl['xyz'][3]>=minpdf1:
            trg_on=1
            trg_off=0
            j=j+1
            j1=0
            #dls.append(dl)
        if (dl['d_pdf']<minpdf and j1>2) or j>40:
            trg_off=1
            trg_on=0
            j=0
            j1=j1+1
            #dls=[]
        if dl['d_pdf']>=minpdf and dl['xyz'][3]>=minpdf1:
            if len(dls)==0:
                dls.append(dl)
                continue
            dl1=dls[-1]
            dt0=UTCDateTime(dl['d_t0'])-UTCDateTime(dl1['d_t0']);dt0=abs(dt0)
            tmp=(dl['xyz'][0]-dl1['xyz'][0])**2+(dl['xyz'][1]-dl1['xyz'][1])**2+\
                (dl['xyz'][2]-dl1['xyz'][2])**2
            dxyz=tmp**0.5
            if dxyz>10 or dt0>3:
                dls.append(dl)
            else:
                dwin=max(dl['mag_dtps'])
                dwin_=max(dl1['mag_dtps'])
                if dwin>10 and dwin<23 and dl['xyz'][3]>0.65:
                #if dl['xyz'][3]>dl1['xyz'][3]:
                    #if dl1['xyz'][3]<0.6:
                     #   dls[-1]=dl
                     #   continue
                    if dwin_<=10 or dwin_>=23:
                        dls[-1]=dl
                        continue
                    #if dl1['mag_median']<dl['mag_median']:
                        #dls[-1]=dl
                    #    continue
                if dl1['xyz'][3]<0.65 and dl['xyz'][3]>dl1['xyz'][3]:
                    dls[-1]=dl
                    continue
                    
                        
                #if dwin>8 and dwin<23:
                    
    print(dls)
    with open('test_dlm_cat.txt','w') as f:
        for i,dl in enumerate(dls):
            mags_s=['[{:.3f},{:.3f}]'.format(i[0],i[1]) for i in dl['mag_mags']]
            dtps_s=['{:.3f}'.format(i) for i in dl['mag_dtps']]
            f.write('{} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {} {:.6f} {:.6f} {} {}\n'.format(
                dl['d_t0'],dl['d_pdf'],dl['xyz'][0],dl['xyz'][1],dl['xyz'][2],dl['xyz'][3],dl['geo'][0],dl['geo'][1],
                dl['d_t0win'],dl['mag_median'],dl['mag_mean'],mags_s,dtps_s))
        
from pytool import read_pkl,read_json,coordinate2geo1,write_sacs,align_mag_imgs
def extr_sac_from_dlm(pklnam='D:/loca_aug/demo/test_dlm_ridgem6.pkl',idx=0,datajs='D:/loca_aug/demo/data_info.json',
                      modeljs='D:/loca_aug/demo/model_info.json',com=2):
    data=read_pkl(pklnam)
    modeljs=read_json(modeljs)
    r=modeljs['xyz_range']
    datajs=read_json(datajs)
    coordinate={'data_org':datajs['data_org'],'model_org':modeljs['model_org']}; rmag=modeljs['mag_range']
    img2d=data[idx]['img2d']
    d_t0win=data[idx]['d_t0win']
    d_pdf=data[idx]['d_pdf']
    geo=data[idx]['geo']
    mag=data[idx]['mag_median']
    wave=data[idx]['wave_test'][:,:,com];print(len(wave[0]))
    mag_mags=data[idx]['mag_mags']
    mag_dtps=data[idx]['mag_dtps']
    mag_img=data[idx]['mag_img']
    d_pdf=data[idx]['d_pdf']
    d_img=data[idx]['d_img'][0,:,0]
    write_sacs(data=[d_img], dt=0.05, nam='pdf')
    write_sacs(data=[i/np.max(i) for i in wave], dt=0.05, nam='wave')
    print(np.array(mag_img[:,:,0]).shape)
    mag_img=align_mag_imgs(mag_imgs=np.array(mag_img[:,:,0]), imags=[int(i/rmag[1]) for i in np.array(mag_mags)[:,0]])['mag_imgs']
    print(np.array(mag_img[:,:]).shape)
    write_sacs(data=mag_img,dt=rmag[1],nam='mag')
    
    x = np.linspace(r[0][0],r[0][0]+r[0][1]*r[0][2],r[0][2])
    y = np.linspace(r[1][0],r[1][0]+r[1][1]*r[1][2],r[1][2])
    #X,Y = np.meshgrid(x,y)
    #print('x',x)
    with open("img2d.txt", "w") as f:
      for i in range(len(x)):
          for j in range(len(y)):
              x1,y1=coordinate2geo1(ex=x[i], ey=y[j], coordinate=coordinate)
              f.write("{:.6f} {:.6f} {:.6f}\n".format(x1,y1,img2d[i][j]+0.0001))
    print(d_t0win,geo,mag_mags,mag_dtps,d_pdf,d_img,np.where(d_img==np.max(d_img)),len(d_img),mag)
    
            
            
        
                
            
    
    
if __name__ == '__main__':
    #stati_catalog()
    extr_sac_from_dlm(idx=2324)
    #extr_sac_from_dlm(idx=2432)
#extr_catalog()

    
