#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:57:38 2019

@author: zxiong
"""
import numpy as np
import matplotlib.pyplot as plt
import moviepy.editor as mpy
from moviepy.video.io.bindings import mplfig_to_npimage
import math
from obspy import UTCDateTime
import gc
import sys
sys.path.append('../moni')
from pytool import read_stns,read_pkl,read_json,coordinate2geo1
from extr_catalog import mag_prep,slt_1event_dlm
arg_fun=sys.argv
modeljs_file=arg_fun[1]
datajs_file=arg_fun[2]
modeljs=read_json(modeljs_file)
datajs=read_json(datajs_file)
print(datajs)
data=read_pkl(datajs['out_nam']+'_dlm.pkl')
out_movie_nam=arg_fun[3]
if len(arg_fun)>4:
    movie_seconds=float(arg_fun[4])
    movie_fps=float(arg_fun[5])
else:
    movie_seconds=600.0
    movie_fps=int(len(data)/movie_seconds+0.5)
    if movie_fps==0:
        movie_fps=1

def plot_1event(data,plt,time=np.array(list(range(600)))*0.05,par={'color':'k-'}):
    nr=len(data)
    #nsmp=len(data[0])
    #data[:,:,2]=data[:,:,2]/np.max(np.abs(data[:,:,2]))
    trc_data=[]
    t=time
    for i in range(nr):
        trc=data[i,:]+i+1
        tmp,=plt.plot(t, trc, par['color'])
        #plt.hold()
        trc_data.append(tmp)
    return trc_data
def stn_text(stns,plt):
    for i in range(len(stns)):
        #print(stns[i])
        plt.text(stns[i][0]-0.08,stns[i][1]+0.02,'{}'.format(stns[i][4]),fontsize=10,color='gray')
def plt_set_vaules(data,trc_data):
    nr=len(trc_data)
    #data[:,:,2]=data[:,:,2]/np.max(np.abs(data[:,:,2]))
    for i in range(nr):
        trc=data[i,0:]
        trc_data[i].set_ydata(trc+i+1)
def prep_mags(magidx,imgs,mags,par={'d_pdf':1.0}):
    magr=[0.0, 0.009765625, 1024];
    imgs1=[np.array([0]*len(imgs[0]))]*len(imgs)
    if par['d_pdf']<0.6:
        return {'img':np.array(imgs1),'dmag':magr[1]}
    for i, idx in enumerate(magidx):
        mag=mags[idx][0]
        img=imgs[idx]
        imag=int((mag-magr[0])/magr[1]+0.5)
        pdf=np.max(img)
        imag_max=np.where(img==pdf)[0][0];
        dt=imag_max-imag
        img=list(img)
        img_len=len(img)
        #if dt<0:
        #    dt=-dt
        #    img_tmp=img[img_len-dt:]+img[0:img_len-dt]
        #elif dt>0:
        #    img_tmp=img[dt:img_len]+img[0:dt]
        #imgs1[idx]=np.array(img_tmp)
         imgs1[idx]=np.array(img)
    return {'img':np.array(imgs1),'dmag':magr[1]}

stn=read_stns(datajs['stations_file'])['stns'];stnxy=np.array([[i[0],i[1]] for i in stn]);
stnnams=[i[4] for i in stn]
print(stn)
#r=[[-2667.0,0.5,96],[-4145.0,0.5,192]]
img2d=data[0]['img2d']
#wins=data[0]['d_t0win']
#fmax=data['fmax'][0]
d_t0win=data[0]['d_t0win']
d_pdf=data[0]['d_pdf']
geo=data[0]['geo']
#mag=sio.loadmat('tmp_mag.mat')
#mag=mag['magm'][0]
mag=data[0]['mag_median']
#magrt=data['magrt'][0]
#print(fmax,wins,mag)
#wave=normtrc_all_event(data['data'])
wave=data[0]['wave_test'];print(len(wave[0]))
mag_mags=data[0]['mag_mags']
mag_dtps=data[0]['mag_dtps']
mag_img=data[0]['mag_img']
d_pdf=data[0]['d_pdf']
d_img=data[0]['d_img']

r=modeljs['xyz_range']
coordinate={'data_org':datajs['data_org'],'model_org':modeljs['model_org']}; center_geo=datajs['data_org']['lonlat']

x = np.linspace(r[0][0],r[0][0]+r[0][1]*r[0][2],r[0][2])
y = np.linspace(r[1][0],r[1][0]+r[1][1]*r[1][2],r[1][2])
X,Y = np.meshgrid(x,y)
X,Y = coordinate2geo1(ex=X, ey=Y, coordinate=coordinate)

fig=plt.figure(figsize=(14, 7))
#plt.title()
fig.subplots_adjust(wspace=0.2,hspace=0)
#title_text=plt.title(wins[0],loc='left')
ax1=fig.add_subplot(1,2,1,label='121')
mags_prep=mag_prep(mag_mags, mag_dtps)
mag_img=prep_mags(mags_prep['mags1_idx'], mag_img[:,:,0], mag_mags,par={'d_pdf':d_pdf})
trc_data_mag=plot_1event(mag_img['img'],ax1,time=np.array(range(1024))*mag_img['dmag'],par={'color':'b-'})
ax1.xaxis.set_ticks_position('top');ax1.tick_params(axis='x',colors='blue'); plt.xlim(0,10.000);plt.xticks(np.arange(0, 10.0, 1.0));plt.ylim(0,14);
plt.xticks(fontsize=15);plt.yticks(range(0,14),());plt.xlabel('Magnitude',fontsize=20,color='blue');
ax1.xaxis.set_label_position('top');
ax1_1=fig.add_subplot(1,2,1,label='121_1',frame_on=False)
wave[:,:,2]=wave[:,:,2]/np.max(wave[:,:,2]);
trc_data=plot_1event(wave[:,:,2],ax1_1,time=np.array(range(len(wave[0])))*0.05)
trc_data_detec=plot_1event(d_img[:,0:600,0]+12,ax1_1,time=np.array(range(600))*0.05,par={'color':'r-'})
#ax1_1.plot([])
plt.xlabel('Time (s)',fontsize=15); plt.xlim(0,30);plt.ylim(0,14)
plt.ylabel('Station Name',fontsize=15)
plt.xticks(fontsize=15);plt.yticks(range(1,14),stnnams+['PDF'],fontsize=15)

ax2_1=fig.add_subplot(1,2,2,label='122_1',frame_on=False)
plt.xlim(0,100);plt.ylim(0,100)
plt.plot(3.0,88+5.5,marker='o', mfc='none',mec='k',mew=2,markersize=2+(2-1)*3,clip_on=False)
plt.text(6.0,87.0+5.5,'M 2',fontsize=12)
plt.plot(3.0,88+6-3.5,marker='o', mfc='none',mec='k',mew=2,markersize=2+(4-1)*3,clip_on=False)
plt.text(6.0,87+6-3.5,'M 4',fontsize=12)
plt.plot(3.0,89+6-9,marker='o', mfc='none',mec='k',mew=2,markersize=2+(6-1)*3,clip_on=False)
plt.text(6.0,87+6-9,'M 6',fontsize=12)
#plt.text(12.9,43.33,'           ',fontsize=20)
#plt.text(12.9,43.33,'           ',fontsize=20)
#plt.scatter(0, 0, 'r.',marker='o',c='',edgecolors='w',s=2+(best_mag[1]-1)*3)
#time_text=plt.text(12.85,42.3,'0')
title_text=plt.text(-10,102,'$T_{win}$: '+d_t0win,fontsize=12)
title_text1=plt.text(5,97,'',fontsize=12)
title_text2=plt.text(5,97,'     ',fontsize=12)
title_text3=plt.text(40,102.5,'      ',fontsize=12)
plt.axis('off')

ax2=fig.add_subplot(1,2,2,label='122')
ax2.patch.set_alpha(0)
img_pcolor=plt.pcolormesh(X,Y,img2d.T)
tmp=plt.colorbar()
tmp.ax.tick_params(labelsize=20)
plt.clim([0,1])
plt.xlabel('Longitude ($\mathregular{^o}$)',fontsize=15)
plt.ylabel('Latitude ($\mathregular{^o}$)',fontsize=15)
plt.axis([center_geo[0]-0.53,center_geo[0]+0.48,center_geo[1]-0.5,center_geo[1]+0.5])
plt.gca().set_aspect(111.19/90.0)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
stn_tr,=plt.plot(stnxy[:,0], stnxy[:,1],color='gray', marker='^',markersize=10,linestyle='None')
#stn_tr,=plt.plot(stn[:,0], stn[:,1],'gray^',markersize=10)
stn_text(stn,plt)
loca_p,=plt.plot(0, 0, 'r*',markersize=15)
loca_p1,=plt.plot(0, 0, 'w.',markersize=10)
loca_p2,=plt.plot(0, 0, 'r.',markersize=13)


duration=movie_seconds
#duration=36
#dt=0.05
t_beg=0
dt=(duration-t_beg)/(len(data)-1)
trg_on=[0,9999,9999]
#best_loc=[0.0,9999,9999,UTCDateTime(),0]
best_loc=[]
best_mag=[0.0,9999]
trg_off=0
best_locamag0=[best_loc,best_mag]
dwin=25.0/0.5*dt
idwin=int(25.0/0.5+0.5)
locamagp=[]
save=[]
pltpcolor=plt.pcolormesh(X,Y,img2d.T,vmin=0.0,vmax=1.0)
print('size pltpcolor',img2d.shape,(img2d.T).shape,X.shape,Y.shape)

def make_frame(t):
    global trg_on,best_loc,best_mag,trg_off,best_locamag0,locamagp,save,data
    it=int(t/dt+0.5)
    #print('s186',t,it,data[1])
    img2d=data[it]['img2d']
    wave=data[it]['wave_test']
    d_t0win=data[it]['d_t0win']
    d_pdf=data[it]['d_pdf']
    l_pdf=data[it]['xyz'][3]
    dep=data[it]['xyz'][2]
    d_t0=data[it]['d_t0']
    torg=data[it]['torg']
    mag=data[it]['mag_median']
    geo=data[it]['geo']
    mag_img=data[it]['mag_img']
    mag_mags=data[it]['mag_mags']
    mag_dtps=data[it]['mag_dtps']
    d_pdf=data[it]['d_pdf']
    d_img=data[it]['d_img']
    
    lon=geo[0];lat=geo[1]
    #plt.pcolormesh(X,Y,img2d[it].T,vmin=0.0,vmax=1.0)
    pltpcolor.set_array((img2d.T).ravel())
    wave[:,:,2]=wave[:,:,2]/np.max(wave[:,:,2]);
    plt_set_vaules(wave[:,:,0],trc_data)
    mags_prep=mag_prep(mag_mags, mag_dtps)
    mag_img=prep_mags(mags_prep['mags1_idx'], mag_img[:,:,0], mag_mags,par={'d_pdf':d_pdf});
    plt_set_vaules(mag_img['img'],trc_data_mag)
    plt_set_vaules(d_img[:,0:600,0]+12,trc_data_detec)
    #plt.plot(stn[:,0], stn[:,1], 'r^')
    #stn_tr.set_data(stn[:,0], stn[:,1])
    #time_text.set_text(str(int(it)))
    win=d_t0win[:21]
    title_text.set_text('$T_{win}$: '+win[:21])
    loca_p.set_data(lon, lat)
    if d_pdf>0.7 and l_pdf>0.6:
        #title_text1.set_text('({:.3f},{:.3f},{:.3f},{:.1f})'.format(round(fmax[it],3),round(lon,3),round(lat,3),round(loca[it][2],3)))
        locastr='({:.3f}'.format(round(lon,3))+'$\mathregular{^o}$, '+'{:.3f}'.format(round(lat,3))+'$\mathregular{^o}$, '+'{:.1f} km)'.format(round(dep))
        #loca_p.set_data(lon, lat)
        #picks,rt=datastream.stalta_wins_(wins=[UTCDateTime(wins[it][:21])],thsh=5.0,idx2t=True)
        title_text3.set_text('    ')
        #title_text2.set_text('                ')
        title_text2.set_text('M {:.1f}'.format(round(mag,1))+', '+locastr)
        #if magrt[it]>0.90:
        #   title_text2.set_text('        '+' M '+'{:.1f}'.format(round(mag[it],1))+', ')
        title_text3.set_text('$T_{org}$: '+'{}'.format(torg)[0:21])
        #print(trg_on[0])
        if len(best_loc)==0:
            best_loc.append(data[it])
        else:
            dt0=UTCDateTime(d_t0)-UTCDateTime(best_loc[-1]['d_t0']);dt0=abs(dt0)
            tmp1=best_loc[-1]['xyz'];tmp2=data[it]['xyz']
            dxyz=(tmp1[0]-tmp2[0])**2.0+(tmp1[1]-tmp2[1])**2.0;dxyz=dxyz**0.5
            if dxyz>10 or dt0>3:
                locamagp,=plt.plot(best_loc[-1]['geo'][0], best_loc[-1]['geo'][1],
                          marker='o', mfc='none',mec='w',mew=2,markersize=2+(best_loc[-1]['mag_median']-1)*3)
                best_loc.append(data[it])
            else:
                best_loc[-1]=slt_1event_dlm(dlm1event=[best_loc[-1],data[it]])['dlm_best']
                    
    gc.collect()
    return mplfig_to_npimage(fig)
def test_tmp():
    global save
    for t in np.arange(0.0, 200.0, 0.05):
        print(t)
        make_frame(t)
    print(save)

#test_tmp()
animation = mpy.VideoClip(make_frame, duration=duration)
#animation.write_gif("tmp.gif", fps=1)
animation.write_videofile(filename=out_movie_nam,fps=movie_fps)


