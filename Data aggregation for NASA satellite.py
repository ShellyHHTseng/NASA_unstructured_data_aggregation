import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyhdf.SD import SD, SDC
from os import listdir
from os.path import isfile, join
from datetime import datetime, timezone
import time
import calendar
import pprint
from mpl_toolkits.basemap import Basemap
import glob
import scipy.ndimage as ndimage
import math
from pyhdf import HDF, VS, V
import pandas as pd
import datetime as dt


t=time.time()

# uncomment the following code if there is output.txt in the directory.
#os.remove("output.txt")

outfile=open('./output.txt', 'a')
#outfile.write('# t    s(t)\n')  # write table header


# you need to change the file path here. 
MO_path='/home/disk/eos12/ktoandy/MODIS_L2/2008/'
MO_name = [os.path.basename(x) for x in glob.glob(MO_path+'MYD06*.hdf')]

CAL_path='/home/disk/sunburn2/hhtseng/CALIPSO_CLay_v4/'
CAL_name = [os.path.basename(x) for x in glob.glob(CAL_path+'CAL*.2008-*.hdf')]

CLO_path1='/home/disk/sunburn2/hhtseng/CloudSat/2B-CWC-RO.P_R04/2008/'
CLO_name = [os.path.basename(x) for x in glob.glob(CLO_path1+'2008*.hdf')]

CLO_path2='/home/disk/sunburn2/hhtseng/CloudSat/2B-GEOPROF.P_R04/2008/'
CLO_name_mask = [os.path.basename(x) for x in glob.glob(CLO_path2+'2008*GEOPROF*.hdf')]
 

###########
##### find the possible collocated CALIPSO and CLOUDSAT based on the filename
###########

o_date=dt.datetime(2008,1,1,0,0,0) ###
second_records_since_calipso=np.zeros(len(CAL_name))
second_records_since_cloudsat=np.zeros(len(CLO_name))

for i in range(0,len(CAL_name)-1): 
    CALt=CAL_name[i]
    yy=CALt[35:39]
    mm=CALt[40:42]
    dd=CALt[43:45]
    hr=CALt[46:48]
    miniute=CALt[49:51]    
    ss=CALt[52:54]
    current_date=dt.datetime(int(yy),int(mm),int(dd),int(hr),int(miniute),int(ss))
    second_records_since_calipso[i]=(current_date-o_date).total_seconds()

for i in range(0,len(CLO_name)):
    CALt=CLO_name[i]
    yy=CALt[0:4]
    
    day_in_year=CALt[4:7]
    zz=dt.datetime(int(yy), 1, 1) + dt.timedelta(int(day_in_year) - 1)        
    
    mm=zz.month
    dd=zz.day
    hr=CALt[7:9]
    miniute=CALt[9:11]    
    ss=CALt[11:13]
    current_date=dt.datetime(int(yy),int(mm),int(dd),int(hr),int(miniute),int(ss))
    second_records_since_cloudsat[i]=(current_date-o_date).total_seconds()
    
cloudsat_index_for_calipso=np.full((len(CAL_name),2),np.nan) 
    
for i in range(0,len(CAL_name)-1): 
    head=second_records_since_calipso[i]
    tail=second_records_since_calipso[i+1]
    candi=[j for j, e in enumerate(second_records_since_cloudsat) if e < head]
    
    cloudsat_index_for_calipso[i,0]=max(candi)
  
    if max(candi)<np.shape(second_records_since_cloudsat)[0]-1 and second_records_since_cloudsat[max(candi)+1]< tail:
        cloudsat_index_for_calipso[i,1]=max(candi)+1

#%%
############
# CALIPSO #
############
zeroarray=np.zeros(shape=(1,10)) 
zeroarray=' '.join(map(str, zeroarray))

#for i in range(6115,len(CAL_name)-1):  ###6115

for i in range(0,869):  ###6115
    print(i)

    pathname=CAL_path+CAL_name[i]
    file = SD(pathname,SDC.READ) #print(file.datasets()) #file.info()

    sds_obj = file.select('Layer_Base_Altitude') # select sds
    CAL_CBH = np.array(sds_obj.get())   # get sds data
    CAL_CBH[CAL_CBH<0]=np.nan
    CAL_Clear=np.nansum(CAL_CBH, axis=1)
    
    sds_obj = file.select('Layer_Top_Altitude') # select sds
    CAL_CTH = np.array(sds_obj.get())   # get sds data
    CAL_CTH[CAL_CTH<0]=np.nan
    
    sds_obj = file.select('Longitude') # select sds
    CAL_Long = np.array(sds_obj.get())  # get sds data
    CAL_Long = CAL_Long[:,1]

    sds_obj = file.select('Latitude') # select sds
    CAL_Lat = np.array(sds_obj.get()) # get sds data
    CAL_Lat = CAL_Lat[:,1]

    sds_obj = file.select('Profile_Time')
    CAL_Time = np.array(sds_obj.get())
    CAL_Time = CAL_Time[:,1]
    
    sds_obj = file.select('Opacity_Flag')
    CAL_Opa = np.array(sds_obj.get())
    CAL_Opa[CAL_Opa==99]=0
    CAL_Opaque=np.sum(CAL_Opa, axis=1)

    sds_obj = file.select('Feature_Optical_Depth_532')
    CAL_OPTDEP = np.array(sds_obj.get()) # get sds data
    CAL_OPTDEP[CAL_OPTDEP==-9999]=0 
    CAL_CTH=CAL_CTH[abs(CAL_Lat)<=30]
    CAL_CBH=CAL_CBH[abs(CAL_Lat)<=30]
    CAL_Time=CAL_Time[abs(CAL_Lat)<=30]
    CAL_Long=CAL_Long[abs(CAL_Lat)<=30]
    CAL_Opaque=CAL_Opaque[abs(CAL_Lat)<=30]
    CAL_Opa=CAL_Opa[abs(CAL_Lat)<=30]
    CAL_OPTDEP=CAL_OPTDEP[abs(CAL_Lat)<=30]
    CAL_Clear=CAL_Clear[abs(CAL_Lat)<=30]

    CAL_Lat=CAL_Lat[abs(CAL_Lat)<=30]

    CALt=CAL_name[i]
    yy=CALt[35:39]
    mm=CALt[40:42]
    dd=CALt[43:45]
    hr=CALt[46:48]
    miniute=CALt[49:51]
    timeymdhrminiute=yy+mm+dd+hr+miniute
    # convert to day of the year (MODIS)
    ss=yy+'.'+mm+'.'+dd
    fmt='%Y.%m.%d'
    dt=(datetime.strptime(ss,fmt)).timetuple()
    yod=str(dt.tm_yday) ###   
    
    if len(yod)==1 and len(str(int(hr)))==1:
        calitimem=yy+'00'+yod+'.0'+str(int(hr))

    elif len(yod)==1 and len(str(int(hr)))==2:
        calitimem=yy+'00'+yod+'.'+str(int(hr))

    elif len(yod)==2 and len(str(int(hr)))==1:
        calitimem=yy+'0'+yod+'.0'+str(int(hr))
        
    elif len(yod)==2 and len(str(int(hr)))==2:
        calitimem=yy+'0'+yod+'.'+str(int(hr))
        
    elif len(yod)==3 and len(str(int(hr)))==1:
        calitimem=yy+yod+'.0'+str(int(hr))
   
    elif len(yod)==3 and len(str(int(hr)))==2:
        calitimem=yy+yod+'.'+str(int(hr))
        
### find MODIS collocate ###
    if i<len(CAL_name)-1 and int(hr)<23: 

        CALt2=CAL_name[i+1]
        yy2=CALt2[35:39]
        mm2=CALt2[40:42]
        dd2=CALt2[43:45]
        hr2=CALt2[46:48]
        miniute2=CALt2[49:51]
        
        # convert to day of the year (MODIS)
        ss=yy2+'.'+mm2+'.'+dd2
        fmt='%Y.%m.%d'
        dt=(datetime.strptime(ss,fmt)).timetuple()
        yod2=str(dt.tm_yday) ###

        if len(yod2)==1 and len(str(int(hr2)))==1:
            calitimem2=yy2+'00'+yod2+'.0'+str(int(hr2))

        elif len(yod2)==1 and len(str(int(hr2)))==2:
            calitimem2=yy2+'00'+yod2+'.'+str(int(hr2))

        elif len(yod2)==2 and len(str(int(hr2)))==1:
            calitimem2=yy2+'0'+yod2+'.0'+str(int(hr2))
        
        elif len(yod2)==2 and len(str(int(hr2)))==2:
            calitimem2=yy2+'0'+yod2+'.'+str(int(hr2))
        
        elif len(yod2)==3 and len(str(int(hr2)))==1:
            calitimem2=yy2+yod2+'.0'+str(int(hr2))
   
        elif len(yod2)==3 and len(str(int(hr2)))==2:
            calitimem2=yy2+yod2+'.'+str(int(hr2)) 

        MODISall_COD=[]
        MODISall_time=[]
        MODISall_lon=[]
        MODISall_lat=[]
        MO_name_coll=[x for x in MO_name if calitimem in x or calitimem2 in x]    
        
        x1=int(hr+miniute)
        x2=int(hr2+miniute2)
        MO_col=[]
        
        if x1==x2:
            continue
        
        x=0 # in order to also include the one just before x1 and the one just after x2
        for q in range(0,len(MO_name_coll)):
            MOMO = MO_name_coll[q]
            MOname = int(MOMO[18:22])
                        
            if (MOname>=x1 and MOname<=x2) and x==0:
                MO_col.append(MO_name_coll[q-1])                
                MO_col.append(MO_name_coll[q])
                x=1
            elif (MOname>=x1 and MOname<=x2) and x==1:
                MO_col.append(MO_name_coll[q])
                xx=q
        if xx+1<=q:
            MO_col.append(MO_name_coll[xx+1])
        
    elif i<len(CAL_name)-1 and int(hr)==23:
        
        CALt2=CAL_name[i+1]
        yy2=CALt2[35:39]
        mm2=CALt2[40:42]
        dd2=CALt2[43:45]
        hr2=CALt2[46:48]
        miniute2=CALt2[49:51]       
        
        # convert to day of the year (MODIS)
        ss=yy2+'.'+mm2+'.'+dd2
        fmt='%Y.%m.%d'
        dt=(datetime.strptime(ss,fmt)).timetuple()
        yod2=str(dt.tm_yday) ###

        if len(yod2)==1 and len(str(int(hr2)))==1:
            calitimem2=yy2+'00'+yod2+'.0'+str(int(hr2))

        elif len(yod2)==1 and len(str(int(hr2)))==2:
            calitimem2=yy2+'00'+yod2+'.'+str(int(hr2))

        elif len(yod2)==2 and len(str(int(hr2)))==1:
            calitimem2=yy2+'0'+yod2+'.0'+str(int(hr2))
        
        elif len(yod2)==2 and len(str(int(hr2)))==2:
            calitimem2=yy2+'0'+yod2+'.'+str(int(hr2))
        
        elif len(yod2)==3 and len(str(int(hr2)))==1:
            calitimem2=yy2+yod2+'.0'+str(int(hr2))
   
        elif len(yod2)==3 and len(str(int(hr2)))==2:
            calitimem2=yy2+yod2+'.'+str(int(hr2))         
        
#MODIS
        MODISall_COD=[]
        MODISall_time=[]
        MODISall_lon=[]
        MODISall_lat=[]
        MO_name_coll=[x for x in MO_name if calitimem in x or calitimem2 in x]    
        
        x1=int(hr+miniute)
        x2=int(hr2+miniute2)
        MO_col=[]
        
        x=0 # in order to also include the one just before x1 and the one just after x2
        for q in range(0,len(MO_name_coll)):
            MOMO = MO_name_coll[q]
            MOname = int(MOMO[18:22])
            
            if (MOname>=x1 or MOname<=x2) and x==0:
                MO_col.append(MO_name_coll[q-1])                
                MO_col.append(MO_name_coll[q])
                x=1
                
            elif (MOname>=x1 or MOname<=x2) and x==1:
                MO_col.append(MO_name_coll[q])
                xx=q
        if xx+1<=q:
            MO_col.append(MO_name_coll[xx+1])        
        
#########   
#CLOUDSAT
#########
    CLOUDSAT_maxtop=[]
    CLOUDSAT_Height=[]
    CLOUDSAT_iwc=[]
    CLOUDSAT_lwc=[]
    CLOUDSAT_lon=[]
    CLOUDSAT_lat=[]
    CLOUDSAT_time=[]
    Cloudsat_col=[]
    cloudsatindex=cloudsat_index_for_calipso[i]
    cloudsatindex=cloudsatindex[cloudsatindex>=0]
    Cloudsat_col=[]
    Cloudsat_col_mask=[]
    
    if len(cloudsatindex)==1:
        ci=int(cloudsatindex[0])
        Cloudsat_col=[x for x in CLO_name if CLO_name[ci] in x] 
        Cloudsat_col_mask=[x for x in CLO_name_mask if CLO_name_mask[ci] in x]
        
    elif len(cloudsatindex)==2:
        for qq in range(2):
            ci=int(cloudsatindex[qq])
            Cloudsat_col.append(CLO_name[ci])         
            Cloudsat_col_mask.append(CLO_name_mask[ci])
    
### Read ####
### Cloudsat ####    
# CWC
    if (Cloudsat_col and Cloudsat_col_mask) and (len(Cloudsat_col)==len(Cloudsat_col_mask)):
        
        for k in range(0,len(Cloudsat_col)): 
            
            pathname=CLO_path1+Cloudsat_col[k]
            file = SD(pathname, SDC.READ)
            
            # Read attributes
            h = HDF.HDF(pathname)
            vs = h.vstart()
            xid = vs.find('Latitude')
            latid = vs.attach(xid)
            latid.setfields('Latitude')
            nrecs, _, _, _, _ = latid.inquire()
            latitude = latid.read(nRec=nrecs)
            latid.detach()
            CLO_latitude = np.array(latitude)
            
            lonid = vs.attach(vs.find('Longitude'))
            lonid.setfields('Longitude')
            nrecs, _, _, _, _ = lonid.inquire()
            longitude = lonid.read(nRec=nrecs)
            lonid.detach()
            CLO_longitude = np.array(longitude)
            
            timeid = vs.attach(vs.find('TAI_start'))
            timeid.setfields('TAI_start')
            nrecs, _, _, _, _ = timeid.inquire()
            tai_time = np.array(timeid.read(nRec=nrecs))
            #units_t = timeid.attr('units').get()
            #longname_t = timeid.attr('long_name').get()
            timeid.detach()
            
            timeid = vs.attach(vs.find('Profile_time'))
            timeid.setfields('Profile_time')
            nrecs, _, _, _, _ = timeid.inquire()
            profile_time = np.array(timeid.read(nRec=nrecs))
            timeid.detach()
            
            CLO_Time=tai_time+profile_time
         
            sds_obj = file.select('Height') # m
            Height = (sds_obj.get()).astype(float)
            Height[Height==-9999]=np.nan       
            sds_obj = file.select('RO_ice_water_content') # mg m^{-3} 
            iwc = (sds_obj.get()).astype(float)
            iwc[iwc<0]=np.nan
            
            sds_obj = file.select('RO_liq_water_content') # mg m^{-3} 
            lwc = (sds_obj.get()).astype(float)
            lwc[lwc<0]=np.nan
            
# GEOPROF
            pathname=CLO_path2+Cloudsat_col_mask[k]
            file = SD(pathname,SDC.READ) 
            #AA=(file.datasets()) #file.info()

            sds_obj = file.select('CPR_Cloud_mask') # m
            cpr = (sds_obj.get()).astype(float)
            Height[cpr<=25]=0
            
            # broadcast latitude to match Height, iwc, lwc
            CLO_latitude2=np.tile(CLO_latitude,(Height.shape[1]))  
            
            CLO_longitude=CLO_longitude[abs(CLO_latitude)<=30]
            CLO_Time=CLO_Time[abs(CLO_latitude)<=30]
            
            Height2=Height[abs(CLO_latitude2)<=30]
            Height=Height2.reshape(int(Height2.shape[0]/CLO_latitude2.shape[1]),CLO_latitude2.shape[1])/1000
            
            iwc2=iwc[abs(CLO_latitude2)<=30]
            iwc=iwc2.reshape(int(iwc2.shape[0]/CLO_latitude2.shape[1]),CLO_latitude2.shape[1])

            lwc2=lwc[abs(CLO_latitude2)<=30]
            lwc=lwc2.reshape(int(lwc2.shape[0]/CLO_latitude2.shape[1]),CLO_latitude2.shape[1])

            CLO_latitude=CLO_latitude[abs(CLO_latitude)<=30]      

            MaxHeight=np.max(Height,axis=1)
            
            CLOUDSAT_lat.append(CLO_latitude)
            CLOUDSAT_lon.append(CLO_longitude)
            CLOUDSAT_time.append(CLO_Time)
            CLOUDSAT_Height.append(Height)
            CLOUDSAT_iwc.append(iwc)
            CLOUDSAT_lwc.append(lwc)
            CLOUDSAT_maxtop.append(MaxHeight)
            
        CLOUDSAT_Height=np.concatenate(CLOUDSAT_Height,axis=0)
        CLOUDSAT_iwc=np.concatenate(CLOUDSAT_iwc,axis=0)
        CLOUDSAT_lwc=np.concatenate(CLOUDSAT_lwc,axis=0)            
                            
    if MO_col:
        for k in range(0, len(MO_col)):
            pathname=MO_path+MO_col[k]
            file = SD(pathname, SDC.READ)
            
            sds_obj = file.select('Cloud_Optical_Thickness') # select sds
            CCOD = sds_obj.get() # get sds data
            
            for key, value in sds_obj.attributes().items():
                if key == 'add_offset':
                    add_offset = value  
                if key == 'scale_factor':
                    scale_factor = value
            MO_COD = (CCOD - add_offset) * scale_factor
            MODISall_COD.append(MO_COD.ravel())
            
            sds_obj = file.select('Scan_Start_Time')
            MO_5kmtime = np.array(sds_obj.get())   # get sds data
            MO_1kmtime = ndimage.zoom(MO_5kmtime, (MO_COD.shape[0]/MO_5kmtime.shape[0],MO_COD.shape[1]/MO_5kmtime.shape[1]))
            MODISall_time.append(MO_1kmtime.ravel())
            
            sds_obj = file.select('Longitude') # select sds
            MO_5kmlon = np.array(sds_obj.get())  # get sds data
            MO_5kmlon_=np.where(MO_5kmlon<=0, MO_5kmlon+360, MO_5kmlon)
            MO_1kmlon = ndimage.zoom(MO_5kmlon_, (MO_COD.shape[0]/MO_5kmlon_.shape[0],MO_COD.shape[1]/MO_5kmlon_.shape[1]))
            MO_1kmlon_=np.where(MO_1kmlon>=180, MO_1kmlon-360, MO_1kmlon)
            MODISall_lon.append(MO_1kmlon_.ravel())
            
            sds_obj = file.select('Latitude') # select sds
            MO_5kmlat = np.array(sds_obj.get())  # get sds data
            MO_1kmlat = ndimage.zoom(MO_5kmlat, (MO_COD.shape[0]/MO_5kmlat.shape[0],MO_COD.shape[1]/MO_5kmlat.shape[1]))
            MODISall_lat.append(MO_1kmlat.ravel())
        
        MODISall_COD=np.hstack(MODISall_COD)
        MODISall_time=np.hstack(MODISall_time)
        MODISall_lon=np.hstack(MODISall_lon)
        MODISall_lat=np.hstack(MODISall_lat)
        
        MODISall_COD=MODISall_COD[abs(MODISall_lat)<=30]
        MODISall_lon=MODISall_lon[abs(MODISall_lat)<=30]
        MODISall_time=MODISall_time[abs(MODISall_lat)<=30]
        MODISall_lat=MODISall_lat[abs(MODISall_lat)<=30]      
                
### 3 cases        
    for j in range(0,CAL_Long.shape[0]):
        
        #clear sky
        if CAL_Clear[j]==0: 
            clear_cloudy_opaque=str(1)
            latitu=str(round(CAL_Lat[j],1))
            longit=str(round(CAL_Long[j],1))
            CT=zeroarray
            CB=zeroarray
            COD=zeroarray
            outfile.write('%s %s %s %s \n' % (clear_cloudy_opaque, latitu, longit, timeymdhrminiute))
            
        # cloudy and all transparent    
        elif CAL_Clear[j]!=0 and CAL_Opaque[j]==0: 
            clear_cloudy_opaque=str(2)
            latitu=str(round(CAL_Lat[j],1))
            longit=str(round(CAL_Long[j],1))
            CT=np.round(CAL_CTH[j].reshape(1,10),2)
            CB=np.round(CAL_CBH[j].reshape(1,10),2)
            COD=np.round(CAL_OPTDEP[j].reshape(1,10),2)
            outfile.write('%s %s %s %s %s %s %s \n' % (clear_cloudy_opaque, latitu, 
                                                       longit,timeymdhrminiute, CT, CB, COD))
            
        # cloudy and opaque    
        elif CAL_Clear[j]!=0 and CAL_Opaque[j]==1: 
            clear_cloudy_opaque=str(3)
            latitu=str(round(CAL_Lat[j],1))
            longit=str(round(CAL_Long[j],1))  
            
            if len(CLOUDSAT_lat)>0:
                    
                CLOUDSAT_lat_=np.hstack(CLOUDSAT_lat)
                CLOUDSAT_lon_=np.hstack(CLOUDSAT_lon)
                CLOUDSAT_maxtop_=np.hstack(CLOUDSAT_maxtop)
                CLOUDSAT_time_=np.hstack(CLOUDSAT_time)
                Dis=np.sqrt((abs(CLOUDSAT_lon_-CAL_Long[j]))**2+(abs(CLOUDSAT_lat_-CAL_Lat[j]))**2)
                Timediff=abs(CLOUDSAT_time_-CAL_Time[j])
                coll_condition = (Dis<=2.5/111) & (Timediff<=30) #### 2.5km   
                
                CLOUDSAT_maxtop_coll=CLOUDSAT_maxtop_[coll_condition]
                
                df=pd.Series(list(CAL_CBH[j]))
                calbaseind=df.last_valid_index()
                
                if any(i>= CAL_CBH[j,calbaseind] for i in CLOUDSAT_maxtop_coll):
                    
                    #CLOUDSAT_Height_=np.hstack(CLOUDSAT_Height)

                    coll_condition_=np.tile(coll_condition,CLOUDSAT_Height.shape[1])
                    coll_condition2=np.reshape(coll_condition_,(CLOUDSAT_Height.shape[1],CLOUDSAT_Height.shape[0])).T     
                    #CLOUDSAT_iwc_=np.hstack(CLOUDSAT_iwc)
                    #CLOUDSAT_lwc_=np.hstack(CLOUDSAT_lwc)
                    
                    CLOUDSAT_Height_coll=CLOUDSAT_Height[coll_condition2]
                    CLOUDSAT_Height_coll_=np.reshape(CLOUDSAT_Height_coll,(int(CLOUDSAT_Height_coll.shape[0]/CLOUDSAT_Height.shape[1])
                    ,CLOUDSAT_Height.shape[1]))
                    
                    CLOUDSAT_iwc_coll=CLOUDSAT_iwc[coll_condition2]
                    CLOUDSAT_iwc_coll_=np.reshape(CLOUDSAT_iwc_coll,(int(CLOUDSAT_iwc_coll.shape[0]/CLOUDSAT_Height.shape[1])
                    ,CLOUDSAT_Height.shape[1]))
                    
                    CLOUDSAT_lwc_coll=CLOUDSAT_lwc[coll_condition2]
                    CLOUDSAT_lwc_coll_=np.reshape(CLOUDSAT_lwc_coll,(int(CLOUDSAT_lwc_coll.shape[0]/CLOUDSAT_Height.shape[1])
                    ,CLOUDSAT_Height.shape[1]))   
                    
                    ### if the highest cloud top that higher than CALIPSO opaque cloud base 
                    ### have more than 1 pixel... find the closest one with CALIPSO pixel
                    if sum(CLOUDSAT_maxtop_coll>CAL_CBH[j,calbaseind])>1: 
                        
                        # filter out the collocate CLOUDSAT pixel that no cloud top is higher than
                        # CALIPSO opaque cloud base
                        oneind=np.where(CLOUDSAT_maxtop_coll>0)[0]
                        Coll_Height=CLOUDSAT_Height_coll_[oneind,:] 
                        Coll_IWC=CLOUDSAT_iwc_coll_[oneind,:]
                        Coll_LWC=CLOUDSAT_lwc_coll_[oneind,:]  
                        CLOUDSAT_maxtop_coll_=CLOUDSAT_maxtop_coll[oneind]
                       
                        ###                             
                        Dis_coll=Dis[coll_condition]
                        Dis_coll=Dis_coll[oneind]
                        mindisindex=np.argmin(Dis_coll)
                        Coll_Height=Coll_Height[mindisindex,:]
                        Coll_IWC=Coll_IWC[mindisindex,:]
                        Coll_LWC=Coll_LWC[mindisindex,:]
                        ###
                        
                        Coll_Height[Coll_Height<=0]=np.nan
                        df=pd.Series(list(Coll_Height))
                        ind=df.first_valid_index()
                        
                        ### MODIS
                        Timediff=abs(MODISall_time-CAL_Time[j])
                        Dis=np.sqrt((abs(MODISall_lon-CAL_Long[j]))**2+(abs(MODISall_lat-CAL_Lat[j]))**2)
                        coll_condition_m = (Timediff<=100) & (Timediff>=62) & (Dis<=1/111) #### 1km
                        MODIS_COD_coll=MODISall_COD[coll_condition_m]                            
                                        
                        if MODIS_COD_coll.shape[0] and max(MODIS_COD_coll)>0:

                            Timed=Timediff[coll_condition_m]
                            Disd=Dis[coll_condition_m]
                            
                            # nearest
                            MODIS_COD=MODIS_COD_coll[MODIS_COD_coll>0]
                            Disd=Disd[MODIS_COD_coll>0]
                            MOCOD=MODIS_COD[Disd==min(Disd)]
                            MOCOD=min(MOCOD)
                            
                            opaque=[i for i, x in enumerate(CAL_Opa[j]) if x]
                            CALIPSO_nonopCOD=sum(CAL_OPTDEP[j,0:opaque[0]])
                            
                            tau_opaque=MOCOD-CALIPSO_nonopCOD
                        
                        ## find first cloud top and cloud base
                        for ij in range(ind,len(Coll_Height)):
                            if np.isnan(Coll_Height[ij])==1:
                                
                                Coll_base=Coll_Height[ij-1]
                                Cloudbaseind=ij-1
                                Coll_top=Coll_Height[ind]
                                Cloudtopind=ind
                                
                                Coll_IWC=Coll_IWC[Cloudtopind:ij]
                                Coll_LWC=Coll_LWC[Cloudtopind:ij]                  
                                ### write out file
                                
                                if Coll_top==Coll_base:
                                    Coll_base=Coll_base-0.12
                                    Coll_top=Coll_top+0.12
                                       
                                    ### write out file
                                    break
                                
                                else:
                                    ### write out file
                                    break
                            
                    ### if the highest cloud top that higher than CALIPSO opaque cloud base 
                    ### has only 1 pixel
                    else:
                        
                        # filter out the collocate CLOUDSAT pixel that no cloud top is higher than
                        # CALIPSO opaque cloud base
                        oneind=next(x[0] for x in enumerate(CLOUDSAT_maxtop_coll) if x[1]>0)
                        #oneind=np.where(CLOUDSAT_maxtop_coll>0)[0]
                        Coll_Height=CLOUDSAT_Height_coll_[oneind,:]
                        Coll_IWC=CLOUDSAT_iwc_coll_[oneind,:]
                        Coll_LWC=CLOUDSAT_lwc_coll_[oneind,:]                            
                        #####
                        
                        Coll_Height[Coll_Height<=0]=np.nan
                        df=pd.Series(list(Coll_Height))
                        ind=df.first_valid_index()
                        
                        ## find first cloud top and cloud base
                        for ij in range(ind,len(Coll_Height)):
                            
                            if np.isnan(Coll_Height[ij])==1:
                                
                                Coll_base=Coll_Height[ij-1]
                                Cloudbaseind=ij-1
                                Coll_top=Coll_Height[ind]
                                Cloudtopind=ind
                                
                                Coll_IWC=Coll_IWC[Cloudtopind:ij]
                                Coll_LWC=Coll_LWC[Cloudtopind:ij]   
                                
                                ### write out file
                                if Coll_top==Coll_base:
                                    Coll_base=Coll_base-0.12
                                    Coll_top=Coll_top+0.12
                                    ### write out file
                                    break
                                else:
                                    ### write out file
                                    break
   
                             
outfile.close()
elapsed=time.time()-t
print(elapsed/3600)