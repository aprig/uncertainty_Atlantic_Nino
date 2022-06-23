
# Author : Arthur Prigent
# Email  : aprigent@geomar.de


###############################################################################################################################
# Load modules 
###############################################################################################################################
dir_data = '/Users/aprigent/Documents/Thesis_GEOMAR/Data/'
from load_librairies import *
from datetime import datetime
from scipy import stats


###############################################################################################################################
# Functions
###############################################################################################################################
def center_map(data):
    if len(data.shape)>3:
        centered_data = ma.concatenate((data[:,0,:,180:360],data[:,0,:,0:180]),axis=2)
    else:
        centered_data = ma.concatenate((data[:,:,180:360],data[:,:,0:180]),axis=2)
    return centered_data

###############################################################################################################################
def center_map_ncep(data):
    if len(data.shape)>3:
        centered_data = ma.concatenate((data[:,0,:,72:],data[:,0,:,:72]),axis=2)
    else:
        centered_data = ma.concatenate((data[:,:,72:],data[:,:,:72]),axis=2)
    return centered_data

###############################################################################################################################

def get_monthly_mean(data):
    monthly_mean_file = np.ones((12,data.shape[1],data.shape[2]))*np.nan
    for i in range(12):
        a = data[np.arange(i,data.shape[0],12),:,:]
        monthly_mean_file[i,:,:] = np.nanmean(a,axis=0)
        
    return monthly_mean_file

###############################################################################################################################
def get_anomaly(data,mean_file):
    anomaly_file = np.ones((data.shape))*np.nan
    n=0
    for j in range(int(data.shape[0]/12)):
        for k in range(12):
            anomaly_file[n,:,:] = np.array(data[n,:,:])- np.array(mean_file[k,:,:])
            n+=1
    return anomaly_file

###############################################################################################################################
def get_std_ano(data):
    std_ano_file = np.ones(((12,180,360)))*np.nan
    for i in range(12):
        d = data[np.arange(i,data.shape[0],12),:,:]
        std_ano_file[i,:,:] = np.nanstd(d,axis=0)
        
    return std_ano_file
  
###############################################################################################################################  
def get_seasonal_mean(data):
    spring_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    summer_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    autumn_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    winter_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    for i in range(1,np.int(data.shape[0]/12)):
            spring_mean[i,:,:] = (data[2+((i)*(12)),:,:] +data[3+(i*12),:,:] +data[4+(i*12),:,:])/3
            summer_mean[i,:,:] = (data[5+((i)*(12)),:,:] +data[6+(i*12),:,:] +data[7+(i*12),:,:])/3
            autumn_mean[i,:,:] = (data[8+((i)*(12)),:,:] +data[9+(i*12),:,:] +data[10+(i*12),:,:])/3
            winter_mean[i,:,:] = (data[11+((i-1)*(12)),:,:] +data[0+(i*12),:,:] +data[1+(i*12),:,:])/3
            
    return spring_mean,summer_mean,autumn_mean,winter_mean    



###############################################################################################################################

def get_all_months(data):
    
    jav_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    feb_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    mar_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    apr_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    may_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    jun_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    jul_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    aug_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    sep_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    oct_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    nov_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan
    dec_mean = np.ones((np.int(data.shape[0]/12),data.shape[1],data.shape[2]))*np.nan

    for i in range(1,np.int(data.shape[0]/12)):
            jav_mean[i,:,:] = data[0+((i)*(12)),:,:]
            feb_mean[i,:,:] = data[1+((i)*(12)),:,:]
            mar_mean[i,:,:] = data[2+((i)*(12)),:,:]
            apr_mean[i,:,:] = data[3+((i)*(12)),:,:]
            may_mean[i,:,:] = data[4+((i)*(12)),:,:]
            jun_mean[i,:,:] = data[5+((i)*(12)),:,:]
            jul_mean[i,:,:] = data[6+((i)*(12)),:,:]
            aug_mean[i,:,:] = data[7+((i)*(12)),:,:]
            sep_mean[i,:,:] = data[8+((i)*(12)),:,:]
            oct_mean[i,:,:] = data[9+((i)*(12)),:,:]
            nov_mean[i,:,:] = data[10+((i)*(12)),:,:]
            dec_mean[i,:,:] = data[11+((i)*(12)),:,:]

            
    return jav_mean,feb_mean,mar_mean,apr_mean,may_mean,jun_mean,jul_mean,aug_mean,sep_mean,oct_mean,nov_mean,dec_mean
###############################################################################################################################
def get_BF1(mean_ano_data_atl3,wind_ano_pointwise):
    
    BF1 = np.zeros((wind_ano_pointwise.shape))*np.nan
    r_BF1 = np.zeros((wind_ano_pointwise.shape))*np.nan
    for j in range(BF1.shape[1]):
        for k in range(BF1.shape[2]):
    
            BF1[:,j,k],_,r_BF1[:,j,k],_,_    = stats.linregress(mean_ano_data_atl3[1:],
                                                  wind_ano_pointwise[1:,j,k])
    
    return BF1, r_BF1

###############################################################################################################################
def plot_BF1(lon,lat,BF1,r_BF1,bounds,ax,minlon,minlat,maxlon,maxlat,cax):
    xlon,ylat = np.meshgrid(lon,lat)
    cmap = mpl.cm.get_cmap('cmo.balance', len(bounds)) 
    m0 = Basemap(projection='merc',
            llcrnrlat=minlat, llcrnrlon=minlon,
            urcrnrlat=maxlat, urcrnrlon=maxlon,
            resolution='c',ax=ax)


    lonLS,latLS     = 18,-33
    scaleunit       = 5000     #scale unit = 10 km 
    ptsize          = 5      # particles size
    fs              = 10      # fontsize
    
    
    x1,y1=m0(-20,-3)
    x2,y2=m0(0,-3)
    ax.plot([x1,x2],[y1,y2],'black',linewidth=2)
    x3,y3=m0(0,-3)
    x4,y4=m0(0,3)
    ax.plot([x3,x4],[y3,y4],'black',linewidth=2)
    x5,y5=m0(-20,-3)
    x6,y6=m0(-20,3)
    ax.plot([x5,x6],[y5,y6],'black',linewidth=2)
    x7,y7=m0(-20,3)
    x8,y8=m0(0,3)
    ax.plot([x7,x8],[y7,y8],'black',linewidth=2)
    
    
    #ax.set_title('January',fontsize=20)
    parallels=np.arange(np.ceil(np.min(lat)),np.ceil(np.max(lat)),10.)
    meridians=np.arange(np.ceil(np.min(lon)),np.ceil(np.max(lon)),20.)
    m0.drawparallels(parallels,labels=[True,False,False,False],fontsize='20',linewidth=0.75)
    m0.drawmeridians(meridians,labels=[False,False,False,True],fontsize='20',linewidth=0.75)

    m0.fillcontinents(color='white',lake_color='aqua')
    m0.drawcoastlines(linewidth=1)
    m0.drawcountries(linewidth=0.25) 
    x, y = m0(xlon,ylat)
    cs0 = m0.contour(x,y,np.nanmean(r_BF1[:,:,:],axis=0),levels=np.arange(-1,1.2,0.3),colors='black')
    ax.clabel(cs0, inline=1, fontsize=13)
    cs0 = m0.contourf(x,y,np.nanmean(BF1[:,:,:],axis=0),levels=bounds,cmap=cmap)
    cbar=plt.colorbar(cs0,cax,orientation='horizontal',extend='both')
    cbar.ax.tick_params(labelsize=20)
    cbar.set_label(r'10$^{-2}$ Pa/$^\circ$C',size=20)





###############################################################################################################################

def get_BF2(mean_ano_data_Watl,z20_ano_pointwise):
    
    BF2 = np.zeros((z20_ano_pointwise.shape))*np.nan
    r_BF2 = np.zeros((z20_ano_pointwise.shape))*np.nan
    for j in range(BF2.shape[1]):
        for k in range(BF2.shape[2]):
    
            BF2[:,j,k],_,r_BF2[:,j,k],_,_    = stats.linregress(mean_ano_data_Watl[1:],
                                                  z20_ano_pointwise[1:,j,k])
    
    return BF2, r_BF2


###############################################################################################################################
def plot_BF2(lon,lat,BF2,r_BF2,bounds,ax,minlon,minlat,maxlon,maxlat,cax):
    xlon,ylat = np.meshgrid(lon,lat)
    cmap = mpl.cm.get_cmap('cmo.balance', len(bounds)) 
    m0 = Basemap(projection='merc',
            llcrnrlat=minlat, llcrnrlon=minlon,
            urcrnrlat=maxlat, urcrnrlon=maxlon,
            resolution='c',ax=ax)


    lonLS,latLS     = 18,-33
    scaleunit       = 5000     #scale unit = 10 km 
    ptsize          = 5      # particles size
    fs              = 10      # fontsize
    
    
    x1,y1=m0(-20,-3)
    x2,y2=m0(0,-3)
    ax.plot([x1,x2],[y1,y2],'black',linewidth=2)
    x3,y3=m0(0,-3)
    x4,y4=m0(0,3)
    ax.plot([x3,x4],[y3,y4],'black',linewidth=2)
    x5,y5=m0(-20,-3)
    x6,y6=m0(-20,3)
    ax.plot([x5,x6],[y5,y6],'black',linewidth=2)
    x7,y7=m0(-20,3)
    x8,y8=m0(0,3)
    ax.plot([x7,x8],[y7,y8],'black',linewidth=2)
    
    
    #ax.set_title('January',fontsize=20)
    parallels=np.arange(np.ceil(np.min(lat)),np.ceil(np.max(lat)),10.)
    meridians=np.arange(np.ceil(np.min(lon)),np.ceil(np.max(lon)),20.)
    m0.drawparallels(parallels,labels=[True,False,False,False],fontsize='20',linewidth=0.75)
    m0.drawmeridians(meridians,labels=[False,False,False,True],fontsize='20',linewidth=0.75)

    m0.fillcontinents(color='white',lake_color='aqua')
    m0.drawcoastlines(linewidth=1)
    m0.drawcountries(linewidth=0.25) 
    x, y = m0(xlon,ylat)
    cs0 = m0.contour(x,y,np.nanmean(r_BF2[:,:,:],axis=0),levels=np.arange(-1,1.2,0.3),colors='black')
    ax.clabel(cs0, inline=1, fontsize=13)
    cs0 = m0.contourf(x,y,np.nanmean(BF2[:,:,:],axis=0),levels=bounds,cmap=cmap)
    cbar=plt.colorbar(cs0,cax,orientation='horizontal',extend='both')
    cbar.ax.tick_params(labelsize=20)
    cbar.set_label(r' m/10$^{-2}$ Pa',size=20)

###############################################################################################################################

def get_BF3(z20_ano_pointwise,SSTa_pointwise):
    
    BF3 = np.zeros((z20_ano_pointwise.shape))*np.nan
    r_BF3 = np.zeros((z20_ano_pointwise.shape))*np.nan
    for j in range(BF3.shape[1]):
        for k in range(BF3.shape[2]):
            if np.isfinite(SSTa_pointwise[1:,j,k]).all()==True:
    
                BF3[:,j,k],_,r_BF3[:,j,k],_,_    = stats.linregress(z20_ano_pointwise[1:,j,k],
                                                  SSTa_pointwise[1:,j,k])
    
    return BF3, r_BF3


###############################################################################################################################
def plot_BF3(lon,lat,BF3,r_BF3,bounds,ax,minlon,minlat,maxlon,maxlat,cax):
    xlon,ylat = np.meshgrid(lon,lat)
    cmap = mpl.cm.get_cmap('cmo.balance', len(bounds)) 
    m0 = Basemap(projection='merc',
            llcrnrlat=minlat, llcrnrlon=minlon,
            urcrnrlat=maxlat, urcrnrlon=maxlon,
            resolution='c',ax=ax)


    lonLS,latLS     = 18,-33
    scaleunit       = 5000     #scale unit = 10 km 
    ptsize          = 5      # particles size
    fs              = 10      # fontsize
    
    
    x1,y1=m0(-20,-3)
    x2,y2=m0(0,-3)
    ax.plot([x1,x2],[y1,y2],'black',linewidth=2)
    x3,y3=m0(0,-3)
    x4,y4=m0(0,3)
    ax.plot([x3,x4],[y3,y4],'black',linewidth=2)
    x5,y5=m0(-20,-3)
    x6,y6=m0(-20,3)
    ax.plot([x5,x6],[y5,y6],'black',linewidth=2)
    x7,y7=m0(-20,3)
    x8,y8=m0(0,3)
    ax.plot([x7,x8],[y7,y8],'black',linewidth=2)
    
    
    #ax.set_title('January',fontsize=20)
    parallels=np.arange(np.ceil(np.min(lat)),np.ceil(np.max(lat)),10.)
    meridians=np.arange(np.ceil(np.min(lon)),np.ceil(np.max(lon)),20.)
    m0.drawparallels(parallels,labels=[True,False,False,False],fontsize='20',linewidth=0.75)
    m0.drawmeridians(meridians,labels=[False,False,False,True],fontsize='20',linewidth=0.75)

    m0.fillcontinents(color='white',lake_color='aqua')
    m0.drawcoastlines(linewidth=1)
    m0.drawcountries(linewidth=0.25) 
    x, y = m0(xlon,ylat)
    cs0 = m0.contour(x,y,np.nanmean(r_BF3[:,:,:],axis=0),levels=np.arange(-1,1.2,0.3),colors='black')
    ax.clabel(cs0, inline=1, fontsize=13)
    cs0 = m0.contourf(x,y,np.nanmean(BF3[:,:,:],axis=0),levels=bounds,cmap=cmap)
    cbar=plt.colorbar(cs0,cax,orientation='horizontal',extend='both')
    cbar.ax.tick_params(labelsize=20)
    cbar.set_label(r' $^\circ$C/10m',size=20)

###############################################################################################################################
def load_ORAS4(detrend =False):
    
    data     = dir_data + 'ORA-S4/'
    file_sst = 'global_SST_ORA_S4_1958_2017_montlhy.nc'
    print('Loading the file :',data+file_sst)
    nc = xr.open_mfdataset(data+file_sst)

    
    SST_tmp= nc.THETAO[:,0,:,:]
    SST_tmp =SST_tmp.sel(TIME = slice(datetime(1958, 1, 1), datetime(2017, 12, 31)))
    tmp = xr.concat([SST_tmp[:,:,180:],SST_tmp[:,:,:180]],dim = 'LON') 
    tmp.coords['LON'] = (tmp.coords['LON'] + 180)%360 - 180
    
    lon = tmp.LON[:]
    lat = tmp.LAT[:]
    time_oras = tmp.TIME[:]
    tmp = tmp.load()
    sst_oras_tmp = tmp

    if detrend ==True:
        print('detrending')
        
        
        for i in range(tmp.shape[1]):
            for j in range(tmp.shape[2]):
                if np.isfinite(tmp[:,i,j]).all()==True:
                    sst_oras_tmp[:,i,j] = detrend_order_n(tmp[:,i,j],1)
                    #sst_oras_tmp[:,i,j] =sc.detrend(tmp[:,i,j],axis=0,type='linear')
        #             sst_oras_tmp[:,i,j] = detrend_order_n(tmp[:,i,j],1)
                    #sst_oras_tmp[:,i,j] = nandetrend(tmp[:,i,j],1)
                    
    elif detrend ==False:
        sst_oras_tmp =    tmp   
        
    X,Y = np.meshgrid(lon,lat) 
    
    sst_oras  = xr.Dataset({'sst': (['time','y','x'],sst_oras_tmp)}
                       ,coords={'time':np.array(time_oras),'lat':(['y','x'],np.array(Y)),'lon':(['y','x'],np.array(X))})
    del sst_oras_tmp,tmp,nc,X,Y
    return sst_oras

###############################################################################################################################
def load_hadi(detrend =False):
    file_sst = 'HadISST_sst.nc'
    nc = xr.open_mfdataset(dir_data+file_sst)
    print('Loading the file :',dir_data+file_sst)

    lon_hadi = nc.longitude
    lat_hadi = nc.latitude
    tmp = nc.sst[:,:,:]
    tmp = tmp.sel(time = slice(datetime(1870, 1, 1), datetime(2017, 12, 31)))
    time_hadi = tmp.time
    tmp = np.array(tmp)

    tmp[tmp<-100] = np.nan

    sst_hadi_tmp = np.ones(tmp.shape)*np.nan
    if detrend ==True: 
        print('detrending')
        for i in range(tmp.shape[1]):
            for j in range(tmp.shape[2]):
                if np.isfinite(tmp[:,i,j]).all()==True:
                #print('yes')
                    sst_hadi_tmp[:,i,j] =sc.detrend(tmp[:,i,j],axis=0,type='linear')
    elif detrend ==False:
        sst_hadi_tmp = tmp
        
    X,Y = np.meshgrid(lon_hadi,lat_hadi) 


    sst_hadi  = xr.Dataset({'sst': (['time','y','x'],sst_hadi_tmp)}
                   ,coords={'time':np.array(time_hadi),'lat':(['y','x'],np.array(Y)),'lon':(['y','x'],np.array(X))})

    del sst_hadi_tmp,tmp,nc,X,Y
    return sst_hadi


###############################################################################################################################
def load_OI(detrend=False):
    
    print('loading file = ',dir_data+'OI_SST/SST_OI_1981_12_to2018_09_mnmean.nc')
    nc = xr.open_mfdataset(dir_data+'OI_SST/SST_OI_1981_12_to2018_09_mnmean.nc')
    tmp = nc.sst[:,:,:]
    tmp = xr.concat([tmp[:,:,180:],tmp[:,:,:180]],dim = 'lon') 
    tmp.coords['lon'] = (tmp.coords['lon'] + 180)%360 - 180
    lon_oiss  = tmp.lon[:]
    #lon_oiss  = np.arange(-179.5,180.5,1)
    lat_oiss  = tmp.lat[:]
    time_oiss = tmp.time[:]
    
    #tmp1 = np.array(xr.concat([tmp[:,:,180:],tmp[:,:,:180]],dim='lon'))

    
    file_sst  = 'lsmask.nc'
    nc        = xr.open_mfdataset(dir_data+'OI_SST/'+file_sst)
    #print(nc.variables.keys())
    mask_tmp      = nc.mask[0,:,:]
    mask = np.array(xr.concat([mask_tmp[:,180:],mask_tmp[:,:180]],dim='lon'))
    tmp2 = np.ones((np.shape(tmp)))* np.nan
    tmp = np.array(tmp)
    for i in range(tmp.shape[0]):
        tmp2[i,:,:] = tmp[i,:,:] * mask 

    tmp2 = tmp2.astype('float')
    tmp2[tmp2 == 0] = np.nan # or use np.nan

    sst_oiss_tmp = np.ones(tmp2.shape)*np.nan


    if detrend ==True:
        print('detrending')
        for i in range(tmp2.shape[1]):
            for j in range(tmp2.shape[2]):
                if np.isfinite(tmp2[:,i,j]).all()==True:
                    #print('yes')
                    sst_oiss_tmp[:,i,j] =sc.detrend(tmp2[:,i,j],axis=0,type='linear')
                    #sst_oiss_tmp[:,i,j] =detrend_order_n(tmp2[:,i,j],1)
    elif detrend ==False:
        sst_oiss_tmp = tmp2

    X,Y = np.meshgrid(lon_oiss,lat_oiss) 

    sst_oiss_tmp = sst_oiss_tmp[:,:,:]
    sst_oiss  = xr.Dataset({'sst': (['time','y','x'],sst_oiss_tmp)}
                       ,coords={'time':np.array(time_oiss),'lat':(['y','x'],np.array(Y)),'lon':(['y','x'],np.array(X))})



    del sst_oiss_tmp,tmp,tmp2,mask_tmp,mask,nc,X,Y
    
    return sst_oiss,lon_oiss,lat_oiss

###############################################################################################################################

def load_u_era_interim(detrend =False):
    
    file_sst = 'ERA_interim_uv_1979_2018_1x1.nc'
    nc = xr.open_mfdataset(dir_data+'ERA_INTERIM/'+file_sst)
    print('Loading the file :',dir_data+file_sst)
    u10 = nc.u10[:]
    u10_tmp = xr.concat([u10[:,:,180:],u10[:,:,:180]],dim = 'longitude') 
    u10_tmp.coords['longitude'] = (u10_tmp.coords['longitude'] + 180)%360 - 180
    lon = u10_tmp.longitude[:]
    lat = u10_tmp.latitude[:]
    time_eraint = u10_tmp.time
    tmp_u = np.ones(u10_tmp.shape)*np.nan
    u10_tmp = np.array(u10_tmp)
    if detrend ==True:
        print('detrending')
        for i in range(u10_tmp.shape[1]):
            for j in range(u10_tmp.shape[2]):
                if np.isfinite(u10_tmp[:,i,j]).all()==True:

                    tmp_u[:,i,j]     = sc.detrend(u10_tmp[:,i,j],axis=0,type='linear')

    else:
        tmp_u = u10_tmp


    X,Y = np.meshgrid(lon,lat) 


    uwnd_eraint  = xr.Dataset({'uwnd': (['time','y','x'],tmp_u)}
                   ,coords={'time':np.array(time_eraint),'lat':(['y','x'],np.array(Y)),'lon':(['y','x'],np.array(X))})

    return uwnd_eraint,lon,lat

###############################################################################################################################

def load_v_era_interim(detrend =False):
    
    file_sst = 'ERA_interim_uv_1979_2018_1x1.nc'
    nc = xr.open_mfdataset(dir_data+'ERA_INTERIM/'+file_sst)
    print('Loading the file :',dir_data+file_sst)
    v10 = nc.v10[:]
    v10_tmp = xr.concat([v10[:,:,180:],v10[:,:,:180]],dim = 'longitude') 
    v10_tmp.coords['longitude'] = (v10_tmp.coords['longitude'] + 180)%360 - 180
    lon = v10_tmp.longitude[:]
    lat = v10_tmp.latitude[:]
    time_eraint = v10_tmp.time
    tmp_v = np.ones(v10_tmp.shape)*np.nan
    v10_tmp = np.array(v10_tmp)
    if detrend ==True:
        print('detrending')
        for i in range(v10_tmp.shape[1]):
            for j in range(v10_tmp.shape[2]):
                if np.isfinite(v10_tmp[:,i,j]).all()==True:

                    tmp_v[:,i,j]     = sc.detrend(v10_tmp[:,i,j],axis=0,type='linear')

    else:
        tmp_v = v10_tmp
        

    X,Y = np.meshgrid(lon,lat) 

    vwnd_eraint  = xr.Dataset({'vwnd': (['time','y','x'],tmp_v)}
                   ,coords={'time':np.array(time_eraint),'lat':(['y','x'],np.array(Y)),'lon':(['y','x'],np.array(X))})


    return vwnd_eraint,lon,lat


###############################################################################################################################


def data_sub(data,lon_min,lon_max,lat_min,lat_max):
    
    '''Define a box between lon_min lon_max lat_min and lat_max and 
    extract the data in the box and drop everything else.
    
    
    Parameters
    ----------
    
    data : xarray_like
    Data to be subdomained. 
    
    lon_min : integer
    Longitude minimum of the subdomain
    
    lon_max : integer
    Longitude maximum of the subdomain
    
    lat_min : integer
    Latitude minimum of the subdomain
    
    lat_max : integer
    Latitude maximum of the subdomain
    
    Returns
    ---------
    
    data_sub : xarray_like
    Subdomain. 
    '''
    
    try:
        data_sub = data.where((  data.lon>=lon_min) & (data.lon<=lon_max) & (data.lat<=lat_max) & (data.lat>=lat_min),
                                                                          drop=True)
    except AttributeError:
        try:
            data_sub = data.where((  data.nav_lon>=lon_min) & (data.nav_lon<=lon_max) & (data.nav_lat<=lat_max) & (data.nav_lat>=lat_min),drop=True)
        except AttributeError:
            try:
                data_sub = data.where((  data.longitude>=lon_min) & (data.longitude<=lon_max) & (data.latitude<=lat_max) & (data.latitude>=lat_min),drop=True)
            except AttributeError:
                try:
                    data_sub = data.where((  data.x>=lon_min) & (data.x<=lon_max) & (data.y<=lat_max) &
                                      (data.y>=lat_min),drop=True)
                except AttributeError:
                    data_sub = data.where((  data.LON>=lon_min) & (data.LON<=lon_max) & (data.LAT<=lat_max) &
                                      (data.LAT>=lat_min),drop=True)
            
    

 
    
    return data_sub




###############################################################################################################################
def ano_norm_t(ds):
    
    '''Compute the anomalies by removing the monthly means. 
    The anomalies are normalized by their corresponding month.
    
    Parameters
    ----------
    
    ds : xarray_like
    Timeserie or 3d field.
    
    Returns
    -----------
    
    ano : xarray_like
    Returns the anomalies of var relative the climatology.
    
    ano_norm : xarray_like
    Returns the anomalies of var relative the climatology normalized by the standard deviation.
    
    '''    
    
    clim     = ds.groupby('time.month').mean('time')
    clim_std = ds.groupby('time.month').std('time')
    ano      = ds.groupby('time.month') - clim
    ano_norm = xr.apply_ufunc(lambda x, m, s: (x - m) / s,
                                    ds.groupby('time.month'),
                                    clim, clim_std)
    
    return ano, ano_norm
#########################################################################################################################################
def ano_norm_tc(ds):
    
    '''Compute the anomalies by removing the monthly means. 
    The anomalies are normalized by their corresponding month.
    
    Parameters
    ----------
    
    ds : xarray_like
    Timeserie or 3d field.
    
    Returns
    -----------
    
    ano : xarray_like
    Returns the anomalies of var relative the climatology.
    
    ano_norm : xarray_like
    Returns the anomalies of var relative the climatology normalized by the standard deviation.
    
    '''
    
    
    clim     = ds.groupby('time_counter.month').mean('time_counter')
    clim_std = ds.groupby('time_counter.month').std('time_counter')
    ano      = ds.groupby('time_counter.month') - clim
    ano_norm = xr.apply_ufunc(lambda x, m, s: (x - m) / s,
                                    ds.groupby('time_counter.month'),
                                    clim, clim_std)
    
    return ano, ano_norm

###############################################################################################################################
###############################################################################################################################
##### COMPUTE THE LATENT HEAT FLUX LIKE IN BENTAMY ET AL 2003 ######

def Ce(u10):
    '''
    Dalton number is estimated from wind speed, air temperature and SST
    '''
    a = -0.146785
    b = -0.292400
    c = -2.206648
    d = 1.6112292
    return 1e-3*(a*np.exp(b*(u10+c)) + (d/u10) + 1)


def l(Ts):
    '''
    Latent heat of evaporation
    '''
    return np.array(4186.8*(597.31 - 0.5625*Ts)) #J/kg

def Tv(T10,qa):

    return T10*(1+0.608*qa)


def rho(P0,Tv):
    '''
    Air density
    '''
    return (100*P0)/(287*Tv) #kg.m-3

def es(Ts):
    Ts=Ts+273.15
    a = -4.928
    b = 23.55
    c = -2937
    return Ts**a *10**(b+(c/Ts))

def qs(es,ps):
    '''
    Saturated surface humidity'''
    return (5/8)*(es /(ps-es)) # kg/kg


def LH(Ts,qa,u10,Ps):
    
    return (l(Ts) * 1.15 * Ce(u10) * u10 * (qs(es(Ts),Ps)*1000 - qa)/1000) #w.m-2


###############################################################################################################################

def nandetrend(x,y):
    ''' Remove the linear trend from the data '''
    
    
    not_nan_ind = ~np.isnan(y)
    m, b, r_val, p_val, std_err = stats.linregress(x[not_nan_ind],np.array(y)[not_nan_ind])
    detrended_y = np.array(y) - m*x
    
    return detrended_y

###############################################################################################################################
###############################################################################################################################

###################################################################################################################################################
def find_benguela_nino_nina(timeserie_ano_ABA_file,treshold):
    
    index_nino_sst = xr.full_like(timeserie_ano_ABA_file,np.nan)
    index_nino_sst.sst[timeserie_ano_ABA_file.to_array()[0,:] >= treshold] = 1
    index_nino_sst.sst[timeserie_ano_ABA_file.to_array()[0,:] <  treshold] = 0
    index_nino_sst_tmp = np.array(index_nino_sst.to_array()[0,:])
    
    index_nina_sst = xr.full_like(timeserie_ano_ABA_file,np.nan)
    index_nina_sst.sst[timeserie_ano_ABA_file.to_array()[0,:] <= -treshold] = 1
    index_nina_sst.sst[timeserie_ano_ABA_file.to_array()[0,:] > - treshold] = 0
    index_nina_sst_tmp = np.array(index_nina_sst.to_array()[0,:])
    
    
    diff_indexes_nino = np.diff(index_nino_sst_tmp[1:-1])
    diff_indexes_nina = np.diff(index_nina_sst_tmp[1:-1])
    
    id_str_nino = []
    id_end_nino = []
    for i in range(len(diff_indexes_nino)):
        if diff_indexes_nino[i]==1.0:
            id_str_nino.append(i)
        elif diff_indexes_nino[i] == -1.0:
            id_end_nino.append(i)


    id_str_nina = []
    id_end_nina = []
    for i in range(len(diff_indexes_nina)):
        if diff_indexes_nina[i]==1.0:
            id_str_nina.append(i)
        elif diff_indexes_nina[i] == -1.0:
            id_end_nina.append(i)
    print('nino str',len(id_str_nino))
    print('nino end',len(id_end_nino))
    try:
        nino_indexes_tmp= np.vstack((id_str_nino,id_end_nino))
        length_events_nino = nino_indexes_tmp[1,:] - nino_indexes_tmp[0,:]
        nino_indexes = np.vstack((nino_indexes_tmp,length_events_nino))
    except ValueError:
        nino_indexes_tmp= np.vstack((id_str_nino,id_end_nino[:-1]))
        length_events_nino = nino_indexes_tmp[1,:] - nino_indexes_tmp[0,:]
        nino_indexes = np.vstack((nino_indexes_tmp,length_events_nino))
        
    try:
        nina_indexes_tmp= np.vstack((id_str_nina,id_end_nina[:]))
        length_events_nina = nina_indexes_tmp[1,:] - nina_indexes_tmp[0,:]
        nina_indexes = np.vstack((nina_indexes_tmp,length_events_nina))
    except ValueError:
        nina_indexes_tmp= np.vstack((id_str_nina,id_end_nina[1:]))
        length_events_nina = nina_indexes_tmp[1,:] - nina_indexes_tmp[0,:]
        nina_indexes = np.vstack((nina_indexes_tmp,length_events_nina))
    
    
    return np.array(nino_indexes), np.array(nina_indexes)





#########################################################################################################################################



def count_events(nino_indexes,nina_indexes):
    ni_events=0
    for i in range(nino_indexes.shape[1]):
        tmp = nino_indexes[2,:]>=2
        if tmp[i]==True:
            ni_events+=1
    na_events = 0    
    for i in range(nina_indexes.shape[1]):
        tmp = nina_indexes[2,:]>=2
        if tmp[i]==True:
            na_events+=1
    return ni_events,na_events



#########################################################################################################################################

def get_peak(data,nino_nina_indexes,nino = False):
    ''' This function find the index of the local maximum between
        the starting months of the events and the end.'''
    
    if nino == True:
        peak_nino_nina_indexes = []
        for i in range(nino_nina_indexes.shape[1]):
            if nino_nina_indexes[2,i]>=3:
                peak_nino_nina_indexes.append(nino_nina_indexes[0,i]+np.argmax(data.to_array()[0,nino_nina_indexes[0,i]:nino_nina_indexes[1,i]]))
    else:
        peak_nino_nina_indexes = []
        for i in range(nino_nina_indexes.shape[1]):
            if nino_nina_indexes[2,i]>=3:
                peak_nino_nina_indexes.append(nino_nina_indexes[0,i]+np.argmin(data.to_array()[0,nino_nina_indexes[0,i]:nino_nina_indexes[1,i]]))


    return np.array(peak_nino_nina_indexes)


#######################################################################################################################

def detrend_order_n(data,order):
    '''
    Remove the linear trend from the data.
    
    Parameters 
    ----------
    data : xarray.data_array
    Input time series.
    
    order : integer
    Order of the detrending.
    order = 1 --> linear detrending
    order = 2 --> quadrtic detrending 
    
    Returns 
    -----------
    data : xarray.data_array
    The detrended input data
    
    References
    -----------
    
    see : https://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
     for more information.
    
    '''
    x = np.arange(0,data.shape[0],1)
    y = data.values
    coeff  = np.polyfit(x,y,order)
    y_poly = np.polyval(coeff,x)
    detrended_data = data.values-y_poly
    
    data.values = detrended_data
    #data['detrend'] = detrended_data
    
    return data

#######################################################################################################################


def draw_proj_rect( minlat, maxlat, minlon, maxlon, m, color):
    
    '''Draw a rectangle with respect to the project of the map (m)
    Please choose the color of the edge of the rectangle. 
    
    Parameters
    ----------
    minlon : integer
    Minimum longitude
    
    maxlon : integer
    Maximum longitude
    
    minlat : integer
    Minimum latitude
    
    maxlat : integer
    Maximum latitude
    
    m : basemap
    Map from basemap already projected.
    
    color : string
    Color of the box to be drawn
    
    lw : integer
    Linewidth of the contour of the box. 
    
    
    Returns
    ----------
    Draw the projected rectangle on the input m map.
    
    '''
    
    lons=np.hstack((np.repeat(minlon,10),\
                  np.linspace(minlon,maxlon, num=10),\
                  np.repeat(maxlon,10),\
                  np.linspace(maxlon,minlon, num=10)))

    lats=np.hstack((np.linspace(minlat,maxlat, num=10),\
                  np.repeat(maxlat,10),\
                  np.linspace(maxlat,minlat, num=10),
                  np.repeat(minlat,10)))

    m.plot(y=lats,x=lons,latlon=True, lw=4, color=color, alpha=0.8)
