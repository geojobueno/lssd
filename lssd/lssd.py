#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:52:57 2025

@author: jobueno
"""

import struct
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

def get_filesize(file):
    with open(file) as f:
        f.seek(0,2)
        s = f.tell()
        f.close()
        return s

class LumiData:
    
    def __init__(self, file='', dtype_print=False, osl_lims=[0,100]):
        self.file = file #file path
        
        #read options
        self.dtype_print = dtype_print
        
        #data declaration
        self.maindataframe, self.measures = self.readbin(self.file)
        #dataframe with all data
        #dictionary with measure of OSL/TL - ID of maindataframe is the kay
                
        
        #processing options/params
        self.osl_lims = osl_lims #s - integral limits to OSL
    
    def __str__(self):
        outstr = f'''file: {self.file}\n{len(self.maindataframe)} elements
samples: {"; ".join(self.maindataframe.sample_name.unique())}'''
        
        return outstr
    
    def __getitem__(self, elid):
        new_dict = {}
        returndata = LumiData()
        returndata.file = self.file
        returndata.osl_lims = self.osl_lims
        
        if isinstance(elid, str):  # Se for um rótulo
            #Usar isso para quando eu quiser buscar um single grain    
            None
        
        elif type(elid) is int:
            returndata.maindataframe = self.maindataframe.loc[elid,:]
            returndata.measures = self.measures[elid]
            return returndata
        
        elif len(elid) > 0:
            for k in elid:
                new_dict[k] = self.measures[k] #copy of some items
            
            returndata.maindataframe = self.maindataframe.loc[elid,:]
            returndata.measures = new_dict
            return returndata
        
    def _get_binrows(self):
        #file = self.file
        fsize = get_filesize(self.file)
        remaining_size = fsize
        l_records = [0]
        acc_record = 0
        c = 0
        
        with open(self.file,'rb') as f:
            while remaining_size > 0:
                f.seek(acc_record+2)
                length_record = struct.unpack('<i',f.read(4))[0]
                #f.seek(length_record)
                acc_record = acc_record + length_record
                remaining_size = remaining_size - length_record
                
                l_records.append(length_record)
                c=c+1
                
        return c,l_records
    
    def _read_row(self, data, pos=0):
        data.seek(pos)
        
        #Header
        version = struct.unpack('<h', data.read(2))[0]
        length_record = struct.unpack('<i',data.read(4))[0]
        length_record_previous = struct.unpack('<i',data.read(4))[0]
        ndatapoints = struct.unpack('<i',data.read(4))[0]
        recordtype=0
        if version==8:
            recordtype = struct.unpack('B',data.read(1))[0]
                   
        
        #Sample characteristics
        run_number = struct.unpack('<h', data.read(2))[0]
        set_number = struct.unpack('<h', data.read(2))[0]
        carousel_pos = struct.unpack('<h', data.read(2))[0]
        grain_number = struct.unpack('<h', data.read(2))[0]
        curve_number = struct.unpack('<h', data.read(2))[0]
        x_pos = struct.unpack('<h', data.read(2))[0]
        y_pos = struct.unpack('<h', data.read(2))[0]
        sample_name_size = struct.unpack('B', data.read(1))[0]
        sample_name = ''
        if sample_name_size>0:
            sample_name = data.read(sample_name_size).decode('latin-1')
        skip_name = data.read(20-sample_name_size)
        
        sample_comment_size = struct.unpack('B', data.read(1))[0]
        sample_comment = ''
        if sample_comment_size >0:
            sample_comment = data.read(sample_comment_size).decode('latin-1')
            # sample_comment = struct.unpack('B', data.read(sample_comment_size))[0]
        skip_comment = data.read(80-sample_comment_size)
        
        #Instrumental and sequences characteristics
        systemid = struct.unpack('<h', data.read(2))[0]
        fname_size = struct.unpack('B', data.read(1))[0]
        fname = data.read(fname_size).decode('latin-1')
        skip_fname = data.read(100-fname_size)    
        username_size = struct.unpack('B', data.read(1))[0]
        username = data.read(username_size).decode('latin-1')
        skip_username=data.read(30-username_size)
        collect_time_s = struct.unpack('B', data.read(1))[0]
        if collect_time_s != 6:
            raise ValueError(f'TIME FORMAT INCORRECT: {collect_time_s} NUMBERS')
        time = data.read(6).decode('ascii') #hh-mm-ss
        collect_date_s = struct.unpack('B', data.read(1))[0]
        if collect_date_s != 6:
            raise ValueError(f'TIME FORMAT INCORRECT: {collect_date_s}')
        date = data.read(6).decode('ascii') #dd-mm-yy
        
        #Analysis
        dtype = struct.unpack('B', data.read(1))[0]
        bleaching_time = struct.unpack('f', data.read(4))[0]
        bleaching_unit = struct.unpack('B', data.read(1))[0]
        norm_factors = struct.unpack('4f', data.read(16)) #norm1, norm2, norm3, BG
        shift_data = struct.unpack('h', data.read(2))[0]
        tag = struct.unpack('B', data.read(1))[0]
        
        norm1,norm2,norm3,bg = norm_factors
        skip_internal = data.read(20)
        
        
        #Mesurements 
        lumi_type = struct.unpack('B', data.read(1))[0]
        lightsource = struct.unpack('B', data.read(1))[0]
        optical_stim_power = struct.unpack('f', data.read(4))[0]
        low = struct.unpack('f', data.read(4))[0]
        high = struct.unpack('f', data.read(4))[0]
        rate = struct.unpack('f', data.read(4))[0]
        sample_temperature = struct.unpack('h', data.read(2))[0]
        measured_temperature = struct.unpack('h', data.read(2))[0]
        preheating_temperature = struct.unpack('f', data.read(4))[0]
        preheating_time = struct.unpack('f', data.read(4))[0]
        
        TOL_delay, TOL_on, TOL_off = struct.unpack('3h', data.read(6)) #delay,on,off
        irradiation_time = struct.unpack('f', data.read(4))[0]
        irradiation_type = struct.unpack('B', data.read(1))[0]
        irradiation_dose_rate = struct.unpack('f', data.read(4))[0]
        irradiation_dose_rate_error = struct.unpack('f', data.read(4))[0]
        time_since_last_irr = struct.unpack('i', data.read(4))[0]
        
        time_unit_pulse = struct.unpack('f', data.read(4))[0]
        ontime_pulse = struct.unpack('i', data.read(4))[0]
        stimulation_period = struct.unpack('i', data.read(4))[0]
        gating_signal = struct.unpack('B', data.read(1))[0] #PMT signal, start, end
        gating_start,gating_end = struct.unpack('2i',data.read(8)) #start,end
        photon_timer = struct.unpack('B', data.read(1))[0]
        PMT_deadtime_corr = struct.unpack('B', data.read(1))[0]
        PMT_deadtime = struct.unpack('f', data.read(4))[0]
        stim_power_100 = struct.unpack('f', data.read(4))[0]
        
        XRF_time = struct.unpack('f', data.read(4))[0] #time (s)
        XRF_V = struct.unpack('f', data.read(4))[0] #voltage (V)
        XRF_A = struct.unpack('i', data.read(4))[0] #current (uA)
        XRF_deadtime = struct.unpack('f', data.read(4))[0]
        
        mk1=mk2=mk3=-1
        extr_start=extr_end=-1
        detectorID=-1
        filters_l=filters_u=-1
        excess_noise_factor=-1
        
        if version==8 or version==7:
            detectorID = struct.unpack('B', data.read(1))[0]
            filters_l,filters_u = struct.unpack('2h', data.read(4)) #lower, upper
            excess_noise_factor = struct.unpack('f', data.read(4))[0]
        if version==8:
            markers = struct.unpack('6f', data.read(24)) #mrk.x,mrk.y
            extr_start, extr_end = struct.unpack('2f', data.read(8)) #start,end
            mk1 = markers[0:2] #xy
            mk2 = markers[2:4] #xy
            mk3 = markers[4:6] #xy
        
            
        
        skippos=42 #version 8
        
        if version==7:
            skippos=15
        elif version==6:
            skippos=24
        
        skip_internal = data.read(skippos)
        
        
        
        #print(data.tell())
        #Get points - depends of recordtype
        
        if recordtype!=128:
            points = np.array(struct.unpack(f'<{ndatapoints}i',data.read(ndatapoints*4)))
            xpoints = np.linspace(low,high,ndatapoints,endpoint=False)
        
        else:
            points = xpoints = None #temporary - adapt it in the future
        
        varlist = [version,length_record,length_record_previous,ndatapoints,recordtype,run_number,
                   set_number,carousel_pos,grain_number,curve_number,x_pos,y_pos,sample_name,sample_comment,systemid,fname,
                   username,time,date,dtype,bleaching_time,bleaching_unit,shift_data,tag,norm1,norm2,norm3,bg,lumi_type,lightsource,
    optical_stim_power,low,high,rate,sample_temperature,measured_temperature,preheating_temperature,
    preheating_time,TOL_delay,TOL_on,TOL_off,irradiation_time,irradiation_type,irradiation_dose_rate,
    irradiation_dose_rate_error,time_since_last_irr,time_unit_pulse,ontime_pulse,stimulation_period,
    gating_signal,gating_start,gating_end,photon_timer,PMT_deadtime_corr,PMT_deadtime,stim_power_100,
    XRF_time,XRF_V,XRF_A,XRF_deadtime,detectorID,filters_l,filters_u,excess_noise_factor,mk1,mk2,mk3,extr_start,extr_end]
        pp = [xpoints,points]
        return varlist, pp

    def readbin(self, file=''):
        
        if len(file)==0:
            return None,None
        
        self.file = file
        fsize = get_filesize(self.file)
        nrows,lengths = self._get_binrows()
        acclens = np.cumsum(lengths)
        
        lumi_dict = {0: "TL",
                     1: "OSL",
                     2: "IRSL",
                     3: "M-IR",
                     4: "M-VIS",
                     5: "TOL",
                     6: "TRPOSL",
                     7: "RIR",
                     8: "RBR",
                     9: "USER",
                     10: "POSL",
                     11: "SGOSL",
                     12: "RL",
                     13: "XRF"}
        
        dtype_dict = {0:'Natural',
                      1:'N+dose',
                      2:'Bleach',
                      3:'Bleach+dose',
                      4:'Natural (Bleach)',
                      5:'N+dose (Bleach)',
                      6:'Dose',
                      7:'Background'}
        
        lightsource_dict = {0:'None',
                       1:'Lamp',
                       2:'IR diodes/IR Laser',
                       3:'Calibration LED',
                       4:'Blue Diodes',
                       5:'White light',
                       6:'Green laser (single grain)',
                       7:'IR laser (single grain)'}
        
        with open(self.file,'rb') as data:
            i = 0
            rows = []
            datapoints={}
            for i in range(nrows):        
                row,datapoints[i] = self._read_row(data,acclens[i])
                row.insert(0, i)
                
                if self.dtype_print==False:
                    row[20] = dtype_dict[row[20]]
                    row[29] = lumi_dict[row[29]]
                    row[30] = lightsource_dict[row[30]]
                
                rows.append(row)
            
            
            colunas = ['id','version','length_record','length_record_previous','ndatapoints','recordtype','run_number','set_number','carousel_pos','grain_number','curve_number','x_pos','y_pos','sample_name','sample_comment','systemid','fname','username','time','date','dtype','bleaching_time','bleaching_unit','shift_data','tag','norm1','norm2','norm3','bg','lumi_type','lightsource','optical_stim_power','low','high','rate','sample_temperature','measured_temperature','preheating_temperature','preheating_time','TOL_delay','TOL_on','TOL_off','irradiation_time','irradiation_type','irradiation_dose_rate','irradiation_dose_rate_error','time_since_last_irr','time_unit_pulse','ontime_pulse','stimulation_period','gating_signal','gating_start','gating_end','photon_timer','PMT_deadtime_corr','PMT_deadtime','stim_power_100','XRF_time','XRF_V','XRF_A','XRF_deadtime','detectorID','filters_l','filters_u','excess_noise_factor','mk1','mk2','mk3','extr_start','extr_end']
            
            df = pd.DataFrame(rows,columns=colunas)
        
        return df, datapoints
    
    def export_data(self):
        #salvar as coisas em csv
        return self.maindataframe, self.measures
    
    def filter_bylumi(self, lumi_type): #filter by lumi type
        del_indexes = self.maindataframe.loc[self.maindataframe['lumi_type']!=lumi_type].index
        print(f'deleting measures: {list(del_indexes)}')
        
        self.maindataframe.drop(del_indexes,inplace=True)
        for key in del_indexes:
            del self.measures[key]
        
        return self
    
    def mergebinx(self, pathfiles:list,update_position=True):
        if type(pathfiles) is not list:
            pathfiles=[pathfiles]
        
        cdf = self.maindataframe 
        mpid = cdf.id.unique().max()+1 #previous max id
        mppos = cdf.carousel_pos.unique().max() #previous max aliq position
        
        newpt = {}
        dfs = [cdf]
        
        for file in pathfiles:
            binx = LumiData(file)
            binxdf = binx.maindataframe
                        
            binxdf.id = binxdf.id+mpid
            binxpt = binx.measures
            
            #binxdf.id=binxdf.id+mpid
            if update_position == True:
                binxdf.carousel_pos = binxdf.carousel_pos+mppos
                
            for pk in range(len(binxdf)):
                #print(pk)
                #ck = pk + mpid
                #print(ck)
                #joao
                ck = binxdf.loc[pk,'id']
                newpt[ck] = binxpt[pk]
            
            
                
            dfs.append(binxdf)
            mpid = binxdf.id.unique().max()+1
            mppos = binxdf.carousel_pos.unique().max()
        
        self.maindataframe=pd.concat(dfs)
        self.maindataframe.reset_index(inplace=True,drop=True)
        self.measures={**self.measures, **newpt}
        
        
        return self
    
    def integrate_OSL(self, osl_lims=None, ids=None, 
                      bgtime=10, background=None, overwrite=False,dt=None): #background = ultimos 10 seg
        if osl_lims==None:
            osl_lims = self.osl_lims
            
        if ids is None:
            sel = (self.maindataframe.lumi_type=='OSL') | (self.maindataframe.lumi_type==1)
            ids = self.maindataframe.loc[sel,'id'] #all aliquots
        

        # adicionar overwrite pra nao sobrescrever as integrais ja feitas 
            
        self.maindataframe['integral_OSL'] = -999.9
        self.maindataframe['background_OSL'] = -999.9
        self.maindataframe['background_OSL'] = -999.9
        for i in ids:
            time_i = self.measures[i][0] #s
            meas_i = self.measures[i][1] #photons
            
            if dt is None:
                dt = time_i[1]-time_i[0]

            time_i = time_i - osl_lims[0]
            cond = (time_i>=0) & (time_i<osl_lims[1])
            
            if background is None:
                condbg = time_i>=time_i[-1]-bgtime
                # print(time_i)
                bgvalue=np.mean(meas_i[condbg])
                # print(time_i[condbg])
                # print(np.mean(meas_i[condbg]))
            else:
                bgvalue = background
                
                
            #time_i = time_i[cond]
            # print(int(osl_lims[0]/dt), int(osl_lims[1]/dt))
            
            
            intg = np.sum(meas_i[cond]-bgvalue)*dt
            #np.trapezoid() for future numpy versions
            self.maindataframe.loc[i,'background_OSL'] = bgvalue
            self.maindataframe.loc[i,'integral_OSL'] = intg

    def plot(self, elid,**kwargs):
        #print(self.maindataframe[self.maindataframe.id==elid]['lumi_type'].item())
        unit = 'T (°C)' if self.maindataframe.loc[self.maindataframe['id']==elid]['lumi_type'].item() not in ['OSL','IRSL','SGOSL','POSL'] else 't (s)'
        default_params = {
            'title': f'Emission plot - {elid}',
            'xlabel': f'{unit}',
            'ylabel': 'Count',
            'color': 'black',
            'ls': '-',
            'lw': 1,
            'grid': True,
            'figsize': (6, 6),
            'show': True,
            'alpha': 1,
            'xlim':(0,np.max(self.measures[elid][0])),
            'ylim':(0,np.max(self.measures[elid][1])+np.max(self.measures[elid][1])*0.05)
        }
        
        plot_params = {**default_params, **kwargs} #plot_params = dict(default_params, **kwargs)
        fig, ax = plt.subplots(1,1,figsize=plot_params['figsize'])
        ax.plot(self.measures[elid][0],self.measures[elid][1],
                ls=plot_params['ls'],lw=plot_params['lw'],color=plot_params['color'])
        
        ax.set_xlabel(plot_params['xlabel'])
        ax.set_ylabel(plot_params['ylabel'])
        ax.set_title(plot_params['title'])
        ax.grid(plot_params['grid'])
        
        ax.set_xlim(plot_params['xlim'])
        ax.set_ylim(plot_params['ylim'])
        
        return fig, ax

# arquivo = '/home/jobueno/Documents/lssd/LxTx_EFC_01_02A_02B_02C.binx'
# binf = LumiData(arquivo)
# # binf.readbin()

# sg = '/home/jobueno/Documents/lssd/XY3SGLaser_Riso CalQz_002.binx'
# sg002 = LumiData().readbin(sg)


# xp = '/home/jobueno/Downloads/Sens_2_Qz_AM1.1_2_3_4_5_6_7_8_240909_R1.binx'
# xops = LumiData(xp)
# xops.integrate_OSL(osl_lims=[0,1],bgtime=10)


