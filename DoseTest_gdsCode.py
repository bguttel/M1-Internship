# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 22:43:48 2016

@author: shintaro
"""
#import gdsCAD as gds
import gdspy as gds
import numpy as np
from jobfile_generator import jobfile

def saveCell2GDS(cell, gdsName):
    """ This function save the given cell to GDS file with the name 'gdsName' """
    layout = gds.GdsLibrary()
    layout.add(cell, overwrite_duplicate=True)
    layout.write_gds(gdsName+'.gds')
    
    
def createIDT3(period=[1, 3],      # [initial period (um), final periof (um) (> initial period)]
              idt_type = 3,             # 1: single, 2: double 3: Split 52
              center_length = 30,       # [um]
              edge_length = 35,         # [um]
              gap_length = 1,           # [um] (if 0 there is no diffraction finger.)
              metal_ratio = 1,        # Ratio set to compensate proximity
#              period_in_time = 40,      # [ns] (period = ceil(period_in_time [ns] * saw_velocity [um/ns] / period [um]))
              total_period = 100,       #Total number of periods
              saw_velocity = 2.77,     # [um/ns]
              cell_name = '',           # name of the cell
              layer = 1,                # layer for idt
              datatype = 0,             # datatype for idt
              pad_layer = 2,            # layer number for bonding pad
              pad_datatype = 0,         # datatype for bonding pad
              inv_GND_layer = 3,        # layer number for inverse of GND plane
              inv_GND_datatype = 0,     # datatype for inverse of GND plane
              ):
    """
    This function create a cell containing IDT with/without bonding pad.
    """
    edge_configs = {
            'single':[1,-1],
            'double':[1,1,-1,-1],
            'split52':[1,1,-1,0,-1],
            'split4':[1,1,1,1,-1,-1,-1,-1],
            'unidir':[-1,-2,-9,1],
            'unidir0':[-1,-2,-9,1,2,9],
            'unidir1':[-2,-9,-1,2,9,1],
            'double3':[-1,-1,3,9],
            'DART':[1,-3,-9,-1],
            'iDART':[1,-1,-3,-9]
            }
    # MEANING OF NUMBERS:
    # 0       ... floating
    # +-1     ... single finger up/down
    # +-2,+-9 ... double-width finger up/down
    # +-3,+-9 ... tripple-width finger up/down; 
    # +-9     ... skip this finger
    pad = 0
    
    if idt_type == 1:
        suffix = 'single'
    elif idt_type == 2:
        suffix = 'double'
    elif idt_type == 3:
        suffix = 'split52'
    elif idt_type == 4:
        suffix = 'split4'
    elif idt_type == 5:
        suffix = 'unidir'
    elif idt_type == 6:
        suffix = 'unidir0'
    elif idt_type == 7:
        suffix = 'unidir1'
    elif idt_type == 8:
        suffix = 'double3'
    elif idt_type == 10:
        suffix = 'DART'
    elif idt_type == 11:
        suffix = 'iDART'
    edge_config = edge_configs[suffix]
    
    if cell_name == '':
        cell_name = 'idt'+str(int(period[0]*1000))+'-'+str(int(period[1]*1000))+'_'+suffix
    idt = gds.Cell(cell_name, exclude_from_current=True)
    
#    total_period = np.ceil(period_in_time * saw_velocity / ((period[0] + period[1])/2))
    total_period = int(40 * saw_velocity / ((period[0] + period[1])/2))
    total_length = total_period * (period[0] + period[1])/2
    noFingerPerPeriod = len(edge_config)
    center_pos = - total_length/2
    x_compensation = (period[1]-period[0])/(2*2*noFingerPerPeriod)/2 # compensate x to keep symmetry about x for chirped IDTs
    
    for i in range(total_period):
        current_period = period[0]+(period[1]-period[0])/(total_period-1) * i
        current_width = current_period /(2*noFingerPerPeriod)
        proximity_width = current_width * metal_ratio
        center_pos += current_period/2
        
        for j, c in enumerate(edge_config):
            x = center_pos - (noFingerPerPeriod - 1 - 2*j)*current_width + x_compensation
            if c > 0:
                sign = +1; connect = 1
            elif c < 0:
                sign = -1; connect = 1
            elif c == 0:
                sign = +1; connect = 0
            y0 = - sign*center_length/2
            y1 = sign*(center_length/2 + connect*edge_length)    
                
            if np.abs(c) == 1:
                tempwidth = proximity_width
            elif np.abs(c) == 2:
                tempwidth = current_width*2.0 * metal_ratio
                x += current_width
            elif np.abs(c) == 3:
                tempwidth = current_width*3.0 * metal_ratio
                x += current_width
                layer = 6

            if np.abs(c) == 9:
                dummy=0 # do nothing
            else:
                p = gds.PolyPath([(x,y0),(x,y1)],width=tempwidth, layer=layer, datatype=datatype)
                idt.add(p)
                layer = 1
              
            # add additional fingers at side    
            if not gap_length == 0:
                y2 = - center_length/2 - edge_length
                y3 = - center_length/2 - gap_length
                if c >= 0:
                    p = gds.PolyPath([(x,y2),(x,y3)],width=tempwidth, layer=layer, datatype=datatype)
                    idt.add(p)
                if c <= 0:
                    p = gds.PolyPath([(x,-y2),(x,-y3)],width=tempwidth, layer=layer, datatype=datatype)
                    idt.add(p)
           
        center_pos += current_period/2
        
    # start creating bondig pad
    overlap = 10        # Keep overlap of 10 [um] between the edge of IDT and pad
    space2GND = 20      # Distance between the IDT and GND plane
    padX_extra = 200    # Extra length for X of signal pad [um]
    padY = 60           # Y size of signal pad [um]
    
    pad = gds.Cell(cell_name+'_pad', exclude_from_current=True)
    
    # Protection of the core part of the IDT
    x0 = - total_length/2 - space2GND
    y0 = - center_length/2 - edge_length + overlap
    p = gds.Rectangle((x0,y0),(-x0,-y0), layer=inv_GND_layer, datatype=inv_GND_datatype)
    pad.add(p)
    
    # pad for signal pad
    x0 = - total_length/2 - space2GND - padX_extra
    y0 = center_length/2 + edge_length - overlap
    x1 = total_length/2 + space2GND
    y1 = y0 + padY
    p = gds.Rectangle((x0,y0),(x1,y1), layer=pad_layer, datatype=pad_datatype)
    pad.add(p)
    
    # Protection of the signal pad
    x0b = x0 - 2*space2GND
    y0b = y0 - 2*space2GND
    x1b = x1 + 2*space2GND
    y1b = y1 + 2*space2GND
    p = gds.Rectangle((x0b,y0b),(x1b,y1b), layer=inv_GND_layer, datatype=inv_GND_datatype)
    pad.add(p)
        
    return (idt, pad)


"""----------------
Utility functions
-------------------"""
def create_IDT_pair(name = 'pair',                      # Name of the cell
                    idt1 = createIDT3(),            # Left IDT
                    idt2 = createIDT3(),            # Right IDT
                    detector = None,                # Detector IDT (Not created if None)
                    distance = 100,                 # [ns]
                    saw_velocity = 2.77,           # [um/ns] SAW velocity
                    ):
    """
    Create a pair of IDTs for VNA measurement.
    - This function is supposed to be used with createIDT3().
    """
    pair = gds.Cell(name, exclude_from_current=True)    # Cell to hold pairs
#    distance_in_length = distance * saw_velocity        # Distance between front of IDTs in [um]
    distance_in_length = 130       # Distance between front of IDTs in [um]
    pattern = [pair, idt1[0], idt1[1], idt2[0], idt2[1]]
    
    # Add left IDT
    bounding_box = idt1[0].get_bounding_box()
    w = bounding_box[1][0] - bounding_box[0][0]
    l = bounding_box[1][1] - bounding_box[0][1]
    
    idt = gds.CellReference(idt1[0], origin=(-(distance_in_length+w)/2, 0))
    pair.add(idt)
    pad = gds.CellReference(idt1[1], origin=(-(distance_in_length+w)/2, 0))
#    pair.add(pad)
    
    # Add right IDT
    bounding_box = idt2[0].get_bounding_box()
    w = bounding_box[1][0] - bounding_box[0][0]
    l = max(l, (bounding_box[1][1] - bounding_box[0][1]))
    
    idt = gds.CellReference(idt2[0], origin=((distance_in_length+w)/2, 0))
    pair.add(idt)
    pad = gds.CellReference(idt2[1], origin=((distance_in_length+w)/2, 0), rotation=180, x_reflection=True)
#    pair.add(pad)
    
    # Get layer for pad and inverse GND from idt1
    #   Here we assume the layer number of pad < the layer number of inverse GND.
    #   Datatype of both is set to 0.
#    layers = pair.get_bounding_box()
    layers = list(idt1[1].get_layers())
    inv_GND_layer = layers[1]
#    pair.add(gds.Rectangle((-distance_in_length/2, -l/2),(distance_in_length/2, l/2),layer=inv_GND_layer, datatype=0))
    
    return pattern

def createMarks(nx = 4,
                ny = 16,
                dx = 3000,
                dy = 600,
                ):
    """
    Create marks at the junction of the lattice
    """
    # Mark for EB
    crmk = gds.Cell('cross', exclude_from_current=True)
    crmk.add(gds.Rectangle((-100,-10),(-10, 10),layer=2, datatype=0))
#    crmk.add(gds.Rectangle((-120, -30),(10, 30),layer=3, datatype=0))
    crmk.add(gds.Rectangle((100,-10),(10,10),layer=2, datatype=0))
#    crmk.add(gds.Rectangle((120,-30),(-10,30),layer=3, datatype=0))
    crmk.add(gds.Rectangle((-10,-100),(10,-10),layer=2, datatype=0))
#    crmk.add(gds.Rectangle((-30,-120),(30, 10),layer=3, datatype=0))
    crmk.add(gds.Rectangle((-10,100),(10,10),layer=2, datatype=0))
#    crmk.add(gds.Rectangle((-30,120),(30,-10),layer=3, datatype=0))
    crmk.add(gds.Rectangle((-4,0),(0,1),layer=2, datatype=0))
    crmk.add(gds.Rectangle((-1,1),(0,4),layer=2, datatype=0))
    crmk.add(gds.Rectangle((4,0),(0,-1),layer=2, datatype=0))
    crmk.add(gds.Rectangle((1,-1),(0,-4),layer=2, datatype=0))
    
    # Make array
    mkar = gds.CellArray(crmk, 2, 2, (2000,2000), origin=(0, 0))
    
    return crmk, mkar
    
"""----------------------------------------------------
Actual pattern for devices
----------------------------------------------------"""
    
 
    
    
def idt_20190407():
    """
    This program is used to create a set of different types of IDTs:
    ---- List of IDTs  -------------
    1. single finger IDTs:  100 nm, 125 nm and 250 nm            A60P60
    2. split 52 IDTs:       100 nm, 125 nm, 275 nm and 550 nm
    3. split 4 IDTs:        125 nm and 250 nm
    

    1, 2 Standard IDT without diffraction fingers [double finger, lambda = 1 um, period = 60, aperature = 30]
    3~64. Split52 IDT (lambda = 2.76 [um]), period = [40, 60, 80, 100, 120] * aperature = [30, 60, 90, 120] * diffraction = [no (10 um), no (20 um), 1*lambda]
    """
    dx = 500   # distance between set of IDTs [um]
    dy = 200    # distance between set of IDTs [um]
    nx = 4
    ny = 10
    space = 1000            # Space along the edge
    idt_distance = 300      # distance between 2 IDTs [ns]
    saw_velocity = 2.77    # [um/ns]
    
    # define top cell
    top = gds.Cell('TOP', exclude_from_current=True)
    
    # Add top cell to pattern list
    pattern = [top]
    
    types1 = np.array( [ 
              [2,2,2,2,2,2,2,2,2,2,2,2,2,2,11,10], # single finger and "special" double finger
              [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,10], # split 52
              [10,10,10,10,10,10,10,10,10,10,10,10,0,10,11,10], # 
              [10,10,10,10,10,10,10,10,10,10,10,10,0,10,2,10]  # split 4 and "special" IDTs
            ])
    
    types2 = np.array( [ 
              [2,2,2,2,2,2,2,2,2,2,2,2,2,2,10,2], # single finger and "special" double finger
              [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2], # split 52
              [11,11,11,11,11,11,11,11,11,11,11,11,11,11,10,2], # 
              [11,11,11,11,11,11,11,11,11,11,11,11,11,11,2,2]  # split 4 and "special" IDTs
            ])
    
    wavelengths = np.array( [ 
              [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,2,2,1,1,1,1],
              [0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.5,0.5,0.75,0.75],
              [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1,1,2,2,2,2],
              [2.0,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,1,1,2,2]
            ]) #[um]
    
    total_periods = np.array( [ 
              [90,90,90,90,90,90,90,90,90,90,110,110,110,110,110,110],
              [90,90,90,90,90,90,90,90,90,90,110,110,90,90,90,90],
              [110,110,110,110,110,110,110,110,90,90,90,90,90,90,90,90],
              [110,110,110,110,90,90,90,90,130,130,130,130,130,130,130,130]
            ]) #Number of periods - transmitter
    
#    total_periods2 = np.array( [ 
#              [90,90,90,90,90,90,90,90,110,110,110,110,110,110,110,110],
#              [130,130,130,130,130,130,130,130,110,110,110,110,90,90,90,90],
#              [110,110,110,110,110,110,110,110,90,90,90,90,90,90,90,90],
#              [110,110,110,110,90,90,90,90,130,130,130,130,130,130,130,130]
#            ]) #Number of periods - reciever
    
    apertures = np.array( [ 
              [30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30],
              [30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30],
              [30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30],
              [30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30]
            ]) #[um]
    
        # -- set jobfile content
    filename = 'GuttelJobfile_DoseTest2'
    path_to_save = 'C:/Users/manip.batm/Desktop/IDTs/Design'
    f = jobfile()
    f.addline('# run nbwrite Takada_Shintaro/Junliang/20190705/GuttelJobfile_0 -1=shintaro:g1  -2=shintaro:g2'.format(filename,path_to_save))
    f.addcl()

    # -- global block
    f.addline('.global')
    f.addel()
    f.addline(['registration',(0,0),(7600*1e3,0)]) #bottom left and bottom right marks
    f.addline(['marktype','car8tmf'])
    f.addline(['focus','auto'])
    f.addel()
    f.addline('.end')
    f.addel()
    
    # -- block default
    f.addcl() 
    f.addline('.block_default')
    f.addel()
    f.addline(['base_dose', 1.0])
    f.addel()
    f.addline('.end')
    f.addel()
    
    # define IDTs
    count = 0
    for i in range(nx):
        vertical_block_count = -1
        for j in range(ny):
            xc = dx/2+i*dx
            yc = dy/2+j*dy

            idt1 = createIDT3(period=[wavelengths[i,j], wavelengths[i,j]],     # [initial period (um), final periof (um) (> initial period)]
                              idt_type = types1[i,j],             # 1: single, 2: double 3: Split 52
                              center_length = apertures[i,j],       # [um]
                              edge_length = 65-apertures[i,j],         # [um]
                              gap_length = 0,           # [um] (if 0 there is no diffraction finger.)
                              metal_ratio = 1,        # Ratio set to compensate proximity
#                              period_in_time = periods_in_time[i,j],      # [ns] (period = ceil(period_in_time [ns] * saw_velocity [um/ns] / period [um]))
                              total_period = total_periods[i,j],    #Number of periods
                              saw_velocity = saw_velocity,     # [um/ns]
                              cell_name = 'didt_%d_%d' % (i,j),           # name of the cell
                              layer = 1,                # layer for idt
                              datatype = 0,             # datatype for idt
                              pad_layer = 5,            # layer number for bonding pad
                              pad_datatype = 0,         # datatype for bonding pad
                              inv_GND_layer = 3,        # layer number for inverse of GND plane
                              inv_GND_datatype = 0,     # datatype for inverse of GND plane
                              )
            
#            if types2[i,j] == 2:
#                total_period2 = 110
#            else:
            total_period2 = total_periods[i,j]
                
            idt2 = createIDT3(period=[wavelengths[i,j], wavelengths[i,j]],     # [initial period (um), final periof (um) (> initial period)]
                                  idt_type = types2[i,j],             # 1: single, 2: double 3: Split 52
                                  center_length = apertures[i,j],       # [um]
                                  edge_length = 65-apertures[i,j],         # [um]
                                  gap_length = 0,           # [um] (if 0 there is no diffraction finger.)
                                  metal_ratio = 1,        # Ratio set to compensate proximity
#                                  period_in_time = periods_in_time[i,j],      # [ns] (period = ceil(period_in_time [ns] * saw_velocity [um/ns] / period [um]))
                                  total_period = total_period2,    #Number of periods
                                  saw_velocity = saw_velocity,     # [um/ns]
                                  cell_name = 'didt_%d_%d_detect' % (i,j),           # name of the cell
                                  layer = 1,                # layer for idt
                                  datatype = 0,             # datatype for idt
                                  pad_layer = 5,            # layer number for bonding pad
                                  pad_datatype = 0,         # datatype for bonding pad
                                  inv_GND_layer = 3,        # layer number for inverse of GND plane
                                  inv_GND_datatype = 0,     # datatype for inverse of GND plane
                                  )
                
            pair = create_IDT_pair(name = 'didt_pair_%d_%d' % (i,j),                      # Name of the cell
                              idt1 = idt1,            # Left IDT
                              idt2 = idt2,            # Right IDT
                              detector = None,                # Detector IDT (Not created if None)
                              distance = idt_distance,                 # [ns]
                              saw_velocity = saw_velocity,           # [um/ns] SAW velocity
                              )
            # Add idt cell to pattern list
            pattern += pair
            idt1 = gds.CellReference(pair[0], origin=(xc, yc))
            top.add(idt1)
#            print('The center position of the (%d,%d) pair cell is (xc=%d,yc=%d)' %(i,j,xc,yc) )
#            print('The center position of the first IDT in the pair is (%d,%d)' %(xc-idt_distance*saw_velocity/2-total_periods[i,j]*wavelengths[i,j]/2,yc))
#            print('The center position of the second IDT in the pair is (%d,%d)' %(xc+idt_distance*saw_velocity/2+total_period2*wavelengths[i,j]/2,yc))
    
                #Block starters code
            block_starts = [0,5,8,11,15]
            block_ends = [4,7,10,14,15]
            x_file_origin = -3800
            y_file_origin = -4000
            if j in block_starts:
                vertical_block_count += 1
                x_origin =  x_file_origin + i*2000
                y_origin =  y_file_origin + vertical_block_count*2000
                f.addline('#'+'-'*30+'Block ('+str(i)+','+str(vertical_block_count)+')'+'-'*30)
                f.addline('.block')
                f.addline(['origin',((x_origin-x_file_origin)*1e3,(y_origin-y_file_origin)*1e3)])
                f.addline(['registration',(0.0, 0.0),(1600*1e3, 0.0)])
                f.addline(['marktype','car8tmf'])
                f.addline(['focus','auto'])
                f.addline(['base_dose',1.0])
                f.addel()
            
            #id's for IDT patterns
            xc_idt1 = (xc-idt_distance*saw_velocity/2-total_periods[i,j]*wavelengths[i,j]/2) #x position of the first IDT in the pair
            xc_idt2 = (xc+idt_distance*saw_velocity/2+total_period2*wavelengths[i,j]/2) #x position of the second IDT in the pair
            
            if types1[i,j] == 2:
                pattern_1 = 'double'
            elif types1[i,j] == 10:
                pattern_1 = 'DART'
            elif types1[i,j] == 11:
                pattern_1 = 'iDART'
                
            if types2[i,j] == 2:
                pattern_2 = 'double'
            elif types2[i,j] == 10:
                pattern_2 = 'DART'
            elif types2[i,j] == 11:
                pattern_2 = 'iDART'
                
            f.addline(['pattern',pattern_1+'_'+str(int(wavelengths[i,j]*1e3)).zfill(4)+'_'+str(total_periods[i,j]),((xc_idt1-x_origin)*1e3,(yc-y_origin)*1e3)])
            f.addline(['pattern',pattern_2+'_'+str(int(wavelengths[i,j]*1e3)).zfill(4)+'_'+str(total_periods[i,j]),((xc_idt2-x_origin)*1e3,(yc-y_origin)*1e3)])
            
            #Block enders code
            if j in block_ends:
                f.addel()
                f.addline('.end')
            
            
    #Pattern Section
    doses = np.array( [
                        [[1,1,1],[1,1,1],[1,1,1],[1,1,1]], #double IDT doses, 4 different wavelengths. 
                        [[1,1,1],[1,1,1],[1,1,1],[1,1,1]], #DART IDT doses, 4 different wavelengths.
                        [[1,1,1],[1,1,1],[1,1,1],[1,1,1]] #iDART IDT doses (same as DART)
                        ], dtype= np.float64) #Number of periods - transmitter
    
    for k,pattern_1 in enumerate(['double','DART','iDART']):
        for i,wavelength in enumerate([0.5,0.75,1.0,2.0]):
            for j,total_period in enumerate([90,110,130]):
                if (pattern_1 == 'iDART' and (wavelength in [0.5,0.75] or total_period == 130)):
                    continue
                f.addcl()
                f.addline('.pattern')
                f.addline(['id',pattern_1+'_'+str(int(wavelength*1e3)).zfill(4)+'_'+str(total_period)])
                f.addline(['filename', 'shintaro/Junliang/20190705/'+pattern_1+'_'+str(int(wavelength*1e3)).zfill(4)+'_'+str(total_period)+'.npf'])
                f.addline(['dose', str(1), doses[k,i,j]])
                f.addel()
                f.addline('.end')
        
                
            
    # -- end of file
    f.addcl()
    f.addline('.write')
    f.addline('current\tauto')
    f.addline('.write')
    
    # -- preview
    f.generatepreview()
#    print(f)

    # -- save file
    if True:
        f.generatefile(filename,path_to_save)  
    
    # GND plane before boolean
#    top.add(gds.Rectangle((-dx*nx/2-space, -dy*ny/2-space),(dx*nx/2+space, dy*ny/2+space), layer=4, datatype=0))
    top.add(gds.Rectangle((0, 0),(2000, 2000), layer=4, datatype=0))
    # Add square at corners for EB positioning
    top.add(gds.Rectangle((-dx*nx/2-space, -dy*ny/2-space),(-dx*nx/2-space+0.05, -dy*ny/2-space+0.05), layer=1, datatype=0))
    top.add(gds.Rectangle((-dx*nx/2-space, dy*ny/2+space),(-dx*nx/2-space+0.05, dy*ny/2+space-0.05), layer=1, datatype=0))
    top.add(gds.Rectangle((dx*nx/2+space, -dy*ny/2-space),(dx*nx/2+space-0.05, -dy*ny/2-space+0.05), layer=1, datatype=0))
    top.add(gds.Rectangle((dx*nx/2+space, dy*ny/2+space),(dx*nx/2+space-0.05, dy*ny/2+space-0.05), layer=1, datatype=0))
    
    # Make EB marks
    crmk, mkar = createMarks(nx = nx,
                            ny = ny,
                            dx = dx,
                            dy = dy,
                            )
    pattern.append(crmk)
    top.add(mkar)
    
    saveCell2GDS(pattern, 'idt_DoseTest_2')
    
if __name__=='__main__':
#    saveCell2GDS(createIDT2(), 'didt125')
#    saveCell2GDS(ebDummy(),'test')
    
#    idt_20181031()
    idt_20190407()