from datetime import datetime, timedelta
import shutil
import numpy as np
import parflow
import sqlite3
import logging
import xarray as xr
from contextlib import closing
from hydrodata.national_mapping.map_wgs84 import ConusMap
from parflow import Run
from parflow.tools.io import read_clm, read_pfb, read_pfb_sequence, write_pfb
import pytz
import os
import pathlib
import subprocess
import sys
import glob
#sys.path.insert(1,'/home/SHARED/data/hydrodata/hydrodata')
import hydrodata.data_catalog.data_access



def get_conus_ij(domain, grid):
    #Eventually add a check if we are getting a shapefile, tif, etc. 
    #for now this can handle huc inputs as string and lat/lon bbox
    if isinstance(domain, str):
        conus_ij = huc_to_ij(domain, grid)
        
    else:
        conus_ij = latlon_to_ij(domain, grid)
        
    
    return conus_ij
        
def huc_to_ij(huc_id, grid):
    huc_len = len(huc_id)
    #****this path to the conus huc tifs will need to be dealt with better...
    conus_hucs = xr.open_dataset(f'/hydrodata/national_mapping/{grid.upper()}/HUC{huc_len}_{grid.upper()}_grid.tif').drop_vars(('x', 'y'))['band_data']
    huc = int(huc_id)
    sel_huc = (conus_hucs == huc).squeeze()
    
    # First find where along the y direction has "valid" cells
    y_mask = (sel_huc.sum(dim='x') > 0).astype(int)
    
    # Then, taking a diff along that dimension let's us see where the boundaries of that mask ar
    diffed_y_mask = y_mask.diff(dim='y')

    # Taking the argmin and argmax get's us the locations of the boundaries
    arr_jmax = np.argmin(diffed_y_mask.values) + 1 #this one is because you want to include this right bound in your slice
    arr_jmin = np.argmax(diffed_y_mask.values) + 1 #because of the point you actually want to indicate from the diff function

    jmin = conus_hucs.shape[1]-arr_jmax 
    jmax = conus_hucs.shape[1]-arr_jmin

    # Do the exact same thing for the x dimension
    diffed_x_mask = (sel_huc.sum(dim='y') > 0).astype(int).diff(dim='x')
    imax = np.argmin(diffed_x_mask.values) + 1
    imin = np.argmax(diffed_x_mask.values) + 1

    ij_bounds = [imin,jmin,imax,jmax]  
  
    return ij_bounds


def latlon_to_ij(latlng_bounds,grid):
    conus_map = ConusMap(grid.lower()) # Creating a ConusMap object
    point0 = conus_map.map_to_grid(latlng_bounds[0][1], latlng_bounds[0][0])  # The ConusMap object is used to transform from lat-lon to i-j indices
    point1 = conus_map.map_to_grid(latlng_bounds[1][1], latlng_bounds[1][0])  # from lat-lon to i-j indexes
    imin, imax = [min(point0[0], point1[0]) , max(point0[0], point1[0])]    # Retrieving xmin, and xmax indices
    jmin, jmax = [min(point0[1], point1[1]) , max(point0[1], point1[1])]    # Retrieving ymin, ymax indices
    
    ij_bounds = [imin,jmin,imax,jmax]
    return ij_bounds

def create_mask_solid(huc_id, grid, write_dir):
    if isinstance(huc_id, str):
        print("Provided HUC ID")
        huc_len = len(huc_id)
        
        conus_hucs = xr.open_dataset(f'/hydrodata/national_mapping/{grid.upper()}/HUC{huc_len}_{grid.upper()}_grid.tif').drop_vars(('x', 'y'))['band_data']
        huc = int(huc_id)
        sel_huc = (conus_hucs == huc).squeeze()

        # First find where along the y direction has "valid" cells
        y_mask = (sel_huc.sum(dim='x') > 0).astype(int)

        # Then, taking a diff along that dimension let's us see where the boundaries of that mask ar
        diffed_y_mask = y_mask.diff(dim='y')

        # Taking the argmin and argmax get's us the locations of the boundaries
        #the location of these boundaries are the ymin/ymax on a numpy NOT parflow grid
        arr_jmax = np.argmin(diffed_y_mask.values) + 1 #this one is because you want to include this right bound in your slice
        arr_jmin = np.argmax(diffed_y_mask.values) + 1 #because of the point you actually want to indicate from the diff function

        #This essentially flips the grid over the x-axis so that we get on the parflow oriented grid
        jmin = conus_hucs.shape[1]-arr_jmax #1888 for conus1
        jmax = conus_hucs.shape[1]-arr_jmin

        # Do the exact same thing for the x dimension
        diffed_x_mask = (sel_huc.sum(dim='y') > 0).astype(int).diff(dim='x')
        imax = np.argmin(diffed_x_mask.values) + 1
        imin = np.argmax(diffed_x_mask.values) + 1

        nj = jmax-jmin
        ni = imax-imin
        
        #checks conus1 / 2 grid and assigns appripriate dz and z_total for making the mask and solid file
        if grid.lower() == "conus1": 
            print("grid is conus1")
            layz = 100
            z_total = str(500)
        else: 
            print("grid is conus2")
            layz = 200
            z_total = str(2000)
        #could look this up in dC and get the information

        #create and write the pfb mask
        mask_clip = np.zeros((1,nj,ni))
        mask_clip[0,:,:] = sel_huc[arr_jmin:arr_jmax,imin:imax] #we need to use numpy iymin / iymax because we are subsetting the tif file
        mask_clip = np.flip(mask_clip,1) #This flip tooks the section we just subset and puts it in the appropriate parflow orientation
        mask_clip = mask_clip.astype(float)
        write_pfb(f'{write_dir}/mask.pfb', mask_clip, dx=1000, dy=1000, dz=layz, dist = False)
        
        print("Wrote mask.pfb")
        
        subprocess.run([
            os.path.join(os.environ["PARFLOW_DIR"], 'bin', 'pfmask-to-pfsol'),
            '--mask', f'{write_dir}/mask.pfb',
            '--pfsol', f'{write_dir}/solidfile.pfsol',
            '--vtk', f'{write_dir}/mask_vtk.vtk',
            '--z-bottom', '0.0',
            '--z-top', z_total
        ], capture_output=True)
        
        print(f"Wrote solidfile.pfsol and mask_vtk.vtk with total z of {z_total} meters")
        
    else:
        print("Unsupported domain boundary, only HUC IDs are supported currently to produce solid files. \nRun ParFlow with a box domain template. ")
    
def subset_static(ij_bounds, dataset, write_dir, var_list=['slope_x','slope_y','pf_indicator','mannings','depth_to_bedrock','pme']): 
    #getting paths and writing subset pfbs for static parameters
    for var in var_list: 
        entry = hydrodata.data_catalog.data_access.get_catalog_entry(dataset = dataset,
                                                                     file_type = "pfb",
                                                                     period = "static",
                                                                     variable = var)
        if entry is not None:
            subset_data = hydrodata.data_catalog.data_access.get_ndarray(entry, grid_bounds = ij_bounds)
            write_pfb(f'{write_dir}/{var}.pfb', subset_data)
            print(f"Wrote {var}.pfb in specified directory.")
        
        else:
            print(f"{var} not found in dataset {dataset}")

def subset_press_init(ij_bounds, dataset, date, write_dir, time_zone = 'UTC'):
    entry = hydrodata.data_catalog.data_access.get_catalog_entry(dataset = dataset,
                                                                 file_type = "pfb",
                                                                 variable = "pressure_head", 
                                                                 period = "hourly")
    
    #getting the correct end day of the previous water year to be the init press
    first_date = datetime.strptime(date, '%Y-%m-%d')
    new_date = first_date - timedelta(hours=1) #assumes time is UTC 0 like CONUS runs, so can remain time unaware and grab the right pressure

    if time_zone != 'UTC':
        print(f"Time zone provided, converting the requested datetime from UTC0 to {time_zone}")
        new_date = new_date.replace(tzinfo=pytz.UTC) #add time awareness as UTC
        new_date = new_date.astimezone(pytz.timezone(time_zone)) #convert to provided timezone
        date_string = new_date.strftime("%Y.%m.%d:%H.%M.%S_UTC0")
        new_date = new_date.strftime("%Y-%m-%d %H:%M:%S")
        
    else:
        date_string = new_date.strftime("%Y.%m.%d:%H.%M.%S_UTC0")
        new_date = new_date.strftime("%Y-%m-%d %H:%M:%S")
    
    subset_data = hydrodata.data_catalog.data_access.get_ndarray(entry, grid_bounds = ij_bounds, start_time=new_date)
    
    if subset_data.size != 0:
        write_pfb(f'{write_dir}/{dataset}_{date_string}_press.pfb', subset_data)

        print(f"Wrote {dataset}_{date_string}_press.pfb in specified directory.")
        
    else:
        print(f"No pressure file found for {new_date} in dataset {dataset}")
            
    
def config_clm(ij_bounds, start, end, dataset, write_dir):  
    file_type_list = ['vegp','vegm', 'drv_clm']

    for file in file_type_list:
        print(file)
        if 'vegp' in file:
            entry = hydrodata.data_catalog.data_access.get_catalog_entry(dataset = dataset,
                                                                        file_type = "vegp",
                                                                        variable = "clm_run",
                                                                        period = "static")
            path = hydrodata.data_catalog.data_access.get_file_paths(entry)
            print(path[0])
            shutil.copyfile(path[0], f'{write_dir}/drv_vegp.dat')
            print(f"copied vegp")
            
        elif 'vegm' in file:
            entry = hydrodata.data_catalog.data_access.get_catalog_entry(dataset = dataset,
                                                                        file_type = "vegm",
                                                                        variable = "clm_run",
                                                                        period = "static")
            path = hydrodata.data_catalog.data_access.get_file_paths(entry)
            print(f"subsetting vegm (this takes awhile)")
            vegm = subset_vegm(path[0],ij_bounds) 
            write_land_cover(vegm, f'{write_dir}/drv_vegm.dat')
            print(f"subset vegm")
            
        elif 'drv' in file:
            entries = hydrodata.data_catalog.data_access.get_catalog_entries(dataset = dataset,
                                                                        file_type = "drv_clm",
                                                                        variable = "clm_run",
                                                                        period = "static")
            path = hydrodata.data_catalog.data_access.get_file_paths(entries[0])
            print(path[0])

            edit_drvclmin(read_path = path[0],write_path = f'{write_dir}/drv_clmin.dat',start = start,end=end)
            print(f"edited drv_clmin")
            
        else:
            print(f"error with file {file}")
    
def subset_vegm(path,ij_bounds):
    
    #read in the target vegm file
    vegm = read_clm(path, type = 'vegm') #returns (j,i,k)
    
    vegm = np.transpose(vegm, (2,0,1)) # transpose to k,j,i
    
    #slice based on the i,j indices
    vegm = vegm[:,(ij_bounds[1]):(ij_bounds[3]),(ij_bounds[0]):(ij_bounds[2])] #slicing on k,j,i
    
    #generate the i, j indices necessary for the vegm file based on the shape of the subset data
    nj, ni = vegm.shape[1:]
    indices = np.indices((nj, ni)) + 1 
    indices = indices[::-1, :, :]
    vegm = np.vstack([indices, vegm])#stack x,y indices on vegm

    #transpose and reshape back into expected 2D vegm file format for the subset
    vegm = vegm.transpose(1,2,0).reshape(-1, 25)
    
    return vegm

def write_land_cover(land_cover_data, out_file):
    """Write the land cover file in vegm format
    Parameters
    ----------
    land_cover_data : ndarray
        formatted vegm data (2d array)
    out_file : str
        path and name to write output
    Returns
    -------
    None
    """
    heading = "x y lat lon sand clay color fractional coverage of grid, by vegetation class (Must/Should Add to " \
              "1.0) "
    col_names = ['', '', '(Deg)', '(Deg)', '(%/100)', '', 'index', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                 '10', '11',
                 '12', '13', '14', '15', '16', '17', '18']
    header = '\n'.join([heading, ' '.join(col_names)])
    np.savetxt(fname=out_file, X=land_cover_data, delimiter=' ', comments='', header=header,
               fmt=['%d'] * 2 +['%.6f'] * 2 + ['%.2f'] * 2 + ['%d'] * 19)


def edit_drvclmin(read_path, write_path=None, start=None, end=None, startcode=2, vegp_name='drv_vegp.dat', vegm_name='drv_vegm.dat'):
    if write_path is not None:
        shutil.copyfile(read_path, write_path)
        lines = open(write_path, 'r').readlines()
    
    else:
        write_path = read_path
        lines = open(read_path, 'r').readlines()
    
    for num, line in enumerate(lines):
        if 'vegtf' in line:
            lines[num] = f'vegtf           {vegm_name}                         Vegetation Tile Specification File\n'
        if 'vegpf' in line:
            lines[num] = f'vegpf           {vegp_name}                         Vegetation Type Parameter\n'
            
        if 'startcode' in line:
                lines[num] = f'startcode       {startcode}                                    1=restart file, 2=defined\n'
        if 'clm_ic' in line:
            lines[num] = f'clm_ic          {startcode}                                    1=restart file, 2=defined\n'
            
    if start is not None:
        
        startdt = datetime.strptime(start, '%Y-%m-%d')
        sd = startdt.strftime('%d')
        sm = startdt.strftime('%m')
        enddt = datetime.strptime(end, '%Y-%m-%d')
        ed = enddt.strftime('%d')
        em = enddt.strftime('%m')

        for num, line in enumerate(lines):
            if 'sda' in line:
                lines[num] = f'sda            {sd}                                    Starting Day\n'
            if 'smo' in line:
                lines[num] = f'smo            {sm}                                    Starting Month\n'
            if 'syr' in line:
                lines[num] = f'syr            {startdt.year}                                  Starting Year\n'

            if 'eda' in line:
                lines[num] = f'eda            {ed}                                    Ending Day\n'
            if 'emo' in line:
                lines[num] = f'emo            {em}                                    Ending Month\n'
            if 'eyr' in line:
                lines[num] = f'eyr            {enddt.year}                                  Ending Year\n'

                
    open(write_path, 'w').writelines(lines)

def subset_forcing(ij_bounds, grid, start, end, dataset, write_dir):
    var_list = ["precipitation","downward_shortwave","downward_longwave","specific_humidity","air_temp","atmospheric_pressure","east_windspeed","north_windspeed"]        # CLM variables requested
    
    for var in var_list:   
        entries = hydrodata.data_catalog.data_access.get_catalog_entry(dataset = dataset,
                                                                         variable = var,
                                                                         grid = grid, 
                                                                         file_type = "pfb",
                                                                         period = "hourly")
        
        print(f'Reading {var} pfb sequence')
        subset_data = hydrodata.data_catalog.data_access.get_ndarray(entries, start_time = start, end_time = end, grid_bounds = ij_bounds)
        print(f'Done reading {var} pfb sequence, starting to write to folder')
        
        paths = hydrodata.data_catalog.data_access.get_file_paths(entries, start_time=start, end_time = end)
        for path in enumerate(paths):
            write_pfb(f'{write_dir}/{os.path.basename(path[1])}', subset_data[path[0],:,:,:])
        
        print(f'finished writing {var} to folder')
    
def edit_runscript_for_subset(ij_bounds, runscript_path, write_dir=None, runname=None, forcing_dir=None):
    if write_dir is None:     
        write_dir = os.path.dirname(runscript_path) 
        
    file_name, file_extension = os.path.splitext(runscript_path)
    file_extension = file_extension[1:]
    #getting the subset ni/nj to update keys
    nj = ij_bounds[3]-ij_bounds[1]
    ni = ij_bounds[2]-ij_bounds[0]
    
    #load in the reference pfidb or yaml specified by the user
    run = Run.from_definition(runscript_path)
    
    if runname is not None: 
        run.set_name(runname)
        print(f"New runname: {runname} provided, a new {file_extension} file will be created")
    else: 
         print(f"No runname provided, old {file_extension} file will be overwritten")
        
    
    run.ComputationalGrid.NY = int(nj)
    run.ComputationalGrid.NX = int(ni) 
    
    print(f"ComputationalGrid.NY set to {nj} and NX to {ni}")
    
    #checks if we're running with clm    
    if forcing_dir is not None:
        print(f"Old path to climate forcing was {run.Solver.CLM.MetFilePath} and has been changed to {forcing_dir} in runscript.")
        
        run.Solver.CLM.MetFilePath = forcing_dir
           
    else:
        print("No forcing directory provided, key not set")
        
    #Checking if we are solid or box  
    domain_type = run.GeomInput.domaininput.InputType

    if domain_type == "SolidFile":
        print(f"GeomInput.domaininput.InputType detected as SolidFile, no additional keys to change for subset")
    
    else:
        print(f"GeomInput.domaininput.InputType detected as Box, updating Geom.domain.Upper.X to {ni * 1000} and Y to {nj * 1000} to match subset")
        
        run.Geom.domain.Upper.X = ni * 1000
        run.Geom.domain.Upper.Y = nj * 1000
        
    
    print(f"Updated runscript written to {write_dir} as detected extension")
    run.write(working_directory=write_dir, file_format=f'{file_extension}') 
    
def copy_static_files(static_input_dir, pf_dir):
    #It's just one line, but it doesn't really seem to fit with the otehr functions... I do think it would be more intuitive to have it in a function but not sure
    os.system('cp -r '+static_input_dir+'*.* '+pf_dir)
    
        
def change_filename_values(runscript_path, write_dir = None, runname=None, slopex=None, slopey=None, solidfile=None, ip=None, indicator=None, depth_to_bedrock=None, mannings=None, evap_trans=None):
    
    file_name, file_extension = os.path.splitext(runscript_path)
    file_extension = file_extension[1:]
    
    if write_dir is None:     
        write_dir = os.path.dirname(runscript_path) 
        print(f"No write directory provided, updated or new {file_extension} file will be written to the runscript path")
    
    run = Run.from_definition(runscript_path)
    
    if runname is not None: 
        run.set_name(runname)
        print(f"New runname: {runname} provided, a new {file_extension} file will be created")
    
    #check which input files are not none, to update the key 
    if slopex is not None:
        run.TopoSlopesX.FileName = slopex
        print(f"X Slopes filename changed to {slopex}")
    if slopey is not None: 
        run.TopoSlopesY.FileName = slopey
        print(f"Y Slopes filename changed to {slopey}")
    if solidfile is not None: 
        run.GeomInput.domaininput.FileName = solidfile
        print(f"Solidfile filename changed to {solidfile}")
    if ip is not None: 
        run.Geom.domain.ICPressure.FileName = ip
        print(f"Initial pressure filename changed to {ip}")
    if indicator is not None:
        run.Geom.indi_input.FileName = indicator
        print(f"Indicator filename changed to {indicator}")
    if depth_to_bedrock is not None: 
        run.Geom.domain.FBz.FileName = depth_to_bedrock
        print(f"Depth to bedrock filename changed to {depth_to_bedrock}")
    if mannings is not None: 
        run.Mannings.FileName = mannings
        print(f"Mannings filename changed to {mannings}")
    if evap_trans is not None: 
        run.Solver.EvapTransFile = evap_trans
        print(f"Evaptrans filename changed to {evap_trans}")
        
    
    print(f"Updated runscript written to {write_dir} as detected file extension")
    run.write(working_directory=write_dir, file_format=f'{file_extension}') 

    
def dist_run(P, Q, runscript_path, pf_run_dir, dist_clim_forcing=False):
    if P != Q:
        print(f"Processor P={P} and Q={Q}, they must be equal.")
  
    run = Run.from_definition(runscript_path)
    
    run.Process.Topology.P = int(P)
    run.Process.Topology.Q = int(Q)
    
    if dist_clim_forcing is True: 
        print("Distributing your climate forcing")
        forcing_dir = run.Solver.CLM.MetFilePath
        for filename_forcing in os.listdir(forcing_dir):
            if filename_forcing[-3:]=='pfb':
                run.dist(f'{forcing_dir}{filename_forcing}')
    else: 
        print("no forcing dir provided, only distributing static inputs")
        
    static_input_paths = list(pathlib.Path(f'{pf_run_dir}').glob('*.pfb'))
    
    for path in static_input_paths:
        input_array = read_pfb(path)
        run.ComputationalGrid.NZ = int(input_array.shape[0])
        run.dist(path)
        print(f"Distributed {os.path.basename(path)} with NZ {int(input_array.shape[0])}")
    
def create_job_script():
    pass

def restart_run(runscript_path):
    pf_dir = os.path.dirname(runscript_path)  
    file_name, file_extension = os.path.splitext(runscript_path)
    file_extension = file_extension[1:]
    run = Run.from_definition(runscript_path)   
    runname = run.get_name()
    
    dump_interval = run.TimingInfo.DumpInterval
    
    path_to_tcl = pf_dir + '/clm_restart.tcl'
    lines = open(path_to_tcl, 'r').readlines()[0]
    istep = [int(i) for i in lines.split() if i.isdigit()][0]
    print(istep)
    
    press_base = runname + ".out.press.{:05d}.pfb"
    print(press_base)
    restart_press = press_base.format(istep-2)
    run.Geom.domain.ICPressure.FileName = restart_press
    print(f"Initial Press filename changed to: {run.Geom.domain.ICPressure.FileName}")
    
    #update pf timing (needs to be done no matter how you're running pf)
    run.TimingInfo.StartCount = istep-1
    run.TimingInfo.StartTime = (istep-1)*dump_interval
    
    print(f"start time is now: {run.TimingInfo.StartTime}")
    print(f"start count is now: {run.TimingInfo.StartCount}")
    
    if run.Solver.LSM == 'CLM':      
        edit_drvclmin(read_path = f'{pf_dir}/drv_clmin.dat', startcode=1)
        
        print(f"Overwrote drv_clmin.dat (changed startcode to 1 (restart file))")
        
    run.write(working_directory=pf_dir, file_format=f'{file_extension}')   
            
if __name__ == "__main__":
    main()
