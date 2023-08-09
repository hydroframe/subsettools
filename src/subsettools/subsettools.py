from datetime import datetime, timedelta
import shutil
import numpy as np
import xarray as xr
from hydrodata.national_mapping.map_wgs84 import ConusMap
from parflow import Run
from parflow.tools.io import read_clm, read_pfb, write_pfb
import pytz
import os
import pathlib
import subprocess
import glob
from hydrodata.data_catalog import data_access

def get_conus_hucs_indices(huc_list, grid):
    huc_len = len(huc_list[0])
    assert all([len(huc) == huc_len for huc in huc_list]), "All huc IDs should have the same length!"
    huc_list = [int(huc) for huc in huc_list]
    entry = data_access.get_catalog_entry(
        dataset = "huc_mapping",
        grid=grid.lower(), 
        file_type = "tiff"
    )

    conus_hucs = data_access.get_ndarray(entry,level = str(huc_len))
    
    sel_hucs = np.isin(conus_hucs, huc_list).squeeze()

    indices_j, indices_i = np.where(sel_hucs > 0)

    return conus_hucs, sel_hucs, indices_j, indices_i

def indices_to_ij(conus_hucs, indices_j, indices_i):
    # # Find the boundaries in the i and j directions
    imin = np.min(indices_i)
    imax = np.max(indices_i) + 1 #right bound inclusive
    arr_jmin = np.min(indices_j)
    arr_jmax = np.max(indices_j) + 1 #right bound inclusive
    
    jmin = conus_hucs.shape[1]-arr_jmax 
    jmax = conus_hucs.shape[1]-arr_jmin

    ij_bounds = [imin,jmin,imax,jmax]  
  
    return ij_bounds

def huc_to_ij(huc_list, grid):
    conus_hucs, sel_hucs, indices_j, indices_i = get_conus_hucs_indices(huc_list, grid)
    
    ij_bounds = indices_to_ij(conus_hucs, indices_j, indices_i)

    return ij_bounds


def latlon_to_ij(latlng_bounds,grid):
    conus_map = ConusMap(grid.lower()) # Creating a ConusMap object
    point0 = conus_map.map_to_grid(latlng_bounds[0][1], latlng_bounds[0][0])  # The ConusMap object is used to transform from lat-lon to i-j indices
    point1 = conus_map.map_to_grid(latlng_bounds[1][1], latlng_bounds[1][0])  # from lat-lon to i-j indexes
    imin, imax = [min(point0[0], point1[0]) , max(point0[0], point1[0])]    # Retrieving xmin, and xmax indices
    jmin, jmax = [min(point0[1], point1[1]) , max(point0[1], point1[1])]    # Retrieving ymin, ymax indices
    
    ij_bounds = [imin,jmin,imax,jmax]
    return ij_bounds


def create_mask_solid(huc_list, grid, write_dir):
    conus_hucs, sel_hucs, indices_j, indices_i = get_conus_hucs_indices(huc_list, grid)
    arr_jmin = np.min(indices_j)
    arr_jmax = np.max(indices_j) + 1 #right bound inclusive
    ij_bounds = indices_to_ij(conus_hucs, indices_j, indices_i)

    nj = ij_bounds[3]-ij_bounds[1]
    ni = ij_bounds[2]-ij_bounds[0]
        
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
    mask_clip[0,:,:] = sel_hucs[arr_jmin:arr_jmax,ij_bounds[0]:ij_bounds[2]] #we need to use numpy iymin / iymax because we are subsetting the tif file
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
    
    
def subset_static(ij_bounds, dataset, write_dir, var_list=['slope_x','slope_y','pf_indicator','mannings','depth_to_bedrock','pme','ss_pressure_head']): 
    #getting paths and writing subset pfbs for static parameters
    for var in var_list: 
        entry = data_access.get_catalog_entry(
            dataset = dataset,
            file_type = "pfb",
            period = "static",
            variable = var
        )
        if entry is not None:
            subset_data = data_access.get_ndarray(entry, grid_bounds = ij_bounds)
            write_pfb(f'{write_dir}/{var}.pfb', subset_data)
            print(f"Wrote {var}.pfb in specified directory.")
        
        else:
            print(f"{var} not found in dataset {dataset}")

def subset_press_init(ij_bounds, dataset, date, write_dir, time_zone = 'UTC'):
    entry = data_access.get_catalog_entry(
        dataset = dataset,
        file_type = "pfb",
        variable = "pressure_head", 
        period = "hourly"
    )
    
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
    
    subset_data = data_access.get_ndarray(entry, grid_bounds = ij_bounds, start_time=new_date)
    out_file_path = f'{write_dir}/{dataset}_{date_string}_press.pfb'
    out_file_name = f'{dataset}_{date_string}_press.pfb'
    if subset_data.size != 0:
        write_pfb(out_file_path, subset_data)
        print(f"Wrote {dataset}_{date_string}_press.pfb in specified directory.")
    else:
        print(f"No pressure file found for {new_date} in dataset {dataset}")
    return out_file_name
            
    
def config_clm(ij_bounds, start, end, dataset, write_dir):  
    file_type_list = ['vegp','vegm', 'drv_clm']
    for file_type in file_type_list:
        print(f"processing {file_type}")    
        entry = data_access.get_catalog_entries(dataset=dataset, file_type=file_type, variable="clm_run", period="static")[0]
        file_path = data_access.get_file_paths(entry)[0]
        print(file_path)
        if file_type == 'vegp':
            shutil.copyfile(file_path, f'{write_dir}/drv_vegp.dat')
            print(f"copied vegp")
        elif file_type == 'vegm':
            land_cover_file = subset_vegm(file_path, ij_bounds) 
            write_land_cover(land_cover_file, f'{write_dir}/drv_vegm.dat')
            print(f"subset vegm")
        elif file_type == 'drv_clm':
            edit_drvclmin(read_path = file_path, write_path = f'{write_dir}/drv_clmin.dat',start = start,end=end)
            print(f"edited drv_clmin")

            
def subset_vegm(path,ij_bounds):
    
    #read in the target vegm file
    vegm = read_clm(path, type = 'vegm')  #returns (j,i,k)
    vegm = np.transpose(vegm, (2,0,1))    #transpose to k,j,i

    imin,jmin,imax,jmax = ij_bounds
    #slice based on the i,j indices
    vegm = vegm[:, jmin:jmax, imin:imax] #slicing on k,j,i
    
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
    return out_file


def edit_drvclmin(read_path, write_path=None, start=None, end=None, startcode=2, vegp_name='drv_vegp.dat', vegm_name='drv_vegm.dat'):
    if write_path is not None:
        shutil.copyfile(read_path, write_path)
        lines = open(write_path, 'r').readlines()
    
    else:
        write_path = read_path
        lines = open(read_path, 'r').readlines()
    
    for i, line in enumerate(lines):
        if 'vegtf' in line:
            lines[i] = f'vegtf           {vegm_name}                         Vegetation Tile Specification File\n'
        elif 'vegpf' in line:
            lines[i] = f'vegpf           {vegp_name}                         Vegetation Type Parameter\n'            
        elif 'startcode' in line:
            lines[i] = f'startcode       {startcode}                                    1=restart file, 2=defined\n'
        elif 'clm_ic' in line:
            lines[i] = f'clm_ic          {startcode}                                    1=restart file, 2=defined\n'
            
    if start is not None:        
        startdt = datetime.strptime(start, '%Y-%m-%d')
        sd = startdt.strftime('%d')
        sm = startdt.strftime('%m')
        enddt = datetime.strptime(end, '%Y-%m-%d')
        ed = enddt.strftime('%d')
        em = enddt.strftime('%m')

        for i, line in enumerate(lines):
            if 'sda' in line:
                lines[i] = f'sda            {sd}                                    Starting Day\n'
            elif 'smo' in line:
                lines[i] = f'smo            {sm}                                    Starting Month\n'
            elif 'syr' in line:
                lines[i] = f'syr            {startdt.year}                                  Starting Year\n'
            elif 'eda' in line:
                lines[i] = f'eda            {ed}                                    Ending Day\n'
            elif 'emo' in line:
                lines[i] = f'emo            {em}                                    Ending Month\n'
            elif 'eyr' in line:
                lines[i] = f'eyr            {enddt.year}                                  Ending Year\n'

    open(write_path, 'w').writelines(lines)
    return write_path


def subset_forcing(ij_bounds, grid, start, end, dataset, write_dir):
    # CLM variables requested
    var_list = [
        "precipitation",
        "downward_shortwave",
        "downward_longwave",
        "specific_humidity",
        "air_temp",
        "atmospheric_pressure",
        "east_windspeed",
        "north_windspeed"
    ]        
    outputs = {}
    
    for var in var_list:   
        entries = data_access.get_catalog_entry(
            dataset = dataset,
            variable = var,
            grid = grid, 
            file_type = "pfb",
            period = "hourly"
        )
        
        print(f'Reading {var} pfb sequence')
        subset_data = data_access.get_ndarray(entries, start_time = start, end_time = end, grid_bounds = ij_bounds)
        print(f'Done reading {var} pfb sequence, starting to write to folder')
        
        paths = data_access.get_file_paths(entries, start_time=start, end_time = end)
        out_paths = [f'{write_dir}/{os.path.basename(p)}' for p in paths]
        for i, op in enumerate(out_paths):
            write_pfb(op, subset_data[i,:,:,:])
        
        print(f'finished writing {var} to folder')
        outputs[var] = out_paths
    return outputs


def edit_runscript_for_subset(ij_bounds, runscript_path, write_dir=None, runname=None, forcing_dir=None):
    if write_dir is None:     
        write_dir = os.path.dirname(runscript_path) 
        
    file_extension = os.path.splitext(runscript_path)[1][1:]

    imin, jmin, imax, jmax = ij_bounds
    #getting the subset ni/nj to update keys
    nj = jmax - jmin
    ni = imax - imin
    
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
    return run.write(working_directory=write_dir, file_format=f'{file_extension}')[0]

    
def copy_static_files(static_input_dir, pf_dir):
    os.system('cp -r '+static_input_dir+'*.* '+pf_dir)
    
        
def change_filename_values(
    runscript_path, 
    write_dir = None, 
    runname=None, 
    slopex=None, 
    slopey=None, 
    solidfile=None, 
    ip=None, 
    indicator=None, 
    depth_to_bedrock=None, 
    mannings=None, 
    evap_trans=None
):
    
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
        run.Solver.EvapTrans.FileName = evap_trans
        print(f"Evaptrans filename changed to {evap_trans}")
        
    
    print(f"Updated runscript written to {write_dir} as detected file extension")
    return run.write(working_directory=write_dir, file_format=f'{file_extension}')[0]

    
def dist_run(P, Q, runscript_path, write_dir, dist_clim_forcing=True):
    assert write_dir is not None

    nz = 0 #starting nz as zero 

    file_name, file_extension = os.path.splitext(runscript_path)
    file_extension = file_extension[1:]
            
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
        
    static_input_paths = list(pathlib.Path(f'{write_dir}').glob('*.pfb'))
    
    for path in static_input_paths:
        input_array = read_pfb(path)
        run.ComputationalGrid.NZ = int(input_array.shape[0])
        
        if int(input_array.shape[0]) > nz: 
            nz = int(input_array.shape[0])
            
        run.dist(path)
        print(f"Distributed {os.path.basename(path)} with NZ {int(input_array.shape[0])}")
        
    run.ComputationalGrid.NZ = nz

    run.write(working_directory=write_dir, file_format=f'{file_extension}') 

    
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
