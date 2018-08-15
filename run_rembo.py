#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Script to run one rembo simulation.
'''
import subprocess as subp 
import sys, os, argparse, shutil, glob, datetime

def run_rembo():
    '''Main subroutine to run rembo.'''

    ### Manage command-line arguments ############################
    parser = argparse.ArgumentParser()
    parser.add_argument('-e','--exe',type=str,default="rembo_test.x",
        help='Define the executable file to use here (rembo_test.x by default)')
    parser.add_argument('-o','--out', type=str, default=None,
        help='Specify the output directory (output/GRID_NAME by default)')
    parser.add_argument('-r','--run',action="store_true",
        help='Run the executable after preparing the job?')
    parser.add_argument('-s','--submit',action="store_true",
        help='Run the executable after preparing the job by submitting to the queue?')
    parser.add_argument('-f','--force',action="store_true",
        help='Execute run_rembo steps without confirmation?')
    parser.add_argument('-w','--wall', type=int, default=12,
        help='Maximum wall time to allow for job (only for jobs submitted to queue)')
    parser.add_argument('par_file', metavar='PAR_FILE', type=str,
        help='Name of the parameter file to be used in Yelmo simulation (with no folder).')
    args = parser.parse_args()

    ### Manage user options ############################
    
    # Copy the executable file to the output directory, or
    # call it from its compiled location?    
    copy_exec  = True 

    par_file   = args.par_file  
    force      = args.force 
    run        = args.run 
    submit     = args.submit 
    wtime      = args.wall 
    executable = args.exe 

    # Submit overrides run 
    if submit: run = True 

    if args.out is None:
        outpath = "output/test"
    else: 
        outpath = args.out 

    ### Start the script ############################

    # Make the job 
    if copy_exec:
        makejob(outpath,par_file,executable,force)
        executable = "./{}".format(executable)
    else:
        makejob(outpath,par_file,executable=None,force=force)
        cwd = os.getcwd()
        executable = "{}/librembo/bin/{}".format(cwd,executable)

    # Run the job if desired 
    if run:

        if submit:
            # Submit job to queue 
            pid = submitjob(outpath,executable,par_file,wtime=wtime) 

        else:
            # Run job in background 
            pid = runjob(outpath,executable,par_file)

    return 

def makejob(outpath,par_file,executable=None,force=False):
    # Define allowed grid names (based on grids.f90 predefined grids)
    
    # Determine the parameter file for constants
    par_fldr   = "par" 
    par_path   = "{}/{}".format(par_fldr,par_file)
    
    # Print info to screen
    print "\n  Output path = {}".format(outpath)
    print "     - parameters = {}".format(par_path) 

    response = ""   # Default response is nothing, check if another is given
    if not force:
        try:
            response = raw_input("\n[Enter to prepare directory]  or  [s to skip] or [ctl-c to exit] ")
            print "\n"
        except:
            print "\n"
            sys.exit()
    
    # Make output directory(-ies) and remove existing nc and nml files
    makedirs(outpath,remove=True)
    
    ## Generate symbolic links to input data folders
    srcname = "input"
    dstname = os.path.join(outpath,srcname)
    srcpath = os.path.abspath(srcname)
    if os.path.islink(dstname): os.unlink(dstname)
    os.symlink(srcpath,dstname)

    # # Generate link to extra data folder for personal data files
    # srcname = "extra_data"
    # dstname = os.path.join(outpath,srcname)
    # srcpath = os.path.abspath(srcname)
    # if os.path.islink(dstname): os.unlink(dstname)
    # os.symlink(srcpath,dstname)

    srcname = "ice_data"
    dstname = os.path.join(outpath,srcname)
    if os.path.islink(dstname): os.unlink(dstname)
    if os.path.islink(srcname):
        linkto = os.readlink(srcname)
        os.symlink(linkto, dstname)
    elif os.path.isdir(srcname):
        srcpath = os.path.abspath(srcname)
        os.symlink(srcpath,dstname)
    else:
        print "Error: path does not exist ",srcname
        sys.exit(2)

    ## Copy executable if given
    if not executable is None:
        shutil.copy("librembo/bin/{}".format(executable),outpath)
    
    ## Copy parameter file(s)
    shutil.copy(par_path,  outpath)

    return 

def runjob(outpath,executable,par_file):
    '''Run a job generated with makejob.'''

    cmd = "cd {} && exec {} {} > {} &".format(outpath,executable,par_file,"out.out")

    print "Running job in background:"
    print cmd 

    #os.system(cmd)
    proc = subp.Popen(cmd,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
    #pid  = proc.pid+1   # This is not necessarily accurate - do not use for anything
    pid = 0

    # Alternative below is supposed to be more safe,
    # and provide the proper pid of the process itself,
    # but doesn't appear to actually work...
    # cmd = ['cd',outpath,'&&','exec',executable,par_file,'>','out.out &']
    # print " ".join(cmd)
    # proc = subp.Popen(cmd,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
    # pid  = proc.pid
    #print "pid = {}".format(pid)

    return pid 

def submitjob(outpath,executable,par_file,wtime=12):
    '''Submit a job to a HPC queue (qsub,sbatch)'''

    # Get info about current system
    username  = os.environ.get('USER')
    usergroup = "anthroia"
    hostname  = os.environ.get('HOSTNAME')

    # Command to be called 
    cmd = "{} {}".format(executable,par_file) 

    # Create the jobscript using current info
    if "cei" in hostname:
        pass # TO DO 
        #script = jobscript_qsub(executable,outpath,username,usergroup,wtime)
    else:
        script = jobscript_slurm(cmd,outpath,username,usergroup,wtime)

    nm_jobscript   = 'job.submit'
    path_jobscript = "{}/{}".format(outpath,nm_jobscript)
    jobfile      = open(path_jobscript,'w').write(script)

    # Copy the job script into output directory for posterity
    if os.path.isfile (path_jobscript):
        print "Created jobscript file: " + nm_jobscript

    # Send the submit command to loadleveler or qsub
    if "cei" in hostname:
        # Change to output directory and submit job
        # stat = command("cd %s && qsub %s" %(outpath,nm_jobscript))
        pass # TO DO
        pid = 0  
    else:
        cmd_job = "cd {} && sbatch {}".format(outpath,nm_jobscript)
        
        #os.system(cmd_job)
        proc = subp.Popen(cmd_job,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
        #pid  = proc.pid+1   # This is not necessarily accurate - do not use for anything
        pid = 0

    return pid 

def makedirs(dirname,remove):
    '''
    Make a directory (including sub-directories),
    but first ensuring that path doesn't already exist
    or some other error prevents the creation.
    '''

    try:
        os.makedirs(dirname)
        print     'Directory created: ', dirname
    except OSError:
        if os.path.isdir(dirname):
            print 'Directory already exists: ', dirname
            if remove:
                for f in glob.glob("{}/*.nc".format(dirname)): 
                    os.remove(f)
                for f in glob.glob("{}/*.nml".format(dirname)):
                    os.remove(f)
                #print 'Removed *.nml and *.nc files.'
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise

    return

def autofolder(params,outfldr0):
    '''Given a list of parameters,
       generate an appropriate folder name.
    '''

    parts = []

    for p in params:
        parts.append( p.short() )

    # Join the parts together, combine with the base output dir
    autofldr = '.'.join(parts)
    outfldr  = outfldr0 + autofldr + '/'

    return outfldr

def jobscript_slurm(cmd,outpath,username,usergroup,wtime):
    '''Definition of the job script'''

    jobname = "rembo" 
    # jobname = "rembo_{}".format(outpath)

# Extra parameter options
##SBATCH --partition=ram_gpu
##SBATCH --mem=50000 

    script = """#! /bin/bash
#SBATCH --qos=short
#SBATCH --job-name={}
#SBATCH --account={}
#SBATCH --mail-user={}@pik-potsdam.de
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --output=./out.out
#SBATCH --error=./out.err
#SBATCH --time={}:00:00

# Run the job
{} 

""".format(jobname,usergroup,username,wtime,cmd)

    return script

if __name__ == "__main__": 

    # Call main run_rembo function...
    run_rembo()



