import sys

from os import system
from os.path import isfile

def background(command,paths=[],work_directory=None,log_file=None):
    '''
    command is a single matlab command.
    paths is a list to add to the matlab path before running.
    if specified, command will execute in work_directory
    '''
    
    # Construct string to move to working directory if specified
    if work_directory is None:
        work_dir_string = ''
    else:
        work_dir_string = "cd('%s');" % work_directory
    
    # Construct strings to extend matlab path if specified
    addpath_string = " ".join("addpath('%s');" % path for path in paths)
    
    # Matlab command string
    matlab_script = " ".join((work_dir_string,addpath_string,command)).strip()
    matlab_cmd = "matlabbg \"%s\"" % matlab_script
    
    return system(matlab_cmd)

def read_commands_from_file(fileName):
    cmdList = open(fileName, 'r').readlines()
    return ';'.join(l.strip() for l in cmdList)

if __name__ == '__main__':

    if len(sys.argv) < 3:
        print ''
        print 'Usage: python matlab.py <paths> "<command>" <workdir>'
        print '\t<paths>: A colon-separated lists of matlab paths.'
        print '\t<command>: Command(s) to execute. If <command> is a file, then the list of commands is read from the file.'
        print '\t<workdir>: Working directory to run from.'
        exit()

    csvPaths = sys.argv[1]
    cmd = sys.argv[2]
    if isfile(cmd):
        print 'Reading commands from file %s' % cmd
        cmdString = read_commands_from_file(cmd)
    else:
        cmdString = cmd

    paths = csvPaths.split(':')

    wdir = None
    if len(sys.argv) == 4:
        wdir = sys.argv[3]

    background(cmdString, paths, wdir)
