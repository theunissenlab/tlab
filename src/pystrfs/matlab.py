import os
import string
import subprocess

class MatlabTemplate:

    def __init__(self):
        self.template = None

    def template_from_file(self, fname):
        if not os.path.isfile(fname):
            print 'Template file does not exist: %s' % fname
            return

        f = open(fname, 'r')
        tLines = f.readlines()
        f.close()

        self.template = ';'.join(l.strip() for l in tLines)


    def fill_template(self, params):

        scriptStr = self.template
        for name,val in params.iteritems():
            nStr = '$%s' % name
            scriptStr = string.replace(scriptStr, nStr, str(val))

        return scriptStr

    def to_matlab_process(self, params):

        cmdStr = self.fill_template(params)
        mp = MatlabProcess(cmdStr)
        return mp


class MatlabProcess:

    def __init__(self, commandString):
        self.cmd_string = commandString

    def get_commands(self):
        cstr = '\"%s\"' % self.cmd_string
        finalCmds = ['matlabbg', cstr]
        return finalCmds

    def run(self):
        finalCmds = self.get_commands()
        fCmdStr = ' '.join(finalCmds)
        print 'Running matlab: %s' % fCmdStr
        p = subprocess.Popen(fCmdStr, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
