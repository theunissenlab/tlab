#
import os, re

f = open('/etc/mysql/debian.cnf')
lines = f.readlines()

client_re = re.match('^\[client\]\s*$',lines[1])
host_re = re.match('^host\s*=\s*localhost\s*$',lines[2])
user_re = re.match('^user\s*=\s*debian-sys-maint\s*$',lines[3])
password_re = re.match('^password\s*=\s*(.*)\s*',lines[4])
socket_re = re.match('^socket\s*=\s*/var/run/mysqld/mysqld.sock\s*$',lines[5])

if all((client_re,host_re,user_re,password_re,socket_re)):
    password = password_re.groups()[0]
    old_password = input('Please enter old password for debian-sys-maint:\n')
    
    mysqladmin_cmd = "mysqladmin -udebian-sys-maint -p%s password %s"\
                        %(old_password,password)
    os.system(mysqladmin_cmd)