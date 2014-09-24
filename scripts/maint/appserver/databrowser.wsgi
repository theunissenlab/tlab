import os, sys
sys.path.append('/var/lib/databrowser')
sys.path.append('/var/lib')

from paste.deploy import loadapp

application = loadapp('config:/var/lib/databrowser/production.ini')
#application = loadapp('config:/var/lib/databrowser/development.ini')
