import sys
sys.path.append('/var/lib')
from tlab.src import *

def application(environ, start_response):
    status = '200 OK'
    output = '\n'.join(p for p in sys.path)

    response_headers = [('Content-type', 'text/plain'),
                        ('Content-Length', str(len(output)))]
    start_response(status, response_headers)

    return [output]