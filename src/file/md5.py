from hashlib import md5

__all__ = ['file_md5hash','check_md5hash']

def file_md5hash(filename):
    '''
    Simple way to get the md5 digest of a file.
    Produces output identical to Linux command-line "md5sum".
    
    Because this reads the entire file into the buffer, it may not be
    the fastest way to do this. One potentially faster method is implemented
    as a third-party Python library called 'fchksum', more info on which at
    http://freshmeat.net/projects/fchksum/
    '''
    
    fobj = open(filename,'rb')
    return md5(fobj.read()).hexdigest().upper()

def check_md5hash(filename,md5hash):
    if file_md5hash(filename).lower() == md5hash.lower():
        return True
    else:
        return False