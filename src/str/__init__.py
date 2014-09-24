def parse(input):
    if input.lower() == 'none':
        return None
    elif input.isdigit():
        return int(input)
    elif input.replace('.','').isdigit():
        return float(input)
    else:
        return input
    
def kv_pairs_string(in_dict):
    pairs = ["'%s',%s" % (k,quote_string(v)) for k,v in in_dict.iteritems()]
    return ','.join(pairs)

def quote_string(value):
    # Wrap strings in quotations
    if isinstance(value,(str,unicode)):
        return "'%s'" % value
    else:
        return value
    
def string_dict(value):
    if isinstance(value,dict):
        return 'struct(%s)' % kv_pairs_string(value)
    else:
        return value
    
def kv_dict_string(in_dict):
    pairs = ["'%s',%s" % (k,string_dict(quote_string(v))) for k,v in in_dict.iteritems()]
    return ','.join(pairs)