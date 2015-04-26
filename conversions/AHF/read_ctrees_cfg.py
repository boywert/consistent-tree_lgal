def read_ctrees_cfg(fname):
    d = {}
    with open(fname, 'r') as f:
        for l in f:
            fields = l.partition('#')[0].partition('=')
            if(fields[1] == '='):
                d[fields[0].strip()] = fields[2].strip().strip('\'"')
    return d
