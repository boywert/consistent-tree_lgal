import numpy as np
import os
import sys
from read_ctrees_cfg import read_ctrees_cfg

cfg = read_ctrees_cfg(sys.argv[1])

out_path = cfg['INBASE']
desc_path = '%s/DESC'%(cfg['INBASE'])
snap_path = '%s/RAW'%(cfg['INBASE'])
snap_list = filter(lambda x: '.AHF_halos' in x, os.listdir(snap_path))
snap_list.sort()
snap_list = snap_list[::-1]
snap_num = range(0, len(snap_list))[::-1]

#ID(0) hostHalo(1) numSubStruct(2) Mvir(3) npart(4) X(5) Y(6) Z(7) VX(8) VY(9)
#VZ(10) Rvir(11) Rmax(12) Rs(13) ... Vmax(16) Vesc(17) sigV(18)

#ID DescID Mass Vmax Vrms Radius Rs Np X Y Z VX VY VZ JX JY JZ Spin

for s, snap in zip(snap_num, snap_list):
    if(s == snap_num[0]):
        get_desc = lambda hid_this: -1
    else:
        desc_file = '%s/desc_%d.bin'%(desc_path, s)
        hid, desc = np.fromfile(desc_file, np.int64).reshape(2, -1)
        if(np.count_nonzero(desc > -1) < 120):
            s_final = s + 1
            break
        hs = np.argsort(hid)
        hid = hid[hs]
        desc = desc[hs]
        def get_desc(hid_this):
            k = np.searchsorted(hid, hid_this)
            if(hid[k] == hid_this): return desc[k]
            return -1
    #
    with open('%s/%s'%(snap_path, snap), 'r') as f, \
            open('%s/out_%d.list'%(out_path, s), 'w') as fo:
        for line in f:
            if(line[0] == '#'): continue
            item = line.split(None, 19)
            desc_this = get_desc(int(item[0]))
            pos = map(lambda x: float(x)*1.e-3, item[5:8])
            rs = float(item[13])
            if(rs > 1e38): rs = 0
            l='%s %d %s %s %s %s %.8g %s %.8g %.8g %.8g %s %s %s 0 0 0 0\n'%(\
                    item[0], desc_this, item[3], item[16], item[18], \
                    item[11], rs, item[4], pos[0], pos[1], pos[2], \
                    item[8], item[9], item[10])
            fo.write(l)
#
snap_data = np.loadtxt('%s/data_snaplist.txt'%(snap_path))
np.savetxt('%s/scales.txt'%(out_path), snap_data[s_final:][:,:2], "%d %f")

