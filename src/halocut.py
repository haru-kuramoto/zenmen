import numpy as np
from tqdm import tqdm
import os, glob

# 傾きを求める
def get_slope(datdir):
    def line2int(line):
        num = line.split(" ")[-1]
        return int(num)
    def get_offset(datdir):
        path = os.path.join(datdir,"qtscan_001.dat")
        f = glob.glob(path)[0]
        f_std = "/Users/kuramoto/ana/data/zenmen_stage_cl/SP8_23Apr/50keV/raster_scan_seg1_Ty+14_Tz-2_p10mm_50keV/qtscan_001.dat"
        data = open(f, "r")
        data_std = open(f_std, "r")
        datalist = data.readlines()
        datalist_std = data_std.readlines()
        diff_sy = line2int(datalist[4]) - line2int(datalist_std[4])
        diff_sz = line2int(datalist[5]) - line2int(datalist_std[5])
        return diff_sy, diff_sz
    path = os.path.join(datdir, "qtscan_*.dat")
    Files = glob.glob(path)
    Files.sort()
    # print(Files)
    print(len(Files))
    Sypls = 5000
    Szpls = 5400
    diff_sy, diff_sz = get_offset(datdir)
    zenmen_pos = []
    for file in Files:
        data = open(file, "r")
        datalist = data.readlines()
        pos_y = line2int(datalist[4]) - diff_sy
        pos_z = line2int(datalist[5]) - diff_sz
        zenmen_pos.append([pos_y/Sypls*(-1), pos_z/Szpls*(-1)])
    Slopes = [pos[1]*(-1)/(pos[0]*(-1)) for pos in zenmen_pos]
    return Slopes

def make_halocutmask(slope, center, shape):
    def linear(x, a, center):
        cen_y, cen_x = center 
        y = a*(x-cen_x) + cen_y
        return y
    mask = np.zeros(shape=shape)
    center = [(shape[0]-1)/2, (shape[1]-1)/2]
    deg = abs(np.arctan(slope))
    reg = abs(150 / (np.sin(np.pi /2 - deg)))
    x = np.arange(0, shape[0], 1)
    lin = linear(x, slope, center)
    for num in range(shape[0]):
        mask[:,num] = [1 if lin[num]-reg <= i <= lin[num]+reg else j for i,j in enumerate(mask[:,num])]
    return mask

def halocut_oneimg(subt, slope):
    shape = subt.shape
    center = [shape[0]/2, shape[1]/2]
    mask = make_halocutmask(slope, center, shape)
    masked_img = mask*subt
    return masked_img

def halocut_main(Subts, datdir, set_tqdm=True):
    Slopes = get_slope(datdir)
    def masking(subt, slope):
        maskedimg = halocut_oneimg(subt, slope)
        return maskedimg
    # Masked = [halocut(subt, slope) for subt,slope in zip(Subts, Slope)]
    Masked = []
    if set_tqdm==True:
        tqdm_intval = 10
        for subt, slope in zip(tqdm(Subts, leave=False,mininterval=tqdm_intval), Slopes):
            masked_img = masking(subt,slope)
            Masked.append(masked_img)
    return Masked
