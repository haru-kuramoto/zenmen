import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob, os, time
from tqdm import tqdm
from scipy import interpolate

def make_mask(center, width, radius, shape):
    x = np.arange(shape[1])
    y = np.arange(shape[0])
    X, Y = np.meshgrid(x, y)
    dist = np.sqrt((X - center[1])**2 + (Y - center[0])**2)
    if width == -1:
        mask = dist <= radius
    else:
        mask = (radius-width/2 <= dist) * (dist <= radius+width/2)
    return mask


def calc_ea_mergedimage(Drctfiles, Reffiles, outdir, segment, energy):
    def fits2arr(fitsfile):
        hdulist = fits.open(fitsfile)
        data = hdulist[0].data
        arr = np.asarray(data).astype(np.float64)
        return arr
    shape = (2048, 2048)
    mergedimg = np.zeros(shape=shape)
    drct_intns = 0
    print("Calcurating direct files ...")
    for fitsfile in Drctfiles:
        arr = fits2arr(fitsfile)
        drct_intns += np.sum(arr)
    print("Calcurating reflected files ...")
    tqdm_intval = 10
    for fitsfile in tqdm(Reffiles, leave=False, mininterval=tqdm_intval):
        arr = fits2arr(fitsfile)
        mergedimg += arr
    outname = "Mergedimage/merged_seg%s_%skeV.fits"%(segment, energy)
    outpath = os.path.join(outdir, outname)
    mergedir = os.path.join(outdir, "Mergedimage")
    if os.path.isdir(mergedir) == False:
        os.makedirs(mergedir)
    new_hdu = fits.PrimaryHDU(mergedimg)
    new_hdu.writeto(outpath, overwrite=True)
    ref_intns = np.sum(mergedimg)
    ea = ref_intns/(drct_intns/len(Drctfiles))
    txtname = os.path.join(outdir,"ea_list.txt")
    outputdata = [segment, energy, str(ea)]
    if os.path.exists(txtname) == False:
        txtfile = open(txtname, "w")
        txtfile.write("# segment, energy, effective area\n")
        print(*outputdata, file=txtfile)
        txtfile.close()
    else:
        txtfile = open(txtname, "a")
        print(*outputdata, file=txtfile)
        txtfile.close()
    return mergedimg, ea


# PSFを何点で描くか
def calc_psf(mergeddata, points):
    print("Calcurating to draw PSF ... ")
    center = (1023.5, 1023.5)
    width = 2
    shape = (2048, 2048)
    PSF_list = []
    R_list = np.linspace(1, 1023.5, points)
    tqdm_intval = 10
    for rad in tqdm(R_list, leave=False, mininterval=tqdm_intval):
        mask = make_mask(center, width, rad, shape)
        masked_data = mergeddata[mask]
        intens = np.sum(masked_data) / len(masked_data)
        PSF_list.append(intens)
    return R_list, PSF_list

def calc_eef(mergeddata, points):
    print("Calcurating to draw EEF ... ")
    center = (1023.5, 1023.5)
    width = -1
    shape = (2048, 2048)
    EEF_list = []
    R_list = np.linspace(1, 1023.5, points)
    tqdm_intval = 10
    for rad in tqdm(R_list, leave=False, mininterval=tqdm_intval):
        mask = make_mask(center, width, rad, shape)
        masked_data = mergeddata[mask]
        intens = np.sum(masked_data)
        EEF_list.append(intens)
    return R_list, EEF_list

if __name__=="__main__":
    t1 = time.time()
    center = (1023.5,1023.5)
    focal_l = 12.761*1000*1000 #um
    pix2min = np.rad2deg(15.48/focal_l) * 60

    # parser = argparse.ArgumentParser()
    # parser.add_argument("-i", "--indir", help="Input dir", dest="indir", type=str, nargs=1)
    # parser.add_argument("-s", "--segment",help="Segment", dest="segment", type=str, nargs=1)
    # parser.add_argument("-e", "--energy",help="Energy", dest="energy", type=str, nargs=1)
    # parser.add_argument("-o", "--outdir", help="Outdir", dest="outdir", type=str, nargs=1)

    # args = parser.parse_args()
    # indir = args.indir[0]
    # segment = args.segment[0]
    # energy = args.energy[0]
    # outdir = args.outdir[0]

    # indir = "/Volumes/HDPH-UT/xl_calibur/work/zenmen/Apr23/CorrectedFits_202304/4CalcHPD"
    indir = "/Volumes/HDPH-UT/xl_calibur/work/zenmen/Jun21/CorrectedFits_202106/4CalcHPD"
    outdir = indir
    path_wc = os.path.join(indir, "**/*.fits")
    print(path_wc)
    Files = glob.glob(path_wc, recursive=True)
    if len(Files) == 0:
        print("Input directory name is wrong!")
        exit()
    Ene_list = ["20", "30", "40", "50", "70"]
    for ene in Ene_list:
        print("energy : %s"%ene)
        File_ene = [i for i in Files if ene in i.split('/')[-2]]
        File_ene.sort()
        if len(File_ene) == 0:
            continue
        R_list = []
        PSF_list_ene = []
        EEF_list_ene = []
        Pix_sum_ene = []
        seg_list = ["1", "2", "3"]
        for segment in seg_list:
            print("segment : %s"%segment)
            Drcts = [i for i in File_ene if i.split('/')[-1].split('_')[0] == "drct" and segment in i.split('/')[-2].split('_')[0]]
            if len(Drcts) == 0:
                continue
            Drcts.sort()
            Refs = [i for i in File_ene if i.split('/')[-1].split('_')[0] == "qtscan"and segment in i.split('/')[-2].split('_')[0]]
            Refs.sort()
            print("Number of direct files : %s"%len(Drcts))
            print("Number of reflected files : %s"%len(Refs))
            mergedimg, ea = calc_ea_mergedimage(Drcts, Refs, outdir, segment, ene)
            print(ea)
            # saving merged image
            imgdir = os.path.join(indir, "MergedFits")
            if os.path.isdir(imgdir) == False:
                os.makedirs(imgdir)
            mergedfile = os.path.join(imgdir, "merged_seg%s_%skeV.fits"%(segment, ene))
            merged_hdu = fits.PrimaryHDU(mergedimg)
            merged_hdu.writeto(mergedfile, overwrite=True)
            # saving normed image
            normedimg = mergedimg / np.amax(mergedimg)
            normedfile = os.path.join(imgdir, "normed_seg%s_%skeV.fits"%(segment, ene))
            normed_hdu = fits.PrimaryHDU(normedimg)
            normed_hdu.writeto(normedfile, overwrite=True)
            print("Seg%s %skeV finished!"%(segment, ene))
        