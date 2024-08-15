import numpy as np
import astropy.io.fits as fits
import argparse
import glob, os, time

import halocut as hc
import imagesubtract as ims
import ionchamber as ic
import opticalvignetcorrect as oc

def calc_ea(Drct, Subt):
    print("Calculating Effect Area ...")
    drct_intns = np.sum(Drct)/len(Drct)
    ref_intns = np.sum(Subt)
    ea = ref_intns / drct_intns
    return ea

def save_as_fits(Drcts, Refs, outdir, segment, energy, expdate, corrtype):
    print("Saving ...")
    output_dir = None
    if corrtype == "OV":
        output_dir = os.path.join(outdir, "CorrectedFits_%s/4CalcRef/seg%s_%skeV/"%(expdate, segment, energy))
    elif corrtype == "HC":
        output_dir = os.path.join(outdir, "CorrectedFits_%s/4CalcHPD/seg%s_%skeV/"%(expdate, segment, energy))  
    if os.path.isdir(output_dir) == False:
        os.makedirs(output_dir)
    for i,dr in enumerate(Drcts):
        num = str(i)
        output_file_name = "drct_%s.fits"%num
        path = os.path.join(output_dir, output_file_name)
        new_hdu = fits.PrimaryHDU(dr)
        new_hdu.writeto(path, overwrite=True)
    for i,ref in enumerate(Refs):
        num = str(i+1)
        output_file_name = "qtscan_%s.fits"%num
        path = os.path.join(output_dir, output_file_name)
        new_hdu = fits.PrimaryHDU(ref)
        new_hdu.writeto(path, overwrite=True)

def calc_all_correction(indir, datdir, outdir ,segment, energy):
    print(indir)
    # ダーク減算
    Subts = ims.get_subtimage(indir)
    ea_nocorr = calc_ea(Subts[0], Subts[1])
    print("EA (nocorr) : %s"%ea_nocorr)
    # IC 補正
    print("Start IC correction ...")
    ic_outdir = os.path.join(outdir, "IC_figure")
    ICCorredDrcts, ICCorredRefs = ic.main_ic_correction(Subts, datdir, segment, energy, ic_outdir)
    ea_ic = calc_ea(ICCorredDrcts, ICCorredRefs)
    print("EA ICcorred : %s"%ea_ic)
    print("End IC correction")
    # Optical vignet correction
    print("Start Optical Vignet correction ...")
    OVCorredRefs = oc.vignet_corr(ICCorredRefs)
    OVCorredDrcts = oc.vignet_corr(ICCorredDrcts)
    txtname = os.path.join(outdir,"ea_list.txt")
    print("End of ea calc")
    ea = calc_ea(OVCorredDrcts, OVCorredRefs)
    outputdata = [segment, energy, str(ea)]
    # save_as_fits(Drcts, Refs, outdir, segment, energy, expdate, corrtype)
    save_as_fits(OVCorredDrcts, OVCorredRefs, outdir, segment, energy, expdate=expdate, corrtype="OV")
    print("End IC correction")
    if os.path.exists(txtname) == False:
        txtfile = open(txtname, "w")
        txtfile.write("# segment, energy, effective area\n")
        print(*outputdata, file=txtfile)
        txtfile.close()
    else:
        txtfile = open(txtname, "a")
        print(*outputdata, file=txtfile)
        txtfile.close()
    print("EA : %s"%ea)
    # Halo cut 
    print("Start Halo cut correction ...")
    HCCorredRefs = hc.halocut_main(OVCorredRefs, datdir, set_tqdm=True)
    print("End Halo Cut correction")
    return OVCorredDrcts, HCCorredRefs



def main(indir, datdir, outdir ,segment, energy, expdate):
    OVCorredDrcts, HCCorredRefs = calc_all_correction(indir, datdir, outdir ,segment, energy)
    save_as_fits(OVCorredDrcts, HCCorredRefs, outdir, segment, energy, expdate=expdate, corrtype="HC")
    ea = calc_ea(OVCorredDrcts, HCCorredRefs)
    print(ea)
"""
python ./AllCorrection.py -i /Users/kuramoto/ana/data/raster_scan_seg2_Ty+56_Tz-16_p10mm_50keV -d /Users/kuramoto/ana/data/zenmen_stage_cl/SP8_23Apr/50keV/raster_scan_seg2_Ty+58_Tz-18_p10mm_50keV -o /Users/kuramoto/ana/work/zenmen_202306 -s 2 -e 50 -ed 202306
"""

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", help="Input dir", dest="indir", type=str, nargs=1)
    parser.add_argument("-d", "--datfile", help="Datfile", dest="datdir", type=str, nargs=1)
    parser.add_argument("-s", "--segment",help="Segment", dest="segment", type=str, nargs=1)
    parser.add_argument("-e", "--energy",help="Energy", dest="energy", type=str, nargs=1)
    parser.add_argument("-o", "--outdir", help="Outdir", dest="outdir", type=str, nargs=1)
    parser.add_argument("-ed", "--expdate", help="Expdate", dest="expdate", type=str, nargs=1)

    # コマンドライン引数を変数に変換
    args = parser.parse_args()
    indir = args.indir[0]
    segment = args.segment[0]
    energy = args.energy[0]
    outdir = args.outdir[0]
    expdate = args.expdate[0]
    datdir = args.datdir[0]

    t1 = time.time()
    main(indir, datdir, outdir ,segment, energy, expdate)
    print(indir)
    print(energy, segment)
    t2 = time.time()
    dt = t2 - t1
    print("Time : %s"%dt)