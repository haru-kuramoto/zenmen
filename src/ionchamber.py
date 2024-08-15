import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.signal import savgol_filter
import datetime
import os, glob

def sortedfiles(Subts, datdir):
    def namelist(file):
        name = file.split('/')[-1].split('_')
        return name

    # def deldir(Files):
    #     no_dir_files = [i.split('/')[-1] for i in Files]
    #     return no_dir_files

    def TimingList(Datfiles):
        Timing = []
        for file in Datfiles:
            data = open(file, 'r')
            datalist = data.readlines()
            date_str = datalist[0]
            date_time = datetime.datetime.strptime(date_str, "Time:  %Y / %m / %d   %H : %M : %S\n")
            Timing.append(date_time.timestamp())
        return Timing
    # Dat 整理
    print("Processing dat files ...")
    path = os.path.join(datdir, "*.dat")
    DatFiles = glob.glob(path)
    DatFiles.sort()
    StageDrct = []
    StageDDark = []
    StageRef = []
    StageTiming = []
    # Datfile の分割
    for i in DatFiles:
        if namelist(i)[0] == "chk":
            continue
        elif namelist(i)[0] == "drct":
            StageDrct.append(i)
        elif namelist(i)[0] == "dark":
            StageDDark.append(i)
        elif namelist(i)[0] == "qtscan":
            StageRef.append(i)
    # Subts と共に並び替え
    Drct_subt, Ref_subt = Subts
    SDrct_tim, Drct_sort = zip(*sorted(zip(TimingList(StageDrct), Drct_subt)))
    SRef_tim, Ref_sort = zip(*sorted(zip(TimingList(StageRef), Ref_subt)))
    StageTiming = np.array([SDrct_tim, SRef_tim])
    Subts_sort = np.array([Drct_sort, Ref_sort])
    return StageTiming, Subts_sort



# まずは与えられたファイルからこのような曲線を得る
def get_ic_corrlist(StageTiming, segment, energy, outdir):
    StageDrct, StageRef = StageTiming
    txtfile = None
    if segment == "1" and energy == "50":
        txtfile = "/Users/kuramoto/ana/data/SP8_23Apr_IC/ic_20230419_1135.txt"
    else:
        txtfile = "/Users/kuramoto/ana/data/SP8_23Apr_IC/ic_20230420_1145.txt"
    timing = np.loadtxt(txtfile, usecols=[0])
    ic_intns = np.loadtxt(txtfile, usecols=[2])
    StageAll = np.concatenate([StageDrct, StageRef])
    Cut_tim = []
    Cut_ic = []
    # 前後何秒とるか
    extra = 30
    # datデータの最初と最後からスキャン間のIC値を取得
    init = np.amin(StageAll)
    end = np.amax(StageAll)
    for i in range(len(timing)):
        if init - extra < timing[i] < end + extra:
            Cut_tim.append(timing[i])
            Cut_ic.append(ic_intns[i])
    print("要素数のチェック")
    print(len(Cut_tim) == len(Cut_ic))
    thresh = None
    seg = int(segment)
    ene = int(energy)
    if seg==1 and ene==20:
        thresh = 105700
    elif seg==2 and ene==20:
        thresh = 108600
    elif seg==3 and ene==20:
        thresh = 110000
    elif seg==1 and ene==30:
        thresh = 74800
    elif seg==2 and ene==30:
        thresh = 75350
    elif seg==3 and ene==30:
        thresh = 74850        
    elif seg==1 and ene==40:
        thresh = 54950
    elif seg==2 and ene==40:
        thresh = 54950
    elif seg==3 and ene==40:
        thresh = 54930
    elif seg==1 and ene==50:
        thresh = 48100
    elif seg==2 and ene==50:
        thresh = 47860
    elif seg==3 and ene==50:
        thresh = 47860
    elif seg==1 and ene==70:
        thresh = 41500
    elif seg==2 and ene==70:
        thresh = 43200
    elif seg==3 and ene==70:
        thresh = 42750
    print(thresh)
    # print(type(thresh))
    Beam_idc = np.where(np.array(Cut_ic) > thresh)
    Timing_beam = [Cut_tim[i] for i in Beam_idc[0]]
    IC_beam = [Cut_ic[i] for i in Beam_idc[0]]
    # X線検出時刻に最も近いIC測定時刻を持ってくる
    def closest_idx(val):
        idx = np.argmin(np.abs(np.array(Timing_beam) - val))
        return idx
    # ダイレクト光
    Closest_drct = [closest_idx(i) for i in StageDrct]
    Closest_drct_tim = [Timing_beam[i] for i in Closest_drct]
    Closest_drct_ic = [IC_beam[i] for i in Closest_drct]
    # 反射光
    Closest_ref = [closest_idx(i) for i in StageRef]
    Closest_ref_tim = [Timing_beam[i] for i in Closest_ref]
    Closest_ref_ic = [IC_beam[i] for i in Closest_ref]
    Drct_ic = []
    for i in range(len(Closest_drct_ic)):
        tim2idx = Closest_drct[i]
        if i == 0 and (Closest_drct_tim[i]-StageDrct[i]) < 0:
            drct = (IC_beam[tim2idx] + IC_beam[tim2idx+1])/2
            Drct_ic.append(drct)
        elif i == 0 and (Closest_drct_tim[i]-StageDrct[i]) > 0:
            drct = IC_beam[tim2idx]
            Drct_ic.append(drct)
        elif i == len(Closest_drct_ic)-1 and (Closest_drct_tim[i]-StageDrct[i]) < 0:
            drct = IC_beam[tim2idx]
            Drct_ic.append(drct)
        elif i == len(Closest_drct_ic)-1 and (Closest_drct_tim[i]-StageDrct[i]) > 0:
            drct = (IC_beam[tim2idx]+IC_beam[tim2idx-1])/2
            Drct_ic.append(drct)
        elif (Closest_drct_tim[i]-StageDrct[i]) > 0:
            drct = (IC_beam[tim2idx] + IC_beam[tim2idx-1])/2
            Drct_ic.append(drct)
        elif (Closest_drct_tim[i]-StageDrct[i]) < 0:
            if IC_beam[tim2idx] == IC_beam[-1]:
                drct = IC_beam[tim2idx]
            else:
                drct = (IC_beam[tim2idx] + IC_beam[tim2idx+1])/2
            Drct_ic.append(drct)
    Ref_ic = []
    check_sort = []
    for i in range(len(Closest_ref_ic)):
        tim2idx = Closest_ref[i]
        if i == 0 and (Closest_ref_tim[i]-StageRef[i]) > 0:
            iccorr = (IC_beam[tim2idx] + IC_beam[tim2idx+1])/2
            Ref_ic.append(iccorr)
            check_sort.append(i)
        elif i == 0 and (Closest_ref_tim[i]-StageRef[i]) < 0:
            iccorr = IC_beam[tim2idx]
            Ref_ic.append(iccorr)
            check_sort.append(i)
        elif i == len(Closest_ref_ic)-1 and (Closest_ref_tim[i]-StageRef[i]) > 0:
            iccorr = (IC_beam[tim2idx-1] + IC_beam[tim2idx])/2
            Ref_ic.append(iccorr)
            check_sort.append(i)
        elif i == len(Closest_ref_ic)-1 and (Closest_ref_tim[i]-StageRef[i]) < 0:
            iccorr = IC_beam[tim2idx]
            Ref_ic.append(iccorr)
            check_sort.append(i)
        elif (Closest_ref_tim[i]-StageRef[i]) > 0:
            iccorr = (IC_beam[tim2idx]+IC_beam[tim2idx-1])/2
            Ref_ic.append(iccorr)
            check_sort.append(i)
        elif (Closest_ref_tim[i]-StageRef[i]) < 0:
            if IC_beam[tim2idx] == IC_beam[-1]:
                iccorr = IC_beam[tim2idx]
            else:
                iccorr = (IC_beam[tim2idx]+IC_beam[tim2idx+1])/2
            Ref_ic.append(iccorr)
            check_sort.append(i) 
    Ref_corrlist = np.array(Ref_ic) / Drct_ic[0] 
    Drct_corrlist = np.array(Drct_ic) / Drct_ic[0] 
    print("-----Check sort when calc reflected ic counts-----")
    print(check_sort)
    print("--------------------------------------------------")
    plt.rcParams["font.size"] = 15
    plt.figure(figsize=(8,6), dpi=300)
    plt.scatter(Timing_beam, IC_beam, c="black", s=5, label="Beam intensity")
    plt.scatter(StageRef, Ref_ic, c="red", s=3, label="IC counts (Reflected)")
    plt.scatter(StageDrct, Drct_ic, c="blue", s=3, label="IC counts (Direct)")
    plt.xlabel("UNIX time")
    plt.ylabel("Counts of ion chamber")
    plt.title("Rasterscan_seg%s_%skeV"%(segment, energy))
    outfile = "beam_fluct_seg%s_%skeV.png"%(segment, energy)
    outpath = os.path.join(outdir,outfile)
    if os.path.isdir(outdir) == False:
        os.makedirs(outdir)
    plt.grid()
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    print("Standard intensity (drct_001_001) : %s"%Drct_ic[0])
    return Drct_corrlist, Ref_corrlist

def main_ic_correction(Subts, datdir, segment, energy, outdir):
    StageTiming, Subts_sort = sortedfiles(Subts, datdir)
    Drct_corrlist, Ref_corrlist = get_ic_corrlist(StageTiming, segment, energy, outdir)
    print(Drct_corrlist)
    # print(len(Subts_sort[0]))
    # print(Ref_corrlist)
    # print(len(Subts_sort[1]))
    CorredDrcts = [Subts_sort[0][i]/Drct_corrlist[i] for i in range(len(Subts_sort[0]))]
    CorredRefs = [Subts_sort[1][i]/Ref_corrlist[i] for i in range(len(Subts_sort[1]))]
    return CorredDrcts, CorredRefs