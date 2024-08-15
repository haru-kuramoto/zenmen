import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import os, glob, time
import queue, threading

def open_worker(Files, q_img):
    for file in Files:
        img = Image.open(file)
        ary = np.asarray(img).astype(np.int64)
        q_img.put([file, ary])
def open_multi_thread(indir, workers=6):
    q_img = queue.Queue()
    Imgs = []
    path = os.path.join(indir, "*.tif")
    filelist = glob.glob(path)
    filelist.sort()
    if len(filelist) == 0:
        print("FIle path will be wrong!")
        exit()
    #ファイル数をスレッド数でわった商
    imgs_per_thread = len(filelist) // workers
    #余り
    quotient = len(filelist) % workers
    #割り切れるだけのファイルをリスト内で分割
    Lists = [filelist[imgs_per_thread*i : imgs_per_thread*(i+1)] for i in range(workers)]
    #割り切れなかったファイルたちを
    quot_list = filelist[imgs_per_thread*workers : imgs_per_thread*workers+quotient]
    #追加
    Lists.append(quot_list)
    #スレッドを立てて処理する
    threads = []
    for filelist in Lists:
        #画像読み込み関数のスレッドを立てる
        t = threading.Thread(target=open_worker, args=(filelist, q_img))
        #デーモン化
        t.setDaemon
        #スレッドの開始
        t.start()
        #スレッドリストに格納
        threads.append(t)
    #あとはファイル名のリスト
    Files = []
    #スレッドの終了を待つ
    for thr in threads:
        thr.join()
    while not q_img.empty():
        temp = q_img.get()
        Files.append(temp[0])
        Imgs.append(temp[1])
    # 順番が揃った Files, Imgs が出てくる
    Files, Imgs = zip(*sorted(zip(Files, Imgs)))
    # print(filelist)
    return Imgs, Files

def filename(path):
    namelist = path.split("/")[-1].split("_")
    return namelist

def get_subtimage(indir):
    print("Opening file ...")
    Imgs, Files = open_multi_thread(indir)
    Drct = []
    Drct_files = []
    DDark = []
    DDark_files = []
    Ref1 = []
    Ref1_files = []
    Ref2 = []
    Ref2_files = []
    RDark = []
    RDark_files = []
    for i in range(len(Files)):
        head = filename(Files[i])[0]
        tail = filename(Files[i])[-1]
        if head == "chk":
            continue
        elif head == "drct":
            Drct.append(Imgs[i])
            Drct_files.append(Files[i])
        elif head == "dark":
            DDark.append(Imgs[i])
            DDark_files.append(Files[i])
        elif "qtscan" in filename[Files[i]] and tail == "001.tif":
            Ref1.append(Imgs[i])
            Ref1_files.append(Files[i])
        elif "qtscan" in filename[Files[i]] and tail == "002.tif":
            Ref2.append(Imgs[i])
            Ref2_files.append(Files[i])
        elif "qtscan" in filename[Files[i]] and tail == "dark.tif":
            Ref1.append(Imgs[i])
            Ref1_files.append(Files[i])
    print(len(Drct))
    Drct_files, Drct = zip(*sorted(zip(Drct_files, Drct)))
    DDark_files, DDark = zip(*sorted(zip(DDark_files, DDark)))
    Ref1_files, Ref1 = zip(*sorted(zip(Ref1_files, Ref1)))
    Ref2_files, Ref2 = zip(*sorted(zip(Ref2_files, Ref2)))
    RDark_files, RDark = zip(*sorted(zip(RDark_files, RDark)))
    # print(Drct_files)
    # print(DDark_files)
    # print(Ref1_files)
    # print(Ref2_files)
    # print(RDark_files)
    # 引き算
    print("Subtract dark ...")
    Drct_subt = [Drct[i]-DDark[i] for i in range(len(Drct))]
    Ref_subt = [(Ref1[i]+Ref2[i])/2-RDark[i] for i in range(len(Ref1))]
    Subts = np.array([Drct_subt, Ref_subt])
    return Subts

if __name__=="__main__":
    t1 = time.time()
    indir = "/Users/kuramoto/ana/data/raster_scan_seg2_Ty+56_Tz-16_p10mm_50keV"
    # Imgs, Files = open_multi_thread(indir)
    Subts = get_subtimage(indir)
    print(Subts.shape)
