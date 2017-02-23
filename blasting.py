from Bio.Blast import NCBIWWW
from Bio import SeqIO
import threading
import os
import time

def query_WWW(db,seq):
    try:
        print("Query",pos,"started.")
        result_handle = NCBIWWW.qblast("blastp", dbase, seq)
        res =  result_handle.read()
        result_handle.close()
        return res
    except Exception:
        return ""

def blast(db="swissprot",seq="sequence.gb"):
    global xml_fd
    global lock
    global tot
    global dic
    xml_fd = dict()
    lock = threading.Lock()
    dic = dict()
    tot = 0
    run_blast(seq,db)
    return xml_fd

def query_WWW(pos):
    global xml_fd
    global lock
    global dic
    global dbase
    global current_num
    try:
        print("Query",pos,"started.")
        result_handle = NCBIWWW.qblast("blastp", dbase, dic[pos].qualifiers["translation"][0], hitlist_size=2,alignments=1)
        lock.acquire()
        xml_fd[dic[pos].qualifiers["locus_tag"][0]]=(result_handle.read())
        print ("Result for Query saved!")
        result_handle.close()
        lock.release()
    except Exception:
        print("Query timed out...")


def run_blast(sequence_path,db):
    global dic
    global tot
    global dbase
    dbase = db
    print("Starting blast on {} database...".format(dbase))
    seq_record = SeqIO.read(sequence_path,"gb")
    threads = []
    for feat in seq_record.features:
        if feat.type=="CDS":
            dic[tot] = feat
            thread = threading.Thread(target = query_WWW, args = (tot,))
            thread.start()
            threads.append(thread)
            time.sleep(5)
            tot = tot+1
    while True:
        flag = True
        for t in threads:
            if (t.isAlive()==True):
                flag = False
            else:
                threads.remove(t)
        if(flag==True):
            break
        print("{} of {} queries left...".format(str(len(threads)),str(tot)))
        time.sleep(20)
    print("Blast done!")

