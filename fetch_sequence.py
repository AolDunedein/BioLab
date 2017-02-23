from Bio import Entrez
from Bio import SeqIO
import os

def fetch_sequence(group):
    print("Fetching sequence from GeneBank for Group {}...".format(group))
    dic = dict()
    dic[1] = (1,270000)
    dic[2] = (270001,505535)
    dic[3] = (505536,753100)
    dic[4] = (753101,998650)
    dic[5] = (998651,1275600)
    dic[6] = (1275601,1533850)
    dic[7] = (1533851,1786150)
    dic[8] = (1786151,2072130)
    dic[9] = (2072131,2327100)
    dic[10] = (2327101,2610800)
    dic[11] = (2610801,2873800)
    dic[12] = (2873801,3124550)
    dic[13] = (3124551,3397754)
    seq_s = dic[group][0]
    seq_e = dic[group][1]
    pre = "sequence"
    ext = ".gb"
    seq_name = pre+ext
    print("start: {}\nstop: {}".format(seq_s,seq_e))

    if (os.path.exists(seq_name)==True):
        i = 0
        while True:
            seq_name=pre+str(i)+ext
            if (os.path.exists(seq_name)==False):
                break
            i = i+1
    Entrez.email = "a70430@alunos.uminho.pt"
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="NC_002942.5", seq_start=seq_s, seq_stop=seq_e)
    seq_record = SeqIO.read(handle, "gb")
    SeqIO.write(seq_record, seq_name, "genbank")
    handle.close()
    print("Done!")
    return seq_name

