from Bio import SeqIO
import anal_prot
from io import StringIO
from Bio import SwissProt
from Bio import ExPASy
from Bio.Blast import NCBIXML
import requests
import re

def create_gene_dict(sequence_name,nr_blast_res,swiss_blast_res):
    dic_res = dict()
    dic_res["default"] = []
    dic_res["default"].append("Locus Tag")
    seq_records = SeqIO.read(sequence_name,"gb")
    for feat in seq_records.features:
        if not (feat.type=="CDS" or feat.type=='gene'):
            continue
        try:
            locus_tag = feat.qualifiers["locus_tag"][0]
        except Exception:
            locus_tag = "-"
        try:
            geneID = feat.qualifiers["db_xref"][0].split(":")[1]
        except Exception:
            geneID = "-"
        try:
            gene_name = feat.qualifiers["gene"][0]
        except Exception:
            gene_name = "-"
        try:
            function = feat.qualifiers["function"][0]
        except Exception:
            function = "-"
        try:
            protein_id = feat.qualifiers["protein_id"][0]
        except Exception:
            protein_id = "-"
        try:
            translation = feat.qualifiers["translation"][0]
        except Exception:
            translation = "-"
        try:
            strand = "Unknown"
            if feat.strand==1:
                strand= "Plus Strand"
            if feat.strand==-1:
                strand= "Minus Strand"
        except Exception:
            strand = "-"
        name, id, locale, status, fmol, bio, function2, leng, ec = "-","-","-","-","-","-","-","-","-"
        try:
            params = {"query": protein_id, "format": "fasta"}
            frecord = requests.get("http://www.uniprot.org/uniprot/", params)
            idfasta = ""
            for record in SeqIO.parse(StringIO(frecord.text), "fasta"):
                aux = record.id
                idfasta = str(aux).split("|")[1]
                name, id, locale, status, fmol, bio, function2, leng, ec = anal_prot.getDataFromProt(idfasta)
        except Exception:
            print("Try going online, no internet connection found.")
        if function2=="-":
            function2 = "Gene coding tRNA"
        local_dict = dict()
        def spec_append(strng,dic):
            if strng not in dic["default"]:
                dic["default"].append(strng)
        local_dict["Status"] = status
        spec_append("Status",dic_res)
        local_dict["Protein Name"] = name
        spec_append("Protein Name",dic_res)
        local_dict["Protein Uniprot ID"] = id
        spec_append("Protein Uniprot ID",dic_res)
        local_dict["Strand"] = strand
        spec_append("Strand",dic_res)
        local_dict["Celular Location"] = locale
        spec_append("Celular Location",dic_res)
        local_dict["GO - Biological Process"] = bio
        spec_append("GO - Biological Process",dic_res)
        local_dict["GO - Molecular function"] = fmol
        spec_append("GO - Molecular function",dic_res)
        local_dict["Gene ID NCBI"] = geneID
        spec_append("Gene ID NCBI",dic_res)
        local_dict["Gene Name"] = gene_name
        spec_append("Gene Name",dic_res)
        #if function!="-":
        local_dict["Function NCBI"] = function
        spec_append("Function NCBI",dic_res)
        #else:
        local_dict["Function UniProt"] = function2
        spec_append("Function UniProt",dic_res)
        spec_append("SwissProt BLAST Function",dic_res)
        local_dict["Protein Accession NCBI"] = protein_id
        spec_append("Protein Accession NCBI",dic_res)
        local_dict["EC Value"] = ec
        spec_append("EC Value",dic_res)
        local_dict["Nr Aminoacidos"] = leng
        spec_append("Nr Aminoacidos",dic_res)
        local_dict["Translation"] = translation
        spec_append("Translation",dic_res)
        E_VALUE_THRESH = 0.05
        try:
            blast_res = swiss_blast_res[locus_tag]
            handle = open("temp.xml","w")
            handle.write(blast_res)
            handle.close()
            handle = open("temp.xml","r")
            blast_record = NCBIXML.read(handle)
            alignment = blast_record.alignments[0]
            hsp = alignment.hsps[0]
            if hsp.expect < E_VALUE_THRESH:
                    handle = ExPASy.get_sprot_raw(alignment.accession)
                    record = SwissProt.read(handle)
                    function3 = ""
                    for cr in record.comments:
                        cr = str(cr).split(":")
                        if cr[0] == "FUNCTION":
                            function3 = cr[1].split("{")[0]
                    if function3 == "":
                        raise ValueError("Not Found")
                    local_dict["SwissProt BLAST Function"] = function3
            else:
                raise ValueError("Not found")
        except Exception:
            local_dict["SwissProt BLAST Function"] = "Not Found"
        dic_res[locus_tag] = local_dict
    return dic_res


# Create CSV Table
def csv_export(genedict):
    sep = ";"
    head = ""
    filecsv = open("table.csv","w")
    for entry in genedict["default"]:
        head += entry+sep
    filecsv.write(head+"\n")
    for entry in genedict:
        dataCSV = ""
        if entry=="default":
            continue
        dataCSV += entry+sep
        local_dic = genedict[entry]
        for key in genedict["default"]:
            if key=="Locus Tag":
                continue
            dataCSV += sanit(str(local_dic[key]))+sep
        filecsv.write(dataCSV + "\n")
    filecsv.close()
    return "table.csv"

def sanit(line):
    res = re.sub('[;]','',line)
    return res

# Create a report with all the info
def create_report(genedict):
    file = open("report.txt", "w")
    file.write("GENERAL GENE REPORT\n")
    for entry in genedict:
        local_dic = genedict[entry]
        if entry == "default":
            continue
        file.write("\n")
        for key in genedict["default"]:
            if key == "Locus Tag":
                continue
            file.write(key+":"+str(local_dic[key]))
    file.close()