import fetch_sequence
import blasting
import dict_csv_gene

# FERRAMENTA CRIADA EM CONJUNÇÃO PELOS GRUPOS 13 e 9

# fetch_sequence
while True:
    inp = input("Insira o seu grupo: ")
    g=9
    try:
        g = int(inp)
        if g<0 or g>13:
            raise ValueError("Not Valid")
    except Exception:
        print("Nao valido, tente outra vez.")
        continue
    break
sequence_name = fetch_sequence.fetch_sequence(g)
#blasting with nr db
#nr_blast_result = blasting2.blast("nr",sequence_name)
swiss_blast_result = blasting.blast("swissprot",sequence_name)
print("Indexing values and creating table...")
# create dictionary for each gene
gene_dict = dict_csv_gene.create_gene_dict(sequence_name,None,swiss_blast_result)
# create csv from dict
csv_filename = dict_csv_gene.csv_export(gene_dict)
print ("Generated CSV file '"+csv_filename+"' successfully!")
report_filename = dict_csv_gene.create_report(gene_dict)
print ("Generated report file '"+report_filename+"' successfully!")


