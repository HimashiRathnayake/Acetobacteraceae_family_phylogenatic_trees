import pandas as pd
import re
from Bio import SeqIO
from prettytable import PrettyTable
from ete3 import Tree

def file_read():
    file_name = "./Data/protein_tables.xlsx"

    xl_file = pd.ExcelFile(file_name)

    dfs = {sheet_name: xl_file.parse(sheet_name) 
            for sheet_name in xl_file.sheet_names}
    return dfs


def get_bacteria_dic():
    bacteria_dict = {
        'NC_017108.1': "Acetobacter pasteurianus", 
        'NC_017150.1': "Acetobacter ghanensis", 
        'NZ_CP022374.1': "Acetobacter oryzifermentans", 
        'NZ_LN609302.1': "Acetobacter ghanensis", 
        'NZ_CP014692.1': "Acetobacter aceti", 
        'NZ_CP015164.1': "Acetobacter ascendens", 
        'NZ_CP011120.1': "Acetobacter oryzifermentans", 
        'NC_017121.1': "Acetobacter pasteurianus", 
        'NZ_CP022699.1': "Acetobacter tropicalis", 
        'NZ_AP018515.1': "Acetobacter orientalis", 
        'NZ_CP014687.1': "Acetobacter persici", 
        'NZ_CP023189.1': "Acetobacter pomorum", 
        'NZ_CP023657.1': "Acetobacter pomorum", 
        'NC_017125.1': "Acetobacter pasteurianus", 
        'NC_017146.1': "Acetobacter pasteurianus", 
        'NZ_AP014881.1': "Acetobacter pasteurianus", 
        'NZ_CP021524.1': "Acetobacter ascendens", 
        'NC_017100.1': "Acetobacter pasteurianus", 
        'NC_017111.1': "Acetobacter pasteurianus", 
        'NZ_CP015168.1': "Acetobacter ascendens", 
        'NZ_LN606600.1': "Acetobacter senegalensis"
    }
    return bacteria_dict

def get_cosmid_proteins_list():
    cosmid_proteins = ["ABC transporter permease", 
                    "LysR family transcriptional regulator",
                    "helix-turn-helix domain-containing protein",
                    "efflux transporter outer membrane subunit"]
    return cosmid_proteins

def get_comman_bacteria_list(cosmid_proteins, dfs):

    accession_ids = list(dfs.keys())

    unique_protein = {}

    for key in accession_ids:
        unique_proteins = list(dfs[key]["Protein name"].unique())
        unique_protein[key] = unique_proteins

    pro = {}

    for protein in cosmid_proteins:
        keys = []
        r2 = re.compile(".*" + protein)
        for key in unique_protein:
            newlist = list(filter(r2.match, unique_protein[key]))
            if (len(newlist))>0:
                keys.append(key)
        pro[protein] = keys

    common_bacteria_set = set(pro[cosmid_proteins[0]]).intersection(set(pro[cosmid_proteins[1]])).intersection(set(pro[cosmid_proteins[2]]))
    return common_bacteria_set

def write_fasta(seqs, fasta_file):
  """Write sequences to a fasta file.

  Parameters
  ----------
  seqs : dict[seq_id] -> seq
      Sequences indexed by sequence id.
  fasta_file : str
      Path to write the sequences to.
  wrap: int
      Number of AA/NT before the line is wrapped.
  """
  with open(fasta_file, 'w') as f:
      for gid, gseq in seqs.items():
          f.write('>{}\n'.format(gid))
          for i in range(0, len(gseq)):
              f.write('{}'.format(gseq[i]))
          f.write('{}'.format("\n")) 

def write_homologous_gene_sequence(cosmid_proteins, common_bacteria_set, dfs):
    homologous_gene_sequence = {}

    data_folder_path = "./Data/ExtractedSequences/"
    table_print = PrettyTable(['Accession ID', 'Protein name', 'Start','End'])
    for p in cosmid_proteins:
        sequence_list = {}
        seq = ''
        fasta_seq = {}
        for s in common_bacteria_set:
            for index, key in enumerate(dfs[s]["Protein name"]):
                if (re.match(".*" + p, key)):
                    start = dfs[s].loc[index]["Start"]
                    stop = dfs[s].loc[index]["Stop"]
                    table_print.add_row([s, p, start, stop])
                    seq_file = "./Data/Sequences/" + s.split(".")[0] + "_1_sequence.fasta"
                    sequence = SeqIO.to_dict(SeqIO.parse(open(seq_file), 'fasta'))[s.split(" ")[0]].seq[start:stop]
                    sequence_list[s] = str(sequence)
                    seq += ">"+s+"\n"
                    seq += str(sequence)+"\n\n"
                    fasta_seq[s] = str(sequence)
                    break
        write_fasta(fasta_seq,data_folder_path+p+'.fasta')
    print(table_print)

def get_robinson_foulds_distance():
    tree_data = ['./Data/Tree Data UPGMA/ABC transporter permease.txt', 
                './Data/Tree Data UPGMA/LysR family transcriptional regulator.txt', 
                './Data/Tree Data UPGMA/helix-turn-helix domain-containing protein.txt', 
                './Data/Tree Data UPGMA/efflux transporter outer membrane subunit.txt']
    return_list = []
    for idx, a in enumerate(tree_data):
        for b in tree_data[idx + 1:]:
            t1 = Tree(a)
            t2 = Tree(b)
            
            rf, max_rf, common_leaves, parts_t1, parts_t2, x4, x5 = t1.robinson_foulds(t2, unrooted_trees=True)
            return_list.append([a.split('/')[-1].split(".")[0], b.split('/')[-1].split(".")[0], rf, max_rf])

    return return_list

if __name__ == "__main__":
    dfs = file_read()
    cosmid_proteins = get_cosmid_proteins_list()

    # Step 1
    print("################ Step 1 ################")
    common_bacteria_set = get_comman_bacteria_list(cosmid_proteins=cosmid_proteins, dfs=dfs)
    print("Common Bacteria set")
    for each in common_bacteria_set:
        print(each)
    print("\n\n")

    # Step 2
    ## downloaded genoms sequence

    # Step 3
    print("################ Step 3 ################")
    write_homologous_gene_sequence(cosmid_proteins=cosmid_proteins, common_bacteria_set=common_bacteria_set, dfs=dfs)

    print("\n\n")

    # Step 4
    ## Phylogenatic trees were generated using the clustal phylogeny and clustal omega softwares

    # Step 5
    print("################ Step 5 ################")
    trees_distances = get_robinson_foulds_distance()
    table_print = PrettyTable(['Tree name - 1', 'Tree name - 2', 'robinson_foulds distance','max robinson_foulds distance'])
    for a,b,rf,max_rf in trees_distances:
        table_print.add_row([a,b,rf,max_rf])
    
    print(table_print)


     
 