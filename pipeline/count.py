import Bio.SeqIO as seqio
keys_list = ['ACRNAN','ANCCAN', 'Acey', 'ACAC', 'Ag', 'Bm', 'BXY', 'nDi', 'Dd', 'EVEC', 'GPLIN', 'GROS', 'HCON', 'Hetgly', 'Hba', 'EN70', 'tig', 'MhA', 'Minc', 'NECAME', 'OVOC', 'PDA', 'PTRK', 'PSAM', 'PPA', 'nRc', 'SC.X', 'L596','SRAE', 'SSTP', 'T4D', 'T01', 'TMUE', 'TTRE']

fasta_path = r'./output/PCs/output_nematodes_1640171400114.fa'
cluster_entry = [f.description for f in seqio.parse(fasta_path, "fasta")]

for key in keys_list:
    count=0
    for value in cluster_entry:
        if key in value:
            count+=1
            print(value)
    print("{}={}".format(key, count))
