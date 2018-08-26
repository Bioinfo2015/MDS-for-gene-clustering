# Ruth_Hashkes
import csv
import numpy as np
from Bio import SeqIO
from sklearn import manifold
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# Question 1:
seqs = []
for record in SeqIO.parse("seqs_to_cluster.fa", "fasta"):
    seqs.append([record.id, record.seq._data])


# Question 2:
def similarity(seq1, seq2):
    """find similarity percentage based on percent of mismatches"""
    counter = 0
    for a in range(len(seq1)):
        if seq1[a] == seq2[a]:
            counter = counter + 1
    similarity_score = (float(counter) / len(seq1)) * 100
    return similarity_score


clusters = {}
already_clustered = []
for i in range(len(seqs)-1):
    if seqs[i][0] in already_clustered:
        continue
    else:
        for j in range((i+1), len(seqs)):
            if seqs[j][0] not in clusters.values() and similarity(seqs[j][1], seqs[i][1]) >= 95:
                if seqs[i][0] not in clusters.keys():
                    clusters[seqs[i][0]] = [seqs[j][0]]
                else:
                    clusters[seqs[i][0]].append(seqs[j][0])
                already_clustered.append(seqs[j][0])


# Question 3:
with open('seq_clusters.tsv', 'w') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(["representative_sequence", "cluster_members"])
    for key in clusters.keys():
        writer.writerow([key, "", "", clusters[key]])

# Question 4:
result_array = np.zeros((100, 100), float)
for i in range(len(seqs)):
    for j in range((i+1), len(seqs)):
        result = 1 - ((similarity(seqs[i][1], seqs[j][1]))/100)
        result_array[i, j] = result
        result_array[j, i] = result


# Question 5:
mds = manifold.MDS(n_components=2, dissimilarity='precomputed')
pos = mds.fit(result_array).embedding_
nmds = manifold.MDS(n_components=2, metric=False, dissimilarity='precomputed')
npos = nmds.fit_transform(result_array, init=pos)

colors = ['#66CDAA', '#00FFFF', '#FF4040', '#FF6103', '#8B7355',
          '#FFD700',  '#FF6A6A', '#EE7AE9',  '#CAE1FF', '#CAFF70']
labels = ['seq_1', 'seq_2', 'seq_3', 'seq_4', 'seq_5',
          'seq_6', 'seq_9', 'seq_14', 'seq_17', 'seq_19']
fig = plt.figure(1)
j = 0
for key in labels:
    x = []
    for i in range(len(seqs)):
        if seqs[i][0] in clusters[key]:
            x.append(i)
    plt.scatter(npos[x, 0], npos[x, 1], color=colors[j], label=key)
    j = j + 1
plt.title('Scaling of 100 Gene Sequences', fontsize=20)
plt.xlabel('NMDS1', fontsize=16)
plt.ylabel('NMDS2', fontsize=16)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.tight_layout()
plt.show()


# Question 6:
with PdfPages('mds_seqs_clusters.pdf') as pdf:
    pdf.savefig(fig)
plt.close()
