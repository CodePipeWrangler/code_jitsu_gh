gen = pd.read_csv('/Users/bjordan/Documents/work/soy_comp/snps_proj_01/blastn.all.nt/out.blastn.Zh13.txt',
                  header=None,names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                  'qend', 'sstart', 'send', 'evalue', 'bitscore'],
                  sep='\t', comment='#')

gen = gen.sort_values(by = ['sseqid','sstart'])
test = gen.length.value_counts()
test2 = test.sort_index()
fig, ax = plt.subplots()
ax.barh(test2.index,test2)
ax.set_title('Glycine Max Zh13')
ax.set_xlabel('frequency')
ax.set_ylabel('base pairs')
plt.show()