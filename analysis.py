import pandas as pd

### Get EAF from Finngen
freq_path = 'finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9.afreq'
snp_path = 'finngen/green/FeiyiWang/all_variants_and_proxies.csv'

freq = pd.read_csv(freq_path, sep='\t')
snp = pd.read_csv(snp_path)

snp = snp.merge(freq[['ID','ALT_FREQS']], left_on='format', right_on='ID')
# len(snp[snp.ALT_FREQS.isna()]) is 0
snp['finn_eaf'] = 1 - snp.ALT_FREQS
