1. generate a txt file with format chrN:12345678-12345678 using chromosome and hg19_position from file all_variants_and_proxies.csv

'''
import pandas as pd
df = pd.read_csv('~/Documents/Materials/iron/all_variants_and_proxies.csv')
df['format'] = 'chr'+df.chromosome.astype(str)+':'+df.hg19_position.astype(str)+'-'+df.hg19_position.astype(str)
df['format'].to_csv('~/Documents/Projects/iron/variants_hg19.txt', index=False)
'''

2. convert hg19 coordinates to hg38 coordinates using liftover tool: https://genome.ucsc.edu/cgi-bin/hgLiftOver

3. rename the file downloaded from the tool as variants_hg19.bed

4. create a list of format that can be recognized in finngen

'''
import pandas as pd
df = pd.read_csv('~/Documents/Materials/iron/all_variants_and_proxies.csv')
df0 = pd.read_csv('~/Documents/Materials/iron/variants_hg38.bed', header=None)
df['hg38_position'] = df0[0].str.extract('\:(\d+)\-')
df['format'] = 'chr'+df.chromosome.astype(str)+'_'+df.hg38_position+'_'+df.ea+'_'+df.nea
df['format_'] = 'chr'+df.chromosome.astype(str)+'_'+df.hg38_position+'_'+df.nea+'_'+df.ea
df.to_csv('~/Documents/Materials/iron/all_variants_and_proxies.csv', index=False)

df.format.to_list()
'''

5. upload the file to green bucket: https://console.cloud.google.com/storage/browser/fg-production-sandbox-6_green

'''
gsutil cp ~/Documents/Materials/iron/all_variants_and_proxies.csv gs://fg-production-sandbox-6_green/FeiyiWang
'''


