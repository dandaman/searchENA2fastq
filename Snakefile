import pandas as pd
import urllib 
import requests
from io import StringIO

query_result_file="ena.txt"
ENA_urls="ena.urls.txt"
fastq_path="fastq"

localrules: all, query_ENA

rule all:
	input:
		dynamic(fastq_path+"/{pair}.fastq.gz")

rule download_fastqs:
	input:
		"urls/{pair}.txt"
	output:
		fastq_path+"/{pair}.fastq.gz"
	shell:
		"wget -P "+fastq_path+"/ -i {input}"

rule query_ENA:
	input:
		ENA_urls
	output:
		temp(dynamic("urls/{pair}.txt"))
	run:
		with open(input[0]) as x: urls = list(map(lambda y: y.rstrip(),x.readlines()))
		frames=list()
		for url in urls:
			response = requests.post(url)
			text=StringIO(response.text)
			d=pd.read_csv(text,sep="\t")
			frames.append(d)
		(s,e,r)=frames
		m=pd.merge(pd.merge(s,e,on=list(set(s.columns) & set(e.columns)),how="inner"),r,on=list(set(e.columns) & set(r.columns)),how="inner")
		m=m[m.fastq_ftp.notnull()]
		m.to_csv(query_result_file,sep="\t",index=False)
		for index,row in m.iterrows():
			counter=1
			for item in row.fastq_ftp.rstrip(";").split(";"):
				filename="urls/{run}_{pair}.txt".format(run=row.run_accession,pair=counter)
				with open(filename, 'w') as f:
					f.write("ftp://%s\n" % item)
				counter+=1
