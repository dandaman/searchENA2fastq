import warnings
from pathlib import Path
import os.path
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
onsuccess:
	d=pd.read_csv(query_result_file,sep="\t")
	c=0
	C=0
	missing=list()
	for i,r in d.iterrows():
		for f in r.fastq_ftp.rstrip(";").split(";"):
			F=os.path.basename(f)
			filename="{fastq_path}/{file}".format(fastq_path=fastq_path,file=F)
			p=Path(filename)
			if p.exists():
				c+=1
			else:
				missing.append("ftp://{url}".format(url=f))
			C+=1
	print("{downloaded}/{total} files downloaded!".format(downloaded=c,total=C))
	if len(missing)>0:
		print("The following URLs are missing:")
		print("\n".join(missing))
