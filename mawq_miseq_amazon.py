# Source activate.sh (of QIIME) before running QIIME on EC2 because non-interactive shell (such as ssh) has $PATH set
# differently compared to interactive (logged in) shell.

import os
import commands
import sys
import re
import subprocess
import logging
from datetime import datetime


def input_file_help():
	print """ 
	Help me please!!
	The input file should be tab-delimited file with .txt extension. The first column should be
	folder name of the MiSeq run. The second column should be the name of the mapping file of the run along 
	with its .txt extension. There should be no trailing white spaces or empty last lines. 
	Following is how a correct file should be:
	
	140401_M01869_0071_000000000-A7YEF	mapping_file_run1.txt
	140407_M01869_0073_000000000-A7WVG	mapping_file_run2.txt
	"""


def input_check(infile):
	""" Checks if input file name is entered correctly """
	if infile == "":
		print "Error: File name not provided!"
		mapfile = raw_input("1) Please provide the full name of the input-file (Type help for input-file format): ")
		return input_check(mapfile)
	elif infile.lower() == "help":
		input_file_help()
		mapfile = raw_input("1) Please provide the full name of the input-file (Type help for input-file format): ")
		return input_check(mapfile)
	else:
		working_folder = commands.getstatusoutput('pwd')[1]
		filelist = os.listdir(working_folder)
		if infile not in filelist:
			print "Error: File doesn't exist!"
			mapfile = raw_input("1) Please provide the full name of the input-file (Type help for input-file format): ")
			return input_check(mapfile)
		else:
			maplist = []
			infl = open(infile, 'rU')
			for line in infl:
				spline = line.strip().split("\t")
				if len(spline) != 2:
					print "Error: File is not in proper format. There's missing data, no tab-seperation and/or extra empty line(s)."
					mapfile = raw_input("1) Please provide the full name of the input-file (Type help for input-file format): ")
					return input_check(mapfile)
				else:
					maplist.append(spline[1])
			return maplist, infile				# Returns list of mapping files along with name of input file


def mapping_check(maplist):
	"""Checks if mapping file name is correct and runs validate_mapping_file.py script """
	for mapfile in maplist:
		with open(os.devnull, "w") as fnull:
			result = subprocess.call(["ls", mapfile], stdout = fnull, stderr = fnull)
			if result != 0:					# Return code is 0 is ls command is successful
				print "Error: One or more of your mapfiles is not present in your current working directory"
				mapfile2 = raw_input("1) Please provide the full name of the input-file (Type help for input-file format): ")
				inp_check = input_check(mapfile2)
				return mapping_check(inp_check[0])
	for mapfile in maplist:
		filename = mapfile.strip().split(".txt")[0]
		os.system("validate_mapping_file.py -m %s -o corrected_%s" % (mapfile.strip(),filename))
		os.system("mv $PWD/corrected_%s/%s_corrected.txt ." % (filename,filename))
	corrected_files = [mapfile.strip().split(".txt")[0]+"_corrected.txt" for mapfile in maplist]
	return corrected_files



def check_value(expression,question,arg):
	""" Function to check if input parameters are correct """
	try:
		if arg == "integer":
			return str(int(expression))
		if arg == "float":
			return str(float(expression))
	except ValueError:
		if expression == "":
			return expression
		else:
			print "Invalid value. Please enter a number or just hit enter for default value."
			checker = raw_input(question)
			return check_value(checker,question,arg)


def log_output(statement):
	""" Logs and prints output messages """
	logging.basicConfig(filename='logging_module_output.txt',level=logging.INFO)
	logging.info(statement)
	print statement


def log_parse(outfile,inputfile):
	output = open(outfile, "w")
	infile = open(inputfile, 'rU')
	for line in infile:
		if line.startswith("INFO:root:"):
			linename = line.strip().split("INFO:root:")
			if linename[1] != '':
				output.write(linename[1]+"\n")
		else:
			output.write(line.strip()+"\n")
	output.close()
	return output


def filtFast(flashread,indexfile):
	""" Shoko's script as a function """
	headerData = open(flashread,"rU").read().strip()
	headers = headerData.split("\n")
	IndexData = open(indexfile,"rU")
	IndexSeqs = {}
	while IndexData:
		headerline=IndexData.readline().split("\n")[0]
		Contents = ''
		if headerline == '':
			break
		for i in range(0,3):
			Contents += IndexData.readline().split("\n")[0] + "\n"
		IndexSeqs[headerline]=Contents
	outdata=''
	for j in headers:
		outdata += j + "\n" + IndexSeqs[j]
	of = open("Index_filtered_ordered.fastq","w")
	of.write(outdata)
	IndexData.close()
	of.close()
	return of

def preproceSteps(seq_data,m_min,read_len,runmap_file):
	""" Unzipping, flashing, pre-processing """
	if m_min == "":
		m_min += "225"
	if read_len == "":
		read_len += "251"
	if seq_data == "":
		seq_data += "/data/MiSeq_16S_data/MiSeqAnalysis"
	log_output("\nRead length for flash: %s" % read_len)
	log_output("Min. overlap for flash: %s" % m_min)
	folders = []
	infile = open(runmap_file, 'rU')
	for line in infile:
		spline = line.strip().split("\t")
		folders.append(spline[0])
	for seqs_id in folders:
		working_folder = commands.getstatusoutput('pwd')[1]
		seq_path = "%s/%s/Data/Intensities/BaseCalls/" % (seq_data,seqs_id)
		os.chdir(seq_path)
		log_output("\n#Step 1: Gunzipping sequence reads files in MiSeqAnalysis folder...")
		os.system("gunzip Undetermined_*")
		log_output("Gunzipping complete!")
		log_output("\n#Step 2:  Assembling R1 and R2 using flash...")			
		os.system("flash -r %s -f 300 -s 30 -m %s -d $PWD/Output_folder_%s/ -q Undetermined_S0_L001_R1_001.fastq Undetermined_S0_L001_R2_001.fastq" % (read_len, m_min,seqs_id))
		log_output("flash complete!")
		os.system("mv -f Output_folder_%s/ %s" % (seqs_id,working_folder))
		os.chdir(working_folder)
		log_output("\n#Step 3: Removing barcode reads from index file that are not in assembled file...")
		os.system("sed -n '1~4'p $PWD/Output_folder_%s/out.extendedFrags.fastq >FLAShReads.txt" % seqs_id)		# Select the headers of all sequences generated. -n flag is for quiet output. '1~4'p means starting from 1, select every 4 lines after it.
		log_output("Barcode removal complete!")
		log_output("\n#Step 4: Extracting those reads from index file and order them the same as flash reads")
		filtFast("FLAShReads.txt","%s/Undetermined_S0_L001_I1_001.fastq" % seq_path)
		log_output("Extraction complete!")
		os.chdir(seq_path)
		log_output("\n#Step 5: Gzipping back the sequence files in MiSeqAnalysis folder...")
		os.system("gzip Undetermined_S0_L001_*")
		os.chdir(working_folder)
		os.system("mv Index_filtered_ordered.fastq Index_filtered_ordered_run_%s.fastq" % seqs_id)
		log_output("Gzip complete!")
	return

def check_ec_permission(ec_key):
	""" Check file permission of EC2 key """
	val = 'ls -l %s' % ec_key
	result = list(subprocess.Popen([val], shell=True, stdout=subprocess.PIPE).stdout)
	if result[0].split("1")[0].strip() == "-r--------":
		print "No need to change file permission."
		return
	else:
		os.system("chmod 400 %s" % ec_key)
		return

def scp_amazon(runmap_file,ec_key,dns,anl_fold):
	""" Copy relevant files and folders to EC2 instance """
	run_map_dict = {}
	infile = open(runmap_file, 'rU')
	for line in infile:
		spline = line.strip().split("\t")
		run_map_dict[spline[0]] = spline[1].strip().split(".txt")[0]+"_corrected"+".txt"			#Run IDs as keys and mapping filenames as values
	mkdir = "mkdir %s/" % anl_fold
	os.system("ssh -i %s ubuntu@%s %s" % (ec_key,dns,mkdir))
	for fold_id in run_map_dict:
		folder = "Output_folder_"+fold_id
		log_output("\n#Step 6: Copying relevant files to the Amazon EC2 instance of QIIME...")
		os.system("scp -r -i %s %s ubuntu@%s:/home/ubuntu/%s" % (ec_key,folder,dns,anl_fold))
		os.system("scp -i %s Index_filtered_ordered_run_%s.fastq ubuntu@%s:/home/ubuntu/%s" % (ec_key,fold_id,dns,anl_fold))
		mapfile = run_map_dict[fold_id]
		os.system("scp -i %s %s ubuntu@%s:/home/ubuntu/%s" % (ec_key,mapfile,dns,anl_fold))
	
	os.system("scp -i %s %s ubuntu@%s:/home/ubuntu/%s" % (ec_key,runmap_file,dns,anl_fold))
	log_output("All files copied to EC2 instace!\n")
	return

def split_library(ec_key,dns,anl_fold,runmap_file,phred,max_bad_run,min_rl_frac,n_chars,barcode,start_seq):
	""" Splits libraries """
	if phred == "":
		phred += "30"
	if max_bad_run == "":
		max_bad_run += "3"
	if min_rl_frac == "":
		min_rl_frac += "0.75"
	if n_chars == "":
		n_chars += "0"
	if barcode == "":
		barcode += "12"
	if start_seq == "":
		start_seq += "0"
	log_output("Phred score: %s"	% phred)
	log_output("Max number of consecutive low quality base calls allowed before truncating a read: %s" % max_bad_run)
	log_output("Min number of consecutive high quality base calls to include a \
read (per single end read) as a fraction of the input read length: %s" % min_rl_frac)
	log_output("Max number of N characters allowed in a sequence to retain it: %s" % n_chars)
	log_output("The type of barcode used: %s" % barcode)
	log_output("The start seq_ids as ascending integers beginning with start_seq_id: %s" % start_seq)
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && mkdir fna_files/'" % (ec_key,dns,anl_fold))
	run_map_dict = {}
	infile = open(runmap_file, 'rU')
	for line in infile:
		spline = line.strip().split("\t")
		run_map_dict[spline[0]] = spline[1].strip().split(".txt")[0]+"_corrected"+".txt"			#Run IDs as keys and mapping filenames as values
	for fold_id in run_map_dict:
		folder = "Output_folder_"+fold_id
		mapfile = run_map_dict[fold_id]
		log_output("\n#Step 7: Splitting libraries using 'split_libraries_fastq.py'...")
		os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && split_libraries_fastq.py -i %s/out.extendedFrags.fastq -m %s \
		-o split_lib_output_%s/ -q %s -r %s -p %s -n %s\
		--rev_comp_barcode -b Index_filtered_ordered_run_%s.fastq \
		--barcode_type %s -s %s'" % (ec_key,dns,anl_fold,folder,mapfile,fold_id,phred,max_bad_run,min_rl_frac,n_chars,fold_id,barcode,start_seq))
		log_output("split_libraries_fastq.py complete!")
		os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && mv split_lib_output_%s/seqs.fna seqs_%s.fna'" % (ec_key,dns,anl_fold,fold_id,fold_id))
		os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && mv seqs_%s.fna fna_files/'" % (ec_key,dns,anl_fold,fold_id))
	return

def open_otus_till_biom(parallel,ref_db,ec_key,dns,anl_fold):
	""" Open OTU picking and other steps """
	if parallel == "":
		parallel += "4"
	if ref_db == "":
		ref_db += "/home/ubuntu/qiime_software/gg_otus-13_8-release/rep_set/97_otus.fasta"
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && cat fna_files/*.fna > fna_files/seqs_cat.fna'" % (ec_key,dns,anl_fold))
	log_output("\n#Step 8: Picking open-references OTUs using 'pick_open_reference_otus.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && pick_open_reference_otus.py -i fna_files/seqs_cat.fna -o open_otus_picked/ -aO %s -r %s'" % (ec_key,dns,anl_fold,parallel,ref_db))
	log_output("OTU picking caused errors, but we'll be able to proceed!")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && cp /home/ubuntu/qiime_software/core_set_aligned.fasta.imputed .'" % (ec_key,dns,anl_fold))
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && mv core_set_aligned.fasta.imputed core_set_aligned_imputed.fasta'" % (ec_key,dns,anl_fold))
	log_output("\n#Step 9.0: Aligning sequences to template using 'parallel_align_seqs_pynast.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && parallel_align_seqs_pynast.py -i open_otus_picked/rep_set.fna -o open_otus_picked/pynast_aligned_seqs \
	-t $PWD/core_set_aligned_imputed.fasta --jobs_to_start %s'" % (ec_key,dns,anl_fold,parallel))
	log_output("parallel_align_seqs_pynast.py complete!")
	log_output("\n#Step 9.1: Making OTU table by filtering alignment to remove sequences that did not align using 'make_otu_table.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && make_otu_table.py -i open_otus_picked/final_otu_map_mc2.txt -o open_otus_picked/otu_table_mc2_no_pynast_failures_w_tax.biom \
	-e open_otus_picked/pynast_aligned_seqs/rep_set_failures.fasta -t open_otus_picked/uclust_assigned_taxonomy/rep_set_tax_assignments.txt'" % (ec_key,dns,anl_fold))
	log_output("make_otu_table.py complete!")
	log_output("\n#Step 9.2: Identifying chimeric sequences using 'parallel_identify_chimeric_seqs.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && parallel_identify_chimeric_seqs.py -i open_otus_picked/pynast_aligned_seqs/rep_set_aligned.fasta -a core_set_aligned_imputed.fasta \
	-m ChimeraSlayer -o chimeraslayer_chimeric_seqs.txt -O %s'" % (ec_key,dns,anl_fold,parallel))
	log_output("parallel_identify_chimeric_seqs.py complete!")
	log_output("\n#Step 10: Filtering chimeric sequences out of the alignment file using 'filter_fasta.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && filter_fasta.py -f open_otus_picked/pynast_aligned_seqs/rep_set_aligned.fasta -o non_chimeric_rep_set_aligned.fasta \
	-s chimeraslayer_chimeric_seqs.txt -n'" % (ec_key,dns,anl_fold))
	log_output("filter_fasta.py complete!")
	log_output("\n#Step 11: Filtering non_chimeric_rep_set_aligned.fasta to remove gaps using 'filter_alignment.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && filter_alignment.py -i non_chimeric_rep_set_aligned.fasta -m /home/ubuntu/qiime_software/lanemask_in_1s_and_0s \
	-o non_chimeric_pynast_filtered/'" % (ec_key,dns,anl_fold))
	log_output("filter_alignment.py complete!")
	log_output("\n#Step 12: Building new phylogenetic tree using 'make_phylogeny.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && make_phylogeny.py -i non_chimeric_pynast_filtered/non_chimeric_rep_set_aligned_pfiltered.fasta \
	-o non_chimeric_rep_set_aligned_pfiltered.tre'" % (ec_key,dns,anl_fold))
	log_output("make_phylogeny.py complete!")
	log_output("\n#Step 13: Filtering chimeric OTUs from the OTU table using 'filter_otus_from_otu_table.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && filter_otus_from_otu_table.py -i open_otus_picked/otu_table_mc2_no_pynast_failures_w_tax.biom \
	-o otu_table_mc2_no_pynast_failures_no_chimeras_w_tax.biom -e chimeraslayer_chimeric_seqs.txt'" % (ec_key,dns,anl_fold))
	log_output("filter_otus_from_otu_table.py complete!")
	log_output("\n#Step 14: Writing biom table summary using 'biom summarize-table'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && biom summarize-table -i otu_table_mc2_no_pynast_failures_no_chimeras_w_tax.biom \
	-o otu_table_mc2_no_pynast_failures_no_chimeras_lowfilter_w_tax_biom_summary_mc2.txt'" % (ec_key,dns,anl_fold))	#0:00:01.602008
	log_output("biom summarize-table complete!")
	return

def summary_view(ec_key,dns,anl_fold,viewtable):
	""" Function to show biom summary table """
	if viewtable.lower() == 'yes':
		os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && less otu_table_mc2_no_pynast_failures_no_chimeras_lowfilter_w_tax_biom_summary_mc2.txt'" % (ec_key,dns,anl_fold))
	elif viewtable.lower() == 'no':
		print "No is not an option!"
		table = raw_input("The summary table of the final OTU table is ready. Type 'yes' to view it. \
Once viewed, you can quit by simply typing q. Are you ready? ")
		return summary_view(ec_key,dns,anl_fold,table)
	else:
		print "I don't understand."
		table = raw_input("The summary table of the final OTU table is ready. Type 'yes' to view it. \
Once viewed, you can quit by simply typing q. Are you ready? ")
		return summary_view(ec_key,dns,anl_fold,table)
	
def rarefaction_check(depth):
	""" Check value of rarefaction depth """
	try:
		return str(int(float(depth)))
	except ValueError:
		if depth == "":
			print "No number of sequences provided to subsample for rarefaction."
			dep = raw_input("1) What is the number of sequences to subsample per sample [-d flag]? (No default): ")
			return rarefaction_check(dep)
		else:
			print "Non-integer value given for number of sequences to subsample for rarefaction."
			dep = raw_input("1) What is the number of sequences to subsample per sample [-d flag]? (No default): ")
			return rarefaction_check(dep)	

def summary_plots(ec_key,dns,anl_fold,depth,merge_metadata):
	""" Create alpha, beta and taxa summary plots """
	log_output("\n#Step 15: Performing single rarefaction on OTU table using 'single_rarefaction.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && single_rarefaction.py -i otu_table_mc2_no_pynast_failures_no_chimeras_w_tax.biom -o single_rarefied_otu_table.biom -d %s'" % (ec_key,dns,anl_fold,depth))
	log_output("single_rarefaction.py complete!")
	log_output("\n#Step 16: Summarizing and plotting taxa using 'summarize_taxa_through_plots.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && summarize_taxa_through_plots.py -o taxa_summary -i single_rarefied_otu_table.biom -m %s'" % (ec_key,dns,anl_fold,merge_metadata))
	log_output("summarize_taxa_through_plots.py complete!")
	log_output("\n#Step 17: Calculating alpha-diversity using 'alpha_rarefaction.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && alpha_rarefaction.py -i single_rarefied_otu_table.biom -o alpha_rarefaction/ -t non_chimeric_rep_set_aligned_pfiltered.tre \
	-m %s --retain_intermediate_files'" % (ec_key,dns,anl_fold,merge_metadata))
	log_output("alpha_rarefaction.py complete!")
	log_output("\n#Step 18: Calculating beta-diversity using 'beta_diversity_through_plots.py'...")
	os.system("ssh -i %s ubuntu@%s 'source $HOME/qiime_software/activate.sh; cd %s && beta_diversity_through_plots.py -i single_rarefied_otu_table.biom -o beta_diversity/ -t non_chimeric_rep_set_aligned_pfiltered.tre \
	-m %s'" % (ec_key,dns,anl_fold,merge_metadata))
	log_output("beta_diversity_through_plots.py complete!")
	return


if __name__ == "__main__":
	print "\n\t\t\t\033[1mWelcome to the Microbiome Analysis through Workflow of QIIME, MAWQ program (pronounced 'mock') brought to you by the Lynch Lab!\033[0m"
	print "\tTo run the script with default parameters, just press enter to each question without entering a value. To \
exit the pipeline at any point in time, press Ctrl+C\n\n"
	try:
		inputfile = raw_input("1) Please provide the full name of the input-file (Type help for input-file format): ")
		checked = input_check(inputfile)
		inputfile = checked[1]
		mapping_check(checked[0])	
		seq_data = raw_input("\n1) What's the path to the MiSeq run folder? (Default: /data/MiSeq_16S_data/MiSeqAnalysis) ")
		print "\nThe following questions are for flash program: \n"
		flash_q1 = "1) What's the minimum overlap length between reads [-m flag]? (Default: 225, if  length of Read 2 > 250) "
		m_min = check_value(raw_input(flash_q1),flash_q1,"integer") 
		flash_q2 = "2) What's the read length [-r flag]? (Default: 251) "
		read_len = check_value(raw_input(flash_q2),flash_q2,"integer")
		print "\nThe following questions are for split_libraries_fastq.py script: \n"
		split_q1 = "1) What's the maximum unacceptable Phred quality score [-q flag]? (Default: 30) "
		phred = check_value(raw_input(split_q1),split_q1,"integer")
		split_q2 = "2) What's the max number of consecutive low quality base calls allowed before truncating a read [-r flag]? (Default: 3) "
		max_bad_run = check_value(raw_input(split_q2),split_q2,"integer")
		split_q3 = "3) What's the min number of consecutive high quality base calls to include a read (per single end read) as a fraction \
of the input read length [-p flag]? (Default: 0.75) "
		min_rl_frac = check_value(raw_input(split_q3),split_q3,"float")
		split_q4 = "4) What's the max number of N characters allowed in a sequence to retain it [-n flag]? (Default: 0) "
		n_chars = check_value(raw_input(split_q4),split_q4,"integer")
		split_q5 = "5) What's the type of barcode used [--barcode_type flag]? (Default: 12) "
		barcode = check_value(raw_input(split_q5),split_q5,"integer")
		split_q6 = "6) What's the start seq_ids as ascending integers beginning with start_seq_id [-s flag]? (Default: 0) "
		start_seq = check_value(raw_input(split_q6),split_q6,"integer")
		print "\nThe following questions are for pick_open_reference_otus.py script: \n"
		otupick_q1 = "1) How many jobs do you wish to run in parallel? (Default: 4) "
		parallel = check_value(raw_input(otupick_q1),otupick_q1,"integer")
		ref_db = raw_input("2) What's the full path to the reference database? \
(Default: /home/ubuntu/qiime_software/gg_otus-13_8-release/rep_set/97_otus.fasta) ")
		print "\nThe following questions are regarding your Amazon EC2 instance\n"
		ec_key = raw_input("1) Please enter the keypair for your EC2 instance: ")	#qiime_wrapper_data_storage.pem
		check_ec_permission(ec_key)
		dns = raw_input("2) Please enter the public DNS for your EC2 instance: ")	#ec2-54-237-163-108.compute-1.amazonaws.com
		anl_fold = raw_input("3) What would you like to call your analysis folder on the EC2 instance? ")
		startTime = datetime.now()
		preproceSteps(seq_data,m_min,read_len,inputfile)
		scp_amazon(inputfile,ec_key,dns,anl_fold)
		split_library(ec_key,dns,anl_fold,inputfile,phred,max_bad_run,min_rl_frac,n_chars,barcode,start_seq)
		open_otus_till_biom(parallel,ref_db,ec_key,dns,anl_fold)
		viewtable = raw_input("The summary table of the final OTU table is ready. Type 'yes' to view it. \
Are you ready? ")
		summary_view(ec_key,dns,anl_fold,viewtable)
		print "\nThe following question is for 'single_rarefaction.py' script: \n"
		depth = raw_input("1) What is the number of sequences to subsample per sample [-d flag]? (No default): ")
		depth = rarefaction_check(depth)
		print "\nThe following question is for 'summarize_taxa_through_plots.py', 'alpha_rarefaction.py', \
and 'beta_diversity_through_plots.py' script: \n"
		merge_metadata = raw_input("1) What is the name of the final mapping data file for all runs [-m flag]? (No default): ")
		merge_metadata_checked = mapping_check([merge_metadata])
		os.system("scp -i %s %s ubuntu@%s:/home/ubuntu/%s" % (ec_key,merge_metadata_checked[0],dns,anl_fold))
		summary_plots(ec_key,dns,anl_fold,depth,merge_metadata_checked[0])	
		log_parse("wrapper_log_file.txt","logging_module_output.txt")
		os.system('rm logging_module_output.txt')
		print "\n"+"Task Completed! Time it took to complete the task: "+ str(datetime.now()-startTime)
	except KeyboardInterrupt:
		print "\n\nThanks for using (or attempting to use) the pipeline. Good-bye!\n"
