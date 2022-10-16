"""
Calculate model M0
"""

import argparse,os,sys
import logging
import csv
import re
import subprocess
from Queue import Queue
from threading import Thread
from ete3 import EvolTree
from ete3.evol.model import Model
from ete3.evol.control import PARAMS
from ete3 import SeqGroup
from Bio import SeqIO


class MyParser(argparse.ArgumentParser):
 	def error(self, message):
 	 	sys.stderr.write('error: %s\n' % message)
 	 	self.print_help()
 	 	sys.exit(2)

parser=MyParser()
#parser = argparse.ArgumentParser()
parser.add_argument('--inputDIR', help='Path of the input directory')
parser.add_argument('--p', help='Number of processors')
parser.add_argument('--log_file', help='Path to the log file')
parser.add_argument('--out_file', help='Path to the out file')

if len(sys.argv)==1:
 	parser.print_help()
 	sys.exit(1)

args = parser.parse_args()

if args.inputDIR:
 	inputDIR = args.inputDIR
 	if inputDIR[-1] != "/":
 	 	inputDIR += "/"

if args.p:
 	num_cores = args.p

if args.log_file:
 	log_file = args.log_file

if args.out_file:
 	out_file = args.out_file


"""
Functions
"""

#Returns the number of sequence of a phylip file using the header
def number_seq_phy(phy_path):
 	with open(phy_path, "r") as phy_file:
 	 	line = phy_file.readline().strip().split()
 	 	return int(line[0])

#Convert sequential phylip files to fasta
def phylip2fasta(path_phylip,path_fasta):
 	#Open the sequential phylip file
 	with open(path_phylip, "r") as phylip_file_path:
 	 	#Parse to a list
 	 	lines = phylip_file_path.readlines()
 	#Open the fasta file and parse
 	with open(path_fasta, "w") as fasta_file_id:
 	 	for i in range(1,len(lines)):
 	 	 	seq = lines[i].split()
 	 	 	fasta_file_id.write(">" + seq[0] + "\n" + seq[1]+ "\n")
 	return path_fasta


#Create a directory
def create_dir(dir_path):
 	if not os.path.exists(dir_path):
 	 	os.makedirs(dir_path)
 	return dir_path


#Implement codeml in parallel
def mutitask(q,files_list,processors):
 	num_threads = int(processors)
 	#Put the list of files
 	for i in range(len(files_list)):
 	 	q.put(files_list[i])

 	#set up the worker threads
 	for i in range(num_threads):
 	 	worker = Thread(target=paml_M0, args=(q,))
 	 	worker.setDaemon(True) 	  #setting threads as "daemon" allows main program to 
 	 	 	 	 	 	 	 	 	#exit eventually even if these dont finish 
 	 	 	 	 	 	 	 	 	#correctly.
 	 	worker.start()
 
 	#now we wait until the queue has been processed
 	q.join()
 	logging.info('All codeml tests are finished.')
 	return True


#Calculate M0 to get the global value of omega using codeML and evoltree.
def paml_M0(q):
 	while not q.empty():
 	 	work = q.get()  
 	 	fasta_file = work[0]
 	 	tree_file = work[1]
 	 	group_dir = work[2]

 	 	#Enter the tree, align and working directory using EvolTree.
 	 	tree = EvolTree(tree_file)
 	 	tree.link_to_alignment(fasta_file,alg_format='fasta')
 	 	tree.workdir = group_dir

 	 	#ClusterID and number of sequences.
 	 	clusterID = os.path.basename(os.path.dirname(work[0]))
 	 	Number_seq = str(len(tree))

 	 	#Test the evolutionary model 0.
 	 	#Creating the directory structure
 	 	dir_M0 = create_dir(group_dir + "M0/") 

 	 	if not os.path.exists(dir_M0 + "out"):
 	 	 	#Setting the initial omega values to 1.
 	 	 	M_0 = Model('M0', tree)
 	 	 	M_0._change_params(dict(list(PARAMS.items())))
 	 	 	M_0.properties['params']['omega'] = 1

 	 	 	#Save the model in tmp.ctl files.
 	 	 	with open(dir_M0 + "tmp.ctl" ,"w") as handle_file:
 	 	 	 	handle_file.write(M_0.get_ctrl_string())

 	 	 	#Save the tree files.
 	 	 	with open(dir_M0 + "tree", "w") as handle_tree:
 	 	 	 	handle_tree.write(tree.write(format=1)+ "\n")

 	 	 	#Save the alignment files.
 	 	 	seqs = SeqGroup(fasta_file,format = "fasta")
 	 	 	seqs.write(format="paml", outfile = dir_M0 + "algn")

 	 	 	cmd = ["codeml", "tmp.ctl"]
 	 	 	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd= dir_M0)
 	 	 	out = p.communicate()

 	 	 	logging.info(clusterID + ": codeml has finished the analysis.")

 	 	else:
 	 	 	logging.info(clusterID + ": Model 0 had previously been run.")

 	 	#Load the paml results.
 	 	tree.link_to_evol_model(dir_M0 + "out", 'M0')

 	 	#Parse the previous results
 	 	Model_0 = tree.get_evol_model('M0')

 	 	#Obtain the global omega.
 	 	Model_0 = tree.get_evol_model('M0')
 	 	w_global = str(Model_0.branches[1]['w'])

 	 	#Calculate dN
 	 	dN=[]
 	 	for value in Model_0.branches:
 	 	 	if 'dN' in Model_0.branches[value]:
 	 	 	 	dN.append(Model_0.branches[value]['dN'])
 	 	tree_length_dN = str(round(sum(dN), 4))

 	 	#Calculate dS
 	 	dS=[]
 	 	for value in Model_0.branches:
 	 	 	if 'dS' in Model_0.branches[value]:
 	 	 	 	dS.append(Model_0.branches[value]['dS'])
 	 	tree_length_dS = str(round(sum(dS), 4))

 	 	#Obtain the lnL of the model 0
 	 	lnL_M0 = str(Model_0.lnL)

 	 	#Store the results in a list
 	 	current_result = [clusterID,Number_seq,lnL_M0,tree_length_dN,tree_length_dS,w_global]
 	 	final_results.append(current_result)
 	 	q.task_done()
 	return True

'''
Main
'''

#Avoid overwrite previous result and log files.
assert not os.path.exists(log_file), "WARNING: The file " + log_file + " already exist. Please rename or move the previous log file to prevent information loss."
assert not os.path.exists(out_file), "WARNING: The file " + out_file + " already exist. Please rename or move the previous results to prevent information loss."

#Create the log file.
logging.basicConfig(filename=log_file, filemode='w', format='%(levelname)s:%(message)s', level=logging.DEBUG)

filecount = 0
q = Queue(maxsize=0)
file_lists_codeml = []
final_results = []

for i in os.listdir(inputDIR):
 	group_dir = inputDIR + i + "/"
 	#Creating the list of tuples with phylip file, tree and working directory for EvolTree.
 	path_new_fasta_file = os.path.join(group_dir,i + "_filtered.fas")
 	path_raxml_bestTree = os.path.join(group_dir,i + "_filtered.fas.raxml.support")
 	three_inputs_ete3 = (path_new_fasta_file,path_raxml_bestTree,group_dir,)
 	file_lists_codeml.append(three_inputs_ete3)
 	filecount = filecount +1


assert filecount > 0, "No phylip files and/or id_names files"

#Running the functions
mutitask(q,file_lists_codeml,num_cores)


#Save the results in branch_site_results file in the results directory creating by POTION. 
with open(out_file,"w") as f:
 	wr = csv.writer(f,delimiter='\t')
 	wr.writerows(final_results)