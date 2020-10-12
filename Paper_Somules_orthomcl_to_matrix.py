'''
Example usage:
 python3 ~alt/python/bin/orthomcl_to_matrix_rewrite.py -s -r

Script to capture as many 1:1s as possible between Schisto genes and Schmidtea genes using orthomcl output and update an input matrix containing
Smps names and data values with Schmidtea names where possible.  Also can return a summary file of the 1:1s found and how they were found.

'''
#----------------IMPORTS----------------

import re
import random
from itertools import islice
from collections import defaultdict
import copy	#Necessary for doing a deep copy of our matrix dictionary
import argparse
import sys
from datetime import datetime



#----------------Random seed----------------

seed=1337


#----------------Make matrix----------------

#Read matrix into dictionary - Smps as keys, string of many counts as values
def make_matrix_dict():

	matrix={}
	with open('/nfs/repository/working_area/SHISTO/V7/cds1_single_cell/matrix_somules_alan.txt','r') as f:
		head=list(islice(f, 1))	#Capture header and move to next line
		for line in f:
			x=line.split()
			y=x[0].replace('"','')
			matrix[y]=x[1:]	#y = Smp, x[1:] = counts

	return matrix, head


#----------------HAPLOTYPES----------------

#All Smps on haplo contigs that we know about
def get_haplos():

	haplo=[]

	with open('/nfs/repository/working_area/SHISTO/V7/haplo_Smps.list', 'r') as f:
		for line in f:
			x=line.split()[0].replace('ID=','').strip()
			haplo.append(x)

	return haplo

#--------------------GET 1:1 lines list---------------------

#Small number of cases where Smp doesn't get a unique dd_smed:
#Smp_346580 ['dd_Smed_v6_13424_0', 'dd_Smed_v6_10765_0']
#Smp_345560 ['dd_Smed_v6_12586_0', 'dd_Smed_v6_6281_0']
#Smp_344030 ['dd_Smed_v6_8914_0', 'dd_Smed_v6_8002_0']
#Smp_337140 ['dd_Smed_v6_3442_0', 'dd_Smed_v6_7430_0']
#Smp_346010 ['dd_Smed_v6_46627_0', 'dd_Smed_v6_9445_0']
#Smp_343910 ['dd_Smed_v6_13660_0', 'dd_Smed_v6_6959_0']
#Smp_337860 ['dd_Smed_v6_5624_0', 'dd_Smed_v6_6147_0']
#Smp_337390 ['dd_Smed_v6_10276_0', 'dd_Smed_v6_2458_0']
#Smp_315740 ['dd_Smed_v6_7564_0', 'dd_Smed_v6_567_0']
#Smp_157710 ['dd_Smed_v6_8460_0', 'dd_Smed_v6_10541_0']
#Smp_126540 ['dd_Smed_v6_6811_0', 'dd_Smed_v6_14517_0']
#Smp_337090 ['dd_Smed_v6_236_0', 'dd_Smed_v6_10589_0']
# We revert to Smps for these cases

#Figure out the obvious 1:1s after taking haplo Smps into account.  Also returns dds + a count for each cluster they occur in as a dictionary
def get_obvious_one_to_ones(haplos, matrix_dict):

	one2ones_dict=defaultdict(list)
	dds_cluster_counts_dict=defaultdict(lambda:0)	#Here lambda initialises the dict values to 0 all in one go
	with open('/nfs/repository/working_area/SHISTO/V7/cds1_single_cell/orthomcl25118_Sm_v7.1.pep_Smed_proteins_planmine/all_orthomcl.out', 'r') as f:
		for line in f:
			Smps=[]
			Smps_minus_haplos=[]
			dds=[]
			x=line.strip().split()
			for i in x:
				if 'Sm_v7.1.pep' in i:
					Smp=i.split('(')[0].split('.')[0]	#This gets just gene root
					assert(len(Smp) == 10)
					Smps.append(Smp)
				elif 'Smed_proteins_planmine.aa' in i:
					ddsmed=('_').join(i.split('(')[0].split('.')[0].split('_')[:5])	#This gets just gene root
					assert(len(ddsmed.split('_')) == 5)
					dds.append(ddsmed)
			for i in set(dds):	#Bug found here - had forgotten to set the list to make entries unique
				dds_cluster_counts_dict[i]+=1
			if len(set(Smps)) == 2 and len(set(dds)) == 1:	#Haplotype cases should only ever be 1 or 2 Smps at gene level.  Haps list is only list of genes on haplo contigs - if we considered lists longer than 2 we might get other cases that are not straightforward haplotypes.  There are only 4 cases in total if we looked at lists longer than 2.
#Example of haplotype case that is now dealt with correctly:
#ORTHOMCL2340(4 genes,2 taxa):	 Smp_175200.1(Sm_v7.1.pep) Smp_175200.2(Sm_v7.1.pep) Smp_319600.1(Sm_v7.1.pep) dd_Smed_v6_6555_0_11833(Smed_proteins_planmine.aa)

				Smp1 = list(set(Smps))[0]
				Smp2 = list(set(Smps))[1]
				if Smp1 or Smp2 in haplos:
					if Smp1 in matrix_dict:		#Here we need to refer to matrix dict to ensure we capture the correct Smp to swap for dd_Smed later
						if Smp1 not in one2ones_dict:
							one2ones_dict[Smp1]=list(set(dds))[0]
					elif Smp2 in matrix_dict:
						if Smp2 not in one2ones_dict:
							one2ones_dict[Smp2]=list(set(dds))[0]
			elif len(set(Smps)) == 1 and len(set(dds)) == 1:
				Smp = list(set(Smps))[0]
				if Smp not in one2ones_dict:
					one2ones_dict[Smp]=list(set(dds))[0]
				else:
					removed=one2ones_dict.pop(Smp)	#Ensures always only 1 dd is added per Smp





	revised_one2ones_dict = copy.deepcopy(one2ones_dict)	#Passing back over to ensure that apparent 1:1s in first pass don't contain multi-cluster dds now we've obtained counts	
	for k, v in one2ones_dict.items():
		if dds_cluster_counts_dict[v] > 1:	#if multicluster ddsmed
			#print(k, v, dds_cluster_counts_dict[v])
			if k in revised_one2ones_dict:
				del revised_one2ones_dict[k]	#Remove the 1:1 from the dictionary


	return revised_one2ones_dict, dds_cluster_counts_dict
				
			
#----------------GET 1:1s by random selection from multidds that don't cross clusters----------------	

def choose_randomly_from_multi_dds(haplos, dds_cluster_counts_dict):

	random_one_to_ones_dict={}
	with open('/nfs/repository/working_area/SHISTO/V7/cds1_single_cell/orthomcl25118_Sm_v7.1.pep_Smed_proteins_planmine/all_orthomcl.out', 'r') as f:
		
		for line in f:
			Smps=[]
			Smps_minus_haplos=[]
			dds=[]
			rand_dd_list=[]
			x=line.strip().split()
			for i in x:
				if 'Sm_v7.1.pep' in i:
					Smp=i.split('(')[0].split('.')[0]	#This gets just gene root
					Smps.append(Smp)
				elif 'Smed_proteins_planmine.aa' in i:
					ddsmed=('_').join(i.split('(')[0].split('.')[0].split('_')[:5])	#This gets just gene root
					dds.append(ddsmed)
			if len(set(Smps)) >= 1:
				for s in Smps:
					if s not in haplos:
						Smps_minus_haplos.append(s)
			if len(set(Smps_minus_haplos)) == 1 and len(set(dds)) > 1:
				for i in set(dds):
					if dds_cluster_counts_dict[i] == 1:
						rand_dd_list.append(i)
				if len(rand_dd_list) >1:	#ie if we have something sensible to add
					selected_dd=random.sample(list(set(rand_dd_list)), 1)[0]
					Smp=list(set(Smps_minus_haplos))[0]
					if Smp not in random_one_to_ones_dict:
						random_one_to_ones_dict[Smp]=selected_dd
					else:
						random_one_to_ones_dict.pop(Smp)	#Ensures always only 1 dd is added per Smp

				
	return random_one_to_ones_dict	#k=Smp, v=1xddsmed	


#----------------Check one2one dictionaries----------------

#Written in response to testing which found 1 edge case where the Smp had been captured once in each dictionary:

#ORTHOMCL681(7 genes,2 taxa):	 Smp_210180.2(Sm_v7.1.pep) Smp_210180.3(Sm_v7.1.pep) dd_Smed_v6_1564_0_1219(Smed_proteins_planmine.aa) dd_Smed_v6_8353_0_1303(Smed_proteins_planmine.aa) dd_Smed_v6_8353_0_2303(Smed_proteins_planmine.aa) dd_Smed_v6_8353_0_3237(Smed_proteins_planmine.aa) dd_Smed_v6_8353_0_4237(Smed_proteins_planmine.aa)
#ORTHOMCL7994(2 genes,2 taxa):	 Smp_210180.1(Sm_v7.1.pep) dd_Smed_v6_9507_0_1507(Smed_proteins_planmine.aa)

def check_dicts(one2ones_dict, random_one_to_ones_dict):

	modified_random_one_to_ones_dict = copy.deepcopy(random_one_to_ones_dict)
	for k, v in random_one_to_ones_dict.items():
		if k in one2ones_dict:
			del modified_random_one_to_ones_dict[k]

	return modified_random_one_to_ones_dict


#----------------Update matrix with 1:1s----------------

def update_matrix_from_one_to_ones(matrix_dict,one2ones_dict):

	matrix_dict2={}		#This is an updated matrix for just simple 1:1s taking haplotypes into consideration

	for k, v in matrix_dict.items():	#7398 lines after header removed

		if k in one2ones_dict:
			ddsmed=one2ones_dict[k]	
			if ddsmed not in matrix_dict2:
				matrix_dict2[ddsmed] = v
			else:
				matrix_dict2[k] = v	#Here ddsmed doesn't have a unique Smp so line stays unchanged with Smp key
		else:
			if k not in matrix_dict2:
				matrix_dict2[k] = v

	return matrix_dict2

#----------------Update matrix with random ddsmeds that don't cross clusters----------------

def update_matrix_with_random_smeds(matrix_dict2,random_one_to_ones_dict):

	matrix_dict3 = copy.deepcopy(matrix_dict2)
	for k, v in matrix_dict2.items():
		if k in random_one_to_ones_dict:	#Considering Smps we want to update....
			ddsmed=random_one_to_ones_dict[k]
			if ddsmed not in matrix_dict3:
				matrix_dict3[ddsmed] = matrix_dict3.pop(k)

	return matrix_dict3

#----------------Scramble random numbers for file naming----------------

def scrambled(orig):

    dest = orig[:]
    random.shuffle(dest)
    return dest

#----------------Make pseudo random string for file naming----------------

def get_now():

	x=[]
	v=str(datetime.now())	#date time string
	chs=[' ', ':', '.', '-']
	for ch in chs:
		if ch in v:
			v=v.replace(ch, '')
	for i in v:
		x.append(i)
	y=scrambled(x)[:9]
	rnd=('').join(y)
	return rnd

#----------------Write matrix----------------

def write_matrix(matrix_dict3, head, pseudorand, matrix_type):

	assert(matrix_type == 'reproducible' or matrix_type == 'variable')
	if matrix_type == 'reproducible':
		outfile='reproducible_seeded_'+str(seed)+'.matrix'
	elif matrix_type == 'variable':
		outfile=pseudorand+'.matrix'
	with open(outfile, 'w') as f:
		h=str(head[0])+'\n'
		f.write(h)	#Write the header
		for k, v in matrix_dict3.items():
			line='"'+k+'"'+'\t'+str("\t".join(v))+'\n'
			f.write(line)

#----------------Write summary----------------

def write_summary(one2ones_dict, modified_random_one_to_ones_dict, pseudorand, summary_type):

	assert(summary_type == 'reproducible' or summary_type == 'variable')
	if summary_type == 'reproducible':
		outfile='reproducible_seeded_'+str(seed)+'_summary.txt'
	elif summary_type == 'variable':
		outfile=pseudorand+'.summary.txt'
	with open(outfile, 'w') as f:
		for k, v in one2ones_dict.items():
			line='Straightforward_1:1s:\t'+k+'\t'+v+'\n'
			f.write(line)
		for k, v in modified_random_one_to_ones_dict.items():
			line='Randomly_selected_1:1s\t'+k+'\t'+v+'\n'
			f.write(line)

#---------------------MAIN------------------------


def main():

	parser = argparse.ArgumentParser(description='Script that takes a matrix of S.mansoni gene names and their counts and using an orthomcl comparison of this genome and Schmidtea attempts to update the Smp names in the matrix with Schmidtea names where direct 1:1s can be found, or chooses a Schmidtea randomly where 1 S.mansoni gene hits several Schmidtea genes that don\'t cross clusters.  The script gives the option to produce reproducible output from the randomisation step or non-reproducible.  Non-reproducible may be desired to test the effect of different randomisation outcomes.  The script also gives the option to produce a summary of 1:1s used.') 

	#positional args
	parser.add_argument('-r','---reproducible', action='store_true', default=False, help='Do you want to get the same result each time you run?')
	parser.add_argument('-s','--summary', action='store_true', default=False, help='Do you want a summary file of 1:1s used?')
	args = parser.parse_args() #gets the arguments	
	if args.reproducible:
		random_seed_we_picked=seed

	(matrix_dict, head)=make_matrix_dict()	#k=Smp, v=counts
	haplos=get_haplos()
	(one2ones_dict, dds_cluster_counts_dict)=get_obvious_one_to_ones(haplos, matrix_dict)
	random_one_to_ones_dict=choose_randomly_from_multi_dds(haplos, dds_cluster_counts_dict)

	modified_random_one_to_ones_dict=check_dicts(one2ones_dict, random_one_to_ones_dict)


	matrix_dict2=update_matrix_from_one_to_ones(matrix_dict,one2ones_dict)
	#At this point, matrix_dict2 has 7398 lines and 3891 ddmseds as keys - ddsmeds are all simple 1:1 cases
	matrix_dict3=update_matrix_with_random_smeds(matrix_dict2,modified_random_one_to_ones_dict)
	#At this point, ddsmeds = 4397/7398 - additional ddsmeds come from random selection of multiple dds where the chosen ddsmed does not span multiple orthomcl clusters

	pseudorand=get_now()	
	if args.summary and args.reproducible:
		write_summary(one2ones_dict, modified_random_one_to_ones_dict, pseudorand, 'reproducible')
		write_matrix(matrix_dict3, head, pseudorand, 'reproducible')
	elif not args.summary and args.reproducible:
		write_matrix(matrix_dict3, head, pseudorand, 'reproducible')
	elif args.summary and not args.reproducible:
		write_summary(one2ones_dict, modified_random_one_to_ones_dict, pseudorand, 'variable')
		write_matrix(matrix_dict3, head, pseudorand, 'variable')
	elif not args.summary and not args.reproducible:
		write_matrix(matrix_dict3, head, pseudorand, 'variable')

	print('\nFINISHED\n')




if __name__=="__main__":
	main()
