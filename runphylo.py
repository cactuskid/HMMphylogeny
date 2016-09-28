from Bio import AlignIO, SeqIO , SearchIO
import Bio.Blast.NCBIStandalone as BlastTools
import numpy as np
import networkx as nx
import os
import glob
from csb.bio.io.hhpred import HHOutputParser
import subprocess , shlex
import pickle
import tempfile
from ete3 import PhyloTree
from modeller import *

from bayes_opt import BayesianOptimization


#use for merging alns
from scipy.spatial import distance as dist
from scipy.cluster import hierarchy as hrc

def save_obj(obj, name ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open( name + '.pkl', 'r') as f:
        return pickle.load(f)

#folder for the starting batch of sequences. this folder should already be created and have your sequence in it
sequences = './finalsequencesNoretro/'


#the following folders will be generated automatically.
#folder where the HHblits results against the uniprot 20 and the HMMs will be kept
modeldir = sequences+'models/'
#folder where the structural prediction results will be kept
structuralpreddir = sequences+'structurepred/'
#final results of the allVall searches and calculated newick trees will be stored here.
resultsdir = sequences+'results/'

#how is your machine configured?
nCPU= 30
#paths to HHblits/HHsearch software
HHsearchpath = 'sudo hhsearch'
hhalignpath = 'sudo hhalign'
hhblitspath = 'sudo hhblits'
hhmakepath = 'sudo hhmake'
AddSSpath = '/usr/lib/hhsuite/scripts/addss.pl'
reformatpath = '/usr/lib/hhsuite/scripts/reformat.pl'
fastmepath = 'sudo fastme'
clustalopath = 'sudo clustalo '
#paths to DBs
UniprotHHDB ='/implacable/sync/uniprot20_2016_02/uniprot20_2016_02' 
PDbHHDB = '/implacable/pdb70/pdb70'


#Turn different steps of the workflow on or off.
prepare_models = False
cleanModels = False

check_structure = False

runphylo = True
runDistmat = True
merge_alns = True

def runHHblits( prot , pargs):
	#Usage: hhblits -i query [options]
	print 'running hhblits on ' + prot
	path , outdir, db , iterations , ncores , runName , SS  = pargs
	name = prot.split('/')[-1]
	outhhr= outdir+name.replace('.fasta','')+runName+".hhr"
	outa3m = outdir+name.replace('.fasta', '.a3m')
	outfasta = outdir+name.replace('.fasta', 'aln.fasta')
	args = path + ' -cpu '+ str(ncores) +' -d ' + db + ' -i ' + prot + ' -oa3m ' + outa3m + ' -ofas ' + outfasta  +' -o '+ outhhr + ' -n ' + str(iterations) + ' -maxmem 15 -b 20 -z 20 -B 2000 -Z 2000 -nodbfilter -mact .09'
	if SS == True:
		 args += ' -ssm 2 -ssw .5 '
	print args

	args = shlex.split(args)
	p = subprocess.call(args )
	return p , [outa3m,outhhr,outfasta]
"""
def runSalignPairwise(prot, pargs):
	# profile-profile alignment using salign

	log.level(1, 0, 1, 1, 1)
	env = environ()

	aln = alignment(env, file='mega_prune.faa', alignment_format='FASTA')

	aln.salign(rr_file='${LIB}/blosum62.sim.mat',
	           gap_penalties_1d=(-500, 0), output='',
	           align_block=15,   # no. of seqs. in first MSA
	           align_what='PROFILE',
	           alignment_type='PAIRWISE',
	           comparison_type='PSSM',  # or 'MAT' (Caution: Method NOT benchmarked
	                                    # for 'MAT')
	           similarity_flag=True,    # The score matrix is not rescaled
	           substitution=True,       # The BLOSUM62 substitution values are
	                                    # multiplied to the corr. coef.
	           #write_weights=True,
	           #output_weights_file='test.mtx', # optional, to write weight matrix
	           smooth_prof_weight=10.0) # For mixing data with priors

	#write out aligned profiles (MSA)
	aln.write(file='salign.ali', alignment_format='PIR')

	# Make a pairwise alignment of two sequences
	aln = alignment(env, file='salign.ali', alignment_format='PIR',
	                align_codes=('12asA', '1b8aA'))
	aln.write(file='salign_pair.ali', alignment_format='PIR')
	aln.write(file='salign_pair.pap', alignment_format='PAP')
	return """


def runAddSS(AddSSpath, prot):
	#perl addss.pl <ali file> [<outfile>] [-fas|-a3m|-clu|-sto]
	#make sure vars are set in config script as specified in the HHblits manual
	if 'a3m' in prot:
		outfile = prot.replace('.a3m' , 'SS.a3m' )
		args = AddSSpath +' ' +prot + ' ' + outfile+ ' -a3m'

	if 'fasta' in prot and 'a3m' not in prot:
		outfile = prot.replace('.fasta' , 'SS.fasta' )
		args = AddSSpath +' ' +prot + ' ' + outfile+ ' -fas'
	
	args = shlex.split(args)
	p= subprocess.call(args )
	return p, [outfile]

def runHHsearch( prot , pargs):
	path , outdir , db, ncores , ssw , mact =pargs
	#Usage: hhblits -i query [options]
	name = prot.split('/')[-1].replace('.hhm','')
	outhhr= outdir+name+".hhr"
	args = path + ' -cpu '+ str(ncores) +' -d ' + db + ' -i ' + prot + ' -o '+ outhhr + ' -local ' + ' -cpu '+ str(ncores) + ' -ssw ' + str(ssw) + ' -ssm 4 -e 1 -p 0 ' + ' -mact ' + str(mact)+ ' -realign'
	print args
	args = shlex.split(args)
	p = subprocess.call(args)
	return p , [outhhr]

def runHHmake(prot , pargs ):
	path, outdir = pargs
	if 'a3m' in prot:
		outhhm = prot.replace('.a3m' , '.hhm')
	if 'fasta' in prot:
		outhhm = prot.replace('.fasta' , '.hhm')
	args = path +' -i ' + prot + ' -o '+ outhhm 
	print args
	args = shlex.split(args)
	p = subprocess.call(args)
	return p , [outhhm]

def runFastme( fastmepath , clusterfile ):
	args =  fastmepath +  ' -i ' + clusterfile + ' -o ' + clusterfile+'_tree.txt' 
	print args
	p = subprocess.call(shlex.split(args) , stdout=subprocess.PIPE )
 	return p,[clusterfile+'_tree.txt' ]

def distmat_to_txt( pdblist , distmat, filedir , name):
	print pdblist
	print distmat.shape
	outstr = str(len(pdblist)) + '\n'
	for i,pdb in enumerate(pdblist):
		#namestr = pdb.replace('.','').replace('_','')[0:20]
		namestr = pdb[0:20]
		outstr += namestr+ ' ' + np.array2string( distmat[i,:], formatter={'float_kind':lambda x: "%.2f" % x} , precision = 8 ).replace('[', '').replace(']', '').replace('\n', '' )  + '\n'
	handle = open(filedir + name + 'fastmemat.txt' , 'w')
	handle.write(outstr)
	handle.close()
	return filedir + name + 'fastmemat.txt' , filedir + name + 'phylipmat.txt'


def HHSearch_parseTo_DMandNX(hhrs ):
	clusternames = []
	for i,hhr in enumerate(hhrs):
		profile = HHOutputParser(alignments=False).parse_file(hhr)
		if profile.query_name not in clusternames:
			clusternames.append(profile.query_name)
	evalDM = np.ones( (len(clusternames),len(clusternames) )) 
	pvalDM = np.ones( (len(clusternames),len(clusternames) )) 
	scoreDM = np.zeros( (len(clusternames),len(clusternames) ))
	SSDM = np.zeros( (len(clusternames),len(clusternames) ))
	probaDM = np.zeros( (len(clusternames),len(clusternames) ))
	lenDM =  np.ones( (len(clusternames),len(clusternames) ))
	NX = nx.Graph()
	for i,hhr in enumerate(hhrs):
		protlist = []
		profile = HHOutputParser(alignments=False).parse_file(hhr)
		for hit in profile:
			DMscore = float(hit.evalue)
			proba = hit.probability
			if 'anchor' not in hit.id and 'anchor' not in profile.query_name:
				i = clusternames.index(hit.id)
				j = clusternames.index(profile.query_name)
				
				if hit.evalue < evalDM[i,j]:
					evalDM[i,j] = hit.evalue
					evalDM[j,i] = evalDM[i,j]

				if hit.pvalue < pvalDM[i,j]:
					pvalDM[i,j] = hit.pvalue
					pvalDM[j,i] = pvalDM[i,j]

				if scoreDM[i,j] < hit.score:
					scoreDM[i,j] = hit.score
					scoreDM[j,i] = scoreDM[i,j]

				if SSDM[i,j] < hit.ss_score:
					SSDM[i,j] = hit.ss_score
					SSDM[j,i] = SSDM[i,j]


				if probaDM[i,j] < hit.probability:
					probaDM[i,j] = hit.probability
					probaDM[j,i] = probaDM[i,j]
				
				#use smallest of the two prots
				if lenDM[i,j] == 1 or lenDM[i,j] > hit.qlength:
					lenDM[i,j] = hit.qlength
					lenDM[j,i] = lenDM[i,j]

			if hit.id != profile.query_name :
				NX.add_edge( hit.id , profile.query_name )
				NX[hit.id][profile.query_name]['score']= hit.score
	return probaDM, evalDM ,pvalDM,  lenDM , scoreDM, SSDM, NX , clusternames


###############################################################algorithm start###################################################################
if prepare_models == True:
	fastas = glob.glob( sequences + '*.fasta')
	for path in [modeldir, structuralpreddir , resultsdir]:
		if os.path.exists(path) == False:
				os.mkdir(path)
	print fastas
	#make models
	models = {}
	fastaModels = {}
	HHblitsArgs = (hhblitspath, modeldir , UniprotHHDB, 3 , nCPU , 'uniprot20', True)
	for fasta in fastas:
		p,output = runHHblits( fasta, HHblitsArgs )
		models[fasta.split('/')[-1].replace('.fasta' , '' )]=output[0]
		fastaModels[fasta.split('/')[-1].replace('.fasta' , '' )] = output[2]
	save_obj(models, modeldir +'baremodels')
	save_obj(fastaModels, modeldir + 'fastaModels')
	
if cleanModels == True:
	baremodels = load_obj(modeldir + 'baremodels')
	print baremodels
	#correct names
	A3mModels = {}
	HMMmodels = {}
	HHmakeargs = (hhmakepath , modeldir)
	for fasta in baremodels:
		p,output = runAddSS(AddSSpath , baremodels[fasta])
		p,output2 = runHHmake(output[0] , HHmakeargs )
		HMMmodels[fasta] = output2[0]
		A3mModels[fasta] = output[0]

	finalmodels = {}
	#clean up names
	for SS in HMMmodels:
		outstr = ''
		with open(HMMmodels[SS] , 'r') as modelfile:
			for line in modelfile:
				already = False
				if 'NAME' in line:
					outstr += 'NAME' + '  ' + SS.split('/')[-1].replace('.fasta' , '' ) + '\n'
				else:
					outstr += line
		handle = open(HMMmodels[SS], 'w')
		handle.write(outstr)
		handle.close()

	save_obj(A3mModels, modeldir + 'A3mModels')
	save_obj(HMMmodels ,modeldir+ 'HMMmodels')

if check_structure == True:
	#check structs
	HMMmodels = load_obj(modeldir+ 'models')
	HHblitsArgs = (hhblitspath, structuralpreddir , PDbHHDB, 1 , nCPU , 'pdb70' , True)
	structmodels= {}
	for SS in HMMmodels:
		p,output  = runHHblits(HMMmodels[SS] , HHblitsArgs)
		structmodels[SS] = output[0]
	save_obj(structmodels,structuralpreddir+ 'structmodels')
	

def writeDB(models,modeldir,db):
	with open(db, 'w') as outdb:
		for model in models.values():
			outdb.write( model + '\n')
	return db


def runHHalign(path , prot1 , prot2 , outfile ):
	#align two hhms... 
	print prot1 
	print prot2
	args = path + ' -i '+ prot1   +' -t ' + prot2 +' -Aa3m '+ outfile + ' -o '+ outfile+'.hhr ' +  '-ssw .3 -ssm 4 -mact 0.0001 '
	print args
	args = shlex.split(args)
	p = subprocess.call(args)
	return p , [outfile]


def runClustalo(clustalopath , fasta1, fasta2, outfile):
	args =  clustalopath + ' --force --p1 ' + fasta1 + ' --p2 ' + fasta2 + ' -o ' + outfile
	print fasta1
	print fasta2
	print args
	args = shlex.split(args)
	p = subprocess.call(args )
	return p,[outfile]

def runReformat(path , infile , out ):
	print infile 
	print out
	args = path + ' a3m fas ' + infile + ' ' +  out 
	
	print args
	args = shlex.split(args)
	p = subprocess.call(args )
	return p,[out]

def reverseReformat(path , infile , out ):
	args = path + ' fas a3m ' + infile + ' ' +  out 
	print args
	args = shlex.split(args)
	p = subprocess.call(args )
	return p,[out]

def merge_alns_hhalign(hhalignpath,hhmakepath, reformatpath , tree, nodedict , probaMAT ,labels , cutoff , filedir = './' ):
	#one round of merging all leaves and returning merged models
	t = PhyloTree( tree, sp_naming_function=None)
	count = 0
	models = {}
	bareA3mModels = {}
	A3mModels = {}
	merged = []
	dist = 10**8
	

	for i, n in enumerate(t.traverse()):
		if len(n.get_leaves()) == 2:
			#check proba cutoff here
			children =  n.get_leaves()
			if probaMAT[labels.index(children[0].name), labels.index(children[1].name) ] > cutoff:
				merged += children
				#merge leaves and assign to parent node.
				file1 = tempfile.NamedTemporaryFile()
				file2 = tempfile.NamedTemporaryFile()
				#p,newaln = runClustalo(clustalopath, runReformat(reformatpath, nodedict[children[0].name], file1.name )[1][0] ,runReformat(reformatpath, nodedict[children[1].name] , file2.name)[1][0]  ,filedir + children[0].name+'_'+ children[1].name + '_merge.fasta' )
				p,newaln = runHHalign(hhalignpath,  nodedict[children[0].name] ,nodedict[children[1].name] , filedir + children[0].name+'_'+ children[1].name + '_merge.a3m' )
				n.name = children[0].name.upper()+'-'+ children[1].name.lower()
				count +=1
				#p,newA3m = reverseReformat( reformatpath , newaln[0] , newaln[0].replace('fasta', 'a3m') )
				bareA3mModels[n.name] = newaln[0]
				p,SS = runAddSS(AddSSpath, newaln[0] )
				A3mModels[n.name] = SS[0]
				models[n.name] = runHHmake( SS[0]  , [hhmakepath, filedir] )[1][0]
				file1.close()
				file2.close()
			
	for i, n in enumerate(t.traverse()):
			if n not in merged and n.name in nodedict:
				A3mModels[n.name] = runAddSS(AddSSpath, nodedict[n.name] ) [1][0]
				bareA3mModels[n.name] = nodedict[n.name]
				models[n.name]= runHHmake(A3mModels[n.name] , [hhmakepath, filedir] )[1][0]
	return models, bareA3mModels , A3mModels


def make_phylofunction(mergeiter, HMMmodels, A3mModels, BareA3mModels, modeldir, resultsdir, HHsearchpath, fastmepath):
	def phylofunction(ssw , mact):
		print 'Running Phylo'
		db = modeldir + 'hhmdb'+ str(mergeiter) +'.pal'
		db = writeDB(HMMmodels,modeldir,db)
		results = {}
		#allvall
		for SS in HMMmodels:
			pargs = (HHsearchpath , resultsdir , db, 30 , ssw , mact  )
			p,outfile = runHHsearch( HMMmodels[SS] , pargs)
			results[SS] = outfile[0]
		save_obj(results ,resultsdir+'results'+str(mergeiter))
		results = load_obj(resultsdir+'results'+str(mergeiter))
		#make distmat
		probaDM, evalDM , pvalDM, lenDM , scoreDM, SSDM, NX , labels = HHSearch_parseTo_DMandNX( results.values())
		#soding distance metric
		HHDM =  1 / ( (scoreDM +  3*SSDM) / lenDM )
		np.fill_diagonal(HHDM , 0 )
		HHDM[HHDM == np.inf ] = 10**5
		save_obj(HHDM, resultsdir+ 'HHDM'+str(mergeiter))
		#theobald distance metric
		Inner = -np.clip(np.log(evalDM) , np.amin(np.log(evalDM) ) , 0)  - .57722
		Inner = np.clip(Inner , 0, np.amax(Inner))
		TheobaldDM = -np.log(Inner)
		TheobaldDM[TheobaldDM == np.inf ] = 10**5
		TheobaldDM[TheobaldDM < 0 ] = 0
		np.fill_diagonal(TheobaldDM , 0 )
		save_obj(TheobaldDM , resultsdir+ 'TheobaldDM'+str(mergeiter))
		#output matrices and calculate Phylogenies
		fastme , phylip = distmat_to_txt( labels , probaDM , resultsdir , 'probaMAT'+str(mergeiter))
		fastme , phylip = distmat_to_txt( labels , np.log(evalDM) , resultsdir , 'evalMAT'+str(mergeiter))
		fastme , phylip = distmat_to_txt( labels , np.log(pvalDM) , resultsdir , 'pvalMAT'+str(mergeiter))
		fastme , phylip = distmat_to_txt( labels , HHDM , resultsdir , 'soddingDISTMAT'+str(mergeiter))
		p,treefile = runFastme(fastmepath,  fastme) 
		save_obj(treefile, resultsdir + 'sodingTree'+str(mergeiter))
		fastme , phylip = distmat_to_txt( labels , TheobaldDM , resultsdir , 'TheboaldDISTMAT'+str(mergeiter))
		p,treefile = runFastme(fastmepath,  fastme) 
		save_obj(treefile, resultsdir + 'theobaldTree'+str(mergeiter))
		save_obj( [probaDM, evalDM , pvalDM, lenDM , scoreDM, SSDM, NX , labels ] , resultsdir + 'DMpack'+str(mergeiter))

		#calculate the result of this phylo run
		print 'DONE'
		Y = np.mean(scoreDM)
		return Y
	return phylofunction


#make a phylogeny function that accepts alignment params and runs the phylo with the presnt configuration


mergers = 4
structmodels = load_obj(structuralpreddir+ 'structmodels')
HMMmodels = load_obj(modeldir+ 'HMMmodels')
A3mModels = load_obj(modeldir + 'A3mModels')
BareA3mModels = load_obj(modeldir + 'BareA3mModels')


for mergeiter in range(mergers):
	print 'merging iteration ' + str(mergeiter)
	#output hhsearch DB
	f = make_phylofunction(mergeiter, HMMmodels, A3mModels, BareA3mModels, modeldir, resultsdir, HHsearchpath, fastmepath)
	#optimize f for alignment parameters using bayesian optimization
	
	#ssw , mact
	bo = BayesianOptimization(f , {'ssw':(0,1) , 'mact':(0,.9999) })
	bo.maximize(init_points=5, n_iter=15, kappa=3.29)
	print(bo.res['max'])
	#run once with best config
	f(bo.res['max']['max_params']['ssw'], bo.res['max']['max_params']['mact'])

	#merge and rerun
	if merge_alns == True:
		[probaDM, evalDM , pvalDM, lenDM , scoreDM, SSDM, NX , labels ] = load_obj(resultsdir + 'DMpack'+str(mergeiter))
		sodingTree = load_obj(resultsdir+ 'sodingTree'+str(mergeiter))
		modeldir = modeldir[:-1]+str( mergeiter ) + '/'
		resultsdir = resultsdir[:-1] +str(mergeiter)+ '/'
		for path in [resultsdir, modeldir]:
			if os.path.exists(path) == False:
					os.mkdir(path)
		HMMmodels , BareA3mModels, A3mModels = merge_alns_hhalign(hhalignpath, hhmakepath , reformatpath , sodingTree[0], BareA3mModels , probaDM, labels , .7,  modeldir)		
		for SS in HMMmodels:
			outstr = ''
			with open(HMMmodels[SS] , 'r') as modelfile:
				for line in modelfile:
					already = False
					if 'NAME' in line:
						outstr += 'NAME' + '  '+SS+ '\n'
					else:
						outstr += line
			handle = open(HMMmodels[SS], 'w')
			handle.write(outstr)
			handle.close()
		save_obj(A3mModels, modeldir + 'A3mModels'+ str(mergeiter) )
		save_obj(BareA3mModels,modeldir + 'BareA3mModels'+str(mergeiter) )
		save_obj(HMMmodels,modeldir + 'HMMmodels'+str(mergeiter) )

		
