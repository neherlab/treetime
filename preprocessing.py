#!/usr/bin/env python
import os
import sys
import subprocess
import multiprocessing

assert len(sys.argv) == 3, "invalid number of arguments"

type = sys.argv[1] # data type
seq = sys.argv[2] # sequence file name
num_threads = multiprocessing.cpu_count()

# path configuration
path = 'data/'
curr_path = os.path.abspath(os.path.curdir)  # absolute current path
abs_path = os.path.join(curr_path, path).replace(' ', '\ ')  # absolute path for data
seq_path = os.path.join(path, seq)  # sequence file relative path
pasta = 'pasta/pasta-tools/run_pasta.py'  # MSA tool
raxml = 'RAxML/raxmlHPC-AVX'  # phylogenetic tree tool

# MSA
subprocess.call('python %s -i %s -d %s --merger=muscle --num-cpus=%d' % (pasta, seq_path, type, num_threads), shell=True)
aln = '{0}pastajob.marker001.{1}.aln'.format(path, seq[:-6])
aln_new = '{0}_aln.fasta'.format(seq_path[:-6])
subprocess.call('mv %s %s' % (aln, aln_new), shell=True)  # rename alignment file

# delete MSA tmp files
subprocess.call('rm -rf {0}pastajob*'.format(path), shell=True)
subprocess.call('rm -rf *.fasta.reduced', shell=True)

# build phylogenetic tree
outfile_suffix = seq[:-6]
cmd = '%s -m PROTGAMMAJTT -s %s -p 12435 -w %s -n %s' % (raxml, aln_new, abs_path, seq[:-6])
subprocess.call(cmd, shell=True)

# rename tree topology (nwk) & info file
tree = '{0}RAxML_bestTree.{1}'.format(path, seq[:-6])
tree_info = '{0}RAxML_info.{1}'.format(path, seq[:-6])
subprocess.call('mv %s %s' % (tree, path+'raxml_tree.nwk'), shell=True)
subprocess.call('mv %s %s' % (tree_info, path+'raxml_tree_info.txt'), shell=True)

# delete RAxML tmp files
subprocess.call('rm -rf {0}RAxML*'.format(path), shell=True)
