import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from treetime import TreeAnc
import treetime

if __name__ == '__main__':
    path = '/Users/yinuo/Desktop/Spring 2020/Computational Genomics/Final project/treetime-ASR/data/'
    #path = '/Users/yinuo/Desktop/Spring 2020/Computational Genomics/Final project/treetime-ASR/data/fp_seqs/'
    
    # result_path = '../results/'
    
    # instantiate treetime
    myTree = TreeAnc(gtr='jtt92', tree=path+'raxml_tree.nwk', aln=path+'fp_seq_aa_aln.fasta', verbose=0)
    
    # test ASRV
    myTree2 = TreeAnc(gtr='jtt92', tree=path+'raxml_tree.nwk', aln=path+'fp_seq_aa_aln.fasta',
                     alpha=0.493483, verbose=0)
  
    
    # test ASRV + Mixture model
    #myTree = TreeAnc(gtr="EX", tree=path+'raxml_tree.nwk', aln=path+'h3n2_aa_aln.fasta',
    #                 asrv=True, alpha_file='raxml_tree_info.txt', verbose=0)
    
    
    myTree.infer_ancestral_sequences(infer_gtr=False, marginal=True)
    myTree2.infer_ancestral_sequences(infer_gtr=False, marginal=True)
    
    
    print(sum(myTree.ancestral_likelihood()))
    print('--------')
    print(sum(myTree2.ancestral_likelihood()))
    
    
    #print(myTree.tree.root.sequence.shape[0])
    
   