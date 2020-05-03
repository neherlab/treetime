import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from treetime import TreeAnc
import treetime

def calc_mutation_rate(t):
    # estimate mutation rate
    node_order = []
    for node in t.tree.get_terminals():
        node_order.append(node)
    
    inconsistent_count = 0
    count = 0
    
    for i in range(len(node_order) - 1):
        for j in range(i + 1, len(node_order)):
            u = node_order[i]
            v = node_order[j]
            for k in range(len(u.sequence)):
                if u.sequence[k] != v.sequence[k]:
                    inconsistent_count += 1
            count += len(u.sequence)
    
    return inconsistent_count / count
    

if __name__ == '__main__':
    path = '../data/'
    #path = '../data/fp_seqs/'
    #path = '../data/treetime_data/h3n2_na/'
    
    
  
    # instantiate treetime
    #myTree = TreeAnc(gtr='jtt92', tree=path+'raxml_tree.nwk', aln=path+'h3n2_aa_aln.fasta', verbose=0)
    #myTree = TreeAnc(gtr='jtt92', tree=path+'phyml_tree.nwk', aln=path+'test_phyml.fasta', mu=5.0, verbose=0)
    
    #myTree = TreeAnc(gtr='t92', tree=path+'h3n2_na_500.nwk',
    #                  aln=path+'h3n2_na_500.fasta', verbose=0)
    #myTree.infer_ancestral_sequences(infer_gtr=False, marginal=True)
    #print(myTree.tree.sequence_LH.mean())

      
    # test ASRV
    #myTree2 = TreeAnc(gtr='jtt92', tree=path+'phyml_tree.nwk', aln=path+'test_phyml.fasta',
    #                  alpha=0.299, verbose=0)
    
    
    #myTree2.infer_ancestral_sequences(infer_gtr=False, marginal=True)
    #print(myTree2.rates)
    #print(myTree2.tree.sequence_LH.mean())
    
    # test ASRV + Mixture model
    pseudo_solvent_accessibility = [np.random.choice(['B', 'E']) for i in range(469)]
    myTree3= TreeAnc(gtr="EX", tree=path+'raxml_tree.nwk', aln=path+'h3n2_aa_aln.fasta',
                    alpha=0.299, struct_propty=pseudo_solvent_accessibility, verbose=0)
    
    myTree3.infer_ancestral_sequences(infer_gtr=False, marginal=True)
    print(myTree3.tree.sequence_LH.mean())
    print(myTree3.data.multiplicity)