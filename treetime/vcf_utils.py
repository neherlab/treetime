import numpy as np
from collections import defaultdict

## Functions to read in and print out VCF files

def read_vcf(vcf_file, ref_file):
    """
    Reads in a vcf/vcf.gz file and associated
    reference sequence fasta (to which the VCF file is mapped).

    Parses mutations, insertions, and deletions and stores them in a nested dict,
    see 'returns' for the dict structure.

    Calls with heterozygous values 0/1, 0/2, etc and no-calls (./.) are
    replaced with Ns at the associated sites.

    Positions are stored to correspond the location in the reference sequence
    in Python (numbering is transformed to start at 0)

    Parameters
    ----------
    vcf_file : string
        Path to the vcf or vcf.gz file to be read in
    ref_file : string
        Path to the fasta reference file to be read in

    Returns
    --------
    compress_seq : nested dict
        In the format: ::

            {
            'reference':'AGCTCGA..A',
            'sequences': { 'seq1':{4:'A', 7:'-'}, 'seq2':{100:'C'} },
            'insertions': { 'seq1':{4:'ATT'}, 'seq3':{1:'TT', 10:'CAG'} },
            'positions': [1,4,7,10,100...]
            }

        references : string
            String of the reference sequence read from the Fasta, to which
            the variable sites are mapped
        sequences : nested dict
            Dict containing sequence names as keys which map to dicts
            that have position as key and the single-base mutation (or deletion)
            as values
        insertions : nested dict
            Dict in the same format as the above, which stores insertions and their
            locations. The first base of the insertion is the same as whatever is
            currently in that position (Ref if no mutation, mutation in 'sequences'
            otherwise), so the current base can be directly replaced by the bases held here.
        positions : list
            Python list of all positions with a mutation, insertion, or deletion.

    """

    #Programming Note:
    # Note on VCF Format
    # -------------------
    # 'Insertion where there are also deletions' (special handling)
    #     Ex:
    #       REF     ALT         Seq1    Seq2
    #       GC      GCC,G       1/1     2/2
    #     Insertions formatted differently - don't know how many bp match
    #     the Ref (unlike simple insert below). Could be mutations, also.
    # 'Deletion'
    #     Ex:
    #       REF     ALT
    #       GC      G
    #     Alt does not have to be 1 bp - any length shorter than Ref.
    # 'Insertion'
    #     Ex:
    #       REF     ALT
    #       A       ATT
    #     First base always matches Ref.
    # 'No indel'
    #     Ex:
    #       REF     ALT
    #       A       G

    #define here, so that all sub-functions can access them
    sequences = defaultdict(dict)
    insertions = defaultdict(dict) #Currently not used, but kept in case of future use.

    #TreeTime handles 2-3 base ambig codes, this will allow that.
    def getAmbigCode(bp1, bp2, bp3=""):
        bps = [bp1,bp2,bp3]
        bps.sort()
        key = "".join(bps)

        return {
            'CT': 'Y',
            'AG': 'R',
            'AT': 'W',
            'CG': 'S',
            'GT': 'K',
            'AC': 'M',
            'AGT': 'D',
            'ACG': 'V',
            'ACT': 'H',
            'CGT': 'B'
        }[key]

    #Parses a 'normal' (not hetero or no-call) call depending if insertion+deletion, insertion,
    #deletion, or single bp subsitution
    def parseCall(snps, ins, pos, ref, alt):

        #Insertion where there are also deletions (special handling)
        if len(ref) > 1 and len(alt)>len(ref):
            for i in range(len(ref)):
                #if the pos doesn't match, store in sequences
                if ref[i] != alt[i]:
                    snps[pos+i] = alt[i] if alt[i] != '.' else 'N' #'.' = no-call
                #if about to run out of ref, store rest:
                if (i+1) >= len(ref):
                    ins[pos+i] = alt[i:]
        #Deletion
        elif len(ref) > 1:
            for i in range(len(ref)):
                #if ref is longer than alt, these are deletion positions
                if i+1 > len(alt):
                    snps[pos+i] = '-'
                #if not, there may be mutations
                else:
                    if ref[i] != alt[i]:
                        snps[pos+i] = alt[i] if alt[i] != '.' else 'N' #'.' = no-call
        #Insertion
        elif len(alt) > 1:
            ins[pos] = alt
        #No indel
        else:
            snps[pos] = alt


    #Parses a 'bad' (hetero or no-call) call depending on what it is
    def parseBadCall(snps, ins, pos, ref, ALT):
        #Deletion
        #   REF     ALT     Seq1    Seq2    Seq3
        #   GCC     G       1/1     0/1     ./.
        # Seq1 (processed by parseCall, above) will become 'G--'
        # Seq2 will become 'GNN'
        # Seq3 will become 'GNN'
        if len(ref) > 1:
            #Deleted part becomes Ns
            if gen[0] == '0' or gen[0] == '.':
                if gen[0] == '0':   #if het, get first bp
                    alt = str(ALT[int(gen[2])-1])
                else: #if no-call, there is no alt, so just put Ns after 1st ref base
                    alt = ref[0]
                for i in range(len(ref)):
                    #if ref is longer than alt, these are deletion positions
                    if i+1 > len(alt):
                        snps[pos+i] = 'N'
                    #if not, there may be mutations
                    else:
                        if ref[i] != alt[i]:
                            snps[pos+i] = alt[i] if alt[i] != '.' else 'N' #'.' = no-call

        #If not deletion, need to know call type
        #if het, see if proposed alt is 1bp mutation
        elif gen[0] == '0':
            alt = str(ALT[int(gen[2])-1])
            if len(alt)==1:
                #alt = getAmbigCode(ref,alt) #if want to allow ambig
                alt = 'N' #if you want to disregard ambig
                snps[pos] = alt
            #else a het-call insertion, so ignore.

        #else it's a no-call; see if all alts have a length of 1
        #(meaning a simple 1bp mutation)
        elif len(ALT)==len("".join(ALT)):
            alt = 'N'
            snps[pos] = alt
        #else a no-call insertion, so ignore.


    #House code is *much* faster than pyvcf because we don't care about all info
    #about coverage, quality, counts, etc, which pyvcf goes to effort to parse
    #(and it's not easy as there's no standard ordering). Custom code can completely
    #ignore all of this.
    import gzip
    from Bio import SeqIO
    import numpy as np

    nsamp = 0
    posLoc = 0
    refLoc = 0
    altLoc = 0
    sampLoc = 9

    #Use different openers depending on whether compressed
    opn = gzip.open if vcf_file.endswith(('.gz', '.GZ')) else open

    with opn(vcf_file, mode='rt') as f:
        for line in f:
            if line[0] != '#':
                #actual data - most common so first in 'if-list'!
                line = line.strip()
                dat = line.split('\t')
                POS = int(dat[posLoc])
                REF = dat[refLoc]
                ALT = dat[altLoc].split(',')
                calls = np.array(dat[sampLoc:])

                #get samples that differ from Ref at this site
                recCalls = {}
                for sname, sa in zip(samps, calls):
                    if ':' in sa: #if proper VCF file (followed by quality/coverage info)
                        gt = sa.split(':')[0]
                    else: #if 'pseudo' VCF file (nextstrain output, or otherwise stripped)
                        gt = sa
                    if gt == '0' or gt == '1': #for haploid calls in VCF
                        gt = '0/0' if gt == '0' else '1/1'

                    #ignore if ref call: '.' or '0/0', depending on VCF
                    if ('/' in gt and gt != '0/0') or ('|' in gt and gt != '0|0'):
                        recCalls[sname] = gt

                #store the position and the alt
                for seq, gen in recCalls.items():
                    ref = REF
                    pos = POS-1     #VCF numbering starts from 1, but Reference seq numbering
                                    #will be from 0 because it's python!
                    #Accepts only calls that are 1/1, 2/2 etc. Rejects hets and no-calls
                    if gen[0] != '0' and gen[2] != '0' and gen[0] != '.' and gen[2] != '.':
                        alt = str(ALT[int(gen[0])-1])   #get the index of the alternate
                        if seq not in sequences.keys():
                            sequences[seq] = {}

                        parseCall(sequences[seq],insertions[seq], pos, ref, alt)

                    #If is heterozygote call (0/1) or no call (./.)
                    else:
                        #alt will differ here depending on het or no-call, must pass original
                        parseBadCall(sequences[seq],insertions[seq], pos, ref, ALT)

            elif line[0] == '#' and line[1] == 'C':
                #header line, get all the information
                header = line.strip().split('\t')
                posLoc = header.index("POS")
                refLoc = header.index('REF')
                altLoc = header.index('ALT')
                sampLoc = header.index('FORMAT')+1
                samps = header[sampLoc:]
                samps = [ x.strip() for x in samps ] #ensure no leading/trailing spaces
                nsamp = len(samps)

            #else you are a comment line, ignore.

    #Gather all variable positions
    positions = set()
    for seq, muts in sequences.items():
        positions.update(muts.keys())

    #One or more seqs are same as ref! (No non-ref calls) So haven't been 'seen' yet
    if nsamp > len(sequences):
        missings = set(samps).difference(sequences.keys())
        for s in missings:
            sequences[s] = {}

    refSeq = SeqIO.read(ref_file, format='fasta')
    refSeq = refSeq.upper() #convert to uppercase to avoid unknown chars later
    refSeqStr = str(refSeq.seq)

    compress_seq = {'reference':refSeqStr,
                    'sequences': sequences,
                    'insertions': insertions,
                    'positions': sorted(positions)}

    return compress_seq


def write_vcf(tree_dict, file_name):#, compress=False):
    """
    Writes out a VCF-style file (which seems to be minimally handleable
    by vcftools and pyvcf) of the alignment. This is created from a dict
    in a similar format to what's created by :py:meth:`treetime.vcf_utils.read_vcf`

    Positions of variable sites are transformed to start at 1 to match
    VCF convention.

    Parameters
    ----------
     tree_dict: nested dict
        A nested dict with keys 'sequence' 'reference' and 'positions',
        as is created by :py:meth:`treetime.TreeAnc.get_tree_dict`

     file_name: str
        File to which the new VCF should be written out. File names ending with
        '.gz' will result in the VCF automatically being gzipped.

    """

#   Programming Logic Note:
#
#    For a sequence like:
#    Pos     1 2 3 4 5 6
#    Ref     A C T T A C
#    Seq1    A C - - - G
#
#    In a dict it is stored:
#    Seq1:{3:'-', 4:'-', 5:'-', 6:'G'}  (Numbering from 1 for simplicity)
#
#    In a VCF it needs to be:
#    POS REF     ALT     Seq1
#    2   CTTA    C       1/1
#    6   C       G       1/1
#
#    If a position is deleted (pos 3), need to get invariable position preceeding it
#
#    However, in alternative case, the base before a deletion is mutant, so need to check
#        that next position isn't a deletion (as otherwise won't be found until after the
#        current single bp mutation is written out)
#
#    When deleted position found, need to gather up all adjacent mutant positions with deletions,
#        but not include adjacent mutant positions that aren't deletions (pos 6)
#
#    Don't run off the 'end' of the position list if deletion is the last thing to be included
#        in the VCF file

    sequences = tree_dict['sequences']
    ref = tree_dict['reference']
    positions = tree_dict['positions']

    def handleDeletions(i, pi, pos, ref, delete, pattern):
        refb = ref[pi]
        if delete: #Need to get the position before
            i-=1    #As we'll next go to this position again
            pi-=1
            pos = pi+1
            refb = ref[pi]
            #re-get pattern
            pattern = []
            for k,v in sequences.items():
                try:
                    pattern.append(sequences[k][pi])
                except KeyError:
                    pattern.append(ref[pi])
            pattern = np.array(pattern)

        sites = []
        sites.append(pattern)

        #Gather all positions affected by deletion - but don't run off end of position list
        while (i+1) < len(positions) and positions[i+1] == pi+1:
            i+=1
            pi = positions[i]
            pattern = []
            for k,v in sequences.items():
                try:
                    pattern.append(sequences[k][pi])
                except KeyError:
                    pattern.append(ref[pi])
            pattern = np.array(pattern)

            #Stops 'greedy' behaviour from adding mutations adjacent to deletions
            if any(pattern == '-'): #if part of deletion, append
                sites.append(pattern)
                refb = refb+ref[pi]
            else: #this is another mutation next to the deletion!
                i-=1    #don't append, break this loop

        #Rotate them into 'calls'
        sites = np.asarray(sites)
        align = np.rot90(sites)
        align = np.flipud(align)

        #Get rid of '-', and put '.' for calls that match ref
        #Only removes trailing '-'. This breaks VCF convension, but the standard
        #VCF way of handling this* is really complicated, and the situation is rare.
        #(*deletions and mutations at the same locations)
        fullpat = []
        for pt in align:
            gp = len(pt)-1
            while pt[gp] == '-':
                pt[gp] = ''
                gp-=1
            pat = "".join(pt)
            if pat == refb:
                fullpat.append('.')
            else:
                fullpat.append(pat)

        pattern = np.array(fullpat)

        return i, pi, pos, refb, pattern


    #prepare the header of the VCF & write out
    header=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+list(sequences.keys())
    with open(file_name, 'w') as the_file:
        the_file.write( "##fileformat=VCFv4.2\n"+
                        "##source=NextStrain\n"+
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        the_file.write("\t".join(header)+"\n")

    vcfWrite = []
    errorPositions = []
    explainedErrors = 0

    #Why so basic? Because we sometimes have to back up a position!
    i=0
    while i < len(positions):
        #Get the 'pattern' of all calls at this position.
        #Look out specifically for current (this pos) or upcoming (next pos) deletions
        #But also distinguish these two, as handled differently.

        pi = positions[i]
        pos = pi+1 #change numbering to match VCF, not python, for output
        refb = ref[pi] #reference base at this position
        delete = False #deletion at this position - need to grab previous base (invariable)
        deleteGroup = False #deletion at next position (mutation at this pos) - do not need to get prev base

        #try/except is much more efficient than 'if' statements for constructing patterns,
        #as on average a 'variable' location will not be variable for any given sequence
        pattern = []
        #pattern2 gets the pattern at next position to check for upcoming deletions
        #it's more efficient to get both here rather than loop through sequences twice!
        pattern2 = []
        for k,v in sequences.items():
            try:
                pattern.append(sequences[k][pi])
            except KeyError:
                pattern.append(ref[pi])

            try:
                pattern2.append(sequences[k][pi+1])
            except KeyError:
                pattern2.append(ref[pi+1])

        pattern = np.array(pattern)
        pattern2 = np.array(pattern2)

        #If a deletion here, need to gather up all bases, and position before
        if any(pattern == '-'):
            if pos != 1:
                deleteGroup = True
                delete = True
            else:
                #If theres a deletion in 1st pos, VCF files do not handle this well.
                #Proceed keeping it as '-' for alt (violates VCF), but warn user to check output.
                #(This is rare)
                print ("WARNING: You have a deletion in the first position of your alignment. VCF format does not handle this well. Please check the output to ensure it is correct.")
        else:
            #If a deletion in next pos, need to gather up all bases
            if any(pattern2 == '-'):
                deleteGroup = True

        #If deletion, treat affected bases as 1 'call':
        if delete or deleteGroup:
            i, pi, pos, refb, pattern = handleDeletions(i, pi, pos, ref, delete, pattern)
        #If no deletion, replace ref with '.', as in VCF format
        else:
            pattern[pattern==refb] = '.'

        #Get the list of ALTs - minus any '.'!
        uniques = np.unique(pattern)
        uniques = uniques[np.where(uniques!='.')]

        #Convert bases to the number that matches the ALT
        j=1
        for u in uniques:
            pattern[np.where(pattern==u)[0]] = str(j)
            j+=1
        #Now convert these calls to #/# (VCF format)
        calls = [ j+"/"+j if j!='.' else '.' for j in pattern ]

        #What if there's no variation at a variable site??
        #This can happen when sites are modified by TreeTime - see below.
        printPos = True
        if len(uniques)==0:
            #If we expect it (it was made constant by TreeTime), it's fine.
            if 'inferred_const_sites' in tree_dict and pi in tree_dict['inferred_const_sites']:
                explainedErrors += 1
                printPos = False #and don't output position to the VCF
            else:
                #If we don't expect, raise an error
                errorPositions.append(str(pi))

        #Write it out - Increment positions by 1 so it's in VCF numbering
        #If no longer variable, and explained, don't write it out
        if printPos:
            output = ["MTB_anc", str(pos), ".", refb, ",".join(uniques), ".", "PASS", ".", "GT"] + calls
            vcfWrite.append("\t".join(output))

        i+=1

    #Note: The number of 'inferred_const_sites' passed back by TreeTime will often be longer
    #than the number of 'site that were made constant' that prints below. This is because given the site:
    # Ref   Alt     Seq
    # G     A       AANAA
    #This will be converted to 'AAAAA' and listed as an 'inferred_const_sites'. However, for VCF
    #purposes, because the site is 'variant' against the ref, it is variant, as expected, and so
    #won't be counted in the below list, which is only sites removed from the VCF.

    if 'inferred_const_sites' in tree_dict and explainedErrors != 0:
        print ( "Sites that were constant except for ambiguous bases were made constant by TreeTime. This happened {} times. These sites are now excluded from the VCF.".format(explainedErrors))

    if len(errorPositions) != 0:
        print ("\n***WARNING: vcf_utils.py"
            "\n{} sites were found that had no alternative bases. If this data has been "
            "run through TreeTime and contains ambiguous bases, try calling get_tree_dict with "
            "var_ambigs=True to see if this clears the error."
            "\n\nAlternative causes:"
            "\n- Not all sequences in your alignment are in the tree (if you are running TreeTime via commandline "
            "this is most likely)"
            "\n- In TreeTime, can be caused by overwriting variants in tips with small branch lengths (debug)"
            "\n\nThese are the positions affected (numbering starts at 0):".format(str(len(errorPositions))))
        print (",".join(errorPositions))

    with open(file_name, 'a') as the_file:
        the_file.write("\n".join(vcfWrite))

    if file_name.endswith(('.gz', '.GZ')):
        import os
        #must temporarily remove .gz ending, or gzip won't zip it!
        os.rename(file_name, file_name[:-3])
        call = ["gzip", file_name[:-3]]
        os.system(" ".join(call))
