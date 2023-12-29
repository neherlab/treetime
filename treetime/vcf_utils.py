import gzip
import numpy as np
from collections import defaultdict
from textwrap import fill
from . import TreeTimeError

## Functions to read in and print out VCF files

def read_vcf(vcf_file, ref_file=None):
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
    ref_file : string, optional
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
    import re
    ALT_CHARS = re.compile(r"^([ACGTNacgtn]+|\*|\.)$") # straight from the VCF 4.3 spec

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
    #     First base _does not_ always match REF, so there may be a SNP as well
    # 'No indel'
    #     Ex:
    #       REF     ALT
    #       A       G

    #define here, so that all sub-functions can access them
    sequences = defaultdict(dict)
    insertions = defaultdict(dict) #Currently not used, but kept in case of future use.
    metadata = {
        'meta_lines': [], # all the VCF meta_lines, i.e. those starting with '##' (they are left unparsed)
        'chrom': None,    # chromosome name (we only allow one)
        'ploidy': None,   # ploidy count -- encoded in how the GT calls are formatted
    }

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
    def parse_homozygous_call(snps, ins, pos, ref, alt):

        # Replace missing allele with N(s). See commentary in test_vcf.py::TestNoCallsOrMissing for more details
        if alt=='*' or alt=='.':
            alt = "N" * len(ref)

        #Insertion where there are also deletions (special handling)
        if len(ref) > 1 and len(alt)>len(ref):
            ## NOTE: the loop below contains a potential double-counting bug. For example,
            ## REF='TCG' ALT='TCAG' (Example 5.1.3 in the VCF 4.2 spec), then when i=2
            ## we'll add both snps[pos+2] = 'A' as well as ins[pos+2] = 'AG'. This has been 
            ## detailed within `test_vcf.py`, as it may also be the expected way to encode
            ## insertions within TreeTime
            for i in range(len(ref)):
                #if the pos doesn't match, store in sequences
                if ref[i] != alt[i]:
                    snps[pos+i] = alt[i]
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
                        snps[pos+i] = alt[i]
        #Insertion + single-base ref
        elif len(alt) > 1:
            ins[pos] = alt
            # If the first base of the allele doesn't match the ref then we _also_ have a mutation
            if ref[0]!=alt[0]:
                snps[pos] = alt[0]
        #No indel
        else:
            snps[pos] = alt


    #Parses a 'bad' (hetero or no-call) call depending on what it is
    #TODO - consider the situation where the alternate allele base(s) is '*' (done for the homozygous case)
    def parse_heterozygous_call(gen, snps, ins, pos, ref, ALT):
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
                            snps[pos+i] = (alt[i] if alt[i] != '.' else 'N') #'.' = no-call

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

    def validate_alt(alt):
        """
        from the VCF 4.3 spec:
        > the ALT field must be a symbolic allele, or a breakend replacement string,
        > or match the regular expression [see ALT_CHARS regex, above]
        Since we only consider GT variation, the symbolic alleles and breakends aren't
        relevant for us.

        The spec also states:
        > Tools processing VCF files are not required to preserve case in the allele String
        
        Return the uppercase allele bases (string) _or_ None if it fails validation
        """
        if ALT_CHARS.match(alt):
            return alt.upper()
        else:
            print(f"WARNING: Encountered invalid allele base(s) {alt!r}. Skipping...")
            return None

    #House code is *much* faster than pyvcf because we don't care about all info
    #about coverage, quality, counts, etc, which pyvcf goes to effort to parse
    #(and it's not easy as there's no standard ordering). Custom code can completely
    #ignore all of this.
    from Bio import SeqIO

    #Use different openers depending on whether compressed
    opn = gzip.open if vcf_file.endswith(('.gz', '.GZ')) else open

    with opn(vcf_file, mode='rt') as f:
        current_block = "meta-information" # The start of VCF files is the meta-information (lines starting with ##)
        header,samps,nsamp=None,None,None

        for line in f:
            if line.startswith("##"):
                if current_block!='meta-information':
                    raise TreeTimeError(f"Malformed VCF file {vcf_file!r} - all the meta-information (lines starting with ##) must appear at the top of the file.")
                metadata['meta_lines'].append(line.strip())
            elif line[0]=='#':
                mandatory_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
                # Note that FORMAT isn't mandatory unless genotype data is present. But every VCF we deal with has genotype data - that's the point!
                header_start = "#"+"\t".join(mandatory_fields)
                if not line.startswith(header_start):
                    raise TreeTimeError(f"Malformed VCF file {vcf_file!r} - the header line must start with {header_start}")
                if current_block!='meta-information':
                    raise TreeTimeError(f"Malformed VCF file {vcf_file!r} - the header must immediately follow the meta-information lines")
                current_block = 'header'
                header = line.strip().split('\t')
                samps = [ x.strip() for x in header[9:] ] #ensure no leading/trailing spaces
                nsamp = len(samps)
                for sample_name in samps:
                    sequences[sample_name] = {}
            else:
                if current_block!='header' and current_block!='data-lines':
                    raise TreeTimeError(f"Malformed VCF file {vcf_file!r} - the data lines must follow the header line")
                current_block='data-lines'
                l = line.strip()
                if not l: # empty line
                    continue
                dat = l.split('\t')
                chrom = str(dat[0])
                if metadata['chrom'] is None:
                    metadata['chrom'] = chrom
                elif metadata['chrom']!=chrom:
                    raise TreeTimeError(f"The VCF file {vcf_file!r} contains multiple chromosomes. TreeTime can not yet handle this.")
                pos = int(dat[1])-1 # Convert VCF 1-based to python 0-based
                REF = dat[3]
                ALT = dat[4].split(',') # List of alternate alleles (strings)
                calls = np.array(dat[9:])
                if len(calls)!=nsamp:
                    raise TreeTimeError(f"Malformed VCF file {vcf_file!r} - the data lines have different number of calls than the number of samples defined in the header")

                # Treetime only parses GT FORMAT calls. Moreover, "the first sub-field must always be the genotype (GT) if it is present"
                # according to the VCF 4.2 spec. Note that if the VCF is for one sample only then the format's a bit different, but this'll
                # raise an error above because the header differs from what we assert. Other variation (CN, BND etc) is not parsed by
                # TreeTime, only GT.
                FORMAT = dat[8]
                if not FORMAT.startswith('GT'):
                    continue

                #get samples that differ from Ref at this site
                for sname, sa in zip(samps, calls):
                    gt = sa.split(':')[0] # May be multiple colon-separated subfields, depending on the line's FORMAT

                    # `gt` is the "index" of the alternate allele for this sample.
                    # The format of `gt` is quite varied:
                    # For haploid genomes, it's simply a single int _or_ "." (meaning a call cannot be made at the locus)
                    # For polyploid genomes, the format is a list of {int, '.'} separated by "/" or "|"
                    # ('/' means unphased, '|' means phased). 
                    # If the index is an int, it's the 1-based (!) lookup index for ALT

                    # Split on the valid separators - if haploid, then we'll get a list len=1
                    gts = gt.split('|') if '|' in gt else gt.split('/')
                    ploidy = len(gts)
                    if metadata['ploidy'] is None:
                        metadata['ploidy'] = ploidy
                    elif metadata['ploidy']!=ploidy:
                        raise TreeTimeError(f"The VCF file {vcf_file!r} had genotype calls of ploidy {metadata['ploidy']} but sample {sname!r} has a genotype of ploidy {ploidy}")

                    # if polyploid, but homozygous, then treat as if haploid
                    if ploidy>1 and len(set(gts))==1:
                        gt = gts[0]

                    if gt.isdigit(): # haploid, and a call has been made (i.e. it's not gt='.')
                        gt = int(gt)
                        if gt==0:
                            continue # reference allele!
                        alt = validate_alt(ALT[gt-1]) # gt is the 1-based lookup, but ALT is 0-indexed
                        if alt:
                            parse_homozygous_call(sequences[sname],insertions[sname], pos, REF, alt)
                        continue

                    if gt=='.': # haploid "call cannot be made" identifier - replace REF with N(s)
                        for i in range(len(REF)):
                            sequences[sname][pos+i] = 'N'
                        continue

                    # ---- heterozygous polyploid call  ----
                    parse_heterozygous_call(gt, sequences[sname],insertions[sname], pos, REF, ALT)

    #Gather all variable positions
    #NOTE: this does not consider positions of insertions
    positions = set()
    for seq, muts in sequences.items():
        positions.update(muts.keys())

    num_insertions = sum([len(list(ins.keys())) for ins in insertions.values()])
    if len(positions)==0 and num_insertions==0:
        raise TreeTimeError(f"VCF file {vcf_file!r} has no data-lines which we could extract genotype information from!")

    #One or more seqs are same as ref! (No non-ref calls) So haven't been 'seen' yet
    if nsamp > len(sequences):
        missings = set(samps).difference(sequences.keys())
        for s in missings:
            sequences[s] = {}

    if ref_file:
        refSeq = SeqIO.read(ref_file, format='fasta')
        refSeq = refSeq.upper() #convert to uppercase to avoid unknown chars later
        refSeqStr = str(refSeq.seq)
    else:
        refSeqStr = None

    compress_seq = {'reference':refSeqStr,
                    'sequences': sequences,
                    'insertions': insertions,
                    'positions': sorted(positions),
                    'metadata': metadata}

    return compress_seq


def write_vcf(tree_dict, file_name, mask=None):#, compress=False):
    """
    Writes out a VCF-style file (which seems to be minimally handleable
    by vcftools and pyvcf) of the alignment. This is created from a dict
    in a similar format to what's created by :py:meth:`treetime.vcf_utils.read_vcf`

    Positions of variable sites are transformed to start at 1 to match
    VCF convention.

    Parameters
    ----------
     tree_dict: nested dict
        A nested dict with required keys of:
        'sequences': maps sampleName -> pos (0-based) -> new base (SNP)
        'reference': string of reference nuc sequence
        'positions': sorted list of 0-based positions with variation in sequences
        And optional keys:
        'inferred_const_sites': list or set, 0-based positions to skip output for.
        This input is often created by :py:meth:`treetime.TreeAnc.get_tree_dict`
        'metadata': dict of information to influence VCF formatting. Only
        the following keys are used:
        'metadata.ploidy': int. Influences how genotype calls are formatted. (default
        of 2 (diploid) used if not provided)
        'metadata.chrom': str. The chromosome name (default of '1' used if not provided)

     file_name: str
        File to which the new VCF should be written out. File names ending with
        '.gz' will result in the VCF automatically being gzipped.

    mask : optional, str of 0 or 1
        Calls at these positions will be skipped
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
    ploidy = tree_dict.get('metadata', {}).get('ploidy', 2)
    chrom_name = tree_dict.get('metadata', {}).get('chrom', '1')
    sample_names = list(sequences.keys())
    inferred_const_sites = set(tree_dict.get('inferred_const_sites', []))

    # For every variable site in sequences, flip the format around so
    # we can have fast lookups later on.
    alleles = {}
    num_samples = len(sample_names)
    for idx, name in enumerate(sample_names):
        for posn, allele in sequences[name].items():
            if posn not in alleles:
                alleles[posn] = np.zeros(num_samples, dtype='U')
            alleles[posn][idx] = allele
    # fill in reference
    for posn,bases in alleles.items():
        bases[bases==''] = ref[posn]

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
            pattern = np.array(pattern).astype('U')

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
            pattern = np.array(pattern).astype('U')

            #Stops 'greedy' behaviour from adding mutations adjacent to deletions
            if any(pattern == '-'): #if part of deletion, append
                sites.append(pattern)
                refb = refb+ref[pi]
            else: #this is another mutation next to the deletion!
                i-=1    #don't append, break this loop

        #Rotate them into 'calls'
        align = np.asarray(sites).T

        #Get rid of '-', and put '0' for calls that match ref
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
                fullpat.append('0')
            else:
                fullpat.append(pat)

        pattern = np.array(fullpat)

        return i, pi, pos, refb, pattern


    #prepare the header of the VCF & write out
    header=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+sample_names

    opn = gzip.open if file_name.endswith(('.gz', '.GZ')) else open
    out_file = opn(file_name, 'w')

    out_file.write( "##fileformat=VCFv4.2\n"+
                        "##source=NextStrain\n"+
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"+
                        f"##contig=<ID={chrom_name}>\n")
    out_file.write("\t".join(header)+"\n")

    vcfWrite = []
    errorPositions = []
    explainedErrors = 0
    mask_skip_count = 0

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

        if mask and mask[pi] == '1':
            mask_skip_count+=1
            i+=1
            continue

        # patterns will be empty if there was no variation in the sequences
        pattern = alleles[pi] if pi in alleles else np.array([]).astype('U')
        pattern2 = alleles[pi+1] if pi+1 in alleles else np.array([]).astype('U')

        #If a deletion here, need to gather up all bases, and position before
        if any(pattern == '-'):
            if pos != 1:
                deleteGroup = True
                delete = True
            else:
                #If theres a deletion in 1st pos, VCF files do not handle this well.
                #Proceed keeping it as '-' for alt (violates VCF), but warn user to check output.
                #(This is rare)
                print(fill("WARNING: You have a deletion in the first position of your"
                           " alignment. VCF format does not handle this well. Please check"
                           " the output to ensure it is correct."))
        else:
            #If a deletion in next pos, need to gather up all bases
            if any(pattern2 == '-'):
                deleteGroup = True

        #If deletion, treat affected bases as 1 'call':
        if delete or deleteGroup:
            if pattern.size==0: # no variation in sequences
                pattern = np.full(num_samples, refb, dtype='U')
            i, pi, pos, refb, pattern = handleDeletions(i, pi, pos, ref, delete, pattern)
        #If no deletion, replace ref with '0' which means the reference base is unchanged
        else:
            pattern[pattern==refb] = '0'

        #Get the list of ALTs - minus any '0' which are unchanged reference sequences!
        uniques = np.unique(pattern)
        uniques = uniques[np.where(uniques!='0')]

        #Convert bases to the number that matches the ALT
        j=1
        for u in uniques:
            pattern[np.where(pattern==u)[0]] = str(j)
            j+=1

        #What if there's no variation at a variable site??
        #This can happen when sites are modified by TreeTime - see below.
        #We don't print to VCF (because no variation!)
        any_variation = len(uniques)!=0
        if not any_variation:
            #If we expect it (it was made constant by TreeTime), it's fine.
            if pi in inferred_const_sites:
                explainedErrors += 1
            else:
                #If we don't expect, raise an error
                errorPositions.append(str(pi))

        #Write it out - Increment positions by 1 so it's in VCF numbering
        #If no longer variable, and explained, don't write it out
        if any_variation:
            calls = [ "/".join([j]*ploidy) for j in pattern ]
            output = [chrom_name, str(pos), ".", refb, ",".join(uniques), ".", "PASS", ".", "GT"] + calls
            vcfWrite.append("\t".join(output))

        i+=1

    #Note: The number of 'inferred_const_sites' passed back by TreeTime will often be longer
    #than the number of 'site that were made constant' that prints below. This is because given the site:
    # Ref   Alt     Seq
    # G     A       AANAA
    #This will be converted to 'AAAAA' and listed as an 'inferred_const_sites'. However, for VCF
    #purposes, because the site is 'variant' against the ref, it is variant, as expected, and so
    #won't be counted in the below list, which is only sites removed from the VCF.

    if mask_skip_count:
        print(f"{mask_skip_count} positions were skipped due to the provided mask")

    if 'inferred_const_sites' in tree_dict and explainedErrors != 0:
        print(fill("Sites that were constant except for ambiguous bases were made" +
                   " constant by TreeTime. This happened {} times. These sites are".format(explainedErrors) +
                   " now excluded from the VCF."))

    if len(errorPositions) != 0:
        print ("\n***WARNING: vcf_utils.py")
        print(fill("\n{} sites were found that had no alternative bases.".format(str(len(errorPositions)))+
                  " If this data has been run through TreeTime and contains ambiguous bases,"
                  " try calling get_tree_dict with var_ambigs=True to see if this clears the error."))
        print(fill("\nAlternative causes:"
                   "\n- Not all sequences in your alignment are in the tree"
                  " (if you are running TreeTime via commandline this is most likely)"
                  "\n- In TreeTime, can be caused by overwriting variants in tips with small branch lengths (debug)"
                  "\n\nThese are the positions affected (numbering starts at 0):"))
        print(fill(", ".join(errorPositions)))

    out_file.write("\n".join(vcfWrite))
    out_file.close()


def process_sparse_alignment(aln, ref, ambiguous_char):
    return process_alignment_dictionary(aln, ref, ambiguous_char)

def process_alignment_dictionary(aln, ref, ambiguous_char):
    """
    prepare the dictionary specifying differences from a reference sequence
    to construct the reduced alignment with variable sites only. NOTE:
        - sites can be constant but different from the reference
        - sites can be constant plus a ambiguous sites

    assigns
    -------
    - self.nonref_positions: at least one sequence is different from ref

    Returns
    -------
    reduced_alignment_const
        reduced alignment accounting for non-variable postitions

    alignment_patterns_const
        dict pattern -> (pos in reduced alignment, list of pos in full alignment)

    variable_positions
        list of variable positions needed to construct remaining

    """

    # number of sequences in alignment
    nseq = len(aln)

    inv_map = defaultdict(list)
    for k,v in aln.items():
        for pos, bs in v.items():
            inv_map[pos].append(bs)

    nonref_positions = np.sort(list(inv_map.keys()))
    constant_up_to_ambiguous = []

    nonref_const = []
    nonref_alleles = []
    ambiguous_const = []
    variable_pos = []
    for pos, bs in inv_map.items(): #loop over positions and patterns
        bases = list(np.unique(bs))
        if len(bs) == nseq: #every sequence is different from reference
            if (len(bases)<=2 and ambiguous_char in bases) or len(bases)==1:
                # all sequences different from reference, but only one state
                # (other than ambiguous_char) in column
                nonref_const.append(pos)
                if len(bases)==1:
                    nonref_alleles.append(bases[0])
                else:
                    nonref_alleles.append([x for x in bases if x!=ambiguous_char][0])

                if ambiguous_char in bases: #keep track of sites 'made constant'
                    constant_up_to_ambiguous.append(pos)
            else:
                # at least two non-reference alleles
                variable_pos.append(pos)
        else: # not every sequence different from reference
            if len(bases)==1 and bases[0]==ambiguous_char:
                ambiguous_const.append(pos)
                constant_up_to_ambiguous.append(pos) #keep track of sites 'made constant'
            else:
                # at least one non ambiguous non-reference allele not in
                # every sequence
                variable_pos.append(pos)

    refMod = np.copy(ref)
    # place constant non reference positions by their respective allele
    refMod[nonref_const] = nonref_alleles
    # mask variable positions
    states = np.unique(refMod)
    refMod[variable_pos] = '.'

    # for each base in the gtr, make constant alignment pattern and
    # assign it to all const positions in the modified reference sequence
    constant_columns = []
    constant_patterns = {}
    for base in states:
        if base==ambiguous_char:
            continue
        p = np.repeat(base, nseq)
        pos = list(np.where(refMod==base)[0])
        #if the alignment doesn't have a const site of this base, don't add! (ex: no '----' site!)
        if len(pos):
            constant_patterns["".join(p.astype('U'))] = [len(constant_columns), pos]
            constant_columns.append(p)

    return {"constant_columns": constant_columns,
            "constant_patterns": constant_patterns,
            "variable_positions": variable_pos,
            "nonref_positions": nonref_positions,
            "constant_up_to_ambiguous": constant_up_to_ambiguous}


