import pytest
from textwrap import dedent
from treetime.vcf_utils import read_vcf, write_vcf
from treetime import TreeTimeError

def create_vcf(dir, content):
    d = dir / 'data'
    d.mkdir()
    fname = d / "no-header.vcf"
    with open(fname, 'w') as fh:
        print(dedent(content), file=fh)
    return str(fname)

class TestMalformedVcf:

    def test_no_header(self, tmp_path):
        with pytest.raises(TreeTimeError):
            read_vcf(create_vcf(tmp_path, """\
                ##fileformat=VCFv4.3
                1\t5\t.\tT\tC\t.\t.\t.\tGT\t1\t1\t1
            """))

    def test_meta_info_after_header(self, tmp_path):
        with pytest.raises(TreeTimeError):
            read_vcf(create_vcf(tmp_path, """\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\tsample_C
                ##fileformat=VCFv4.3
                1\t5\t.\tT\tC\t.\t.\t.\tGT\t1\t1\t1
            """))
    
    def test_inconsistent_sample_numbers(self, tmp_path):
        with pytest.raises(TreeTimeError):
            read_vcf(create_vcf(tmp_path, """\
                ##fileformat=VCFv4.3
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\tsample_C
                1\t5\t.\tT\tC\t.\t.\t.\tGT\t1\t1\t1
                1\t6\t.\tT\tC\t.\t.\t.\tGT\t1\t1
            """))

    def test_no_data_lines(self, tmp_path):
        with pytest.raises(TreeTimeError):
            read_vcf(create_vcf(tmp_path, """\
                ##fileformat=VCFv4.3
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\tsample_C
            """))

    def test_consistent_ploidy(self, tmp_path):
        """
        The spec heavily implies that ALT alleles must be supplied for each copy
        of the chromosome, even for no-calls: "If a call cannot be made for a
        sample at a given locus, '.' must be specified for each missing allele
        in the GT field (for example './.' for a diploid genotype and '.' for
        haploid genotype)."
        """
        with pytest.raises(TreeTimeError):
            read_vcf(create_vcf(tmp_path, """\
                ##fileformat=VCFv4.3
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B
                1\t5\t.\tT\tC\t.\t.\t.\tGT\t1\t1\t1/1
            """))

class TestMultipleChromosomes:
    """
    These are valid VCFs but TreeTime cannot yet parse them
    """
    def test_multiple_chromosomes(self, tmp_path):
        with pytest.raises(TreeTimeError):
            read_vcf(create_vcf(tmp_path, """\
                ##fileformat=VCFv4.3
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\tsample_C
                1\t5\t.\tT\tC\t.\t.\t.\tGT\t1\t1\t1
                2\t6\t.\tT\tC\t.\t.\t.\tGT\t1\t1\t1
            """))

def zero_based(idx):
    """Useful as a form of code commentary. Convert idx from 1-based to 0-based"""
    return idx-1

def write_data(tmp_path_factory, sample_names, data_lines, reference_seq, meta_lines=None, reference_name='reference_name'):
    """helper function to create and write a VCF and FASTA file. Returns a tuple of filenames"""
    if meta_lines:
        vcf = "\n".join(meta_lines)+"\n"
    else:
        vcf = dedent("""\
            ##fileformat=VCFv4.3
            ##contig=<ID=1,length=50>
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        """)
    vcf += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n"
    for line in data_lines:
        vcf += "\t".join(line) + "\n"
    reference=dedent(f"""\
        >{reference_name}
        {reference_seq}
    """)
    dir = tmp_path_factory.mktemp("data")
    vcf_fn = dir / "snps.vcf"
    ref_fn = dir / "reference.fasta"
    with open(vcf_fn, 'w') as fh:
        print(vcf, file=fh)
    with open(ref_fn, 'w') as fh:
        print(reference, file=fh)
    return (str(vcf_fn), str(ref_fn))


class TestSimpleVcf:
    @pytest.fixture(scope="class")
    def data(self, tmp_path_factory):
        """
        Creates the VCF and reference files (in a temporary location) to be used by
        the tests in this class.
        Returns their filenames, as well as information about their contents which
        we can use as comparisons in tests.
        NOTE: this function is only run once - _not_ once per test.
        """
        sample_names = ["sample_A", "sample_B", "sample_C"]
        data_lines = [
            ["1", "5",  ".", "T", "C",   ".", ".", ".", "GT", "1", "1", "1"],
            ["1", "7",  ".", "T", "G",   ".", ".", ".", "GT", "1", "1", "0"],
            ["1", "14", ".", "C", "T",   ".", ".", ".", "GT", "1", "1", "0"],
            ["1", "18", ".", "C", "T",   ".", ".", ".", "GT", "0", "0", "1"],
            ["1", "28", ".", "A", "N",   ".", ".", ".", "GT", "0", "0", "1"],
            ["1", "29", ".", "A", "N",   ".", ".", ".", "GT", "0", "0", "1"],
            ["1", "30", ".", "A", "N",   ".", ".", ".", "GT", "0", "0", "1"],
            ["1", "33", ".", "A", "C,G", ".", ".", ".", "GT", "1", "2", "0"],
            ["1", "39", ".", "C", "T",   ".", ".", ".", "GT", "1", "0", "0"],
            ["1", "42", ".", "G", "A",   ".", ".", ".", "GT", "0", "1", "0"],
            ["1", "43", ".", "A", "T",   ".", ".", ".", "GT", "1", "1", "0"],
        ]
        positions = [int(line[1]) for line in data_lines]
        reference_seq = 'AAAAAAAAAATGCCCTGCGGGTAAAAAAAAAAAAACTACTTGACCATAAA'
        filenames = write_data(tmp_path_factory, sample_names, data_lines, reference_seq)
        return {
            "filenames": filenames,
            "positions": positions,
            "data_lines": data_lines,
            "sample_names": sample_names,
            "reference_seq": reference_seq,
        }

    def test_mutation_structure(self, data):
        vcf_data = read_vcf(*data['filenames'])
        # The structure of insertions & sequences is a dict a key for each sample, irregardless of whether
        # any insertions/mutations are defined for that sample.
        # samples = {'sample_A', 'sample_B', 'sample_C'}
        assert(set(data['sample_names'])==set(vcf_data['sequences'].keys()))
        assert(set(data['sample_names'])==set(vcf_data['insertions'].keys()))

    def test_mutation_parsing(self, data):
        """
        Test that the correct SNPs (ALTs) are extracted from the VCF for the three samples defined in the VCF
        """
        vcf_data = read_vcf(*data['filenames'])
        def summarise(sample):
            # uses 1-based position so we can easily compare with the VCF file
            return ",".join([f"{pos+1}{alt}" for pos,alt in vcf_data['sequences'][sample].items()])

        assert("5C,7G,14T,33C,39T,43T"==summarise('sample_A'))
        assert("5C,7G,14T,33G,42A,43T"==summarise('sample_B'))
        assert("5C,18T,28N,29N,30N"==summarise('sample_C'))

    def test_no_insertions_parsed(self, data):
        vcf_data = read_vcf(*data['filenames'])
        for sample, ins in vcf_data['insertions'].items():
            assert(ins=={})

    def test_reference_sequence(self, data):
        vcf_data = read_vcf(*data['filenames'])
        assert(data['reference_seq']==vcf_data['reference'])

    def test_positions_with_mutations(self, data):
        vcf_data = read_vcf(*data['filenames'])
        assert(data['positions'] == [pos+1 for pos in vcf_data['positions']]) 


class TestNoCallsOrMissing:
    """
    Tests a few overlapping concepts:
    - The genotype field (i.e. within the column of a sample name) may be "." which results in a
      no call allele. So we call this as SNP(s) changing the REF to N(s). The 4.2 spec refers
      to this as "a call cannot be made for a sample at a given locus". Spec 4.3 defines this
      as "If a call cannot be made for a sample at a given locus, '.' must be specified for each
      missing allele in the GT field (for example ./. for a diploid genotype and . for haploid
      genotype)"
    - The ALT allele may be "*". v4.2. defines this as "missing due to a upstream deletion (sic)"
      and 4.3 defines this as "allele missing due to overlapping deletion". Either way we encode
      this similarly to the no-call above, i.e. changing the REF to N(s). Note that the allele
      must be the one-character "*" - this can't appear within other bases.
    - The ALT allele may be ".". This was introduced in version 4.3 as "a MISSING value (no variant)"
      So we treat this the same as "*" above. Again, this must be the entire allele, '.' can't appear
      as part of other bases. Note that while this doesn't exist in v4.2, it is commonly found in
      v4.2 VCF files, e.g. those produced by `snp_sites`.
    """

    @pytest.fixture(scope="class")
    def data(self, tmp_path_factory):
        reference="ATCGA"
        sample_names = ["sample_A", "sample_B", "sample_C", "sample_D"]
        data_lines = [
            ## No call alleles
            ["1", "2",  ".", "T", "A",    ".", ".", ".", "GT", ".", "1", "0", "0"], # sample A has T2N
            ["1", "4",  ".", "GA", "C",   ".", ".", ".", "GT", ".", "1", "0", "0"], # sample A has GA -> NN
            ## star alleles & dot alleles indicate missing
            ["1", "3",  ".", "C", "*,G",  ".", ".", ".", "GT", "1", "2", "0", "0"], # sample A has C->N
            ["1", "3",  ".", "C", "T,.",  ".", ".", ".", "GT", "0", "0", "1", "2"], # sample D has C->N
            ["1", "4",  ".", "GA", ".,*", ".", ".", ".", "GT", "0", "0", "1", "2"], # both samples C & D have G->N and A->N
        ]
        filenames = write_data(tmp_path_factory, sample_names, data_lines, reference)
        return {"filenames": filenames}

    def test_no_call_allele(self, data):
        vcf_data = read_vcf(*data['filenames'])
        assert(vcf_data['sequences']['sample_A'][zero_based(2)]=='N')
        assert(vcf_data['sequences']['sample_A'][zero_based(4)]=='N')
        assert(vcf_data['sequences']['sample_A'][zero_based(5)]=='N')
        assert(vcf_data['sequences']['sample_B'][zero_based(4)]=='C') 
        assert(vcf_data['sequences']['sample_B'][zero_based(5)]=='-')

    def test_star_allele(self, data):
        vcf_data = read_vcf(*data['filenames'])
        assert(vcf_data['sequences']['sample_A'][zero_based(3)]=='N')
        assert(vcf_data['sequences']['sample_B'][zero_based(3)]=='G')

        assert(vcf_data['sequences']['sample_D'][zero_based(4)]=='N')
        assert(vcf_data['sequences']['sample_D'][zero_based(5)]=='N')

    def test_dot_allele(self, data):
        vcf_data = read_vcf(*data['filenames'])
        assert(vcf_data['sequences']['sample_C'][zero_based(3)]=='T')
        assert(vcf_data['sequences']['sample_D'][zero_based(3)]=='N')

        assert(vcf_data['sequences']['sample_C'][zero_based(4)]=='N')
        assert(vcf_data['sequences']['sample_C'][zero_based(5)]=='N')


    def test_malformed_star_allele(self, tmp_path):
        """* must be an entire allele, not together with other bases"""
        with pytest.raises(TreeTimeError):
            read_vcf(create_vcf(tmp_path, """\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\tsample_C
                ##fileformat=VCFv4.3
                1\t5\t.\tT\tC*\t.\t.\t.\tGT\t1\t1\t1
            """))

    def test_malformed_dot_allele(self, tmp_path):
        """. must be an entire allele, not together with other bases"""
        with pytest.raises(TreeTimeError):
            read_vcf(create_vcf(tmp_path, """\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\tsample_C
                ##fileformat=VCFv4.3
                1\t5\t.\tT\tC.,G\t.\t.\t.\tGT\t1\t1\t1
            """))


class TestVcfSpecExamples:
    """
    Test the examples provided in the VCF specs. Note that the examples
    aren't provided as per-sample VCF files, so there's a bit of interpretation
    required to create the example data. Note that the actual `test_<name>`
    methods are added dynamically after the class definition.
    """
    def create(self, tmp_path, input):
        lines = ["##fileformat=VCFv4.2", # Actual version may differ, but we don't parse this so it doesn't matter (yet)
            "##contig=<ID=1,length=50>",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(input['samples']),
        ]
        for d in input['data']:
            lines.append(f"{d['chrom']}\t{d['pos']}\t.\t{d['ref']}\t{d['alt']}\t.\t.\t.\tGT\t" + "\t".join([str(x) for x in d['values']]))
        return create_vcf(tmp_path, "\n".join(lines))

    @staticmethod
    def add_test(example_data):
        def test(self, tmp_path):
            vcf = read_vcf(self.create(tmp_path, example_data))
            for sample_name, sample_data in example_data['expect_sequences'].items():
                assert(sample_name in vcf['sequences'])
                for snp in sample_data:
                    assert(len(snp)==2)
                    assert(vcf['sequences'][sample_name][zero_based(int(snp[0]))] == snp[1])
                assert(len(vcf['sequences'][sample_name].keys()) == len(sample_data))
            for sample_name, sample_data in example_data['expect_insertions'].items():
                assert(sample_name in vcf['insertions'])
                for ins in sample_data:
                    assert(vcf['insertions'][sample_name][zero_based(int(ins[0]))] == ins[1:])
                assert(len(vcf['insertions'][sample_name].keys()) == len(sample_data))
        return test

"""
-------- comment re: TreeTime's encoding of insertions --------

Give a reference base of 'C' at 1-based position 2, and an
insertion _after_ this of 'A', TreeTime encodes this as an insertion
of 'CA' at (0-based) position 1. This was unexpected to me,
I would have expected an insertion of just 'A'. However I've written
the tests to conform to TreeTime's existing behaviour (version 0.11.1)

---------------------------------------------------------------
"""

vcf_spec_data = {
    'version_4_2': {
        '5_1_1': {
            'samples': ['A', 'B', 'C'],
            'data': [
                {'chrom': 20, 'pos': 3, 'ref': 'C', 'alt': 'G', 'values': [1, 0, 0]},
                {'chrom': 20, 'pos': 2, 'ref': 'TC', 'alt': 'T,TCA', 'values': [0, 1, 2]}
            ],
            # Note that with TC->TCA, there are no sequence changes (T->T, C->C) but one insertion ("A")
            'expect_sequences': {'A': ['3G'], 'B': ['3-']},
            'expect_insertions': {'C': ['3CA']} # see comment above. Actually only a single base "A" insertion
        },
        '5_1_2': {
            'samples': ['A', 'B'],
            'data': [
                {'chrom': 20, 'pos': 2, 'ref': 'TC', 'alt': 'TG,T', 'values': [1,2]}
            ],
            'expect_sequences': {'A': ['3G'], 'B': ['3-']},
            'expect_insertions': {}
        },
        '5_1_3': {
            'samples': ['A', 'B', 'C'],
            ## This example is quite hard to understand for sample A. The alignment provided indicates that
            ## TCG -> T-G is encoded as "ref: TCG, alt: TG". But without aligning the alt to the ref, the only
            ## sane interpretation of this is 2 changes: C->G and G->-. This is what the spec means by
            ## "the molecular equivalence explicitly listed above in the per-base alignment is discarded so the
            ## actual placement of equivalent g isn't retained" (also explained in example 5.2.4)
            ## Similarly, for sample C, the actual event is described as "A base is inserted wrt the reference
            ## sequence" but the way the VCF file reads we are going to have a SNP (G->A) + a insertion (G)
            'data': [
                {'chrom': 20, 'pos': 2, 'ref': 'TCG', 'alt': 'TG,T,TCAG', 'values': [1,2,3]}
            ],
            'expect_sequences': {'A': ['3G', '4-'], 'B': ['3-', '4-'], 'C': ['4A']},
            'expect_insertions': {'C': ['4AG']} # See comment above
        },
        '5_2_1': {
            'samples': ['A'],
            'data': [
                {'chrom': 20, 'pos': 3, 'ref': 'C', 'alt': 'T', 'values': [1]}
            ],
            'expect_sequences': {'A': ['3T']},
            'expect_insertions': {}
        },
        '5_2_2': {
            'samples': ['A'],
            'data': [
                {'chrom': 20, 'pos': 3, 'ref': 'C', 'alt': 'CTAG', 'values': [1]}
            ],
            'expect_sequences': {},
            'expect_insertions': {'A': ['3CTAG']} # See comment above
        },
        '5_2_3': {
            'samples': ['A'],
            'data': [
                {'chrom': 20, 'pos': 2, 'ref': 'TCG', 'alt': 'T', 'values': [1]}
            ],
            'expect_sequences': {'A': ['3-', '4-']},
            'expect_insertions': {}
        },
    }
}
# dynamically create a test for each example in the above data (by adding methods to the class)
for spec_version, spec_data in vcf_spec_data.items():
    for example_key, example_data in spec_data.items():
        setattr(TestVcfSpecExamples, f"test_{spec_version}_example_{example_key}", TestVcfSpecExamples.add_test(example_data))


class TestMutationAndInsertion:
    """
    Tests the situation where a reference base is mutated _and_ there's an insertion
    """

    @pytest.fixture(scope="class")
    def data(self, tmp_path_factory):
        reference="ATCGA"
        sample_names = ["sample_A", "sample_B"]
        data_lines = [
            ["1", "2",  ".", "T", "GA",    ".", ".", ".", "GT", "1", "0"], # sample A has both a C->G mutation _and_ a subsequent "A" insertion
            ["1", "3",  ".", "CGA", "NTTTA,CCGT",   ".", ".", ".", "GT", "1", "2"], # complex! Both samples have multiple mutations + an insertion
        ]
        filenames = write_data(tmp_path_factory, sample_names, data_lines, reference)
        return {"filenames": filenames}

    def test_single_ref_base_mutation_and_insertion(self, data):
        # This case was missed in treetime 0.11.1
        vcf_data = read_vcf(*data['filenames'])
        assert(vcf_data['sequences']['sample_A'][zero_based(2)]=='G')
        assert(vcf_data['insertions']['sample_A'][zero_based(2)]=='GA') # see comment above re: insertion encoding

    def test_multi_ref_base_mutations_and_insertion(self, data):
        # This case was missed in treetime 0.11.1
        vcf_data = read_vcf(*data['filenames'])
        assert(vcf_data['sequences']['sample_A'][zero_based(3)]=='N')
        assert(vcf_data['sequences']['sample_A'][zero_based(4)]=='T')
        assert(vcf_data['sequences']['sample_A'][zero_based(5)]=='T')
        assert(vcf_data['insertions']['sample_A'][zero_based(5)]=='TTA') # see comment above re: insertion encoding

        assert(vcf_data['sequences']['sample_B'][zero_based(4)]=='C')
        assert(vcf_data['sequences']['sample_B'][zero_based(5)]=='G')
        assert(vcf_data['insertions']['sample_B'][zero_based(5)]=='GT') # see comment above re: insertion encoding

class TestMetadataParsing:
        
    def test_simple_haploid(self, tmp_path):
        data = read_vcf(create_vcf(tmp_path, """\
            ##fileformat=VCFv4.3
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\tsample_C
            1\t5\t.\tT\tC\t.\t.\t.\tGT\t1\t1\t1
            1\t6\t.\tT\tC\t.\t.\t.\tGT\t1\t1\t1
        """))
        assert(data['metadata']['chrom']=='1')
        assert(data['metadata']['ploidy']==1)
        assert(data['metadata']['meta_lines']==['##fileformat=VCFv4.3'])

    def test_simple_diploid(self, tmp_path):
        data = read_vcf(create_vcf(tmp_path, """\
            ##fileformat=VCFv4.3
            ##contig=<ID=1,length=50>
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A
            foo\t5\t.\tT\tC\t.\t.\t.\tGT\t1/1
            foo\t6\t.\tT\tC\t.\t.\t.\tGT\t.|.
        """))
        assert(data['metadata']['chrom']=='foo')
        assert(data['metadata']['ploidy']==2)
        assert(data['metadata']['meta_lines']==['##fileformat=VCFv4.3', '##contig=<ID=1,length=50>', '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'])



    
def roundtrip(tmp_path, sample_names, data_lines, reference, meta_lines, pass_metadata=True):
    """
    Write the provided data to a (temporary) VCF file.
    Then read it via `read_vcf` and store this as `vcf_a`.
    Then write this data via `write_vcf` and read it back, storing as `vcf_b`.
    Returns a tuple of (vcf_a, vcf_b)
    """

    dir = tmp_path / 'data'
    dir.mkdir()
    vcf_filename_a = str(dir / "a.vcf")
    vcf_filename_b = str(dir / "b.vcf")
    reference_filename = str(dir / 'reference.fasta')

    with open(str(dir / 'reference.fasta'), 'w') as fh:
        print(dedent(f"""\
            >reference_name
            {reference}
        """), file=fh)

    with open(vcf_filename_a, 'w') as fh:
        vcf_lines = meta_lines[:] + \
            ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names)] + \
            ["\t".join(data) for data in data_lines]
        print("\n".join(vcf_lines), file=fh)

    vcf_a = read_vcf(vcf_filename_a, reference_filename)

    ## take this data and write it out as a VCF
    tree_dict = {'sequences': vcf_a['sequences'], 'reference': vcf_a['reference'], 'positions': vcf_a['positions']}
    if pass_metadata:
        tree_dict['metadata'] = vcf_a['metadata']
    write_vcf(tree_dict, vcf_filename_b)

    ## then read in this (just-created) VCF
    vcf_b = read_vcf(vcf_filename_b, reference_filename)

    ## uncomment the following & run pytest with "-s" to see the contents of the VCF file being written by TreeTime
    # print("_________________________________")
    # with open(vcf_filename_b) as fh:
    #     for line in fh:
    #         print(line, end='')
    # print("\n_________________________________")

    return (vcf_a, vcf_b)


class TestWriting:
    """
    Write a simple VCF file out (created by hand), parse it with `read_vcf` (called "vcf_a")
    then write that data out and read it back in (called "vcf_b"). vcf_a should equal vcf_b
    for the data we care about.
    """

    def test_basic_roundtripping_vcf(self, tmp_path):
        [vcf_a, vcf_b] = roundtrip(tmp_path,
            reference = "ATCGACC",
            meta_lines = [
                "##fileformat=VCFv4.3",
                "##contig=<ID=1,length=7>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            ],
            sample_names = ["sample_A", "sample_B"],
            data_lines = [
                ["foo", "2",  ".", "T", "G",    ".", ".", ".", "GT", "1", "0"],
                ["foo", "3",  ".", "C", "GA,*", ".", ".", ".", "GT", "1", "2"], # "A" insertion will be lost on roundtrip!
                ["foo", "5",  ".", "A", "TC",   ".", ".", ".", "GT", ".", "1"],
                ["foo", "6",  ".", "CC", "C",   ".", ".", ".", "GT", "1", "0"], # del of (1-based) pos 7 in sample A
            ],
            pass_metadata=False
        )

        ## reference is certainly the same, the same FASTA file is being used for both read_vcf commands, but check anyway
        assert(vcf_a['reference']==vcf_a['reference'])
        ## insertions, if there are any, _won't_ be the same because `write_vcf` doesn't read insertions even if they're provided
        ## as input. Note that vcf_a does have an insertion!

        ## The meta-lines are different (as expected) and since we are not passing the metadata information to `write_vcf`
        ## we get back the defaults, which differ from the input
        assert(vcf_a['metadata']['chrom'] == "foo" and vcf_b['metadata']['chrom'] == "1")
        assert(vcf_a['metadata']['ploidy'] == 1 and vcf_b['metadata']['ploidy'] == 2)

        ## positions should be the same (positions are only reflective of data in `sequences`, so the ignoring of insertions is ok)
        assert(vcf_a['positions']==vcf_a['positions'])
        ## check sequences are the same
        assert(vcf_a['sequences']==vcf_b['sequences'])

    def test_basic_roundtripping_vcf_with_metadata(self, tmp_path):
        [vcf_a, vcf_b] = roundtrip(tmp_path,
            reference = "ATCGACC",
            meta_lines = [
                "##fileformat=VCFv4.3",
                "##contig=<ID=1,length=7>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            ],
            sample_names = ["sample_A", "sample_B"],
            data_lines = [
                ["1", "2",  ".", "T", "G",    ".", ".", ".", "GT", "1", "0"],
                ["1", "3",  ".", "C", "GA,*", ".", ".", ".", "GT", "1", "2"], # "A" insertion will be lost on roundtrip!
                ["1", "5",  ".", "A", "TC",   ".", ".", ".", "GT", ".", "1"],
                ["1", "6",  ".", "CC", "C",   ".", ".", ".", "GT", "1", "0"], # del of (1-based) pos 7 in sample A
            ],
            pass_metadata=True
        )

        ## reference is certainly the same, the same FASTA file is being used for both read_vcf commands, but check anyway
        assert(vcf_a['reference']==vcf_a['reference'])
        ## insertions, if there are any, _won't_ be the same because `write_vcf` doesn't read insertions even if they're provided
        ## as input. Note that vcf_a does have an insertion!

        ## The meta-lines are different (as expected) but the chrom + ploidy should be the same
        assert(vcf_a['metadata']['chrom'] == vcf_b['metadata']['chrom'])
        assert(vcf_a['metadata']['ploidy'] == vcf_b['metadata']['ploidy'])

        ## positions should be the same (positions are only reflective of data in `sequences`, so the ignoring of insertions is ok)
        assert(vcf_a['positions']==vcf_a['positions'])
        ## check sequences are the same
        assert(vcf_a['sequences']==vcf_b['sequences'])
