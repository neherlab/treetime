import pytest
from textwrap import dedent
from treetime.vcf_utils import read_vcf
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

def zero_based(idx):
    """Useful as a form of code commentary. Convert idx from 1-based to 0-based"""
    return idx-1

def write_data(tmp_path_factory, sample_names, data_lines, reference_seq, meta_lines=None, reference_name='reference_name'):
    """helper function to create and write a VCF and FASTA file. Returns a tuple of filenames"""
    if meta_lines:
        vcf = "\n".join(meta_lines)
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


class TestNoCallAllele:
    @pytest.fixture(scope="class")
    def data(self, tmp_path_factory):
        sample_names = ["sample_A", "sample_B"]
        data_lines = [
            ["1", "2",  ".", "T", "A",   ".", ".", ".", "GT", ".", "1"], # No-call allele -> make a SNP of T2N in sample_A
            ["1", "4",  ".", "C", "G",   ".", ".", ".", "GT", "1", "1"],
        ]
        filenames = write_data(tmp_path_factory, sample_names, data_lines, "ATGC")
        return {"filenames": filenames}

    def test_no_call_allele(self, data):
        vcf_data = read_vcf(*data['filenames'])
        assert(vcf_data['sequences']['sample_A'][zero_based(2)]=='N')
