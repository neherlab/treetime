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
