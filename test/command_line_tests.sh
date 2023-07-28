all_tests=0

treetime homoplasy --aln treetime_examples/data/h3n2_na/h3n2_na_20.fasta --tree treetime_examples/data/h3n2_na/h3n2_na_20.nwk
retval="$?"
if [ "$retval" == 0 ]; then
	echo "homoplasy_scanning ok"
else
	echo "homoplasy_scanning failed $retval"
	((all_tests++))
fi

treetime ancestral --aln treetime_examples/data/h3n2_na/h3n2_na_20.phylip --tree treetime_examples/data/h3n2_na/h3n2_na_20.nwk
retval="$?"
if [ "$retval" == 0 ]; then
	echo "ancestral_reconstruction ok"
else
	((all_tests++))
	echo "ancestral_reconstruction failed $retval"
fi

treetime clock --tree treetime_examples/data/h3n2_na/h3n2_na_20.nex --dates treetime_examples/data/h3n2_na/h3n2_na_20.metadata.csv --sequence-length 1400
retval="$?"
if [ "$retval" == 0 ]; then
	echo "temporal_signal ok"
else
	((all_tests++))
	echo "temporal_signal failed $retval"
fi

treetime --aln treetime_examples/data/h3n2_na/h3n2_na_20.fasta --tree treetime_examples/data/h3n2_na/h3n2_na_20.nwk --dates treetime_examples/data/h3n2_na/h3n2_na_20.metadata.csv --confidence --covariation
retval="$?"
if [ "$retval" == 0 ]; then
	echo "timetree_inference ok"
else
	((all_tests++))
	echo "timetree_inference failed $retval"
fi

treetime mugration --tree treetime_examples/data/zika/zika.nwk --states treetime_examples/data/zika/zika.metadata.csv --weights treetime_examples/data/zika/zika.country_weights.csv --attribute country
retval="$?"
if [ "$retval" == 0 ]; then
	echo "mugration ok"
else
	((all_tests++))
	echo "mugration failed $retval"
fi

treetime --aln treetime_examples/data/tb/lee_2015.vcf.gz --vcf-reference treetime_examples/data/tb/tb_ref.fasta --tree treetime_examples/data/tb/lee_2015.nwk --clock-rate 1e-7 --dates treetime_examples/data/tb/lee_2015.metadata.tsv
retval="$?"
if [ "$retval" == 0 ]; then
	echo "timetree_inference on vcf data ok"
else
	((all_tests++))
	echo "timetree_inference on vcf data failed $retval"
fi

treetime --tree treetime_examples/data/ebola/ebola.nwk --dates treetime_examples/data/ebola/ebola.metadata.csv --aln treetime_examples/data/ebola/ebola.fasta  --coalescent skyline --gen-per-year 100
retval="$?"
if [ "$retval" == 0 ]; then
	echo "skyline approximation ok"
else
	((all_tests++))
	echo "skyline approximation failed $retval"
fi

# From https://github.com/neherlab/treetime/issues/250 
treetime --tree treetime_examples/data/ebola/ebola.nwk --dates treetime_examples/data/ebola/ebola.metadata.csv --sequence-length 1000
retval="$?"
if [ "$retval" == 0 ]; then
	echo "sequence length only ok"
else
	((all_tests++))
	echo "sequence length only failed $retval"
fi

if [ "$all_tests" == 0 ];then
	echo "All tests passed"
	exit 0
else
	exit "$all_tests"
fi

