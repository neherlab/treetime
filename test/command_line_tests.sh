all_tests=0

treetime homoplasy --aln ../data/H3N2_NA_allyears_NA.20.fasta --tree ../data/H3N2_NA_allyears_NA.20.nwk
retval="$?"
if [ "$retval" == 0 ]; then
	echo "homoplasy_scanning ok"
else
	echo "homoplasy_scanning failed $retval"
	((all_tests++))
fi

treetime ancestral --aln ../data/H3N2_NA_allyears_NA.20.phylip --tree ../data/H3N2_NA_allyears_NA.20.nwk
retval="$?"
if [ "$retval" == 0 ]; then
	echo "ancestral_reconstruction ok"
else
	((all_tests++))
	echo "ancestral_reconstruction failed $retval"
fi

treetime clock --tree ../data/H3N2_NA_allyears_NA.20.nex --dates ../data/H3N2_NA_allyears_NA.20.metadata.csv --sequence-length 1400
retval="$?"
if [ "$retval" == 0 ]; then
	echo "temporal_signal ok"
else
	((all_tests++))
	echo "temporal_signal failed $retval"
fi

treetime tree --aln ../data/H3N2_NA_allyears_NA.20.fasta --tree ../data/H3N2_NA_allyears_NA.20.nwk --dates ../data/H3N2_NA_allyears_NA.20.metadata.csv
retval="$?"
if [ "$retval" == 0 ]; then
	echo "timetree_inference ok"
else
	((all_tests++))
	echo "timetree_inference failed $retval"
fi

treetime mugration --tree ../data/Zika_tree.newick --states ../data/Zika_metadata.csv --weights ../data/Zika_country_weights.csv --attribute country
retval="$?"
if [ "$retval" == 0 ]; then
	echo "mugration ok"
else
	((all_tests++))
	echo "mugration failed $retval"
fi

if [ "$all_tests" == 0 ];then
	echo "All tests passed"
	exit 0
else
	exit "$all_tests"
fi

