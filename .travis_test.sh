cd test
git clone https://github.com/neherlab/treetime_examples.git
bash command_line_tests.sh
OUT=$?
if [ "$OUT" != 0 ]; then
	exit 1
fi
nosetests test_treetime.py
