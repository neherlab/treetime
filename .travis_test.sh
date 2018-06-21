cd test
bash command_line_tests.sh
OUT=$?
if [ "$OUT" != 0 ]; then
	exit 1
fi
nosetests test_treetime.py
