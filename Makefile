test:
	if python --version 2>&1 | grep -q 2.7; then \
	  python -munittest discover -v; \
	else \
	  $(MAKE) test_individually; \
	fi
	
# ugly kludge for python 2.6 that doesn't have "-munittest discover"
test_individually:
	export PYTHONPATH=$$PWD:$$PYTHONPATH; \
	for i in tests/test_*.py; do \
	  python $$i -v; done
