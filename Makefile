ifdef FES_DIR
uptide/fes2012.so: fes/fes2012.c fes/fes2012.pyx
	python setup.py build_ext --with-fes -i
TEST_REQUIREMENTS = uptide/fes2012.so
export FES_DATA = $(FES_DIR)/data/
else
TEST_REQUIREMENTS =
endif
test: $(TEST_REQUIREMENTS)
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
