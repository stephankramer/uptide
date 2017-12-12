ifdef FES_DIR
fes/fes2012.c: fes/fes2012.pyx
	cd fes; cython fes2012.pyx
uptide/fes2012.so: fes/fes2012.c
	python setup.py build_ext --with-fes -i
TEST_REQUIREMENTS = uptide/fes2012.so
export FES_DATA = $(FES_DIR)/data/
else
TEST_REQUIREMENTS =
endif

test: $(TEST_REQUIREMENTS)
	py.test tests/

lint:
	@echo "    Linting uptide codebase"
	@python -m flake8 uptide
	@echo "    Linting uptide test suite"
	@python -m flake8 tests
