test:
	python -mpytest tests/

lint:
	@echo "    Linting uptide codebase"
	@python -m flake8 uptide
	@echo "    Linting uptide test suite"
	@python -m flake8 tests
