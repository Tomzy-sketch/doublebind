# Unix and macOS shortcut. On Windows, run the python commands in README.md directly.
.PHONY: test data verify variants numbers

test:
	python -m pytest

data:
	python -m doublebind data

verify:
	python -m doublebind verify

variants:
	python -m doublebind variants gyra

numbers:
	python -m doublebind numbers
