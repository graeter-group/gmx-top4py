.PHONY: Makefile setup-docs preview docs clear-docs watch

test:
	pytest tests --cov=gmxtop --cov-report=term-missing

format:
	black src tests

venv:
	uv lock
	uv sync
	source .venv/bin/activate

failing:
	echo None

package:
	rm -r dist
	python -m build 

testpypi:
	python -m twine upload --repository testpypi dist/*

pypi:
	python -m twine upload dist/*

# Write explicitly all dependencies to a requirements_explicit.txt file
# requirements:
# 	source .venv/bin/activate
# 	uv export --format requirements.txt > requirements_explicit.txt

