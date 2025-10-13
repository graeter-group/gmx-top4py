.PHONY: Makefile setup-docs preview docs clear-docs watch

test:
	pytest tests --cov=gmx_top4py --cov-report=term-missing

format:
	black src tests

venv:
	uv lock
	uv sync
	source .venv/bin/activate

# Write explicitly all dependencies to a requirements_explicit.txt file
# requirements:
# 	source .venv/bin/activate
# 	uv export --format requirements.txt > requirements_explicit.txt

