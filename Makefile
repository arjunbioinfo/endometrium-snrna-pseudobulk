.PHONY: test run clean docker-build docker-run

test:
	pytest -v

run:
	snakemake --cores 2

clean:
	rm -rf results figures .snakemake .pytest_cache

docker-build:
	docker build -t endometrium-snrna-pseudobulk .

docker-run:
	docker run --rm -v $(PWD)/data:/app/data -v $(PWD)/results:/app/results -v $(PWD)/figures:/app/figures endometrium-snrna-pseudobulk
