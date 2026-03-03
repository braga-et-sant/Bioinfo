# Bioinfo Lab pipeline - Exercise 6
# Usage:
#   make            # runs everything
#   make clean      # removes generated files

ACC=NC_001416

FASTA=lambda_phage.fasta
PYOUT=lambda_phage_py_rev.fasta
EMBOUT=lambda_phage_rev.fasta

.PHONY: all download pyrev emboss validate clean

all: validate

download: $(FASTA)

$(FASTA):
	efetch -db nuccore -id $(ACC) -format fasta > $(FASTA)

pyrev: $(PYOUT)

$(PYOUT): $(FASTA) revcomp.py
	python revcomp.py $(FASTA) > $(PYOUT)

emboss: $(EMBOUT)

$(EMBOUT): $(FASTA)
	revseq -sequence $(FASTA) -outseq $(EMBOUT) -notag

validate: download pyrev emboss
	@py_md5=$$(grep -v '^>' $(PYOUT)  | tr -d '\n\r' | md5sum | awk '{print $$1}'); \
	em_md5=$$(grep -v '^>' $(EMBOUT) | tr -d '\n\r' | md5sum | awk '{print $$1}'); \
	if [ "$$py_md5" = "$$em_md5" ]; then \
		echo "Validation Passed"; \
	else \
		echo "Validation Failed"; \
		echo "python md5: $$py_md5"; \
		echo "emboss  md5: $$em_md5"; \
		exit 1; \
	fi

clean:
	rm -f $(FASTA) $(PYOUT) $(EMBOUT)