TOOL := $(realpath ../../WHTreeMapper)
QUERY := query.faa
OUTDIR := out

.PHONY: test clean

test: clean
	$(TOOL) -q $(QUERY) -o $(OUTDIR) -n 1
	@echo ""
	@echo "=== detected.tsv ==="
	@cat $(OUTDIR)/detected.tsv
	@echo ""
	@echo "=== test passed ==="

clean:
	rm -rf $(OUTDIR)
