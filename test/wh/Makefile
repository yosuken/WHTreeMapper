TOOL := $(realpath ../../wh_tree_mapper)
QUERY := query.faa
OUTDIR := out

.PHONY: test clean

test: clean
	$(TOOL) --prot -i $(QUERY) -o $(OUTDIR) -n 1
	@echo ""
	@echo "=== detected.tsv ==="
	@cat $(OUTDIR)/detected.tsv
	@echo ""
	@echo "=== test passed ==="

clean:
	rm -rf $(OUTDIR)
