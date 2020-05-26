MATH=$(shell readlink `which math`)

GRAPHS_SRC=$(wildcard Graph_*.m)
GRAPHS_PDF=$(patsubst %.m,%.pdf,$(wildcard Graph_*.m))

ALLGRAPHS=all_graphs.pdf

all: $(ALLGRAPHS)

Graph_%.pdf: Graph_%.m
	@echo "Generating graph $<..."
	$(MATH) -run < $< > /dev/null 

$(ALLGRAPHS): $(sort $(GRAPHS_PDF))
	pdfunite $^ $@

clean:
	rm -f *.pdf
