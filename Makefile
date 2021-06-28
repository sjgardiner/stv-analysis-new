all: analyzer respmat

respmat: respmat.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

analyzer: analyzer.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

.PHONY: clean

clean:
	$(RM) respmat analyzer
