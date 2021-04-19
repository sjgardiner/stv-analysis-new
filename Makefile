respmat:
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ respmat.C

.PHONY: clean

clean:
	$(RM) respmat
