rmm:
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ test_rmm.C

.PHONY: clean

clean:
	$(RM) rmm
