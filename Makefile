selector:
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ selector_test.C

.PHONY: clean

clean:
	$(RM) selector
