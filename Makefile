build:
	mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Debug && make -j$(nproc)

test: build
	cd build && ctest --output-on-failure

clean:
	rm -rf build

