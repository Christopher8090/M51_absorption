LLIB = /usr/lib/gcc-lib/i386-linux/2.95.4

NUrad: ./NUrad.f ../../dust/NUdust.f
	ifort -o NUrad -O2 -wp_ipo -w ./NUrad.f ../../dust/NUdust.f -L$(LLIB)

clean:
	rm NUrad
