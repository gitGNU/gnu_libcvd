x11_objs=

ifeq ($(shell $(CC) make/test.c -Dtest_x11 -o tmp -I$(X11_I) -L$(X11_L) -lX11 -lGL 2>/dev/null;echo $$?),0)
	x11_objs=cvd_src/videodisplay.o
	options+=videodisplay
	TESTLIB+=-L$(X11_L) -lX11 -lGL
endif
