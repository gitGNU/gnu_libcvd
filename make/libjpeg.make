jpeg_objs=
ifeq ($(shell $(CC) make/test.c -Dtest_jpeg -o tmp -ljpeg 2>/dev/null;echo $$?),0)
	jpeg_objs=pnm_src/jpeg.o
	options += jpeg
	images+=jpeg
	TESTLIB+=-ljpeg
endif

