include ../makefile.include
plot_DOS_interp.x: plot_DOS_interp.f90
	ifort -I$(IOTK_DIR) constant.f90 upfpsp_module.f90 transform_to_upper.f90 get_ALI.f gaussj.f d_matrix.f90\
		scan_key_words.f90 read_key_words.f90 readusp_head.f90 read_upf_v2.f90 \
		plot_DOS_interp.f90 $(IOTK_DIR)libiotk.a -lmkl_intel_lp64 \
		-lmkl_sequential -lmkl_core -qopenmp -o $@

#	ifort -I/opt/iotk/src constant.f90 upfpsp_module.f90 transform_to_upper.f90 \
#		scan_key_words.f90 read_key_words.f90 readusp_head.f90 read_upf_v2.f90 \
#		plot_DOS_interp.f90 $(IOTK_DIR)libiotk.a -lmkl_intel_lp64 \
#		-lmkl_sequential -lmkl_core -g -traceback  -check all -qopenmp -o $@
clean:
	rm -f *.x *.mod
