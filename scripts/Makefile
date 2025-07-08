.PHONY: params

params: params-setup params-abdlop params-quad params-many

params-setup: 
	if [ -d "lattice-estimator" ]; then \
		cp -r lattice-estimator/estimator estimator; \
		rm -rf lattice-estimator;\
	fi\
	

params-abdlop: params-setup
	sage rf_abdlop_param.sage ../params/params1.sage > ../tests/rf-abdlop-params1.h
	sage rf_abdlop_param.sage ../params/params2.sage > ../tests/rf-abdlop-params2.h
	sage rf_abdlop_param.sage ../params/params3.sage > ../tests/rf-abdlop-params3.h
	sage rf_abdlop_param.sage ../params/params4.sage > ../tests/rf-abdlop-params4.h
	sage rf_abdlop_param.sage ../params/params5.sage > ../tests/rf-abdlop-params5.h

	rm rf_abdlop_param.sage.py

params-quad: params-setup
	sage rf_quad_param.sage ../params/params1.sage > ../tests/rf-quad-params1.h
	sage rf_quad_param.sage ../params/params2.sage > ../tests/rf-quad-params2.h
	sage rf_quad_param.sage ../params/params3.sage > ../tests/rf-quad-params3.h
	sage rf_quad_param.sage ../params/params4.sage > ../tests/rf-quad-params4.h
	sage rf_quad_param.sage ../params/params5.sage > ../tests/rf-quad-params5.h

	rm rf_quad_param.sage.py

params-many: params-setup
	sage rf_quad_eval_param.sage ../params/params1.sage > ../tests/rf-quad-eval-params1.h
	sage rf_quad_eval_param.sage ../params/params2.sage > ../tests/rf-quad-eval-params2.h
	sage rf_quad_eval_param.sage ../params/params3.sage > ../tests/rf-quad-eval-params3.h
	sage rf_quad_eval_param.sage ../params/params4.sage > ../tests/rf-quad-eval-params4.h
	sage rf_quad_eval_param.sage ../params/params5.sage > ../tests/rf-quad-eval-params5.h
	rm rf_quad_eval_param.sage.py