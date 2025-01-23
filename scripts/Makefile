.PHONY: params

params: params-setup params-abdlop params-quad params-many

params-setup: 
	if [ -d "lattice-estimator" ]; then \
		cp -r lattice-estimator/estimator estimator; \
		rm -rf lattice-estimator;\
	fi\
	

params-abdlop: params-setup
	sage modified-abdlop-codegen.sage ../params/params1.sage > ../tests/modified-abdlop-params1.h
	sage modified-abdlop-codegen.sage ../params/params2.sage > ../tests/modified-abdlop-params2.h
	sage modified-abdlop-codegen.sage ../params/params3.sage > ../tests/modified-abdlop-params3.h
	sage modified-abdlop-codegen.sage ../params/params4.sage > ../tests/modified-abdlop-params4.h
	sage modified-abdlop-codegen.sage ../params/params5.sage > ../tests/modified-abdlop-params5.h

	rm modified-abdlop-codegen.sage.py

params-quad: params-setup
	sage modified-quad-codegen.sage ../params/params1.sage > ../tests/modified-quad-params1.h
	sage modified-quad-codegen.sage ../params/params2.sage > ../tests/modified-quad-params2.h
	sage modified-quad-codegen.sage ../params/params3.sage > ../tests/modified-quad-params3.h
	sage modified-quad-codegen.sage ../params/params4.sage > ../tests/modified-quad-params4.h
	sage modified-quad-codegen.sage ../params/params5.sage > ../tests/modified-quad-params5.h

	rm modified-quad-codegen.sage.py

params-many: params-setup
	sage modified-quad-eval-codegen.sage ../params/params1.sage > ../tests/modified-quad-eval-params1.h
	sage modified-quad-eval-codegen.sage ../params/params2.sage > ../tests/modified-quad-eval-params2.h
	sage modified-quad-eval-codegen.sage ../params/params3.sage > ../tests/modified-quad-eval-params3.h
	sage modified-quad-eval-codegen.sage ../params/params4.sage > ../tests/modified-quad-eval-params4.h
	sage modified-quad-eval-codegen.sage ../params/params5.sage > ../tests/modified-quad-eval-params5.h
	rm modified-quad-eval-codegen.sage.py