
// Copyright 2019 Adam Campbell, Seth Hall, Andrew Ensor
// Copyright 2019 High Performance Computing Research Laboratory, Auckland University of Technology (AUT)

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DEGRIDDER_H_
#define DEGRIDDER_H_

	#ifndef C
		#define C 299792458.0
	#endif

	typedef struct Config {
		int grid_size;
		double cell_size;
		bool right_ascension;
		double frequency_hz;
		int kernel_size;
		int oversampling;
		double uv_scale;
		int num_visibilities;
		char *grid_real_source_file;
		char *grid_imag_source_file;
		char *kernel_real_source_file;
		char *kernel_imag_source_file;
		char *visibility_source_file;
		char *visibility_dest_file;
	} Config;
	
	typedef struct Visibility {
		double u;
		double v;
		double w;
	} Visibility;

	typedef struct Complex {
		double real;
		double imag;
	} Complex;

void init_config(Config *config);

bool load_grid(Config *config, Complex *grid);

bool load_visibilities(Config *config, Visibility **vis_uvw, Complex **vis_intensities);

void save_visibilities(Config *config, Visibility *vis_uvw, Complex *vis_intensity);

void execute_degridding(Config *config, Complex *grid, Visibility *vis_uvw, Complex *vis_intensities, Complex *kernel, int num_visibilities);

bool load_kernel(Config *config, Complex *kernel);

Complex complex_multiply(Complex z1, Complex z2);

void clean_up(Complex **grid, Visibility **visibilities, Complex **vis_intensities, Complex **kernel);

void unit_test_init_config(Config *config);

double unit_test_generate_approximate_visibilities(void);

#endif // DEGRIDDER_H_

#ifdef __cplusplus
}
#endif
