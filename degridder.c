
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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "degridder.h"

void init_config(Config *config)
{
	// Single dimension of grid (dirty residual image)
	config->grid_size = 16384;
	
	config->right_ascension = true;
	
	config->cell_size = 6.39708380288950e-6;
	
	config->frequency_hz = 100e6;
	
	// Single dimension of basic convolution kernel
	config->kernel_size = 9;
	
	// Kernel oversampling factor
	config->oversampling = 4; // Oxford configuration
	
	// Used to convert visibility uvw coordinates into grid coordinates
	config->uv_scale = config->grid_size * config->cell_size;
	
	// Number of visibilities to process
	config->num_visibilities = 1;
	
	// File location to load grid
	config->grid_real_source_file = "../data/grid_real.csv";
	config->grid_imag_source_file = "../data/grid_imag.csv";

	// File location to load pre-calculated w-projection kernel
	config->kernel_real_source_file = "../data/wproj_kernel_real.csv";
	config->kernel_imag_source_file = "../data/wproj_kernel_imag.csv";

	// File location to load visibility uvw coordinates
	config->visibility_source_file = "../data/el82-70.txt";
	
	// File location to store extracted visibilities
	config->visibility_dest_file = "../data/visibility_dest_file.csv";
}

void execute_degridding(Config *config, Complex *grid, Visibility *vis_uvw, Complex *vis_intensities, Complex *kernel, int num_visibilities)
{
	double uv_scale = config->uv_scale;
	int grid_size = config->grid_size;
	int half_grid_size = grid_size / 2;
	int kernel_size = config->kernel_size;
	int half_kernel_size = (kernel_size - 1) / 2; // 4
	int oversampling = config->oversampling;
	
	Visibility current_vis;
	int vis_grid_u_center = 0;
	int vis_grid_v_center = 0;
	
	int grid_u_start = 0;
	int grid_v_start = 0;
	int grid_u_end = 0;
	int grid_v_end = 0;
	int grid_index = 0;
	
	int kernel_u = 0;
	int kernel_v = 0;
	int kernel_u_offset = 0;
	int kernel_v_offset = 0;
	int kernel_index = 0;
	
	Complex current_grid_point;
	Complex current_kernel_point;
	Complex grid_kernel_product;
	Complex predicted_visibility;
	
	for(int vis_index = 0; vis_index < num_visibilities; ++vis_index)
	{
		predicted_visibility = (Complex) {.real = 0.0, .imag = 0.0};
		current_vis = vis_uvw[vis_index];
		// Calculate the central grid point of current visibility
		vis_grid_u_center = round(current_vis.u * uv_scale) + half_grid_size;
		vis_grid_v_center = round(current_vis.v * uv_scale) + half_grid_size;
		// Calculate the starting indices for the convolution kernel
		grid_u_start = vis_grid_u_center - half_kernel_size;
		grid_u_end = vis_grid_u_center + half_kernel_size;
		grid_v_start = vis_grid_v_center - half_kernel_size;
		grid_v_end = vis_grid_v_center + half_kernel_size;
		
		kernel_u_offset = (int) round((current_vis.u - (int) current_vis.u) * oversampling);
		kernel_v_offset = (int) round((current_vis.v - (int) current_vis.v) * oversampling);
		//printf("U offset: %d, V offset: %d\n", kernel_u_offset, kernel_v_offset);
		kernel_u = -half_kernel_size * oversampling + kernel_u_offset;
		kernel_v = -half_kernel_size * oversampling + kernel_v_offset;
		
		// Iterate over grid, extracting convolved values
		for(int grid_v = grid_v_start; grid_v <= grid_v_end; ++grid_v, kernel_v += oversampling)
		{	
			for(int grid_u = grid_u_start; grid_u <= grid_u_end; ++grid_u, kernel_u += oversampling)
			{
				// Get grid point
				grid_index = grid_v * grid_size + grid_u;
				current_grid_point = grid[grid_index];				
				
				// Get kernel sample
				kernel_index = abs(kernel_v * half_kernel_size) * 5 + abs(kernel_u);
				
				//printf("Grid R/C: %d, %d = %f+%f\n", grid_v, grid_u, current_grid_point.real, current_grid_point.imag);
				current_kernel_point = kernel[kernel_index];
				
				// printf("Colonel V/U (index): %d, %d, %d\n", kernel_v, kernel_u, kernel_index);
				// Calculate complex product
				grid_kernel_product = complex_multiply(current_grid_point, current_kernel_point);
				// Add complex product to predicted visibility
				predicted_visibility.real += grid_kernel_product.real;
				predicted_visibility.imag += grid_kernel_product.imag;  
			}
			
			// Reset kernel index
			kernel_u = -half_kernel_size * oversampling + kernel_u_offset;
		}
		
		//printf("Predicted Vis: %f+%f\n", predicted_visibility.real, predicted_visibility.imag);
		vis_intensities[vis_index] = predicted_visibility;
	}
}

void save_visibilities(Config *config, Visibility *vis_uvw, Complex *vis_intensity)
{
	FILE *vis_file = fopen(config->visibility_dest_file, "w");
	
	if(vis_file == NULL)
	{
		printf("Unable to open file...\n");
		return; // unsuccessfully saved visibility data
	}
	
	// Define the number of processed visibilities
	fprintf(vis_file, "%d\n", config->num_visibilities);
	
	double meters_to_wavelengths = config->frequency_hz / C;
	Visibility current_vis;
	Complex current_intensity;
	
	for(int vis_index = 0; vis_index < config->num_visibilities; ++vis_index)
	{
		current_vis = vis_uvw[vis_index];
		current_intensity = vis_intensity[vis_index];
		
		current_vis.u /= meters_to_wavelengths;
		current_vis.v /= meters_to_wavelengths;
		current_vis.w /= meters_to_wavelengths;
		
		if(config->right_ascension)
			current_vis.u *= -1.0;
		
		// u, v, w, vis(real), vis(imag), weighting
		fprintf(vis_file, "%f %f %f %f %f %f\n", 
			current_vis.u,
			current_vis.v,
			current_vis.w,
			current_intensity.real,
			current_intensity.imag,
			1.0); // static weight (for now)
	}
	
	fclose(vis_file);
}

bool load_kernel(Config *config, Complex *kernel)
{
	FILE *kernel_real_file = fopen(config->kernel_real_source_file, "r");
	FILE *kernel_imag_file = fopen(config->kernel_imag_source_file, "r");
	
	if(kernel_real_file == NULL || kernel_imag_file == NULL)
	{
		printf("Unable to open kernel source files...\n");
		if(kernel_real_file != NULL) fclose(kernel_real_file);
		if(kernel_imag_file != NULL) fclose(kernel_imag_file);
		return false; // unsuccessfully loaded data
	}
	
	int kernel_size = config->kernel_size;
	int half_kernel_size = (kernel_size / 2) + 1;
	int oversampling = config->oversampling;
	int half_kernel_oversampled = half_kernel_size * oversampling;
	int kernel_index = 0;
	double kernel_real = 0.0;
	double kernel_imag = 0.0;
	
	for(int row_index = 0; row_index < half_kernel_oversampled; ++row_index)
	{
		for(int col_index = 0; col_index < half_kernel_oversampled; ++col_index)
		{
			fscanf(kernel_real_file, "%lf ", &kernel_real);
			fscanf(kernel_imag_file, "%lf ", &kernel_imag);
			kernel_index = row_index * half_kernel_oversampled + col_index;
			kernel[kernel_index] = (Complex) {.real = kernel_real, .imag = kernel_imag};
		}
	}
	
	fclose(kernel_real_file);
	fclose(kernel_real_file);
	fclose(kernel_imag_file);
	return true;
}

bool load_grid(Config *config, Complex *grid)
{
	FILE *grid_real_file = fopen(config->grid_real_source_file, "r");
	FILE *grid_imag_file = fopen(config->grid_imag_source_file, "r");
	
	if(grid_real_file == NULL || grid_imag_file == NULL)
	{
		printf("Unable to open grid files...\n");
		if(grid_real_file != NULL) fclose(grid_real_file);
		if(grid_imag_file != NULL) fclose(grid_imag_file);
		return false; // unsuccessfully loaded data
	}
	
	int grid_size = config->grid_size;
	int grid_index = 0;
	double grid_real = 0.0;
	double grid_imag = 0.0;
	
	for(int row_index = 0; row_index < grid_size; ++row_index)
	{
		for(int col_index = 0; col_index < grid_size; ++col_index)
		{
			fscanf(grid_real_file, "%lf ", &grid_real);
			fscanf(grid_imag_file, "%lf ", &grid_imag);
			grid_index = row_index * grid_size + col_index;
			grid[grid_index] = (Complex) {.real = grid_real, .imag = grid_imag};
		}
	}
	
	fclose(grid_real_file);
	fclose(grid_imag_file);
	return true;
}

bool load_visibilities(Config *config, Visibility **vis_uvw, Complex **vis_intensities)
{
	// Attempt to open visibility source file
	FILE *vis_file = fopen(config->visibility_source_file, "r");
	if(vis_file == NULL)
	{
		printf("Unable to open visibility file...\n");
		return false; // unsuccessfully loaded data
	}
	
	// Configure number of visibilities from file
	int num_visibilities = 0;
	fscanf(vis_file, "%d", &num_visibilities);
	num_visibilities = 30;
	config->num_visibilities = num_visibilities;
	
	// Allocate memory for incoming visibilities
	*vis_uvw = calloc(num_visibilities, sizeof(Visibility));
	*vis_intensities = calloc(num_visibilities, sizeof(Complex));
	if(!(*vis_uvw) || !(*vis_intensities))
	{
		printf("Unable to allocate memory...\n");
		return false;
	}
	
	// Load visibility uvw coordinates into memory
	double vis_u = 0.0;
	double vis_v = 0.0;
	double vis_w = 0.0;
	double vis_real = 0.0;
	double vis_imag = 0.0;
	double vis_weight = 0.0;
	double wavelength_to_meters = config->frequency_hz / C;
	for(int vis_index = 0; vis_index < num_visibilities; ++vis_index)
	{
		// Discard vis(real), vis(imag), and weighting (for now)
		fscanf(vis_file, "%lf %lf %lf %lf %lf %lf\n", &vis_u, &vis_v,
			&vis_w, &vis_real, &vis_imag, &vis_weight);
			
		(*vis_uvw)[vis_index] = (Visibility) {
			.u = vis_u * wavelength_to_meters,
			.v = vis_v * wavelength_to_meters,
			.w = 0.0 //vis_w * wavelength_to_meters
		};

		(*vis_intensities)[vis_index] = (Complex) {
			.real = vis_real,
			.imag = vis_imag
		};
		
		if(config->right_ascension)
			(*vis_uvw)[vis_index].u *= -1.0;
		
		printf("%f, %f, %f, %f, %f, %f\n", vis_u, vis_v,
			vis_w, vis_real, vis_imag, vis_weight);
		// printf("%f, %f, %f\n", (*vis_uvw)[vis_index].u, (*vis_uvw)[vis_index].v, (*vis_uvw)[vis_index].w);
	}
		
	// Clean up
	fclose(vis_file);
	return true;
}

Complex complex_multiply(Complex z1, Complex z2)
{
	Complex z3;
	z3.real = z1.real * z2.real - z1.imag * z2.imag;
	z3.imag = z1.imag * z2.real + z1.real * z2.imag;
	return z3;
}

void clean_up(Complex **grid, Visibility **vis_uvw, Complex **vis_intensities, Complex **kernel)
{
	if(*grid) 			 free(*grid);
	if(*vis_uvw) 	 	 free(*vis_uvw);
	if(*vis_intensities) free(*vis_intensities);
	if(*kernel) 		 free(*kernel);
}

/***************************************
*      UNIT TESTING FUNCTIONALITY      *
***************************************/

void unit_test_init_config(Config *config)
{
	config->grid_size = 16384;
	config->right_ascension = true;
	config->cell_size = 6.39708380288950e-6;
	config->frequency_hz = 100e6;
	config->kernel_size = 9;
	config->oversampling = 4;
	config->uv_scale = config->grid_size * config->cell_size;
	config->num_visibilities = 1;
	config->grid_real_source_file = "../data/grid_real.csv";
	config->grid_imag_source_file = "../data/grid_imag.csv";
	config->kernel_real_source_file = "../data/wproj_kernel_real.csv";
	config->kernel_imag_source_file = "../data/wproj_kernel_imag.csv";
	config->visibility_source_file = "../data/el82-70.txt";
	config->visibility_dest_file = "../data/visibility_dest_file.csv";
}

double unit_test_generate_approximate_visibilities(void)
{
	// used to invalidate the unit test
	double error = DBL_MAX;

	Config config;
	unit_test_init_config(&config);

	// Prepare required memory
	Complex *grid = (Complex*) calloc(config.grid_size * config.grid_size, sizeof(Complex));
	size_t kernel_size = pow(((config.kernel_size / 2) + 1) * config.oversampling, 2.0);
	Complex *kernel = (Complex*) calloc(kernel_size, sizeof(Complex));
	
	// Evaluate memory allocation success
	if(!grid || !kernel)
	{
		clean_up(&grid, NULL, NULL, &kernel);
		return error;
	}
	
	// Load in w-projection kernel for w == 0
	bool loaded_kernel = load_kernel(&config, kernel);
	if(!loaded_kernel)
	{
		clean_up(&grid, NULL, NULL, &kernel);
		return error;
	}
	
	bool loaded_grid = load_grid(&config, grid);
	Visibility *vis_uvw = NULL;
	Complex *vis_intensities = NULL;
	bool loaded_vis = load_visibilities(&config, &vis_uvw, &vis_intensities);
	if(!loaded_grid || !loaded_vis || !vis_uvw)
	{
		clean_up(&grid, &vis_uvw, &vis_intensities, &kernel);
		return error;
	}

	Visibility test_visibility_uvw;
	Complex test_visibility;
	Visibility approx_visibility_uvw[1]; // testing one at a time
	Complex approx_visibility[1];

	double difference = 0.0;

	for(int vis_index = 0; vis_index < config.num_visibilities; ++vis_index)
	{
		test_visibility_uvw = vis_uvw[vis_index];
		test_visibility = vis_intensities[vis_index];

		approx_visibility_uvw[0] = (Visibility) {
			.u = test_visibility_uvw.u,
			.v = test_visibility_uvw.v,
			.w = test_visibility_uvw.w,
		};

		approx_visibility[0] = (Complex) {
			.real = 0.0,
			.imag = 0.0
		};
		
		execute_degridding(&config, grid, approx_visibility_uvw, approx_visibility, kernel, 1);

		double current_difference = sqrt(pow(approx_visibility[0].real - test_visibility.real, 2.0)
	  		+ pow(approx_visibility[0].imag - test_visibility.imag, 2.0));

		if(current_difference > difference)
			difference = current_difference;
	}

	return difference;	
}