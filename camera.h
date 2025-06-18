#ifndef CAMERA_H
#define CAMERA_H

#include "hittable.h"
#include "material.h"
#include "pdf.h"
#include <vector>
#include <omp.h>

#include <chrono>

class camera {
public:
	double aspect_ratio = 1.0;			// Ratio of image width over height
	double vfov = 90;					// Vertical view angle (field of view)
	int    image_width = 100;			// Rendered image width in pixel count
	int    samples_per_pixel = 10;		// Count of random samples for each pixel
	int    max_depth = 10;				// Maximum number of ray bounces into scene
	color  background;					// Scene background color

	point3 lookfrom = point3(0, 0, 0);	// Point camera is looking from
	point3 lookat = point3(0, 0, -1);	// Point camera is looking at
	vec3   vup = vec3(0, 1, 0);			// Camera-relative "up" direction

	double defocus_angle = 0;			// Variation angle of rays through each pixel
	double focus_dist = 10;				// Distance from camera lookfrom point to plane of perfect focus


	void render(const hittable& world, const hittable& lights) {
		initialize();
		auto start = std::chrono::steady_clock::now();
		std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";
		std::vector<std::vector<color>> img(image_width, std::vector<color>(image_height, color(0, 0, 0)));
		#pragma omp parallel shared(img)
		{
		int id = omp_get_thread_num();
		for (int j = 0; j < image_height; j++) {
			if (id == 0) std::clog << "\rProgress: " << (j * 100 / image_height) << '%' << std::flush;
			#pragma omp for
			for (int i = 0; i < image_width; i++) {
				color pixel_color(0, 0, 0);
				for (int s_j = 0; s_j < sqrt_spp; s_j++) {
					for (int s_i = 0; s_i < sqrt_spp; s_i++) {
						ray r = get_ray(i, j, s_i, s_j);
						pixel_color += ray_color(r, max_depth, world, lights);
					}
				}
				img[i][j] = pixel_color;
			}
		}
		}
		for (int j = 0; j < image_height; j++) {
			//std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
			for (int i = 0; i < image_width; i++) {
				write_color(std::cout, pixel_samples_scale * img[i][j]);
			}
		}
		auto end = std::chrono::steady_clock::now();
		auto diff = duration_cast<std::chrono::milliseconds>(end - start);
		std::clog << "\rDone.                        \n";
		std::clog << "\r" << diff.count() / 1000 << "s \n";
	}

private:
	int    image_height;   // Rendered image height
	double pixel_samples_scale;  // Color scale factor for a sum of pixel samples
	int    sqrt_spp;             // Square root of number of samples per pixel
	double recip_sqrt_spp;       // 1 / sqrt_spp
	point3 center;         // Camera center
	point3 pixel00_loc;    // Location of pixel 0, 0
	vec3   pixel_delta_u;  // Offset to pixel to the right
	vec3   pixel_delta_v;  // Offset to pixel below
	vec3   u, v, w;		   // Camera frame basis vectors
	vec3   defocus_disk_u; // Defocus disk horizontal radius
	vec3   defocus_disk_v; // Defocus disk vertical radius

	void initialize() {
		// Calculate the image height, and ensure that it's at least 1.
		image_height = int(image_width / aspect_ratio);
		image_height = (image_height < 1) ? 1 : image_height;

		sqrt_spp = int(std::sqrt(samples_per_pixel));
		pixel_samples_scale = 1.0 / (sqrt_spp * sqrt_spp);
		recip_sqrt_spp = 1.0 / sqrt_spp;

		center = lookfrom;

		// Viewport dimensions
		auto theta = degrees_to_radians(vfov);
		auto h = std::tan(theta / 2);
		auto viewport_height = 2 * h * focus_dist;
		auto viewport_width = viewport_height * (double(image_width) / image_height);


		// Calculate camera basis vector uvw
		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);

		// Calculate the vectors across the horizontal and down the vertical viewport edges.
		auto viewport_u = viewport_width * u;	// Vector across viewport horizontal edge
		auto viewport_v = viewport_height * -v;	// Vector down viewport vertical edge

		// Calculate the horizontal and vertical delta vectors from pixel to pixel.
		pixel_delta_u = viewport_u / image_width;
		pixel_delta_v = viewport_v / image_height;

		// Calculate the location of the upper left pixel.
		auto viewport_upper_left = center - (focus_dist * w) - viewport_u / 2 - viewport_v / 2;
		pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

		// Calculate defocus disk base vectors
		auto defocus_radius = focus_dist * std::tan(degrees_to_radians(defocus_angle/2.0));
		defocus_disk_u = u * defocus_radius;
		defocus_disk_v = v * defocus_radius;
	}

	ray get_ray(int i, int j, int s_i, int s_j) {
		// Construct a camera ray originating from the defocus disk and directed at a randomly
		// sampled point around the pixel location i, j.
		auto offset = sample_square_stratified(s_i, s_j);

		auto pixel_sample = pixel00_loc 
						  + ((i + offset.x()) * pixel_delta_u) 
						  + ((j + offset.y()) * pixel_delta_v);

		auto ray_origin = (defocus_angle<=0) ? center : defocus_disk_sample();
		auto ray_direction = pixel_sample - ray_origin;

		return ray(ray_origin, ray_direction);
	}

	vec3 sample_square() const {
		// Returns the vector to a random point in the [-.5,-.5]-[+.5,+.5] unit square.
		return vec3(random_double() - 0.5, random_double() - 0.5, 0.0);
	}

	vec3 sample_square_stratified(int s_i, int s_j) const {
		// Returns the vector to a random point in the square sub-pixel specified by grid
		// indices s_i and s_j, for an idealized unit square pixel [-.5,-.5] to [+.5,+.5].

		auto px = ((s_i + random_double()) * recip_sqrt_spp) - 0.5;
		auto py = ((s_j + random_double()) * recip_sqrt_spp) - 0.5;

		return vec3(px, py, 0);
	}

	vec3 defocus_disk_sample() const {
		// Returns a random point in the camera defocus disk.
		auto p = random_in_unit_disk();
		return center + (p[0] * defocus_disk_u + p[1] * defocus_disk_v);
	}

	color ray_color(const ray& r, double depth, const hittable& world, const hittable& lights) {
		if (depth <= 0) {
			return color(0, 0, 0);
		}

		hit_record rec;
		if (!world.hit(r, interval(0.001, infinity), rec)) {
			return background;
		}

		scatter_record srec;
		color emitted_light = rec.mat->emitted(r, rec, rec.u, rec.v, rec.p);

		if (!rec.mat->scatter(r, rec, srec))
			return emitted_light;

		if (srec.skip_pdf) {
			return srec.attenuation * ray_color(srec.skip_pdf_ray, depth - 1, world, lights);
		}

		auto light_pdf_ptr = make_shared<hittable_pdf>(lights, rec.p);
		mixture_pdf mixed_pdf(light_pdf_ptr, srec.pdf_ptr);

		ray scattered = ray(rec.p, mixed_pdf.generate());
		auto pdf_value = mixed_pdf.value(scattered.direction());

		double scatter_pdf = rec.mat->scattering_pdf(r, rec, scattered);

		color scattered_light 
			= srec.attenuation * scatter_pdf * ray_color(scattered, depth - 1, world, lights) / pdf_value;
		return scattered_light + emitted_light;

		//vec3 unit_direction = unit_vector(r.direction());
		//auto a = 0.5 * (unit_direction.y() + 1.0);
		//return (1.0 - a) * color(1.0, 1.0, 1.0) + a * color(0.5, 0.7, 1.0);
	}
};

#endif
