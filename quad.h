#ifndef QUAD_H
#define QUAD_H

#include "hittable.h"

class quad : public hittable {
public:
	quad(const point3& Q, const vec3& u, const vec3& v, shared_ptr<material> mat) 
		: Q(Q), u(u), v(v), mat(mat) 
	{
		set_bounding_box();
		auto n = cross(u, v);
		normal = unit_vector(n);
		D = dot(normal, Q);
		w = n / dot(n, n);

		area = n.length();
	}
	virtual void set_bounding_box() {
		auto box1 = aabb(Q, Q + u + v);
		auto box2 = aabb(Q + u, Q + v);
		bbox = aabb(box1, box2);
	}

	aabb bounding_box() const override { return bbox; }

	bool hit(const ray& r, interval ray_t, hit_record& rec) const override {
		auto ndotd = dot(normal, r.direction());

		// No hit if the ray is parallel to the plane.
		if (std::fabs(ndotd) < 1e-8) return false;

		auto t = (D - dot(normal, r.origin())) / ndotd;

		// Return false if the hit point parameter t is outside the ray interval.
		if (!ray_t.contains(t)) return false;

		auto I = r.at(t);
		auto p = I - Q;
		double a = dot(w, cross(p, v));
		double b = dot(w, cross(u, p));

		// Return false if intersection is not interior to the primitive.
		if (!is_interior(a, b, rec)) return false;
		
		rec.t = t;
		rec.p = I;
		rec.mat = mat;
		rec.set_face_normal(r, normal);

		return true; 
	}

	virtual bool is_interior(double a, double b, hit_record& rec) const {
		interval unit_interval = interval(0, 1);
		if (!unit_interval.contains(a) || !unit_interval.contains(b)) return false;
		rec.u = a;
		rec.v = b;
		return true;
	}

	double pdf_value(const point3& origin, const vec3& direction) const override {
		hit_record rec;
		if (!this->hit(ray(origin, direction), interval(0.001, infinity), rec))
			return 0;

		auto distance_squared = rec.t * rec.t * direction.length_squared();
		auto cosine = std::fabs(dot(direction, rec.normal) / direction.length());

		return distance_squared / (cosine * area);
	}

	vec3 random(const point3& origin) const override {
		auto p = Q + (random_double() * u) + (random_double() * v);
		return p - origin;
	}

private:
	point3 Q;
	vec3 u, v;
	vec3 w; 
	shared_ptr<material> mat;
	aabb bbox;
	vec3 normal;
	double D;
	double area;
};

inline shared_ptr<hittable_list> box(const point3& a, const point3& b, shared_ptr<material> mat)
{
	// Returns the 3D box (six sides) that contains the two opposite vertices a & b.

	auto sides = make_shared<hittable_list>();

	// Construct the two opposite vertices with the minimum and maximum coordinates.
	auto min = point3(std::fmin(a.x(), b.x()), std::fmin(a.y(), b.y()), std::fmin(a.z(), b.z()));
	auto max = point3(std::fmax(a.x(), b.x()), std::fmax(a.y(), b.y()), std::fmax(a.z(), b.z()));

	auto dx = vec3(max.x() - min.x(), 0, 0);
	auto dy = vec3(0, max.y() - min.y(), 0);
	auto dz = vec3(0, 0, max.z() - min.z());

	sides->add(make_shared<quad>(point3(min.x(), min.y(), max.z()), dx, dy, mat)); // front
	sides->add(make_shared<quad>(point3(max.x(), min.y(), max.z()), -dz, dy, mat)); // right
	sides->add(make_shared<quad>(point3(max.x(), min.y(), min.z()), -dx, dy, mat)); // back
	sides->add(make_shared<quad>(point3(min.x(), min.y(), min.z()), dz, dy, mat)); // left
	sides->add(make_shared<quad>(point3(min.x(), max.y(), max.z()), dx, -dz, mat)); // top
	sides->add(make_shared<quad>(point3(min.x(), min.y(), min.z()), dx, dz, mat)); // bottom

	return sides;
}

#endif // !QUAD_H
