#ifndef ONB_H
#define ONB_H 

class onb {
public:
	onb(const vec3& n) {
		basis[2] = unit_vector(n);
		vec3 a = (std::fabs(basis[2].x()) > 0.9) ? vec3(0, 1, 0) : vec3(1, 0, 0);
		basis[1] = unit_vector(cross(basis[2], a));
		basis[0] = cross(basis[2], basis[1]);
	}

	const vec3& u() const { return basis[0]; }
	const vec3& v() const { return basis[1]; }
	const vec3& w() const { return basis[2]; }

	vec3 transform(const vec3& v) const {
		return (v[0] * basis[0]) + (v[1] * basis[1]) + (v[2] * basis[2]);
	}
private:
	vec3 basis[3];
};

#endif // !ONB_H

