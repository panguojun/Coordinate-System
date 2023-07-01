struct crd_curve 
{

    std::vector<coord3> crdlist;



    // Default constructor

    crd_curve() {}



    // Constructor with coordinates

    crd_curve(const std::vector<coord3>& crds) {

        crdlist = crds;

    }



    // Add a coordinate to the curve

    void add(const coord3& crd) {

        crdlist.push_back(crd);

    }



    // Concatenate two curves

    crd_curve operator + (const crd_curve& other) const {

        crd_curve result = *this;

        for (const auto& crd : other.crdlist) {

            result.add(crd);

        }

        return result;

    }



    // Subtract another curve from this curve

    crd_curve operator - (const crd_curve& other) const {

        crd_curve result = *this;

        result.crdlist.erase(result.crdlist.end() - other.crdlist.size(), result.crdlist.end());

        return result;

    }



    // Multiply this curve with another curve

    crd_curve operator * (const crd_curve& other) const {

        crd_curve result;

        for (const auto& crd1 : crdlist) {

            for (const auto& crd2 : other.crdlist) {

                result.add(crd1 * crd2);

            }

        }

        return result;

    }



    // Divide this curve by another curve

    crd_curve operator / (const crd_curve& other) const {

        crd_curve result;

        for (const auto& crd1 : crdlist) {

            for (const auto& crd2 : other.crdlist) {

                result.add(crd1 / crd2);

            }

        }

        return result;

    }



    // Get coordinate at a given parameter t

    vec3 getcrd(float t) const {

        if (crdlist.empty()) {

            return vec3::ZERO;

        }

        if (t <= 0) {

            return crdlist.front().pos();

        }

        if (t >= 1) {

            return crdlist.back().pos();

        }

        float step = 1.0 / (crdlist.size() - 1);

        int index = t / step;

        float sub_t = (t - index * step) / step;

        return crdlist[index].pos() * (1 - sub_t) + crdlist[index + 1].pos() * sub_t;

    }



    // Convert the curve to a polyline with given step size

    std::vector<vec3> to_poly(float step) const {

        std::vector<vec3> polyline;

        for (float t = 0; t <= 1; t += step) {

            polyline.push_back(getcrd(t));

        }

        return polyline;

    }



    // Project the curve to a given coordinate

    void project_to(const coord3& coord) {

        for (auto& crd : crdlist) {

            crd = crd / coord;

        }

    }



    // Project the curve to another curve by dividing the coordinates

    void project_to(const crd_curve& curve) {

        for (auto& crd : crdlist) {

            crd = crd / curve.crdlist.front();

        }

    }

};
