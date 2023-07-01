/**   coordnate surface
*/
struct crd_surface
{
  crd_curve curve1;
  crd_curve curve2;
  // Default constructor
  crd_surface() {}
  
  // Constructor with curves
  crd_surface(const crd_curve& c1, const crd_curve& c2) {
      curve1 = c1;
      curve2 = c2;
  }
  
  // Project a curve onto the surface
  crd_curve project_curve(const crd_curve& curve) {
      crd_curve result;
      result.crdlist.reserve(curve.crdlist.size());
  
      for (const auto& crd : curve.crdlist) {
          result.add(project_coord(crd));
      }
  
      return result;
  }
  
  // Project a coordinate onto the surface
  coord3 project_coord(const coord3& crd) {
      coord3 result;
  
      // Perform projection calculation here
  
      return result;
  }
};
