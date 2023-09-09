/**    【插值单元对象】
*      通过乘法可以计算插值
*      使用插值单元可以代替线性函数
*/
template<typename T>
struct unit
{
  static const T ZERO;
  static const T ONE;
  T a, b;
  T operator *(T t)
  {
    return lerp(a, b, t);
  }
};
template<typename T>
const T unit<T>::ZERO = T(0);
template<typename T>
const T unit<T>::ONE = T(1);
