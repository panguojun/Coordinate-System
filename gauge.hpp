/**
*							规范变换
*				向量可以跟某个相位（单位复数z）结合，称为"相位投影":
*						P(v,z) = v * exp(z)
*				然后在某个频率n上投影，称为"频率投影":
*					G(v,z,n) = (P(v,z))^n = (v*exp(z))^n
*				最后在所有频率求和：
*					SUM(N, G(v,z,n)) = SUM(N, (v*exp(z))^n)
* 
*				规范变换本质上是相位变换，然后合成不同频率总体效应
*/
scope gauge_math
{
	// 向量向某个相位"投影"变换（注意是单位复数!)
	compx phase(const compx& v, const compx& z)
	{
		compx nz = exp(z);
		return v * nz;
	}
	quaternion phase(const quaternion& v, const quaternion& q)
	{
		quaternion nq = exp(q);
		return v * nq;
	}

	// 各个频率求和
	quaternion freqsum(const quaternion& v, int N)
	{
		quaternion sum;
		for (int i = 1; i <= N; i++)
		{
			quaternion q = v ^ i;
			sum.x += q.x;
			sum.y += q.y;
			sum.z += q.z;
			sum.w += q.w;
		}
		return sum;
	}
	vec3 test(crvec v)
	{
		quat q(1.0f,0,1,0);
		quat nv = phase(quat(0, v), q);
		nv = freqsum(nv, 8);
		nv.normalize();
		PRINTV3(nv.xyz())
	}
}