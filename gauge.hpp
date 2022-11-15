/**
*							�淶�任
*				�������Ը�ĳ����λ����λ����z����ϣ���Ϊ"��λͶӰ":
*						P(v,z) = v * exp(z)
*				Ȼ����ĳ��Ƶ��n��ͶӰ����Ϊ"Ƶ��ͶӰ":
*					G(v,z,n) = (P(v,z))^n = (v*exp(z))^n
*				���������Ƶ����ͣ�
*					SUM(N, G(v,z,n)) = SUM(N, (v*exp(z))^n)
* 
*				�淶�任����������λ�任��Ȼ��ϳɲ�ͬƵ������ЧӦ
*/
scope gauge_math
{
	// ������ĳ����λ"ͶӰ"�任��ע���ǵ�λ����!)
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

	// ����Ƶ�����
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