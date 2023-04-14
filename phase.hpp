/*********************************************************************
*					����λ / ��λ����ϵ��
*
*				��һ�ֶ�����ʱ��ά���ϵ�����ϵ
*				������ϵ�������Ƶ���ʽ
*/

// *******************************************************************
//  th
// C    Phase
// *******************************************************************
struct phase
{
	static const phase ONE;
	real angle = 0;

	phase() {}
	phase(real _angle)
	{
		angle = _angle;
	}
	void rot(real _angle)
	{
		angle = _angle;
	}
	phase operator + (const phase& c) const
	{
		phase rc;
		rc.angle = angle + c.angle;
		return rc;
	}
	phase operator - (const phase& c) const
	{
		phase rc;
		rc.angle = angle - c.angle;
		return rc;
	}
	friend vec3 operator * (crvec p, const phase& c)
	{
		number_math::complex cx;
		cx.from_ra(1, c.angle);
		return p * cx.x;
	}
	phase operator * (const phase& c) const
	{
		phase rc;
		rc.angle = angle + c.angle;
		return rc;
	}
	friend crvec operator / (crvec p, const phase& c)
	{
		number_math::complex cx(1, -c.angle);
		return p * cx.x;
	}
	phase operator / (const phase& c) const
	{
		phase rc;
		rc.angle = angle - c.angle;
		return rc;
	}
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	phase reversed() const
	{
		return ONE / (*this);
	}
	static phase gradcoord(const phase& c1, const phase& c2)
	{
		return c1.reversed() * c2;
	}
	real getreal()
	{
		return cos(angle);
	}
	void dump() const
	{
		PRINTV(angle);
	}
};
const phase phase::ONE = phase();