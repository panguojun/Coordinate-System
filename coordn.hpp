/*********************************************************************
*				      【坐标系】
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
* 	坐标系类是我单独封装，用于简化坐标变换，衍生出许多算法，能解决一些
* 	坐标系变换相关的问题。
* 	坐标系的运算跟李群很相似。
*	坐标系由三个部分组成：C = M(位置） + S（缩放） * R（旋转）
*
*  *  *  *  *  *  *  *  *  *  详解  *  *  *  *  *  *  *  *  *  *  *  *
*	坐标系变换分为投影（/), 平移（^), 还原（*）三个步骤，
*	坐标系本体符号 C，坐标系在单位距离上的变换梯度可以写成：
*			G = C2 / C1 - I,
*	除法分为两种（左，右）：
*			oper(/)  = C1 * C2^-1
*			oper(\)  = C1^-1 * C2
*	在曲面上可以进行微分几何运算：
*	定义一个内禀坐标系(假设它是平直空间，向量可以随意移动而不变)下V,在弯
*	曲坐标系下观察V，不同点上V是不同的，故而坐标系跟位置有关，取相邻两点
*	（1),(2)点处有向量V1,V2，对应坐标系C1,C2，那么：
*			V = V1 * C1 = V2 * C2 =>
*			V2 = V1 * C1 / C2, 令 G12 = C1 / C2 =>
*			V2 = V1 * G12
*
*	可以使用坐标系计算空间曲率，在u,v坐标系下黎曼曲率张量为：
*			Ruv = Gu*Gv - Gv*Gu - G[u,v]
*			其中：Gu = C2 / C1
*			     W = [u, v] : 李括号运算
*			     G[u,v] = Gu^Wu * Gv^Wv
*/

//#define	Parallel_Projection		 // 非正交坐标系下平行投影
// *******************************************************************
//  |_
// UC     ND Rotation Coordinate System
// *******************************************************************
struct ucoordn {
    static const ucoord ZERO;
    static const ucoord ONE;

    vector<vecn> dirs;

    // 乘法：在坐标系下定义一个向量，或者向量向父空间还原
    friend vecn operator * (const vecn& p, const ucoord& c)
    {
        vecn ret(p.dim());
        for (int i = 0; i < c.dirs.size(); i++) {
            ret += c.dirs[i] * p[i];
        }
        return ret;
    }

    ucoord operator * (const ucoord& c) const
    {// Cchild * Cparent * ...
        ucoord rc;
        rc.dirs.resize(dirs.size());
        for (int i = 0; i < dirs.size(); i++) {
            rc.dirs[i] = dirs[i] * c.dirs[i];
        }
        return rc;
    }

    // 除法：向量向坐标系投影
    friend vecn operator / (const vecn& v, const ucoord& c)
    {
#ifdef Parallel_Projection
        {// 对于非正交情况
            vecn proj(v.dim());
            for (int i = 0; i < c.dirs.size(); i++) {
                vecn axis = c.dirs[i];
                proj += pl_prj(v, axis, c.dirs - axis) * axis;
            }
            proj /= c.dirs.norms();
            return proj;
        }
#endif
        vecn res(v.dim());
        for (int i = 0; i < c.dirs.size(); i++) {
            res[i] = vecn::dot(v, c.dirs[i]);
        }
        return res;
    }
};
