// ****************************************************************************
// 管道流体的渲染版本
// ****************************************************************************
namespace PIPE_FLOW_FIELD_RENDER
{
	float unit_scale = 2.0f;
	float border_deta = 1;

	struct center_t{vec3 p;float r;};
	using CENTERLIST = std::vector<center_t>;
	std::vector<CENTERLIST> centerlinestack;

	boundingboxi aabb;
	vec3 cur_center1, cur_center2;
	real cur_radius = 0;

	// -----------------------------------------------
	real field2tris(submesh& sm, int i, int j, int k, real scale, std::function<real(vec3& p)> field)
	{
		GRIDCELL grid;
		grid.p[0] = vec3(i, j, k);
		grid.p[1] = vec3(i + 1, j, k);
		grid.p[2] = vec3(i + 1, j, k + 1);
		grid.p[3] = vec3(i, j, k + 1);
		grid.p[4] = vec3(i, j + 1, k);
		grid.p[5] = vec3(i + 1, j + 1, k);
		grid.p[6] = vec3(i + 1, j + 1, k + 1);
		grid.p[7] = vec3(i, j + 1, k + 1);

		grid.val[0] = field(grid.p[0]);
		grid.val[1] = field(grid.p[1]);
		grid.val[2] = field(grid.p[2]);
		grid.val[3] = field(grid.p[3]);
		grid.val[4] = field(grid.p[4]);
		grid.val[5] = field(grid.p[5]);
		grid.val[6] = field(grid.p[6]);
		grid.val[7] = field(grid.p[7]);

		TRIANGLE tri[5];
		int num = Polygonise(grid, 1.0, tri);
		num = _MIN(num, 5);

		gsearchcomvertex = true;
		gcommonvertex = true;

		if (cur_center1 == cur_center2)
		{
			for (int i = 0; i < num; i++)
			{
				vertex v1, v2, v3;
				v1.p = tri[i].p[0] * scale;
				v2.p = tri[i].p[1] * scale;
				v3.p = tri[i].p[2] * scale;

				v1.n = (v1.p - cur_center1).normcopy();
				v2.n = (v2.p - cur_center1).normcopy();
				v3.n = (v3.p - cur_center1).normcopy();

				binvnorm = 1;
				triang(v1, v2, v3);
			}
		}
		else
		{
			vec3 vc12 = cur_center2 - cur_center1; vc12.norm();
			for (int i = 0; i < num; i++)
			{
				vertex v1, v2, v3;
				v1.p = tri[i].p[0] * scale;
				v2.p = tri[i].p[1] * scale;
				v3.p = tri[i].p[2] * scale;

				{
					v1.n = (v1.p - cur_center1);
					v1.n = v1.n - vc12 * v1.n.dot(vc12); v1.n.norm();

					v2.n = (v2.p - cur_center1);
					v2.n = v2.n - vc12 * v2.n.dot(vc12); v2.n.norm();

					v3.n = (v3.p - cur_center1);
					v3.n = v3.n - vc12 * v3.n.dot(vc12); v3.n.norm();
				}

				binvnorm = 1;
				triang(v1, v2, v3);
			}
		}
		return grid.val[0];
	}
	// -----------------------------------------------
	void updateaabb(const point3_t& p)
	{
		if (p.x < aabb.a.x + border_deta)
			aabb.a.x = p.x - border_deta;
		if (p.y < aabb.a.y + border_deta)
			aabb.a.y = p.y - border_deta;
		if (p.z < aabb.a.z + border_deta)
			aabb.a.z = p.z - border_deta;

		if (p.x > aabb.b.x - border_deta)
			aabb.b.x = p.x + border_deta;
		if (p.y > aabb.b.y - border_deta)
			aabb.b.y = p.y + border_deta;
		if (p.z > aabb.b.z - border_deta)
			aabb.b.z = p.z + border_deta;
	}
	// -----------------------------------------------
	real getdis_onpoly(const CENTERLIST& poly, crvec p, real minddis)
	{
		int currentpoint = -1;
		real alpha = 0.0f;

		if(poly.size() == 1)
		{
			real ddis = (p - poly[0].p).sqrlen() / (poly[0].r * poly[0].r);
			if (ddis < minddis)
			{
				currentpoint = 0;
				minddis = ddis;
			}
		}
		else
		{
			for (int i = 1; i < poly.size(); i++)
			{
				crvec p1 = poly[i - 1].p;
				crvec p2 = poly[i].p;
				{
					vec3 v12 = (p2 - p1);
					real d12 = v12.len();
					v12 /= d12;
					vec3 vp = p - p1;
					real dot = vp.dot(v12);
					if (dot > 0 && dot <= d12)
					{
						real radius = blend(poly[i - 1].r, poly[i].r, dot / d12);
						real ddis = (vp - v12 * dot).sqrlen() / (radius * radius);
						if (ddis < minddis)
						{
							alpha = dot / d12;
							currentpoint = i;
							cur_center1 = poly[currentpoint - 1].p;
							cur_center2 = poly[currentpoint].p;
							cur_radius = radius;
							minddis = ddis;
						}
					}
				}
				{
					real radius = poly[i - 1].r;
					real ddis = (p - p1).sqrlen() / (radius * radius);
					if (ddis < minddis)
					{
						alpha = 0;
						currentpoint = i - 1;
						cur_radius = radius;
						cur_center1 = cur_center2 = poly[currentpoint].p;
						minddis = ddis;
					}
				}
				{
					real radius = poly[i].r;
					real ddis = (p - p2).sqrlen() / (radius * radius);
					if (ddis < minddis)
					{
						alpha = 0;
						currentpoint = i;
						cur_radius = radius;
						cur_center1 = cur_center2 = poly[currentpoint].p;
						minddis = ddis;
					}
				}
			}
		}

		return minddis;
	}
	// -----------------------------------------------
	void _render(submesh& sm)
	{
		for (int i = aabb.a.x; i < aabb.b.x; i++)
			for (int j = aabb.a.y; j < aabb.b.y; j++)
				for (int k = aabb.a.z; k < aabb.b.z; k++)
				{
					field2tris(sm, i, j, k, unit_scale / 5.0f,
						[i, j, k](crvec p)->float
						{
							float mdd = 1e10;
							for (auto& it : centerlinestack)
							{
								float dd = getdis_onpoly(it, p, mdd);
								if (dd < mdd)
									mdd = dd;
							}
							return (mdd);
						}
					);
				}
	}
	// -----------------------------------------------
	void add_centerpoint(crvec p, float r, int line)
	{
		r /= unit_scale;
		if (r > border_deta)
			border_deta = r + 1;

		vec3 tp = p / unit_scale;
		updateaabb(point3_t(tp.x, tp.y, tp.z));

		while (line >= centerlinestack.size())
			centerlinestack.push_back(CENTERLIST());
		centerlinestack[line].push_back({ tp, r });

	}
	// -----------------------------------------------
	void render_pipe(submesh& sm)
	{
		if (centerlinestack.empty())
			return;

		PRINTVEC3((aabb.b - aabb.a));

		_render(sm);
	}
	// -----------------------------------------------
	void clear()
	{
		centerlinestack.clear();
	}

	// -----------------------------------------------
	#include "artwork/cellflow.hpp"
	// -----------------------------------------------
	void test()
	{
		setup();

		render_pipe(SUBMESH);

		clear();
	}
};

void flow()
{
	PIPE_FLOW_FIELD_RENDER::test();
}
