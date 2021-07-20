// ****************************************************************************
// 管道流体的渲染版本
// ****************************************************************************
namespace PIPE_FLOW_FIELD_RENDER
{
	float mesh_scale = 0.5f;
	int border_deta = 1;

	struct center_t
	{
		vec3 p;
		float r;
	};
	using CENTERLIST = std::vector<center_t>;
	std::vector<CENTERLIST> centerlinestack;
	boundingboxi aabb;
	vec3 cur_centerpoint;
	real cur_radius = 0;
	// -----------------------------------------------
	inline real field2tris(submesh& sm, int i, int j, int k, real scale, std::function<real(vec3& p)> field)
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
		for (int i = 0; i < num; i++)
		{
			vertex v1, v2, v3;
			v1.p = tri[i].p[0] * scale;
			v2.p = tri[i].p[1] * scale;
			v3.p = tri[i].p[2] * scale;

			v1.n = (v1.p - cur_centerpoint).normcopy();
			v2.n = (v2.p - cur_centerpoint).normcopy();
			v3.n = (v3.p - cur_centerpoint).normcopy();

			binvnorm = 1;
			//bautonorm = false;
			triang(v1, v2, v3);
		}
		return grid.val[0];
	}
	// -----------------------------------------------
	inline void updateaabb(const point3_t& p)
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
	inline real getdis_onpoly(const CENTERLIST& poly, crvec p, real minddis)
	{
		int currentpoint = -1;
		real alpha = 0.0f;

		if(poly.size() == 1)
		{
			real ddis = (p - poly[0].p).sqrlen();
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
						real ddis = (vp - v12 * dot).sqrlen();
						if (ddis < minddis)
						{
							alpha = dot / d12;
							currentpoint = i;
							minddis = ddis;
						}
					}
				}
				{
					real ddis = (p - p1).sqrlen();
					if (ddis < minddis)
					{
						currentpoint = i - 1;
						minddis = ddis;
					}
				}
				{
					real ddis = (p - p2).sqrlen();
					if (ddis < minddis)
					{
						currentpoint = i;
						minddis = ddis;
					}
				}
			}
		}
		if (currentpoint >= 0)
		{
			if (alpha == 0)
			{
				cur_centerpoint = poly[currentpoint].p;
				cur_radius = poly[currentpoint].r;
			}
			else
			{
				cur_centerpoint = blend(poly[currentpoint - 1].p, poly[currentpoint].p, alpha);
				cur_radius = blend(poly[currentpoint - 1].r, poly[currentpoint].r, alpha);
			}
		}
		return minddis;
	}
	// -----------------------------------------------
	inline void _render(submesh& sm)
	{
		for (int i = aabb.a.x; i < aabb.b.x; i++)
			for (int j = aabb.a.y; j < aabb.b.y; j++)
				for (int k = aabb.a.z; k < aabb.b.z; k++)
				{
					field2tris(sm, i, j, k, mesh_scale / 10.0f,
						[i, j, k](crvec p)->float
						{
							float dd = 1e5;
							for (auto& it : centerlinestack)
							{
								dd = getdis_onpoly(it, p, dd);
							}
							return dd / cur_radius;
						}
					);
				}
	}
	// -----------------------------------------------
	inline void add_centerpoint(crvec p, float r, int line)
	{
		r /= mesh_scale;
		if (r > border_deta)
			border_deta = r;

		vec3 tp = p / mesh_scale;
		updateaabb(point3_t(tp.x, tp.y, tp.z));

		while (line >= centerlinestack.size())
			centerlinestack.push_back(CENTERLIST());
		centerlinestack[line].push_back({ tp, r });

	}
	// -----------------------------------------------
	inline void render_pipe(submesh& sm)
	{
		if (centerlinestack.empty())
			return;

		PRINTVEC3((aabb.b - aabb.a));

		_render(sm);
	}
	// -----------------------------------------------
	inline void clear()
	{
		centerlinestack.clear();
	}
};
