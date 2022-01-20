/*****************************************************************************
		管道流体的渲染版本（21-11月优化过）
	构型内核带一个脚本接口，通过脚本来建构形态，另外可以加入UI操作
*****************************************************************************/
namespace PIPE_FLOW_FIELD_RENDER
{
	const float unit_scale = 1.0f;  // 每个方格的尺寸
	float border_deta = 1;		// 边界宽度

	struct center_t{vec3 p;float r;};
	using CENTERLIST = std::vector<center_t>;
	std::vector<CENTERLIST> centerlinestack;

	boundingboxi aabb;
	vec3 cur_center1, cur_center2;
	real cur_radius = 1;
	int cur_side = -1;

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

		grid.val[0] = field(grid.p[0]); if (grid.val[0] < 0) return 0; // 优化！
		grid.val[1] = field(grid.p[1]); if (grid.val[1] < 0) return 0; // 优化！
		grid.val[2] = field(grid.p[2]); if (grid.val[2] < 0) return 0; // 优化！
		grid.val[3] = field(grid.p[3]); if (grid.val[3] < 0) return 0; // 优化！
		grid.val[4] = field(grid.p[4]); if (grid.val[4] < 0) return 0; // 优化！
		grid.val[5] = field(grid.p[5]); if (grid.val[5] < 0) return 0; // 优化！
		grid.val[6] = field(grid.p[6]); if (grid.val[6] < 0) return 0; // 优化！
		grid.val[7] = field(grid.p[7]); if (grid.val[7] < 0) return 0; // 优化！


		/*if (rand() % 10 == 0)
		{
			color = 0xFFFF0000;
			pyramid(grid.p[0], scale * 5);
		}*/

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
	inline bool getdis_onpoly(const CENTERLIST& poly, crvec p, real& minddis)
	{
		int currentpoint = -1;

		if(poly.size() == 1)
		{
			real ddis = (p - poly[0].p).sqrlen() / (poly[0].r * poly[0].r);
			if (ddis < minddis)
			{
				currentpoint = 0;
				minddis = ddis;
			}
			if (ddis < 1.0)
				return true;
		}
		else
		{
			for (int i = 1; i < poly.size(); i++)
			{
				crvec p1 = poly[i - 1].p;
				crvec p2 = poly[i].p;
				{
					vec3 v12 = (p2 - p1);
					real d12 = v12.len(); v12 /= d12;
					vec3 vp = p - p1;
					real dot = vp.dot(v12);
					if (dot < 0 && dot > -2)
					{// 前端
						real radius = poly[i - 1].r;
						float ddis = _MAX(dot * dot, vp.sqrlen() - radius * radius);

						ddis = 1 + ddis / (radius * radius);
						if (ddis < minddis)
						{
							cur_radius = radius;
							minddis = ddis;
							cur_side = 1;
						}
					}
					else if (dot > d12 && dot < d12 + 2)
					{// 后端
						real radius = poly[i].r;
						float ddis = (dot - d12);
						ddis = _MAX(ddis * ddis, (p - p2).sqrlen() - radius * radius);

						ddis = 1 + ddis / (radius * radius);
						if (ddis < minddis)
						{
							cur_radius = radius;
							minddis = ddis;
							cur_side = 2;
						}
					}
					else 
					if(dot >= 0 && dot <= d12)
					{
						real alpha = dot / d12;
						real radius = blend(poly[i - 1].r, poly[i].r, alpha);
						float ddis = (vp - v12 * dot).sqrlen();
						
						float rrad = radius * radius;

						if(dot > 2 && dot < d12 - 2)
						{// 优化！
							if ((rrad - ddis) > radius * 4)
							{
								minddis = -1;
								return true; 
							}
						}
						ddis = ddis / rrad;
						{// 前后边界
							float dotdeta = _MIN(dot, d12 - dot);
							dotdeta = 0.95f - dotdeta * dotdeta / rrad;

							ddis = _MAX(ddis, dotdeta);
						}

						if (ddis < minddis)
						{// inside the pipe
							currentpoint = i;
							cur_center1 = poly[currentpoint - 1].p;
							cur_center2 = poly[currentpoint].p;
							cur_radius = radius;
							minddis = ddis;

							if (ddis < 1.0f)
								return true;
						}
					}
				}
			}
		}
		return false;
	}
	// -----------------------------------------------
	inline void _render(submesh& sm)
	{
		int tk = GetTickCount();
		
		for (int i = aabb.a.x; i < aabb.b.x; i++)
			for (int j = aabb.a.y; j < aabb.b.y; j++)
				for (int k = aabb.a.z; k < aabb.b.z; k++)
				{
					field2tris(sm, i, j, k, unit_scale / 5.0f,
						[i, j, k](crvec p)->float
						{
							cur_side = -1;
							float mdd = 1e10;
							for (auto& it : centerlinestack)
							{
								if (getdis_onpoly(it, p, mdd))
									return mdd;
							}
							if (cur_radius > 1.65f &&
								(mdd -1.0f) * (cur_radius * cur_radius) > cur_radius * 4) return -1;// 优化！

							//if (cur_side > 0)
							//	return 2;
							return mdd;
						}
					);
				}

		tk = GetTickCount() - tk;
		PRINTV(tk);
	}
	// -----------------------------------------------
	inline void add_centerpoint(crvec p, float r, int line)
	{
		r /= unit_scale;
		if (r < 1.8f)
			border_deta = r + 2;
		else
			border_deta = r + 1;

		vec3 tp = p / unit_scale;
		updateaabb(point3_t(tp.x + 0.5f, tp.y + 0.5f, tp.z + 0.5f));

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
	
	// -----------------------------------------------
	// PHG
	// -----------------------------------------------
	static real centerpoint(RealPHG::code& cd, int args)
	{
		real x = PHG_PARAM(1);
		real y = PHG_PARAM(2);
		real z = PHG_PARAM(3);

		real r = PHG_PARAM(4);
		real group = PHG_PARAM(5);

		PIPE_FLOW_FIELD_RENDER::add_centerpoint(vec3(x, y, z), r, group);
		return 0;
	}
	// -----------------------------------------------
	static real draw(RealPHG::code& cd, int args)
	{
		render_pipe(SUBMESH);

		clear();
	}
	// -----------------------------------------------
	inline void reg()
	{
		PRINT("reg");
		RealPHG::register_api("centerpoint", centerpoint);
		RealPHG::register_api("draw", draw);
	}
	
};
