//
//  MCFDriver.cpp
//
//  Christopher Batty, Fang Da 2014
//
//

#include "MCFDriver.h"

using namespace LosTopos;

bool MCFDriver::step(SurfTrack * st, double dt, bool bbwall)
{
    std::vector<Vec3d> v;
    
    // adaptive time stepping
    double dt_left = dt;
    double dt_largest_possible = 0;

    // set newpositions to position
    st->set_all_newpositions(st->get_positions());
    
    int substep_count = 0;
    while (dt_left > 0 && substep_count < 100)
    {
        evaluateV(st, v);
        
        // determine max dt
        dt_largest_possible = determineMaxDt(st, v);
        std::cout << "dt_left = " << dt_left << " max dt = " << dt_largest_possible << std::endl;
        if (dt_largest_possible >= dt_left)
        {
            dt_largest_possible = dt_left; // clamp dt_largest_possible

            // Setting dt_left to 0 manually as otherwise floating point precision
            // problems from doing dt_left -= dt_largest_possible can cause
            // dt_left to be a very small number near zero leading to some issues.
            dt_left = 0;
        } else
        {
            dt_left -= dt_largest_possible;
        }

        for (size_t i = 0; i < st->m_mesh.nv(); i++)
        {
            if (st->m_mesh.vertex_is_deleted(i))
                continue;

            Vec3d newpos = st->get_newposition(i);

            Vec3d m(1, 1, 1);
            if (bbwall)
            {
                const int onwall = onBBWall(newpos);
                m[0] = ((onwall & (1 << 0)) || (onwall & (1 << 3)) ? std::numeric_limits<double>::infinity() : 1);
                m[1] = ((onwall & (1 << 1)) || (onwall & (1 << 4)) ? std::numeric_limits<double>::infinity() : 1);
                m[2] = ((onwall & (1 << 2)) || (onwall & (1 << 5)) ? std::numeric_limits<double>::infinity() : 1);
            }

            newpos += Vec3d(v[i][0] / m[0], v[i][1] / m[1], v[i][2] / m[2]) * dt_largest_possible;
            assert(newpos == newpos);
            if (bbwall) newpos = Vec3d(std::min(1.0, std::max(0.0, newpos[0])), std::min(1.0, std::max(0.0, newpos[1])), std::min(1.0, std::max(0.0, newpos[2])));
            st->set_newposition(i, newpos);
        }
        
        substep_count++;
    }
    
    if (dt_left > 0)
    {
        std::cout << "[MCFDriver::step] Warning: time step incomplete after maximum number of substeps." << std::endl;
        return false;
    } else
    {
        std::cout << "[MCFDriver::step] Notice: time step finished after " << substep_count << " substeps." << std::endl;
        return true;
    }
}

double MCFDriver::determineMaxDt(const SurfTrack * st, const std::vector<Vec3d> & v)
{
    double global_max_dt = 1.0;
    for (size_t i = 0; i < st->m_mesh.ne(); i++)
    {
       if (st->m_mesh.edge_is_deleted(i))
          continue;
        const size_t v0 = st->m_mesh.m_edges[i][0];
        const size_t v1 = st->m_mesh.m_edges[i][1];
        const Vec3d edge = st->get_newposition(v1) - st->get_newposition(v0);
        const Vec3d rv = v[v1] - v[v0];

        const double alpha = 0.2; // allow the edge to be shrunk by at most 80% in one time substep
        const double a = dot(rv, rv);
        if (a < 1e-6) continue; // both endpoints have same velocity vectors
        const double b = 2 * dot(rv, edge);
        const double c = (1 - alpha*alpha) * dot(edge, edge);
        const double discriminant = b*b - 4*a*c;

        // If discriminant < 0 (no solutions) or discriminant = 0 (one solution),
        // then there is no constraint on t. Else if discriminant > 0 (two solutions),
        // and the smaller solution is nonnegative, then that solution is an
        // upper bound on t. Note that in this case both solutions are nonnegative
        // or nonpositive.
        double max_dt = std::numeric_limits<double>::infinity();
        if (discriminant > 1e-6)
        {
            const double smaller_root = (-b - std::sqrt(discriminant)) / (2*a);
            max_dt = smaller_root < 0 ? std::numeric_limits<double>::infinity() : smaller_root;
        }

        global_max_dt = std::min(global_max_dt, max_dt);
    }
    
    return global_max_dt;
}

void MCFDriver::evaluateV(const LosTopos::SurfTrack * st, std::vector<LosTopos::Vec3d> & v)
{
    const double COEF = 1;
    
    v.resize(st->m_mesh.nv(), Vec3d(0, 0, 0));
    std::fill(v.begin(), v.end(), Vec3d(0, 0, 0));

    for (size_t i = 0; i < st->m_mesh.nt(); i++)
    {
       if (st->m_mesh.triangle_is_deleted(i))
          continue;

        const Vec3d p0 = st->get_newposition(st->m_mesh.m_tris[i][0]);
        const Vec3d p1 = st->get_newposition(st->m_mesh.m_tris[i][1]);
        const Vec3d p2 = st->get_newposition(st->m_mesh.m_tris[i][2]);
        
        const Vec3d v01 = p1 - p0;
        const Vec3d v20 = p0 - p2;
        
        const Vec3d A = cross(v01, -v20);
        const double Anorm = mag(A);
        const Vec3d mul = A / Anorm;
        if ((Anorm == 0) || (mul != mul)) // cross product may sometimes be 0 due to degenerate triangles
            continue;
        
        const Vec3d p2part = COEF * cross(v01, mul);
        const Vec3d p1part = COEF * cross(mul, -v20);
        const Vec3d p0part = -(p1part + p2part);

        v[st->m_mesh.m_tris[i][0]] += p0part;
        v[st->m_mesh.m_tris[i][1]] += p1part;
        v[st->m_mesh.m_tris[i][2]] += p2part;
    }
}

int MCFDriver::onBBWall(const Vec3d & pos)
{
    static const double WALL_THRESHOLD = 1e-6;
    
    int walls = 0;
    if (pos[0] < 0 + WALL_THRESHOLD) walls |= (1 << 0);
    if (pos[1] < 0 + WALL_THRESHOLD) walls |= (1 << 1);
    if (pos[2] < 0 + WALL_THRESHOLD) walls |= (1 << 2);
    if (pos[0] > 1 - WALL_THRESHOLD) walls |= (1 << 3);
    if (pos[1] > 1 - WALL_THRESHOLD) walls |= (1 << 4);
    if (pos[2] > 1 - WALL_THRESHOLD) walls |= (1 << 5);
    
    return walls;
}
