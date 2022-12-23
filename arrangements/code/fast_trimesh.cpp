/*****************************************************************************************
 *              MIT License                                                              *
 *                                                                                       *
 * Copyright (c) 2022 G. Cherchi, M. Livesu, R. Scateni, M. Attene and F. Pellacini      *
 *                                                                                       *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this  *
 * software and associated documentation files (the "Software"), to deal in the Software *
 * without restriction, including without limitation the rights to use, copy, modify,    *
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to    *
 * permit persons to whom the Software is furnished to do so, subject to the following   *
 * conditions:                                                                           *
 *                                                                                       *
 * The above copyright notice and this permission notice shall be included in all copies *
 * or substantial portions of the Software.                                              *
 *                                                                                       *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,   *
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT    *
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION     *
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE        *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                *
 *                                                                                       *
 * Authors:                                                                              *
 *      Gianmarco Cherchi (g.cherchi@unica.it)                                           *
 *      https://people.unica.it/gianmarcocherchi/                                        *
 *                                                                                       *
 *      Marco Livesu (marco.livesu@ge.imati.cnr.it)                                      *
 *      http://pers.ge.imati.cnr.it/livesu/                                              *
 *                                                                                       *
 *      Riccardo Scateni (riccardo@unica.it)                                             *
 *      https://people.unica.it/riccardoscateni/                                         *
 *                                                                                       *
 *      Marco Attene (marco.attene@ge.imati.cnr.it)                                      *
 *      https://www.cnr.it/en/people/marco.attene/                                       *
 *                                                                                       *
 *      Fabio Pellacini (fabio.pellacini@uniroma1.it)                                    *
 *      https://pellacini.di.uniroma1.it                                                 *
 *                                                                                       *
 * ***************************************************************************************/

#include "fast_trimesh.h"

#include "utils.h"

inline FastTrimesh::FastTrimesh(const genericPoint *tv0, const genericPoint *tv1, const genericPoint *tv2, const uint *tv_id, const Plane &ref_p)
{
    addVert(tv0, tv_id[0]);
    addVert(tv1, tv_id[1]);
    addVert(tv2, tv_id[2]);
    addTri(0, 1, 2);

    triangle_plane = ref_p;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//根据输入的外部顶点与三角形信息，建立自己的顶点，边，三角形结构,这里的边不是半边
inline FastTrimesh::FastTrimesh(const std::vector<genericPoint *> &in_verts, const std::vector<uint> &in_tris, bool parallel)
{
    if(parallel)
    {
        // add vertices
        vertices.resize(in_verts.size());
        v2e.resize(in_verts.size());
        for(uint v_id = 0; v_id < (uint)in_verts.size(); v_id++)
        {
            vertices[v_id] = iVtx(in_verts[v_id], 0);
        }

        // add triangles
        triangles.resize(in_tris.size() / 3);
        for(uint t_id = 0; t_id < in_tris.size() / 3; t_id++)
        {
            auto tv0_id = in_tris[3 * t_id], tv1_id = in_tris[3 * t_id + 1], tv2_id = in_tris[3 * t_id + 2];
            triangles[t_id] = {tv0_id, tv1_id, tv2_id, 0};
        }

        // we want to remove duplicates while keeping track of their provenance
        // build edges
        std::vector<std::array<uint, 2>> sorted_edges(in_tris.size());
        for(uint t_id = 0; t_id < in_tris.size() / 3; t_id++)
        {
            auto tv0_id = in_tris[3 * t_id], tv1_id = in_tris[3 * t_id + 1], tv2_id = in_tris[3 * t_id + 2];
            sorted_edges[t_id * 3 + 0] = {std::min(tv0_id, tv1_id), std::max(tv0_id, tv1_id)};
            sorted_edges[t_id * 3 + 1] = {std::min(tv1_id, tv2_id), std::max(tv1_id, tv2_id)};
            sorted_edges[t_id * 3 + 2] = {std::min(tv2_id, tv0_id), std::max(tv2_id, tv0_id)};
        }
        tbb::parallel_sort(sorted_edges.begin(), sorted_edges.end());
        sorted_edges.erase(std::unique(sorted_edges.begin(), sorted_edges.end()), sorted_edges.end());

        // fix edges
        edges.resize(sorted_edges.size());
        for(uint e_id = 0; e_id < sorted_edges.size(); e_id++)
        {
            auto [v0, v1] = sorted_edges[e_id];
            edges[e_id] = {v0, v1, false};
        }

        // fix adjacencies
        v2e.resize(vertices.size());
        for(uint e_id = 0; e_id < sorted_edges.size(); e_id++)
        {
            auto [v0, v1] = sorted_edges[e_id];
            v2e[v0].push_back(e_id);
            v2e[v1].push_back(e_id);
        }

        // fix adjacencies
        e2t.resize(sorted_edges.size());
        std::vector<tbb::spin_mutex> mutexes(sorted_edges.size());
        tbb::parallel_for((uint)0, (uint)triangles.size(), [&](uint t_id) {
            auto tv0_id = triangles[t_id].v[0], tv1_id = triangles[t_id].v[1], tv2_id = triangles[t_id].v[2];
            auto e0_id = (uint)(std::lower_bound(sorted_edges.begin(), sorted_edges.end(), std::array<uint, 2>{std::min(tv0_id, tv1_id), std::max(tv0_id, tv1_id)}) - sorted_edges.begin());
            auto e1_id = (uint)(std::lower_bound(sorted_edges.begin(), sorted_edges.end(), std::array<uint, 2>{std::min(tv1_id, tv2_id), std::max(tv1_id, tv2_id)}) - sorted_edges.begin());
            auto e2_id = (uint)(std::lower_bound(sorted_edges.begin(), sorted_edges.end(), std::array<uint, 2>{std::min(tv2_id, tv0_id), std::max(tv2_id, tv0_id)}) - sorted_edges.begin());
            {
                std::lock_guard<tbb::spin_mutex> lock(mutexes[e0_id]);
                e2t[e0_id].push_back(t_id);
            }
            {
                std::lock_guard<tbb::spin_mutex> lock(mutexes[e1_id]);
                e2t[e1_id].push_back(t_id);
            }
            {
                std::lock_guard<tbb::spin_mutex> lock(mutexes[e2_id]);
                e2t[e2_id].push_back(t_id);
            }
        });
    }
    else
    {
        vertices.reserve(in_verts.size());
        edges.reserve(in_verts.size() / 2);
        triangles.reserve(in_tris.size() / 3);
        v2e.reserve(in_verts.size());

        for(auto &v : in_verts) addVert(v);

        for(uint t_id = 0; t_id < in_tris.size() / 3; t_id++)
            addTri(in_tris[3 * t_id], in_tris[3 * t_id + 1], in_tris[3 * t_id + 2]);
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void FastTrimesh::preAllocateSpace(uint estimated_num_verts)
{
    vertices.reserve(estimated_num_verts);
    rev_vtx_map.reserve(estimated_num_verts);
    edges.reserve(estimated_num_verts / 2);
    triangles.reserve(estimated_num_verts / 3);
    v2e.reserve(estimated_num_verts);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//重置所有三角形的nodeid
inline void FastTrimesh::resetTrianglesInfo()
{
    for(uint t_id = 0; t_id < numTris(); t_id++)
        triangles[t_id].info = 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint FastTrimesh::numVerts() const
{
    return static_cast<uint>(vertices.size());
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint FastTrimesh::numEdges() const
{
    return static_cast<uint>(edges.size());
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint FastTrimesh::numTris() const
{
    return static_cast<uint>(triangles.size());
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline Plane FastTrimesh::refPlane() const
{
    return triangle_plane;
}

/************************************************************************************************
 *          VERTICES
 * *********************************************************************************************/
//返回顶点v_id的point指针
inline const genericPoint* FastTrimesh::vert(uint v_id) const
{
    assert(v_id < vertices.size() && "vtx id out of range");
    return vertices[v_id].p;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回new_v_id的重复顶点的id，一个顶点只能有一个重复顶点
inline uint FastTrimesh::vertOrigID(uint new_v_id) const
{
    assert(new_v_id < vertices.size() && "vtx id out of range");
    return vertices[new_v_id].info;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回orig_v_id的重复顶点的id,一个顶点只能有一个重复顶点
inline  uint FastTrimesh::vertNewID(uint orig_v_id) const
{
    auto it = rev_vtx_map.find(orig_v_id);
    assert(it != rev_vtx_map.end() && "vtx id not found in reverse map");

    return it->second;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//有多少条边共用了顶点v_id
inline uint FastTrimesh::vertValence(uint v_id) const
{
    assert(v_id < vertices.size() && "vtx id out of range");
    return static_cast<uint>(v2e[v_id].size());
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回一个数组，共用了顶点v_id的所有边的ids
inline const fmvector<uint> &FastTrimesh::adjV2E(uint v_id) const
{
    assert(v_id < vertices.size() && "vtx id out of range");
    return v2e[v_id];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回一个数组，里面是共用顶点v_id的所有三角形的id
inline fmvector<uint> FastTrimesh::adjV2T(uint v_id) const
{
    assert(v_id < vertices.size() && "vtx id out of range");
    fmvector<uint> v2t;

    for(uint e_id : v2e[v_id])
    {
        for(uint t_id : e2t[e_id])
            v2t.push_back(t_id);
    }
    remove_duplicates(v2t);

    return v2t;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void FastTrimesh::resetVerticesInfo()
{
    for(auto &v : vertices)
        v.info = 0;
}
//设置顶点的额外的一个uint信息，算是顶点属性吧
//用于存储origin顶点的id（重复顶点）
inline void FastTrimesh::setVertInfo(const uint v_id, const uint info)
{
    assert(v_id < vertices.size() && "vtx id out of range");
    vertices[v_id].info = info;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回顶点v_id的info（通常是跟它重复的origin顶点的id)
inline uint FastTrimesh::vertInfo(const uint v_id) const
{
    assert(v_id < vertices.size() && "vtx id out of range");
    return vertices[v_id].info;
}

/************************************************************************************************
 *          EDGES
 * *********************************************************************************************/
//返回一个pair<int,int>,是边e_id的两个顶点的id
inline const std::pair<uint, uint> &FastTrimesh::edge(uint e_id) const
{
    assert(e_id < edges.size() && "edge id out of range");
    return edges[e_id].v;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回边e_id的第0个或者第1个顶点的id
inline uint FastTrimesh::edgeVertID(uint e_id, uint off) const
{
    assert(e_id < edges.size() && "edge id out of range");
    if(off == 0) return edges[e_id].v.first;
    return edges[e_id].v.second;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//根据两个顶点id，返回这两个顶点组成的边的id,类里面有顶点-边的map
inline int FastTrimesh::edgeID(uint ev0_id, uint ev1_id) const
{
    assert(ev0_id != ev1_id && "edge with equal endpoints");
    assert((ev0_id < vertices.size() && ev1_id < vertices.size()) && "vtx id out of range");

    for(uint e_id : v2e[ev0_id])
    {
        if(edgeContainsVert(e_id, ev0_id) && edgeContainsVert(e_id, ev1_id))
            return static_cast<int>(e_id);
    }

    return -1;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//判断边是否被访问过
inline bool FastTrimesh::edgeIsConstr(uint e_id) const
{
    assert(e_id < edges.size() && "edge id out of range");
    return edges[e_id].constr;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//设置边被访问了
inline void FastTrimesh::setEdgeConstr(uint e_id)
{
    assert(e_id < edges.size() && "edge id out of range");
    edges[e_id].constr = true;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回三角形t_id中顶点v_id对面的那条边的id
inline uint FastTrimesh::edgeOppToVert(uint t_id, uint v_id) const
{
    assert(t_id < triangles.size() && "tri id out of range");
    assert(triContainsVert(t_id, v_id) && "tri doesn't contain vtx");

    int e_id = -1;
    if(triangles[t_id].v[0] == v_id) e_id = edgeID(triangles[t_id].v[1], triangles[t_id].v[2]); else
    if(triangles[t_id].v[1] == v_id) e_id = edgeID(triangles[t_id].v[0], triangles[t_id].v[2]); else
    if(triangles[t_id].v[2] == v_id) e_id = edgeID(triangles[t_id].v[0], triangles[t_id].v[1]);

    assert(e_id != -1 && "opposite edge not found in tri");
    return static_cast<uint>(e_id);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//边e_id是不是有1个三角形共用它
inline bool FastTrimesh::edgeIsBoundary(uint e_id) const
{
    assert(e_id < edges.size() && "edge id out of range");
    return e2t[e_id].size() == 1;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//边e_id是不是有两个三角形共用它
inline bool FastTrimesh::edgeIsManifold(uint e_id) const
{
    assert(e_id < edges.size() && "edge id out of range");
    return e2t[e_id].size() == 2;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回共用了边e_id的三角形，最多两个
inline const fmvector<uint> &FastTrimesh::adjE2T(uint e_id) const
{
    assert(e_id < edges.size() && "edge id out of range");
    return e2t[e_id];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//设置边被访问了
inline void FastTrimesh::edgeSetVisited(uint e_id, const bool &vis)
{
    assert(e_id < edges.size() && "edge id out of range");
    edges[e_id].constr = vis;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//判断边是否被访问过
inline bool FastTrimesh::edgeIsVisited(uint e_id) const
{
    assert(e_id < edges.size() && "edge id out of range");
    return edges[e_id].constr;
}


/************************************************************************************************
 *          TRIANGLES
 * *********************************************************************************************/
//返回三角形t_id的三个顶点的id的数组的地址
inline const uint *FastTrimesh::tri(uint t_id) const
{
    assert(t_id < triangles.size() && "tri id out of range");
    return triangles[t_id].v;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//根据三个顶点的id得到三角形的id，如果不存在返回-1
inline int FastTrimesh::triID(uint tv0_id, uint tv1_id, uint tv2_id) const
{
    assert((tv0_id < vertices.size() && tv1_id < vertices.size() && tv2_id < vertices.size()) && "vtx id out of range");

    int e_id = edgeID(tv0_id, tv1_id);
    if(e_id == -1) return -1;

    for(uint t_id : e2t[static_cast<uint>(e_id)])
        if(triContainsVert(t_id, tv2_id))
            return static_cast<int>(t_id);

    return -1;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回三角形t_id的第off（0，1，2）个顶点的id
inline uint FastTrimesh::triVertID(uint t_id, uint off) const
{
    assert(t_id < triangles.size() && "tri id out of range");
    return triangles[t_id].v[off];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回三角形t_id的第off（0，1，2）个顶点的顶点指针
inline const genericPoint *FastTrimesh::triVert(uint t_id, uint off) const
{
    assert(t_id < triangles.size() && "tri id out of range");
    return vertices[triangles[t_id].v[off]].p;
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回三角形t_id的第off（0，1，2）条边的id
inline int FastTrimesh::triEdgeID(uint t_id, uint off) const
{
    assert(t_id < triangles.size() && "tri id out of range");
    return edgeID(triVertID(t_id, off), triVertID(t_id, ((off +1) %3)));
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//三角形的nodeid存储在三角形的info上
inline uint FastTrimesh::triNodeID(uint t_id) const
{
    assert(t_id < triangles.size() && "tri id out of range");
    return triangles[t_id].info;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//三角形的nodeid存储在三角形的info上
inline void FastTrimesh::setTriNodeID(uint t_id, uint n_id)
{
    assert(t_id < triangles.size() && "tri id out of range");
    triangles[t_id].info = n_id;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
* 找到三角形t_id中除了v0_id,v1_id两个点外的第三个点
*/
inline uint FastTrimesh::triVertOppositeTo(uint t_id, uint v0_id, uint v1_id) const
{
    assert(t_id < triangles.size() && "tri id out of range");
    assert(v0_id != v1_id && "verts are equal");
    assert((triContainsVert(t_id, v0_id) && triContainsVert(t_id, v0_id)) && "tri dosn't contain vtx");

    for(uint off = 0; off < 3; off++)
    {
        uint v_id = triangles[t_id].v[off];
        if(v_id != v0_id && v_id != v1_id) return v_id;
    }

    assert(false && "This should not happen");
    return 0; // warning killer
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回于三角形t_id共用同一条边e_id的另外一个三角形，如果没有返回-1
inline int FastTrimesh::triOppToEdge(uint e_id, uint t_id) const
{
    assert(e_id < edges.size() && "edge id out of range");
    assert(t_id < triangles.size() && "tri id out of range");
    assert(e2t[e_id].size() <= 2 && "no opposite tri found");

    if(e2t[e_id].size() == 1) return -1; // boundary edge

    if(e2t[e_id][0] == t_id)
        return static_cast<int>(e2t[e_id][1]);
    else if(e2t[e_id][1] == t_id)
        return static_cast<int>(e2t[e_id][0]);

    assert(false && "no opposite tri found");
    return -1; // warning killer
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回三角形t_id的三条边的id
inline fmvector<uint> FastTrimesh::adjT2E(uint t_id) const
{
    assert(t_id < triangles.size() && "tri id out of range");

    return {static_cast<uint>(triEdgeID(t_id, 0)),
            static_cast<uint>(triEdgeID(t_id, 1)),
            static_cast<uint>(triEdgeID(t_id, 2))};
}
//返回所有的三角形的边id
inline std::vector<std::array<uint, 3>> FastTrimesh::adjT2EAll(bool parallel) const {
    if(parallel) {
        std::vector<std::array<uint, 3>> adjT2E(triangles.size());
        tbb::parallel_for((uint)0, (uint)triangles.size(), [this, &adjT2E](uint t_id) {
            adjT2E[t_id] = {static_cast<uint>(triEdgeID(t_id, 0)),
                            static_cast<uint>(triEdgeID(t_id, 1)),
                            static_cast<uint>(triEdgeID(t_id, 2))};
        });
        return adjT2E;
    } else {
        std::vector<std::array<uint, 3>> adjT2E(triangles.size());
        for(uint t_id = 0; t_id < (uint)triangles.size(); t_id++) {
            adjT2E[t_id] = {static_cast<uint>(triEdgeID(t_id, 0)),
                            static_cast<uint>(triEdgeID(t_id, 1)),
                            static_cast<uint>(triEdgeID(t_id, 2))};
        }
        return adjT2E;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回与三角形t_id相邻的所有三角形的id
inline fmvector<uint> FastTrimesh::adjT2T(uint t_id) const
{
    assert(t_id < triangles.size() && "tri id out of range");
    fmvector<uint> res;

    for(uint e_id : adjT2E(t_id))
    {
        for(uint nbr_t : adjE2T(e_id))
        {
            if(nbr_t != t_id)
                res.push_back(nbr_t);
        }
    }

    return res;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//根据两个顶点的前后关系，判断这个三角形是不是逆时针的
//逆时针的顶点顺序是0，1，2
inline bool FastTrimesh::triVertsAreCCW(uint t_id, uint curr_v_id, uint prev_v_id) const
{
    uint prev_off = triVertOffset(t_id, prev_v_id);
    uint curr_off = triVertOffset(t_id, curr_v_id);
    if(curr_off == ((prev_off +1) %3)) return true;
    return false;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//判断t_id这个三角形在它的投影平面上是顺时针还是逆时针
inline int FastTrimesh::triOrientation(uint t_id) const
{
    assert(t_id < triangles.size() && "tri id out of range");
    switch (triangle_plane)
    {
        case XY: return genericPoint::orient2Dxy(*triVert(t_id, 0), *triVert(t_id, 1), *triVert(t_id, 2));
        case YZ: return genericPoint::orient2Dyz(*triVert(t_id, 0), *triVert(t_id, 1), *triVert(t_id, 2));
        case ZX: return genericPoint::orient2Dzx(*triVert(t_id, 0), *triVert(t_id, 1), *triVert(t_id, 2));
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//判断顶点v_id是否在三角形t_id中
inline bool FastTrimesh::triContainsVert(uint t_id, uint v_id) const
{
    assert(t_id < triangles.size() && "tri id out of range");

    if(triangles[t_id].v[0] == v_id) return true;
    if(triangles[t_id].v[1] == v_id) return true;
    if(triangles[t_id].v[2] == v_id) return true;
    return false;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//如果顶点v_id在三角形t_id中，返回它是三角形中的第几个点（0，1，2）
inline uint FastTrimesh::triVertOffset(uint t_id, uint v_id) const
{
    for(uint off = 0; off < 3; off++)
        if(triangles[t_id].v[off] == v_id) return off;

    assert(false && "This should not happen");
    return 0; // warning killer
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//拿到三角形t_id的额外的uint
inline uint FastTrimesh::triInfo(uint t_id) const
{
    assert(t_id < triangles.size() && "tri id out of range");
    return triangles[t_id].info;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//为三角形t_id设置一个额外的uint，每个三角形都可以设置一个
inline void FastTrimesh::setTriInfo(uint t_id, uint val)
{
    assert(t_id < triangles.size() && "tri id out of range");
    triangles[t_id].info = val;
}

/************************************************************************************************
 *          MESH MANIPULATION
 * **********************************************************************************************/
//添加顶点，这个顶点是某个顶点orig_v_id的重复
inline uint FastTrimesh::addVert(const genericPoint *v, uint orig_v_id)
{
    uint v_id = static_cast<uint>(vertices.size());
    vertices.emplace_back(v, orig_v_id);

    v2e.emplace_back();
    rev_vtx_map[orig_v_id] = v_id;

    return v_id;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//添加顶点，但是这个顶点的index设置为0了
inline void FastTrimesh::addVert(const genericPoint *v)
{
    vertices.emplace_back(v, 0);
    v2e.emplace_back();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//根据三个顶点id构建一个三角形，同时构建相应的边，跟新e2t索引
//如果三角形已经存在直接返回它的id
inline uint FastTrimesh::addTri(uint tv0_id, uint tv1_id, uint tv2_id)
{
    assert((tv0_id < vertices.size() && tv1_id < vertices.size() && tv2_id < vertices.size()) && "vtx id out of range");
    assert((tv0_id != tv1_id && tv0_id != tv2_id && tv1_id != tv2_id) && "degenerate triangle");

    int t_id = triID(tv0_id, tv1_id, tv2_id);
    if(t_id != -1) return static_cast<uint>(t_id);

    t_id = static_cast<int>(triangles.size());

    triangles.emplace_back(tv0_id, tv1_id, tv2_id, 0);

    // adding missing edges
    int e0_id = addEdge(tv0_id, tv1_id);
    int e1_id = addEdge(tv1_id, tv2_id);
    int e2_id = addEdge(tv2_id, tv0_id);

    e2t[static_cast<uint>(e0_id)].push_back(static_cast<uint>(t_id));
    e2t[static_cast<uint>(e1_id)].push_back(static_cast<uint>(t_id));
    e2t[static_cast<uint>(e2_id)].push_back(static_cast<uint>(t_id));

    return static_cast<uint>(t_id);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//删除边e_id，其实是删除两个三角形（或者1个）
inline void FastTrimesh::removeEdge(uint e_id)
{
    assert(e_id < edges.size() && "edge id out of range");
    removeTris(e2t[e_id]);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//先找到三角形的三条边，删除每条边e2t里面的t_id，
//遍历这三条边，每条边取出两个顶点，如果这个边已经没有三角形了，
//在v2e里面删除这条边，
//移除每一条边（边占用的空间）。
//最后移除这个三角形（三角形占用的空间）
inline void FastTrimesh::removeTri(uint t_id)
{
    assert(t_id < triangles.size() && "tri id out of range");

    uint e0_id = static_cast<uint>(triEdgeID(t_id, 0));
    uint e1_id = static_cast<uint>(triEdgeID(t_id, 1));
    uint e2_id = static_cast<uint>(triEdgeID(t_id, 2));
    removeFromVec(e2t[e0_id], t_id);
    removeFromVec(e2t[e1_id], t_id);
    removeFromVec(e2t[e2_id], t_id);

    int len = 0;
    std::array< uint, 3 > dangling_edges; // higher ids first
    if(e2t[e0_id].empty()) dangling_edges[len++] = e0_id;
    if(e2t[e1_id].empty()) dangling_edges[len++] = e1_id;
    if(e2t[e2_id].empty()) dangling_edges[len++] = e2_id;
    std::sort(dangling_edges.begin(), dangling_edges.begin()+len, std::greater<uint>());

    for(auto i = 0; i < len; i++)
    {
        uint e_id = dangling_edges[i];
        uint v0_id = edges[e_id].v.first;
        uint v1_id = edges[e_id].v.second;
        if (v1_id > v0_id) std::swap(v0_id, v1_id);
        removeFromVec(v2e[v0_id], e_id);
        removeFromVec(v2e[v1_id], e_id);
        removeEdgeUnref(e_id);
    }

    removeTriUnref(t_id);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//先给t_ids数组排序，再移除每一个三角形
inline void FastTrimesh::removeTris(const std::vector<uint> &t_ids)
{
    std::vector<uint> tmp_t_ids = t_ids;
    std::sort(tmp_t_ids.rbegin(), tmp_t_ids.rend());
    for(uint &t_id : tmp_t_ids) removeTri(t_id);
}
//先给t_ids数组排序，再移除每一个三角形
inline void FastTrimesh::removeTris(const fmvector<uint> &t_ids)
{
    fmvector<uint> tmp_t_ids = t_ids;
    std::sort(tmp_t_ids.rbegin(), tmp_t_ids.rend());
    for(uint &t_id : tmp_t_ids) removeTri(t_id);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//把一条边打断为两条边，把e_id这条边，用v_id与e_id上两个顶点组成的两条边替换
inline void FastTrimesh::splitEdge(const uint  &e_id, uint v_id)
{
    assert(e_id < edges.size() && "edge id out of range");

    uint ev0_id = edges[e_id].v.first;
    uint ev1_id = edges[e_id].v.second;

    for(uint t_id : e2t[e_id])
    {
        uint v_opp = triVertOppositeTo(t_id, ev0_id, ev1_id);
        if(triVertsAreCCW(t_id, ev0_id, ev1_id)) std::swap(ev0_id, ev1_id);

        addTri(v_opp, ev0_id, v_id);
        addTri(v_opp, v_id, ev1_id);
    }

    removeTris(e2t[e_id]);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//把一条边打断为两条边，把e_id这条边，用v_id与e_id上两个顶点组成的两条边替换
inline void FastTrimesh::splitEdge(const uint  &e_id, uint v_id, Tree &tree)
{
    assert(e_id < edges.size() && "edge id out of range");

    uint ev0_id = edges[e_id].v.first;
    uint ev1_id = edges[e_id].v.second;

    for(uint t_id : e2t[e_id])
    {
        uint v_opp = triVertOppositeTo(t_id, ev0_id, ev1_id);
        if(triVertsAreCCW(t_id, ev0_id, ev1_id)) std::swap(ev0_id, ev1_id);

        uint t0_id = addTri(v_opp, ev0_id, v_id);
        uint t1_id = addTri(v_opp, v_id, ev1_id);

        uint n0_id = tree.addNode(v_opp, ev0_id, v_id);
        uint n1_id = tree.addNode(v_opp, v_id, ev1_id);

        uint node_id = triNodeID(t_id);
        tree.addChildren(node_id, n0_id, n1_id);

        setTriNodeID(t0_id, n0_id);
        setTriNodeID(t1_id, n1_id);
    }

    removeTris(e2t[e_id]);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//将三角形t_id划分为三个子三角形，删掉这个三角形t_id
inline void FastTrimesh::splitTri(uint t_id, uint v_id)
{
    assert(t_id < triangles.size() && "tri id out of range");
    assert(v_id < vertices.size() && "vtx id out of range");

    addTri(triVertID(t_id, 0), triVertID(t_id, 1), v_id);
    addTri(triVertID(t_id, 1), triVertID(t_id, 2), v_id);
    addTri(triVertID(t_id, 2), triVertID(t_id, 0), v_id);

    removeTri(t_id);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//将三角形t_id划分为三个子三角形，删掉这个三角形t_id
inline void FastTrimesh::splitTri(uint t_id, uint v_id, Tree &tree)
{
    assert(t_id < triangles.size() && "tri id out of range");
    assert(v_id < vertices.size() && "vtx id out of range");

    uint node_id = triNodeID(t_id);

    uint t0_id = addTri(triVertID(t_id, 0), triVertID(t_id, 1), v_id);
    uint t1_id = addTri(triVertID(t_id, 1), triVertID(t_id, 2), v_id);
    uint t2_id = addTri(triVertID(t_id, 2), triVertID(t_id, 0), v_id);

    uint n0_id = tree.addNode(triVertID(t_id, 0), triVertID(t_id, 1), v_id);
    uint n1_id = tree.addNode(triVertID(t_id, 1), triVertID(t_id, 2), v_id);
    uint n2_id = tree.addNode(triVertID(t_id, 2), triVertID(t_id, 0), v_id);
    tree.addChildren(node_id, n0_id, n1_id, n2_id);

    // the triangle label contains the node position
    setTriNodeID(t0_id, n0_id);
    setTriNodeID(t1_id, n1_id);
    setTriNodeID(t2_id, n2_id);

    removeTri(t_id);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//翻转三角形，交换两个顶点的位置而已，顺时针变逆时针或者相反
inline void FastTrimesh::flipTri(uint t_id)
{
    assert(t_id < triangles.size() && "tri id out of range");

    uint tmp = triangles[t_id].v[0];
    triangles[t_id].v[0] = triangles[t_id].v[2];
    triangles[t_id].v[2] = tmp;
}

/***********************************************************************************************
 *          PRIVATE METHODS
 * ********************************************************************************************/
//根据两个顶点id，创建一条边，如果边已经存在直接返回它的id
//会完善其他的顶点-边索引，e2t数组会新增一个成员
//返回新增的边的id=edges.size();
inline int FastTrimesh::addEdge(uint ev0_id, uint ev1_id)
{
    int e_id = edgeID(ev0_id, ev1_id);
    if(e_id != -1) return e_id;

    e_id = static_cast<int>(edges.size());

    edges.emplace_back(ev0_id, ev1_id, false);

    e2t.emplace_back();
    v2e[ev0_id].push_back(static_cast<uint>(e_id));
    v2e[ev1_id].push_back(static_cast<uint>(e_id));

    return e_id;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//判断边e_id是否包含v_id这个顶点
inline bool FastTrimesh::edgeContainsVert(uint e_id, uint v_id) const
{
    if(edges[e_id].v.first == v_id) return true;
    if(edges[e_id].v.second == v_id) return true;
    return false;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//删除数组vec中所有值为elem的元素
inline void FastTrimesh::removeFromVec(fmvector<uint> &vec, uint elem)
{
    vec.erase(std::remove(vec.begin(), vec.end(), elem), vec.end());
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//交换两个三角形在各种数组中的位置
inline void FastTrimesh::triSwitch(uint t0_id, uint t1_id)
{
    if(t0_id == t1_id) return;

    std::swap(triangles[t0_id], triangles[t1_id]);

    std::array<int, 6> edges_to_update{triEdgeID(t0_id, 0), triEdgeID(t0_id, 1), triEdgeID(t0_id, 2), triEdgeID(t1_id, 0), triEdgeID(t1_id, 1), triEdgeID(t1_id, 2)};
    auto len = (int)remove_duplicates(edges_to_update);

    for(int i = 0; i < len; i++)
    {
        uint e_id = edges_to_update[i];
        if(e_id == -1) continue; // the triangle t_id could contain an already removed edge

        for(uint &t_id : e2t[static_cast<uint>(e_id)])
        {
            if(t_id == t0_id)  t_id = t1_id; else
            if(t_id == t1_id)  t_id = t0_id;
        }
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//交换两条边在各种数组中的位置，边所在的数组的索引要交换，边对应的三角形数组里面也要交换
//边上的顶点，顶点对应的边的数组里面也要交换。
inline void FastTrimesh::edgeSwitch(uint e0_id, const uint e1_id)
{
    if(e0_id == e1_id) return;

    std::swap(edges[e0_id], edges[e1_id]);
    std::swap(e2t[e0_id], e2t[e1_id]);

    std::array<uint, 4> verts_to_update{edges[e0_id].v.first, 
      edges[e0_id].v.second, 
      edges[e1_id].v.first, 
      edges[e1_id].v.second};
    auto len = (int)remove_duplicates(verts_to_update);

    for(int i = 0; i < len; i++)
    {
        uint v_id = verts_to_update[i];
        for(uint &e_id : v2e[v_id])
        {
            if(e_id == e0_id)  e_id = e1_id; else
            if(e_id == e1_id)  e_id = e0_id;
        }
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//通过与最后一条边交换的方式，移除一条边
inline void FastTrimesh::removeEdgeUnref(uint e_id)
{
    e2t[e_id].clear();
    edgeSwitch(e_id, static_cast<uint>(edges.size() -1));
    edges.pop_back();
    e2t.pop_back();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//通过与最后一个三角形交换的方式移除一个三角形
inline void FastTrimesh::removeTriUnref(uint t_id)
{
    triSwitch(t_id, static_cast<uint>(triangles.size() -1));
    triangles.pop_back();
}


