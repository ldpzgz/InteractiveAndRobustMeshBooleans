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


#include "aux_structure.h"
#include "utils.h"
//init all kinds of array size in auxiliary,insert vertex-index  into v_map
inline void AuxiliaryStructure::initFromTriangleSoup(TriangleSoup &ts)
{
    num_original_vtx = ts.numVerts();
    num_original_tris = ts.numTris();

    coplanar_tris.resize(ts.numTris());

    tri2pts.resize(ts.numTris());
    edge2pts.resize(ts.numEdges());
    tri2segs.resize(ts.numTris());
    tri_has_intersections.resize(ts.numTris(), false);

    num_intersections = 0;
    num_tpi = 0;

    for(uint v_id = 0; v_id < ts.numVerts(); v_id++)
    {
        v_map.insert({ts.vert(v_id), v_id});
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline std::vector<std::pair<uint, uint> > &AuxiliaryStructure::intersectionList()
{
    return intersection_list;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline const std::vector<std::pair<uint, uint> > &AuxiliaryStructure::intersectionList() const
{
    return intersection_list;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//把一个点插入到三角形t_id里面
inline bool AuxiliaryStructure::addVertexInTriangle(uint t_id, uint v_id)
{
    assert(t_id < tri2pts.size());
    auto& points = tri2pts[t_id];
    if(contains(points, v_id)) return false;
    if(points.empty()) points.reserve(8);
    points.push_back(v_id);
    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//把一个点插入到e_id这条边里面
inline bool AuxiliaryStructure::addVertexInEdge(uint e_id, uint v_id)
{
    assert(e_id < edge2pts.size());
    auto& points = edge2pts[e_id];
    if(contains(points, v_id)) return false;
    if(points.empty()) points.reserve(8);
    points.push_back(v_id);
    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//把一条线段插入到三角形t_id里面，线段用两个点的索引来表示
inline bool AuxiliaryStructure::addSegmentInTriangle(uint t_id, const UIPair &seg)
{
    assert(t_id < tri2segs.size());
    UIPair key_seg = uniquePair(seg);
    auto& segments = tri2segs[t_id];
    if(contains(segments, key_seg)) return false;
    if(segments.empty()) segments.reserve(8);
    segments.push_back(key_seg);
    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//把一条线段跟两个三角形关联起来，表面这个线段是这两个三角形的交线，
inline void AuxiliaryStructure::addTrianglesInSegment(const UIPair &seg, uint tA_id, uint tB_id)
{
    UIPair key_seg = uniquePair(seg);
    auto& tris = seg2tris[key_seg];
    if(tris.empty()) tris.reserve(4);
    if(tA_id == tB_id)
    {
        if(!contains(tris, tA_id)) tris.push_back(tA_id);
    }
    else
    {
        if(!contains(tris, tA_id)) tris.push_back(tA_id);
        if(!contains(tris, tB_id)) tris.push_back(tB_id);
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//把一条线段一分为二，更新新线段的seg2tris信息。
inline void AuxiliaryStructure::splitSegmentInSubSegments(uint orig_v0, uint orig_v1, uint midpoint)
{
    auto& tris = segmentTrianglesList(std::make_pair(orig_v0, orig_v1));
    UIPair sub_seg0 = uniquePair(std::make_pair(orig_v0, midpoint));
    UIPair sub_seg1 = uniquePair(std::make_pair(midpoint, orig_v1));
    seg2tris[sub_seg0] = tris;
    seg2tris[sub_seg1] = tris;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//每一个三角形都记录哪些三角形跟他共面
inline void AuxiliaryStructure::addCoplanarTriangles(uint ta, uint tb)
{
    assert(ta != tb);
    assert(ta < coplanar_tris.size() && tb < coplanar_tris.size());

    if(coplanar_tris[ta].empty()) coplanar_tris[ta].reserve(8);
    if(coplanar_tris[tb].empty()) coplanar_tris[tb].reserve(8);
    coplanar_tris[ta].push_back(tb);
    coplanar_tris[tb].push_back(ta);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//获取跟t_id这个三角形共面的其他三角形
inline const auxvector<uint> &AuxiliaryStructure::coplanarTriangles(uint t_id) const
{
    assert(t_id < coplanar_tris.size());
    return coplanar_tris[t_id];
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//t_id这个三角形是否有跟其他三角形共面，
inline bool AuxiliaryStructure::triangleHasCoplanars(uint t_id) const
{
    assert(t_id < coplanar_tris.size());
    return (coplanar_tris[t_id].size() > 0);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//设置三角形有与其他三角形相交
inline void AuxiliaryStructure::setTriangleHasIntersections(uint t_id)
{
    assert(t_id < tri_has_intersections.size());
    tri_has_intersections[t_id] = true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回三角形t_id是否有与其他三角形相交
inline bool AuxiliaryStructure::triangleHasIntersections(uint t_id) const
{
    assert(t_id < tri_has_intersections.size());
    return tri_has_intersections[t_id];
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//获取插入到三角形里面的所有点
inline const auxvector<uint> &AuxiliaryStructure::trianglePointsList(uint t_id) const
{
    assert(t_id < tri2pts.size());
    return tri2pts[t_id];
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//获取插入到e_id边里面的所有点
inline const auxvector<uint> &AuxiliaryStructure::edgePointsList(uint e_id) const
{
    assert(e_id < edge2pts.size());
    return edge2pts[e_id];
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回插入到三角形t_id里面的所有线段
inline const auxvector<UIPair> &AuxiliaryStructure::triangleSegmentsList(uint t_id) const
{
    assert(t_id < tri2segs.size());
    return tri2segs[t_id];
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//返回线段seg所属的所有三角形 
inline const auxvector<uint> &AuxiliaryStructure::segmentTrianglesList(const UIPair &seg) const
{
    UIPair key_seg = uniquePair(seg);

    auto res = seg2tris.find(key_seg);
    assert(res != seg2tris.end());

    return res->second;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//把pair(点,点的位置) 插入到map里面
inline std::pair<uint, bool> AuxiliaryStructure::addVertexInSortedList(const genericPoint *v, uint pos)
{
    auto ins = v_map.insert({v, pos});

    return std::make_pair(ins.first->second, // the position of v (pos if first time, or the previous saved position otherwise)
                          ins.second);       // the result of the insert operation /true or false)
}

//把pair(polygon,pos) 存储到pockets_map里面如果之前还没有这个pair返回-1，否则返回pos
inline int AuxiliaryStructure::addVisitedPolygonPocket(const std::vector<uint> &polygon, uint pos)
{
    auto poly_it = pockets_map.insert(std::make_pair(polygon, pos));

    if(poly_it.second) return -1; // polygon not present yet

    return static_cast<int>(poly_it.first->second);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//按uip.first < uip.second，否则交换顺序
inline UIPair AuxiliaryStructure::uniquePair(const UIPair &uip) const
{
    if(uip.first < uip.second) return  uip;
    return std::make_pair(uip.second, uip.first);
}