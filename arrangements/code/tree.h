/*****************************************************************************************
 *              MIT License                                                              *
 *                                                                                       *
 * Copyright (c) 2020 Gianmarco Cherchi, Marco Livesu, Riccardo Scateni e Marco Attene   *
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
 *      https://www.gianmarcocherchi.com                                                 *
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
 * ***************************************************************************************/

#ifndef TREE_H
#define TREE_H

#include "predicates/indirect_predicates.h"
#include <vector>

typedef unsigned int uint;
//一个node代表一个三角形，
//node里面是3个顶点的id，以及3个子节点（一个三角形最多被划分为3个子三角形）的id
struct Node
{
    Node(){}

    Node(const uint &_v0, const uint &_v1, const uint &_v2) :v0(_v0), v1(_v1), v2(_v2)
    {
        children[0] = -1;
        children[1] = -1;
        children[2] = -1;
    }

    uint v0, v1, v2;
    int children[3];
};
//一个树代表一个三角形，这个三角形可能被三角化为多个三角形，用树型结构来组织管理
//树里面有一个Node数组，a node装的是一个三角形的三个顶点的索引，还有三个子节点(在node数组里面的index)
//一个三角形只有3个子节点，说明一个三角形最多只被划分为3个三角形
class Tree
{
    public:

        inline Tree() {}

        inline Tree(const uint &size)
        {
            nodes.reserve(size);
        }
        //用三个顶点的index 构造一个节点，添加到节点数组的末尾，返回节点id
        inline uint addNode(const uint &v0, const uint &v1, const uint &v2)
        {
            nodes.emplace_back(v0, v1, v2);
            return static_cast<uint>(nodes.size()) -1;
        }

        inline const Node &getNode(const uint &node_id) const
        {
            assert(node_id < nodes.size() && "out fo range node id");
            return nodes[node_id];
        }
        //为node_id节点添加两个子节点，(node_id必须还没有子节点）
        inline void addChildren(const uint &node_id, const uint &c0, const uint &c1)
        {
            assert(node_id < nodes.size() && "out fo range node id");
            assert(nodes[node_id].children[0] == -1 && "assigning no empty children list");

            nodes[node_id].children[0] = static_cast<int>(c0);
            nodes[node_id].children[1] = static_cast<int>(c1);
        }
        //为node_id节点添加三个子节点，(node_id必须还没有子节点）
        inline void addChildren(const uint &node_id, const uint &c0, const uint &c1, const uint &c2)
        {
            assert(node_id < nodes.size() && "out fo range node id");
            assert(nodes[node_id].children[0] == -1 && "assigning no empty children list");

            nodes[node_id].children[0] = static_cast<int>(c0);
            nodes[node_id].children[1] = static_cast<int>(c1);
            nodes[node_id].children[2] = static_cast<int>(c2);
        }


    private:

        std::vector<Node> nodes;
};



#endif // TREE_H
