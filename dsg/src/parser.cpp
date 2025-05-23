#include "parser.h"
#include "graph.h"
#include <cstdio>
#include <cstdlib>

namespace {
    inline char GET_CHAR(FILE *f)
    {
        const int maxn = 131072;
        static char buf[maxn];
        static size_t p1 = 0, p2 = 0;
        static FILE *current = nullptr;

        if (f != current)
        { // Reset buffer if new file
            current = f;
            p1 = p2 = 0;
        }

        if (p1 == p2)
        {
            p2 = fread(buf, 1, maxn, f);
            p1 = 0;
        }
        return p1 == p2 ? EOF : buf[p1++];
    }

    inline int getInt(FILE *f)
    {
        int res = 0;
        char c = GET_CHAR(f);
        while (c < '0')
            c = GET_CHAR(f);
        while (c >= '0')
        {
            res = res * 10 + (c - '0');
            c = GET_CHAR(f);
        }
        return res;
    }
} // namespace

Graph Parser::read_graph_from_structured_file(const std::string &filename)
{
    FILE *f = fopen(filename.c_str(), "r");
    if (!f)
    {
        perror(("Failed to open file: " + filename).c_str());
        exit(1);
    }

    int n = getInt(f);
    int m = getInt(f);
    Graph graph(n, m);
    for (int i = 0; i < m; i++)
    {
        int u = getInt(f);
        int v = getInt(f);
        graph.add_edge(u, v);
    }

    fclose(f);
    return graph;
}
