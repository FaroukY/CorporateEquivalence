#include <math.h>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>

#include "incremental.h"

typedef long long int llint;

int readInt()
{
    int res = 0;
    char c = getchar_unlocked();
    while (c < '0' || c > '9')
    {
        c = getchar_unlocked();
    }

    while (c >= '0')
    {
        res = res * 10 + (c - '0');
        c = getchar_unlocked();
    }
    return res;
}

static char *
getNextWord(char *line, char *word)
{
    int wordlength = 0;

    while (((*line) == ' ') || ((*line) == '\t') || ((*line) == '\n') || ((*line) == '\0'))
    {
        ++line;
    }

    while (((*line) != ' ') && ((*line) != '\t') && ((*line) != '\n') && ((*line) != '\0'))
    {
        word[wordlength] = (*line);
        ++wordlength;
        ++line;
    }

    word[wordlength] = '\0';
    return line;
}

double
fmax(double a, double b)
{
    return (a > b) ? a : b;
}

double
fmin(double a, double b)
{
    return (a < b) ? a : b;
}

int max(int a, int b)
{
    return (a > b) ? a : b;
}

int min(int a, int b)
{
    return (a < b) ? a : b;
}

long long
lmax(long long a, long long b)
{
    return (a > b) ? a : b;
}

long long
lmin(long long a, long long b)
{
    return (a < b) ? a : b;
}

long long
ceil_div(long long a, long long b)
{
    return (a + b - 1) / b;
}

static double
timer(void)
{
    struct rusage r;

    getrusage(0, &r);
    return (double)(r.ru_utime.tv_sec + r.ru_utime.tv_usec / (double)1000000);
}

struct node;

typedef struct arc
{
    long long flow;
    long long capacity;
    long long intercept;
    long long slope;
    int from;
    int to;
    char direction;
} Arc;

typedef struct node
{
    struct node *parent;
    struct node *childList;
    struct node *nextScan;
    struct node *next;
    struct node *prev;
    Arc **outOfTree;
    Arc *arcToParent;
    long long excess;
    int numOutOfTree;
    int nextArc;
    int visited;
    int numAdjacent;
    int label;
    int breakpoint;
} Node;

typedef struct root
{
    Node *start;
    Node *end;
} Root;

// edge of the original graph
typedef struct edge
{
    int u;
    int v;
    int w;
} Edge;

//---------------  Global variables ------------------
static long long APP_VAL = 10000LL;
static int DECIMALS = 4;
static int n = 0;
static int m = 0;
static int numNodes = 0;
static int numArcs = 0;
static int source = 0;
static int sink = 0;
static int numParams = 0;
static int lastInternalArc = 0;
static int weightedEdges = 0;
static int weightedNodes = 0;

static double initTime = 0.0f;
static double injectedLambda = 0.0f;

long long currLambda = 0;

static int highestStrongLabel = 1;

static Node *adjacencyList = NULL;
static Root *strongRoots = NULL;
static int *labelCount = NULL;
static Arc *arcList = NULL;
static Edge *edges = NULL;

static double TOL = 1e-5;
static double totDegreeSourceSet = 0.0;
static double totWeightSourceSet = 0.0;

static int *degrees = NULL;
static int *adjacents = NULL;
static int *value = NULL;
static int *weights = NULL;
static char *inSourceSet = NULL;
static char *bestSourceSet = NULL;
static int *edgeEnd = NULL;
static int *firstEdge = NULL;

static FILE *dumpSourceSetFile = NULL;
//-----------------------------------------------------

#ifdef STATS
static llint numPushes = 0;
static int numMergers = 0;
static int numRelabels = 0;
static int numGaps = 0;
static llint numArcScans = 0;
#endif

static long long
computeNumerator(char *sourceSet)
{
    int i;
    long long res = 0;
    for (i = 0; i < lastInternalArc; ++i)
    {
        int from = arcList[i].from;
        int to = arcList[i].to;

        long long weight = arcList[i].capacity;

        res += weight * (sourceSet[from - 2] & sourceSet[to - 2]);
    }
#ifdef DIRECTED_CASE
    return res;
#else
    return res / 2;
#endif
}

static void
dumpSourceSet(FILE *out)
{
    int i;
    int sourceSetSize = 0;
    for (i = 2; i < numNodes; i++)
    {
        if (bestSourceSet[i - 2])
        {
            ++sourceSetSize;
            fprintf(out, "%d\n", i - 1);
        }
    }
}

static long long
computeSourceWeight(char *sourceSet)
{
    int i;
    long long res = 0;
    int sourceSetSize = 0;
    for (i = 2; i < numNodes; ++i)
    {

        if (sourceSet[i - 2])
        {
            res += weights[i - 2];
            ++sourceSetSize;
        }
        // if(weights[i-2] < 0 ) printf("c Error: weight[%d] is negative = %lf\n", i, weights[i-2]);
    }

    return res;
}

static void
copySourceSet(void)
{
    int i;
    for (i = 0; i < numNodes - 2; i++)
    {
        bestSourceSet[i] = inSourceSet[i];
    }
}

static long long
computeArcCapacity(const Arc *ac, long long param)
{
    if (ac->to == sink)
    {

        return lmax(ac->intercept - param * (ac->slope), 0);
    } else if (ac->from == source)
    {
        return lmax(param * (ac->slope) - ac->intercept, 0);
    }
    return ac->intercept;
}

static void
initializeNode(Node *nd, const int n)
{
    nd->label = 0;
    nd->excess = 0;
    nd->parent = NULL;
    nd->childList = NULL;
    nd->nextScan = NULL;
    nd->nextArc = 0;
    nd->numOutOfTree = 0;
    nd->arcToParent = NULL;
    nd->next = NULL;
    nd->prev = NULL;
    nd->visited = 0;
    nd->numAdjacent = 0;
    nd->outOfTree = NULL;
    nd->breakpoint = (numParams + 1);
}

static void
initializeRoot(Root *rt)
{
    rt->start = (Node *)malloc(sizeof(Node));
    rt->end = (Node *)malloc(sizeof(Node));

    if ((rt->start == NULL) || (rt->end == NULL))
    {
        printf("%s Line %d: Out of memory\n", __FILE__, __LINE__);
        exit(1);
    }

    initializeNode(rt->start, 0);
    initializeNode(rt->end, 0);

    rt->start->next = rt->end;
    rt->end->prev = rt->start;
}

static void
freeRoot(Root *rt)
{
    free(rt->start);
    rt->start = NULL;

    free(rt->end);
    rt->end = NULL;
}

static void
liftAll(Node *rootNode, const int theparam)
{
    Node *temp, *current = rootNode;

    current->nextScan = current->childList;

    --labelCount[current->label];
    current->label = numNodes;
    current->breakpoint = (theparam + 1);

    inSourceSet[(int)(current - adjacencyList) - 2] = 0;

    for (; (current); current = current->parent)
    {
        while (current->nextScan)
        {
            temp = current->nextScan;
            current->nextScan = current->nextScan->next;
            current = temp;
            current->nextScan = current->childList;

            inSourceSet[(int)(current - adjacencyList) - 2] = 0;
            --labelCount[current->label];
            current->label = numNodes;
            current->breakpoint = (theparam + 1);
        }
    }
}

static void
addToStrongBucket(Node *newRoot, Node *rootEnd)
{
    newRoot->next = rootEnd;
    newRoot->prev = rootEnd->prev;
    rootEnd->prev = newRoot;
    newRoot->prev->next = newRoot;
}

static void
createOutOfTree(Node *nd)
{
    if (nd->numAdjacent)
    {
        if ((nd->outOfTree = (Arc **)malloc(nd->numAdjacent * sizeof(Arc *))) == NULL)
        {
            printf("%s Line %d: Out of memory\n", __FILE__, __LINE__);
            exit(1);
        }
    }
}

static void
initializeArc(Arc *ac)
{
    int i;

    ac->from = 0;
    ac->to = 0;
    ac->capacity = 0;
    ac->flow = 0;
    ac->direction = 1;
    ac->intercept = 0;
    ac->slope = 0;
    // ac->capacities = NULL;
}

static void
addOutOfTreeNode(Node *n, Arc *out)
{
    n->outOfTree[n->numOutOfTree] = out;
    ++n->numOutOfTree;
}

int cmp_deg(const void *lhs, const void *rhs)
{
    int a = *((int *)lhs);
    int b = *((int *)rhs);

    return (firstEdge[b + 1] - firstEdge[b]) - (firstEdge[a + 1] - firstEdge[a]);
}

int cmp_edge(const void *lhs, const void *rhs)
{
    Edge *el = (Edge *)lhs;
    Edge *er = (Edge *)rhs;

    if (el->u < er->u)
        return -1;
    if (el->u > er->u)
        return 1;

    if (el->v < er->v)
        return -1;
    if (el->v > er->v)
        return 1;

    throw std::runtime_error("Error: edges are equal");
}

#ifdef USE_HEURISTIC
static long long
computeHeuristicLambda()
{

    long long numerator = 0;
    long long denominator = 0;

    int i;
    int u;
    int v;
    int e;

    int *order = NULL;
    long long bestLambda = 0;
    if ((order = (int *)malloc((numNodes - 2) * sizeof(int))) == NULL)
    {
        printf("%s, %d: could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }

    for (i = 0; i < numNodes - 2; i++)
    {
        inSourceSet[u] = 0;
        order[i] = i;
    }

    qsort(order, numNodes - 2, sizeof(int), cmp_deg);

    for (i = 0; i < (numNodes - 2); ++i)
    {
        u = order[i];

        if (inSourceSet[u])
            continue;
        inSourceSet[u] = 1;

        denominator += weights[u];

        for (e = firstEdge[u]; edges[e].u == u; ++e)

        {

            v = edges[e].v; // edgeEnd[e];
            if (!inSourceSet[v])
                continue;

            numerator += edges[e].w * APP_VAL;
        }

        long long thisLambda = ceil_div(numerator, denominator);

        bestLambda = lmax(thisLambda, bestLambda);
    }
    free(order);

    return bestLambda;
}
#endif

static void
readGraphFile(void)
{
    int i;

    n = 0;
    m = 0;
    currLambda = 0;
#ifdef FAST_IO
    n = readInt();
    m = readInt();
#else
    scanf("%d %d", &n, &m);
#endif

    printf("c n=%d m=%d\n", n, m);
    source = 0;
    sink = 1;

    numNodes = n + 2;
    numArcs = m + 2 * n;

    adjacents = NULL;

    if ((adjacents = (int *)malloc(numNodes * sizeof(int))) == NULL)
    {
        printf("%s, %d: could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }

    if ((weights = (int *)malloc(n * sizeof(int))) == NULL)
    {
        printf("%s, %d: could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }

    if ((degrees = (int *)malloc(n * sizeof(int))) == NULL)
    {
        printf("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }
#ifdef USE_HEURISTIC

#ifdef DIRECTED_CASE
    int totalEdges = 2 * m;
#else
    int totalEdges = m;
#endif

    if ((edges = (Edge *)malloc((totalEdges + 1) * sizeof(Edge))) == NULL)
    {
        printf("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }

    if ((firstEdge = (int *)malloc((n + 1) * sizeof(int))) == NULL)
    {
        printf("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }
#endif

    if ((inSourceSet = (char *)malloc(n * sizeof(char))) == NULL)
    {
        printf("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }

    if ((bestSourceSet = (char *)malloc(n * sizeof(char))) == NULL)
    {
        printf("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }

    if ((adjacencyList = (Node *)malloc(numNodes * sizeof(Node))) == NULL)
    {
        printf("%s, %d: could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }

    if ((strongRoots = (Root *)malloc(numNodes * sizeof(Root))) == NULL)
    {
        printf("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }

    if ((labelCount = (int *)malloc(numNodes * sizeof(int))) == NULL)
    {
        printf("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }

    if ((arcList = (Arc *)malloc(numArcs * sizeof(Arc))) == NULL)
    {
        printf("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
        exit(1);
    }

    for (i = 0; i < n; ++i)
    {

        if (weightedNodes)
        {
            int weight = 0;
#ifdef FAST_IO
            weight = readInt();
#else
            scanf("%d", &weight);
#endif
            weights[i] = weight;
        } else
        {

            weights[i] = 1;
        }
        // we consider everything to be in the source set at start
#ifndef USE_HEURISTIC
        inSourceSet[i] = 1;
#else
        inSourceSet[i] = 0;
#endif

        degrees[i] = 0;
    }

    // Initialize lambda as

    for (i = 0; i < numNodes; ++i)
    {
        initializeRoot(&strongRoots[i]);
        initializeNode(&adjacencyList[i], i);
        adjacents[i] = 0;
        labelCount[i] = 0;
    }

    for (i = 0; i < numArcs; ++i)
    {
        initializeArc(&arcList[i]);
    }

    int currEdge = 0;
    int currArc = 0;
    Arc *ac = NULL;
    // add arc for each edge
    for (i = 0; i < m; i++)
    {

        int u = 0;
        int v = 0;
        int w = 1;
#ifdef FAST_IO
        u = readInt();
        v = readInt();
        if (weightedEdges)
        {
            w = readInt();
        }
#else
        scanf("%d %d", &u, &v);
        if (weightedEdges)
        {
            scanf("%d", &w);
        }
#endif

        degrees[u - 1] += w;
#ifdef USE_HEURISTIC
        edges[currEdge].u = u - 1;
        edges[currEdge].v = v - 1;
        edges[currEdge++].w = w;
#ifdef DIRECTED_CASE

        edges[currEdge].v = u - 1;
        edges[currEdge].u = v - 1;
        edges[currEdge++].w = w;
#endif
#endif
        // we invert the orientation of the network
        // to make source adjacent monotone increasing
        // and sink adjacent monotone decreasing
        int from = v + 1;
        int to = u + 1;

        ac = &arcList[currArc++];

        ac->from = from;
        ac->to = to;
        ac->intercept = w * APP_VAL;
        ac->capacity = w * APP_VAL;

        ++adjacents[ac->from];
        ++adjacents[ac->to];
    }
    lastInternalArc = currArc;

    // dummy edge as end
#ifdef USE_HEURISTIC
    edges[currEdge].u = n + 10;
    firstEdge[n] = totalEdges;
    qsort(edges, totalEdges, sizeof(Edge), cmp_edge);

    for (i = totalEdges - 1; i >= 0; i--)
    {
        int u = edges[i].u;
        firstEdge[u] = i;
    }

#endif
    // add arcs (source, node) (node, sink)
    for (i = 0; i < n; i++)
    {
        ac = &arcList[currArc++];

        int from = source;
        int to = i + 2;
        int weight = weights[i];
        int degree = degrees[i];
        ac->from = from;
        ac->to = to;
        ac->intercept = degree * APP_VAL;

#ifdef DIRECTED_CASE
        ac->slope = weight;
#else
        ac->slope = 2 * weight;
#endif

        ++adjacents[ac->from];
        ++adjacents[ac->to];
        ac = &arcList[currArc++];

        from = i + 2;
        to = sink;
        weight = weights[i];
        degree = degrees[i];
        ac->from = from;
        ac->to = to;
        ac->intercept = degree * APP_VAL;
#ifdef DIRECTED_CASE
        ac->slope = weight;
#else
        ac->slope = 2 * weight;
#endif

        ++adjacents[ac->from];
        ++adjacents[ac->to];
    }

    for (i = 0; i < numNodes; ++i)
    {
        adjacencyList[i].numAdjacent = adjacents[i];
    }

#ifdef USE_HEURISTIC

    currLambda = computeHeuristicLambda();

#else
    currLambda = computeNumerator(inSourceSet) / computeSourceWeight(inSourceSet);
#endif

    currLambda = max(currLambda, (long long)(injectedLambda * APP_VAL));

    printf("c Initial lambda = %.*lf\n", DECIMALS, (double)currLambda / APP_VAL);

    for (i = 0; i < numArcs; ++i)
    {
        ac = &arcList[i];
        ac->capacity = computeArcCapacity(ac, currLambda);
    }

    for (i = 0; i < numNodes; ++i)
    {
        createOutOfTree(&adjacencyList[i]);
    }

    for (i = 0; i < numArcs; i++)
    {
        int to = arcList[i].to;
        int from = arcList[i].from;
        long long capacity = arcList[i].capacity;

        if (!((source == to) || (sink == from) || (from == to)))
        {
            if ((source == from) && (to == sink))
            {
                arcList[i].flow = capacity;
            } else if (from == source)
            {
                addOutOfTreeNode(&adjacencyList[from], &arcList[i]);
            } else if (to == sink)
            {
                addOutOfTreeNode(&adjacencyList[to], &arcList[i]);
            } else
            {
                addOutOfTreeNode(&adjacencyList[from], &arcList[i]);
            }
        }
    }

    if (adjacents != NULL)
    {
        free(adjacents);
    }
    adjacents = NULL;

#ifdef USE_HEURISTIC
    if (edges != NULL)
    {
        free(edges);
    }

    if (firstEdge != NULL)
    {
        free(firstEdge);
    }
#endif
}

static void
simpleInitialization(void)
{
    int i, size;
    Arc *tempArc;

    size = adjacencyList[source].numOutOfTree;
    for (i = 0; i < size; ++i)
    {
        tempArc = adjacencyList[source].outOfTree[i];
        tempArc->flow = tempArc->capacity;
        adjacencyList[tempArc->to].excess += tempArc->capacity;
    }

    size = adjacencyList[sink].numOutOfTree;
    for (i = 0; i < size; ++i)
    {
        tempArc = adjacencyList[sink].outOfTree[i];
        tempArc->flow = tempArc->capacity;
        adjacencyList[tempArc->from].excess -= tempArc->capacity;
    }

    adjacencyList[source].excess = 0;
    adjacencyList[sink].excess = 0;

    for (i = 0; i < numNodes; ++i)
    {
        if (adjacencyList[i].excess > 0)
        {
            adjacencyList[i].label = 1;
            ++labelCount[1];

            addToStrongBucket(&adjacencyList[i], strongRoots[1].end);
        }
    }

    adjacencyList[source].label = numNodes;
    adjacencyList[source].breakpoint = 0;
    adjacencyList[sink].label = 0;
    adjacencyList[sink].breakpoint = (numParams + 2);
    labelCount[0] = (numNodes - 2) - labelCount[1];
}

static inline int
addRelationship(Node *newParent, Node *child)
{
    child->parent = newParent;
    child->next = newParent->childList;

    // mi cambio

    if (newParent->childList != NULL)
    {
        newParent->childList->prev = child;
    }

    // end mi cambio

    newParent->childList = child;

    return 0;
}

static inline void
breakRelationship(Node *oldParent, Node *child)
{
    Node *current;

    child->parent = NULL;

    if (oldParent->childList == child)
    {
        oldParent->childList = child->next;
        child->next = NULL;
        return;
    }

    // for (current = oldParent->childList; (current->next != child); current = current->next);
    current = child->prev;

    current->next = child->next;
    if (child->next != NULL)
        child->next->prev = current;
    child->next = NULL;
}

static void
merge(Node *parent, Node *child, Arc *newArc)
{
    Arc *oldArc;
    Node *current = child, *oldParent, *newParent = parent;

#ifdef STATS
    ++numMergers;
#endif

    while (current->parent)
    {
        oldArc = current->arcToParent;
        current->arcToParent = newArc;
        oldParent = current->parent;
        breakRelationship(oldParent, current);
        addRelationship(newParent, current);
        newParent = current;
        current = oldParent;
        newArc = oldArc;
        newArc->direction = 1 - newArc->direction;
    }

    current->arcToParent = newArc;
    addRelationship(newParent, current);
}

static inline void
pushUpward(Arc *currentArc, Node *child, Node *parent, const long long resCap)
{
#ifdef STATS
    ++numPushes;
#endif

    if (resCap >= child->excess)
    {
        parent->excess += child->excess;
        currentArc->flow += child->excess;
        child->excess = 0;
        return;
    }

    currentArc->direction = 0;
    parent->excess += resCap;
    child->excess -= resCap;
    currentArc->flow = currentArc->capacity;
    parent->outOfTree[parent->numOutOfTree] = currentArc;
    ++parent->numOutOfTree;
    breakRelationship(parent, child);

    addToStrongBucket(child, strongRoots[child->label].end);
}

static inline void
pushDownward(Arc *currentArc, Node *child, Node *parent, long long flow)
{
#ifdef STATS
    ++numPushes;
#endif

    if (flow >= child->excess)
    {
        parent->excess += child->excess;
        currentArc->flow -= child->excess;
        child->excess = 0;
        return;
    }

    currentArc->direction = 1;
    child->excess -= flow;
    parent->excess += flow;
    currentArc->flow = 0;
    parent->outOfTree[parent->numOutOfTree] = currentArc;
    ++parent->numOutOfTree;
    breakRelationship(parent, child);

    addToStrongBucket(child, strongRoots[child->label].end);
}

static void
pushExcess(Node *strongRoot)
{
    Node *current, *parent;
    Arc *arcToParent;

    for (current = strongRoot; (current->excess && current->parent); current = parent)
    {
        parent = current->parent;
        arcToParent = current->arcToParent;
        if (arcToParent->direction)
        {
            pushUpward(arcToParent, current, parent, (arcToParent->capacity - arcToParent->flow));
        } else
        {
            pushDownward(arcToParent, current, parent, arcToParent->flow);
        }
    }

    if (current->excess > 0)
    {
        if (!current->next)
        {
            addToStrongBucket(current, strongRoots[current->label].end);
        }
    }
}

static Arc *
findWeakNode(Node *strongNode, Node **weakNode)
{
    int i, size;
    Arc *out;

    size = strongNode->numOutOfTree;

    for (i = strongNode->nextArc; i < size; ++i)
    {

#ifdef STATS
        ++numArcScans;
#endif

        if (adjacencyList[strongNode->outOfTree[i]->to].label == (highestStrongLabel - 1))
        {
            strongNode->nextArc = i;
            out = strongNode->outOfTree[i];
            (*weakNode) = &adjacencyList[out->to];
            --strongNode->numOutOfTree;
            strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
            return (out);
        } else if (adjacencyList[strongNode->outOfTree[i]->from].label == (highestStrongLabel - 1))
        {
            strongNode->nextArc = i;
            out = strongNode->outOfTree[i];
            (*weakNode) = &adjacencyList[out->from];
            --strongNode->numOutOfTree;
            strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
            return (out);
        }
    }

    strongNode->nextArc = strongNode->numOutOfTree;

    return NULL;
}

static void
checkChildren(Node *curNode)
{
    for (; (curNode->nextScan); curNode->nextScan = curNode->nextScan->next)
    {
        if (curNode->nextScan->label == curNode->label)
        {
            return;
        }
    }

    --labelCount[curNode->label];
    ++curNode->label;
    ++labelCount[curNode->label];

#ifdef STATS
    ++numRelabels;
#endif

    curNode->nextArc = 0;
}

static void
processRoot(Node *strongRoot)
{
    Node *temp, *strongNode = strongRoot, *weakNode;
    Arc *out;

    strongRoot->nextScan = strongRoot->childList;

    if ((out = findWeakNode(strongRoot, &weakNode)))
    {
        merge(weakNode, strongNode, out);
        pushExcess(strongRoot);
        return;
    }

    checkChildren(strongRoot);

    while (strongNode)
    {
        while (strongNode->nextScan)
        {
            temp = strongNode->nextScan;
            strongNode->nextScan = strongNode->nextScan->next;
            strongNode = temp;
            strongNode->nextScan = strongNode->childList;

            if ((out = findWeakNode(strongNode, &weakNode)))
            {
                merge(weakNode, strongNode, out);
                pushExcess(strongRoot);
                return;
            }

            checkChildren(strongNode);
        }

        if ((strongNode = strongNode->parent))
        {
            checkChildren(strongNode);
        }
    }

    addToStrongBucket(strongRoot, strongRoots[strongRoot->label].end);

    ++highestStrongLabel;
}

static Node *
getHighestStrongRoot(const int theparam)
{
    int i;
    Node *strongRoot;

    for (i = highestStrongLabel; i > 0; --i)
    {
        if (strongRoots[i].start->next != strongRoots[i].end)
        {
            highestStrongLabel = i;
            if (labelCount[i - 1])
            {
                strongRoot = strongRoots[i].start->next;
                strongRoot->next->prev = strongRoot->prev;
                strongRoot->prev->next = strongRoot->next;
                strongRoot->next = NULL;
                return strongRoot;
            }

            while (strongRoots[i].start->next != strongRoots[i].end)
            {

#ifdef STATS
                ++numGaps;
#endif
                strongRoot = strongRoots[i].start->next;
                strongRoot->next->prev = strongRoot->prev;
                strongRoot->prev->next = strongRoot->next;
                liftAll(strongRoot, theparam);
            }
        }
    }

    if (strongRoots[0].start->next == strongRoots[0].end)
    {
        return NULL;
    }

    while (strongRoots[0].start->next != strongRoots[0].end)
    {
        strongRoot = strongRoots[0].start->next;
        strongRoot->next->prev = strongRoot->prev;
        strongRoot->prev->next = strongRoot->next;

        strongRoot->label = 1;
        --labelCount[0];
        ++labelCount[1];

#ifdef STATS
        ++numRelabels;
#endif

        addToStrongBucket(strongRoot, strongRoots[strongRoot->label].end);
    }

    highestStrongLabel = 1;

    strongRoot = strongRoots[1].start->next;
    strongRoot->next->prev = strongRoot->prev;
    strongRoot->prev->next = strongRoot->next;
    strongRoot->next = NULL;

    return strongRoot;
}

static void
updateCapacities(const int theparam)
{
    int i, size;
    long long delta;
    Arc *tempArc;
    Node *tempNode;

    long long param = currLambda; /// lambdaStart + ((lambdaEnd-lambdaStart) *  ((double)theparam/(numParams-1)));
    size = adjacencyList[source].numOutOfTree;
    for (i = 0; i < size; ++i)
    {
        tempArc = adjacencyList[source].outOfTree[i];
        // delta = (tempArc->capacities[theparam] - tempArc->capacity);
        delta = (computeArcCapacity(tempArc, param) - tempArc->capacity);
        // delta = overestimate(delta);
        if (delta < 0)
        {
            printf("c Error on source-adjacent arc (%d, %d): capacity decreases by %lld at parameter %d.\n",
                   tempArc->from,
                   tempArc->to,
                   (-delta),
                   (theparam + 1));
            exit(0);
        }

        tempArc->capacity += delta;
        tempArc->flow += delta;
        adjacencyList[tempArc->to].excess += delta;
        // tempArc->capacity = overestimate(tempArc->capacity);
        // tempArc->flow = overestimate(tempArc->flow);
        // adjacencyList[tempArc->to].excess = overestimate(adjacencyList[tempArc->to].excess);

        if ((adjacencyList[tempArc->to].label < numNodes) && (adjacencyList[tempArc->to].excess > 0))
        {
            pushExcess(&adjacencyList[tempArc->to]);
        }
    }

    size = adjacencyList[sink].numOutOfTree;
    for (i = 0; i < size; ++i)
    {
        tempArc = adjacencyList[sink].outOfTree[i];
        // delta = (tempArc->capacities[theparam] - tempArc->capacity);
        delta = (computeArcCapacity(tempArc, param) - tempArc->capacity);
        // delta = overestimate(delta);
        if (delta > 0)
        {
            printf("c Error on sink-adjacent arc (%d, %d): capacity %lld increases to %lld at parameter %d.\n",
                   tempArc->from,
                   tempArc->to,
                   tempArc->capacity,
                   tempArc->capacity + delta,
                   (theparam + 1));
            exit(0);
        }

        tempArc->capacity += delta;
        tempArc->flow += delta;
        adjacencyList[tempArc->from].excess -= delta;

        // tempArc->capacity = overestimate(tempArc->capacity);
        // tempArc->flow = overestimate(tempArc->flow);
        // adjacencyList[tempArc->from].excess = overestimate(adjacencyList[tempArc->from].excess);
        if ((adjacencyList[tempArc->from].label < numNodes) && (adjacencyList[tempArc->from].excess > 0))
        {
            pushExcess(&adjacencyList[tempArc->from]);
        }
    }

    highestStrongLabel = (numNodes - 1);
}

static long long
computeMinCut(void)
{
    int i;
    long long mincut = 0;

    for (i = 0; i < numArcs; ++i)
    {
        if ((adjacencyList[arcList[i].from].label >= numNodes) && (adjacencyList[arcList[i].to].label < numNodes))
        {
            mincut += arcList[i].capacity;
        }
    }
    return mincut;
}

static void
incrementalCut(Stats &stats)
{
    Node *strongRoot;

    double thetime;
    int theparam = 0;
    thetime = timer();
    stats.start_timer();

    // currLambda = overestimate(currLambda); //ceil( currLambda * APP_VAL ) / APP_VAL;
    long long initLambda = currLambda;
    long long bestLambda = currLambda;
    copySourceSet();
    int iteration = 0;
    printf("c Iteration %d\n", iteration++);
    printf("c Computing mincut for lambda = %.*lf\n", DECIMALS, (double)currLambda / APP_VAL);
    while ((strongRoot = getHighestStrongRoot(theparam)))
    {
        processRoot(strongRoot);
    }

    long long css = computeNumerator(inSourceSet);
    long long qs = computeSourceWeight(inSourceSet);
    // currLambda = css / qs;

    // printf("c C(S,S)/q(S) = %lld/%lld = %lld\n",css, qs, ceil_div(css, qs));
    long long val = css - currLambda * qs;

    printf("C(S,S) - lambda * q(S) = %lf\n", (double)val / APP_VAL);
    stats.pause_timer();

    auto best_density = (double)bestLambda / APP_VAL;
    stats.push(iteration, best_density, 0.0 /* load_norm not applicable */);

    while (val > 0)
    {
        stats.start_timer();
        // compute new lambda

        printf("c Iteration %d\n", iteration++);
        currLambda = ceil_div(css, qs);

        if (currLambda > bestLambda)
        {
            bestLambda = currLambda;
            copySourceSet();
        }

#ifdef ANYTIME
        double totalTime = timer() - thetime;
        printf("c %d,%d,%.4lf,%.4lf,%.*lf,%.*lf,%d\n", n, m, initTime, totalTime, DECIMALS, (double)initLambda / APP_VAL, DECIMALS, (double)bestLambda / APP_VAL, iteration);
#endif

        printf("c Updating capacities for lambda = %.*lf\n", DECIMALS, (double)currLambda / APP_VAL);
        updateCapacities(theparam);

        while ((strongRoot = getHighestStrongRoot(theparam)))
        {
            processRoot(strongRoot);
        }
        css = computeNumerator(inSourceSet);
        qs = computeSourceWeight(inSourceSet);
        // mincut = computeMinCut();

        val = css - currLambda * qs;
        printf("C(S,S) - lambda * q(S) = %lf\n", (double)val / APP_VAL);

        // Store stats
        stats.pause_timer();

        const auto best_density = (double)bestLambda / APP_VAL;
        stats.push(iteration, best_density, 0.0 /* load_norm not applicable */);
    }

    if (val < -TOL)
    {
        printf("c ERROR: value is negative\n");
    }

    // css = computeNumerator();
    // qs = computeSourceWeight();
    if (qs != 0)
    {
        currLambda = css / qs;
        bestLambda = fmax(currLambda, bestLambda);
    }
    double totalTime = timer() - thetime;

    css = computeNumerator(bestSourceSet);
    qs = computeSourceWeight(bestSourceSet);
    best_density = (double)bestLambda / APP_VAL;
    // printf("c %d,%d,%.4lf,%.4lf,%.*lf,%.*lf,%d\n", n, m, initTime, totalTime, DECIMALS, (double)initLambda/APP_VAL, DECIMALS, (double)bestLambda/APP_VAL,iteration);
    printf("\n");
    printf("c Optimal solution found: density* = %.*lf\n", DECIMALS, best_density);
    printf("c C(S, S) = %lld, q(S) = %lld\n", css / APP_VAL, qs);
    if (dumpSourceSetFile != NULL)
    {
        dumpSourceSet(dumpSourceSetFile);
        fclose(dumpSourceSetFile);
    } else
    {
        printf("\n");
        printf("c To dump the list of nodes in the densest subgraph into\n");
        printf("c a file, use option \"--dumpDensest FILE\"\n");
    }
}

static void
checkOptimality(void)
{
    int i, check = 1;
    long long mincut = 0, *excess;

    excess = (long long *)malloc(numNodes * sizeof(long long));
    if (!excess)
    {
        printf("%s Line %d: Out of memory\n", __FILE__, __LINE__);
        exit(1);
    }

    for (i = 0; i < numNodes; ++i)
    {
        excess[i] = 0;
    }

    for (i = 0; i < numArcs; ++i)
    {
        if ((adjacencyList[arcList[i].from].label >= numNodes) && (adjacencyList[arcList[i].to].label < numNodes))
        {
            mincut += arcList[i].capacity;
        }

        if ((arcList[i].flow > arcList[i].capacity) || (arcList[i].flow < 0))
        {
            check = 0;
            printf("c Capacity constraint violated on arc (%d, %d)\n",
                   arcList[i].from,
                   arcList[i].to);
        }
        excess[arcList[i].from] -= arcList[i].flow;
        excess[arcList[i].to] += arcList[i].flow;
    }

    for (i = 0; i < numNodes; i++)
    {
        if ((i != (source)) && (i != (sink)))
        {
            if (excess[i])
            {
                check = 0;
                printf("c Flow balance constraint violated in node %d. Excess = %lld\n",
                       i + 1,
                       excess[i]);
            }
        }
    }

    if (check)
    {
        printf("c\nc Solution checks as feasible.\n");
    }

    check = 1;

    if (excess[sink] != mincut)
    {
        check = 0;
        printf("c Flow is not optimal - max flow does not equal min cut!\nc\n");
    }

    if (check)
    {
        printf("c\nc Solution checks as optimal.\nc \n");
        printf("s Max Flow            : %lld\n", mincut);
    }

    free(excess);
    excess = NULL;
}

static void
quickSort(Arc **arr, const int first, const int last)
{
    int i, j, left = first, right = last, x1, x2, x3, mid, pivot, pivotval;
    Arc *swap;

    if ((right - left) <= 5)
    { // Bubble sort if 5 elements or less
        for (i = right; (i > left); --i)
        {
            swap = NULL;
            for (j = left; j < i; ++j)
            {
                if (arr[j]->flow < arr[j + 1]->flow)
                {
                    swap = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = swap;
                }
            }

            if (!swap)
            {
                return;
            }
        }

        return;
    }

    mid = (first + last) / 2;

    x1 = arr[first]->flow;
    x2 = arr[mid]->flow;
    x3 = arr[last]->flow;

    pivot = mid;

    if (x1 <= x2)
    {
        if (x2 > x3)
        {
            pivot = left;

            if (x1 <= x3)
            {
                pivot = right;
            }
        }
    } else
    {
        if (x2 <= x3)
        {
            pivot = right;

            if (x1 <= x3)
            {
                pivot = left;
            }
        }
    }

    pivotval = arr[pivot]->flow;

    swap = arr[first];
    arr[first] = arr[pivot];
    arr[pivot] = swap;

    left = (first + 1);

    while (left < right)
    {
        if (arr[left]->flow < pivotval)
        {
            swap = arr[left];
            arr[left] = arr[right];
            arr[right] = swap;
            --right;
        } else
        {
            ++left;
        }
    }

    swap = arr[first];
    arr[first] = arr[left];
    arr[left] = swap;

    if (first < (left - 1))
    {
        quickSort(arr, first, (left - 1));
    }

    if ((left + 1) < last)
    {
        quickSort(arr, (left + 1), last);
    }
}

static void
sort(Node *current)
{
    if (current->numOutOfTree > 1)
    {
        quickSort(current->outOfTree, 0, (current->numOutOfTree - 1));
    }
}

static void
minisort(Node *current)
{
    Arc *temp = current->outOfTree[current->nextArc];
    int i, size = current->numOutOfTree, tempflow = temp->flow;

    for (i = current->nextArc + 1; ((i < size) && (tempflow < current->outOfTree[i]->flow)); ++i)
    {
        current->outOfTree[i - 1] = current->outOfTree[i];
    }
    current->outOfTree[i - 1] = temp;
}

static void
decompose(Node *excessNode, const int source, int *iteration)
{
    Node *current = excessNode;
    Arc *tempArc;
    long long bottleneck = excessNode->excess;

    for (; (((int)(current - adjacencyList)) != source) && (current->visited < (*iteration));
         current = &adjacencyList[tempArc->from])
    {
        current->visited = (*iteration);
        tempArc = current->outOfTree[current->nextArc];

        if (tempArc->flow < bottleneck)
        {
            bottleneck = tempArc->flow;
        }
    }

    if ((int)((int)(current - adjacencyList)) == source)
    {
        excessNode->excess -= bottleneck;
        current = excessNode;

        while ((int)((int)(current - adjacencyList)) != source)
        {
            tempArc = current->outOfTree[current->nextArc];
            tempArc->flow -= bottleneck;

            if (tempArc->flow)
            {
                minisort(current);
            } else
            {
                ++current->nextArc;
            }
            current = &adjacencyList[tempArc->from];
        }
        return;
    }

    ++(*iteration);

    bottleneck = current->outOfTree[current->nextArc]->flow;

    while (current->visited < (*iteration))
    {
        current->visited = (*iteration);
        tempArc = current->outOfTree[current->nextArc];

        if (tempArc->flow < bottleneck)
        {
            bottleneck = tempArc->flow;
        }
        current = &adjacencyList[tempArc->from];
    }

    ++(*iteration);

    while (current->visited < (*iteration))
    {
        current->visited = (*iteration);

        tempArc = current->outOfTree[current->nextArc];
        tempArc->flow -= bottleneck;

        if (tempArc->flow)
        {
            minisort(current);
            current = &adjacencyList[tempArc->from];
        } else
        {
            ++current->nextArc;
            current = &adjacencyList[tempArc->from];
        }
    }
}

static void
recoverFlow(void)
{
    int i, j, iteration = 1;
    Arc *tempArc;
    Node *tempNode;

    for (i = 0; i < adjacencyList[sink].numOutOfTree; ++i)
    {
        tempArc = adjacencyList[sink].outOfTree[i];
        if (adjacencyList[tempArc->from].excess < 0)
        {
            tempArc->flow -= (long long)(-1 * adjacencyList[tempArc->from].excess);
            adjacencyList[tempArc->from].excess = 0;
        }
    }

    for (i = 0; i < adjacencyList[source].numOutOfTree; ++i)
    {
        tempArc = adjacencyList[source].outOfTree[i];
        addOutOfTreeNode(&adjacencyList[tempArc->to], tempArc);
    }

    adjacencyList[source].excess = 0;
    adjacencyList[sink].excess = 0;

    for (i = 0; i < numNodes; ++i)
    {
        tempNode = &adjacencyList[i];

        if ((i == (source)) || (i == (sink)))
        {
            continue;
        }

        if (tempNode->label >= numNodes)
        {
            tempNode->nextArc = 0;
            if ((tempNode->parent) && (tempNode->arcToParent->flow))
            {
                addOutOfTreeNode(&adjacencyList[tempNode->arcToParent->to], tempNode->arcToParent);
            }

            for (j = 0; j < tempNode->numOutOfTree; ++j)
            {
                if (!tempNode->outOfTree[j]->flow)
                {
                    --tempNode->numOutOfTree;
                    tempNode->outOfTree[j] = tempNode->outOfTree[tempNode->numOutOfTree];
                    --j;
                }
            }

            sort(tempNode);
        }
    }

    for (i = 0; i < numNodes; ++i)
    {
        tempNode = &adjacencyList[i];
        while (tempNode->excess > 0)
        {
            ++iteration;
            decompose(tempNode, source, &iteration);
        }
    }
}

static void
displayBreakpoints(void)
{
    int i;
    for (i = 0; i < numNodes; ++i)
    {
        printf("n %d %d\n", (i + 1), adjacencyList[i].breakpoint);
    }
}

static void
freeMemory(void)
{
    int i;

    for (i = 0; i < numNodes; ++i)
    {
        freeRoot(&strongRoots[i]);
    }

    free(strongRoots);

    for (i = 0; i < numNodes; ++i)
    {
        if (adjacencyList[i].outOfTree)
        {
            free(adjacencyList[i].outOfTree);
        }
    }

    if (adjacencyList != NULL)
    {
        free(adjacencyList);
    }

    if (labelCount != NULL)
    {
        free(labelCount);
    }

    if (arcList != NULL)
    {
        free(arcList);
    }

    if (degrees != NULL)
    {
        free(degrees);
    }
    if (weights != NULL)
    {
        free(weights);
    }
    if (inSourceSet != NULL)
    {
        free(inSourceSet);
    }
    if (bestSourceSet != NULL)
    {
        free(bestSourceSet);
    }
}

void printHelp(int argc, char *argv[])
{
    printf("Usage : %s [OPTIONS]                                        \n\n", argv[0]);
    printf("Options:\n");

    printf("  --help                 print this help message                \n");
    printf("  --accuracy N           set the accuracy to N decimal places   \n");
    printf("                         note that a big N may induce errors on \n");
    printf("                         graphs that are big or have big weights\n");
    printf("                         by default, 4 decimal places           \n");
    printf("  --startLambda L        start incremental procedure with L     \n");
    printf("  --weightedEdges        enable if input has weights on edges   \n");
    printf("  --weightedNodes        enable if input has weights on nodes   \n");
    printf("  --dumpDensest FILE     dump list of nodes in optimal solution \n");
    printf("                         to FILE                                \n");

    printf("\n\nNote that this program read from stdin\n");
    printf("An example command is \'%s < mygraph.txt\' to read from file\n", argv[0]);
    printf("\n");
    printf("The format of the input is as follows: \n");
    printf("\n");
    printf("N M\n");
    printf("[WEIGHT_NODE_1]\n");
    printf("[WEIGHT_NODE_2]\n");
    printf("...\n");
    printf("[WEIGHT_NODE_N]\n");
    printf("FROM_EDGE_1 TO_EDGE_1 [WEIGHT_EDGE_1]\n");
    printf("FROM_EDGE_2 TO_EDGE_2 [WEIGHT_EDGE_2]\n");
    printf("...\n");
    printf("FROM_EDGE_M TO_EDGE_M [WEIGHT_EDGE_M]\n");
    printf("\n");
    printf("Remember that for weighted graphs, options --weightedEdges \n");
    printf("and/or --weightedNodes must be used.  \n");
}

void parseParameters(int argc, char *argv[])
{
    int i = 0;
    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--help") == 0)
        {
            printHelp(argc, argv);
            exit(0);
        } else if (strcmp(argv[i], "--accuracy") == 0)
        {
            int accuracy = atoi(argv[++i]);
            DECIMALS = accuracy;
            int z = 0;
            APP_VAL = 1;
            for (z = 0; z < accuracy; z++)
            {
                APP_VAL *= 10LL;
            }
        } else if (strcmp(argv[i], "--startLambda") == 0)
        {
            injectedLambda = atof(argv[++i]);
        } else if (strcmp(argv[i], "--weightedEdges") == 0)
        {
            weightedEdges = 1;
        } else if (strcmp(argv[i], "--weightedNodes") == 0)
        {
            weightedNodes = 1;
        } else if (strcmp(argv[i], "--dumpDensest") == 0)
        {
            dumpSourceSetFile = fopen(argv[++i], "w");
        } else
        {
            printf("%s: invalid argument %s\n", argv[0], argv[1]);
            printf("Try \'%s --help\' for more information.\n", argv[0]);
            exit(-1);
        }
    }
}

int run_incremental(int argc, char **argv, Stats &stats)
{
    parseParameters(argc, argv);
    printf("c Incremental cut procedure for maximum subgraph density\n");

    double thetime = timer();
    readGraphFile();

#ifdef PROGRESS
    printf("c Finished reading file.\n");
    fflush(stdout);
#endif

    stats.start_timer();
    simpleInitialization();
    stats.pause_timer();
    stats.push(0, 1.0 * m / n, 0.0 /* load_norm not applicable */);

#ifdef PROGRESS
    printf("c Finished initialization.\n");
    fflush(stdout);
#endif

    initTime = timer() - thetime;
    incrementalCut(stats);

#ifdef PROGRESS
    printf("c Finished phase 1.\n");
    fflush(stdout);
#endif

#ifdef RECOVER_FLOW
    recoverFlow();
    checkOptimality();
#endif

    // printf ("c Number of nodes     : %d\n", numNodes);
    // printf ("c Number of arcs      : %d\n", numArcs);
#ifdef STATS
    printf("c Number of arc scans : %lld\n", numArcScans);
    printf("c Number of mergers   : %d\n", numMergers);
    printf("c Number of pushes    : %lld\n", numPushes);
    printf("c Number of relabels  : %d\n", numRelabels);
    printf("c Number of gaps      : %d\n", numGaps);
#endif

#ifdef BREAKPOINTS
    displayBreakpoints();
#endif

    // Commented because it segfaults :(
    // freeMemory();

    return 0;
}