package triangleTwoClub;

import robustTwoClub.graph.Edge;
import robustTwoClub.graph.RtcGraph;
import robustTwoClub.graph.Triangle;
import triangleTwoClub.graph.T2CGraph;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class ConditionChecker {

    private ITriangleTwoCliqueAlgorithm algorithm;

    public ConditionChecker(ITriangleTwoCliqueAlgorithm algorithm){
        this.algorithm = algorithm;
    }

    /**
     * Two vertices are compatible if the distance between them is at most 2. Uses conflict graph if possible.
     * @param v a fist vertex id
     * @param w a second vertex id
     * @return true if v and w are compatible
     */
    public boolean verticesAreCompatible(int v, int w){
        if(v == w)
            return true;

        boolean compatible;
        if(algorithm.getConflictGraph() != null)
            compatible = !algorithm.getConflictGraph().adjacent(v, w);
        else
            compatible =  verticesAreCompatible(algorithm.getGraph(), v, w);

        if(!algorithm.isEdgeConditionEnabled() || !compatible)
            return compatible;

        // Edge L>1 addition:
        // vw \in E => v, w incompatible if |N(w)\cap N(v)| < \ell
        // vw \notin E => there must be vertex u with |N(w)\cap N(u)|\ge \ell  and |N(u)\cap N(v)|\ge \ell
        int l = algorithm.getL();
        if(algorithm.getGraph().adjacent(v, w)){
            return countCommonNeighbors(v, w, l) == l;
        }
        for(int u : algorithm.getGraph().getNeighbors(v)){
            if(!algorithm.getGraph().getNeighbors(w).contains(u)) continue;
            if(countCommonNeighbors(v, u, l) == l && countCommonNeighbors(w, u, l) == l){
                return true;
            }
        }

        return false;

    }

    /**
     * Two vertices are compatible if the distance between them is at most 2.
     * @param graph the graph v and w are part of
     * @param v a fist vertex id
     * @param w a second vertex id
     * @return true if v and w are compatible
     */
    public static boolean verticesAreCompatible(RtcGraph graph, int v, int w){
        if(v == w) return true;
        if(graph == null) return false;

        if(graph.getNeighbors(v).contains(w)) //
            return true;

        HashSet<Integer> small;
        int toFind;
        if(graph.getNeighbors(v).size() < graph.getNeighbors(w).size()){
            small = graph.getNeighbors(v);
            toFind = w;
        }else{
            small = graph.getNeighbors(w);
            toFind = v;
        }
        return small.stream().anyMatch( u -> graph.getNeighbors(u).contains(toFind));
    }

    /**
     * Two vertices are compatible if the distance between them is at most 2.
     * @param neighborhood the graph v and w are part of
     * @param v a fist vertex id
     * @param w a second vertex id
     * @return true if v and w are compatible
     */
    public static boolean verticesAreCompatible(HashSet<Integer> neighborhood,
                                                HashMap<Integer, HashSet<Integer>> adjacency, int v, int w){
        if(v == w) return true;
        if(neighborhood == null || adjacency == null) return false;

        if(adjacency.get(v).contains(w))
            return true;

        HashSet<Integer> small;
        int toFind;
        if(adjacency.get(v).size() < adjacency.get(w).size()){
            small = adjacency.get(v);
            toFind = w;
        }else{
            small = adjacency.get(w);
            toFind = v;
        }
        return small.stream().anyMatch( u -> adjacency.get(u).contains(toFind));
    }

    /**
     * Counts the number of triangles.
     * @param graph the graph
     * @param triangleMap the triangle map of the graph, if null then the triangle map gets computed
     * @param v the vertex id
     * @return the number of vertex triangles v is part of in the graph
     */

    public static boolean enoughVertexTriangles(RtcGraph graph, HashMap<Integer, HashSet<Triangle>> triangleMap,
                                                int v, int l){
        // compute the triangle map if necessary
        if(triangleMap == null)
            triangleMap = Triangle.getTriangleHashMap(graph);

        // check if there are triangles of v at all
        if(!graph.contains(v) || !triangleMap.containsKey(v))
            return false;

        if(graph.degree(v)*(graph.degree(v)-1)/2<l)
            return false;


        int n = 0;
        for(Triangle t : triangleMap.get(v)){
            if(t.exists(graph)) // only count triangles that currently exist
                n++;
            if(n >= l)
                return true;
        }
        return false;
    }


    public static int countVertexTrianglesFromMap(RtcGraph graph, HashMap<Integer, HashSet<Triangle>> triangleMap,
                                                  int v, int max){
        // compute the triangle map if necessary
        if(triangleMap == null)
            triangleMap = Triangle.getTriangleHashMap(graph);

        // check if there are triangles of v at all
        if(!graph.contains(v) || !triangleMap.containsKey(v))
            return 0;

        int n = 0;
        for(Triangle t : triangleMap.get(v)){
            if(t.exists(graph)) // only count triangles that currently exist
                n++;
            if(n >= max)
                break;
        }
        return n;
    }

    /**
     * Counts the triangles of a vertex.
     * @param v the vertex id
     * @return the number of triangles v is part of
     */
    public int countTriangles(int v){
        if(!algorithm.getGraph().contains(v) || !algorithm.getTriangles().containsKey(v))
            return 0;

        int c = 0;
        for(Triangle t : algorithm.getTriangles().get(v))
            if(t.exists(algorithm.getGraph()))
                c++;
        return c;
    }

    /**
     * Counts the vertices that are compatible to a given one in the working graph. Uses the conflict graph if possible.
     * @param v the vertex to count the compatibles from
     * @param max the maximum number to count to
     * @return the number of vertices v is compatible to
     */
    public int countCompatibles(int v, int max, boolean countV){
        if(algorithm.getConflictGraph() != null){
            return algorithm.getConflictGraph().size() - algorithm.getConflictGraph().getNeighbors(v).size();
        }
        else
            return countCompatibles(algorithm.getGraph(), v, max, countV);
    }

    /**
     * Counts the vertices that are compatible to a given one.
     * @param graph the graph containing the vertices
     * @param v the vertex to count the compatibles from
     * @param max the maximum number to count to
     * @return the number of vertices v is compatible to
     */
    public static int countCompatibles(RtcGraph graph, int v, int max, boolean countV){
        if(graph == null || !graph.contains(v))
            return 0;


        int n = graph.size();
        if(max > (countV ? n : n - 1))
            max = (countV ? n : n - 1);

        int compatibles = countV ? 1 : 0;
        for(int w : graph.getVertices()){
            n--;
            if(v == w)
                continue;
            if(verticesAreCompatible(graph, v, w)){
                if(++compatibles >= max)
                    break;
            }
            if (compatibles+n < max){
                return compatibles;
            }
        }

        return compatibles;
    }

    /**
     * Counts the number of triangles an edge is part of in the working graph. Uses the precomputed triangles of the
     * algorithm if possible.
     * @param v the first vertex the edge is adjacent to
     * @param w the second vertex the edge is adjacent to
     * @param max the maximum number of triangles to count to
     * @return the number of triangles the edge (v, w) is part of
     */
    public int countEdgeTriangles(int v, int w, int max){
        if(algorithm.getTriangles() == null)
            return countEdgeTriangles(algorithm.getGraph(), v, w, max);

        if(algorithm.getTriangles().get(v) == null)
            return 0;

        int c = 0;
        for(Triangle t : algorithm.getTriangles().get(v)){
            if(t.contains(w) && t.exists(algorithm.getGraph())){
                if(++c >= max)
                    break;
            }
        }
        return c;
    }

    /**
     * Counts the number of triangles an edge is part of in a graph.
     * @param graph the graph
     * @param v the first vertex the edge is adjacent to
     * @param w the second vertex the edge is adjacent to
     * @param max the maximum number of triangles to count to
     * @return the number of triangles the edge (v, w) is part of
     */
    public static int countEdgeTriangles(RtcGraph graph, int v, int w, int max){
        if(graph == null || v == w || !graph.contains(v) || graph.contains(w) || max <= 0)
            return 0;

        int c = 0;
        for(int u : graph.getNeighbors(v)){
            if(u != w && graph.getNeighbors(w).contains(u)){
                if(++c >= max)
                    break;
            }
        }
        return c;
    }

    /**
     * Counts edge triangles of a given edge.
     * @param triangles the map of triangles
     * @param v the first vertex the edge is adjacent to
     * @param w the second vertex the edge is adjacent to
     * @param max the maximum number of triangles to count to
     * @return the number of triangles the edge (v, w) is part of
     */
    public static int countEdgeTriangles(HashMap<Integer, HashSet<Triangle>> triangles, int v, int w, int max){
        int c = 0;
        if(!triangles.containsKey(v))
            return 0;
        for(Triangle t : triangles.get(v)){
            if(t.contains(w))
                if(++c >= max)
                    break;
        }
        return c;
    }

    /**
     * Counts the number of triangles an edge is part of in a graph using the precomputed triangles.
     * @param graph the graph
     * @param triangles the map of precomputed triangles
     * @param v the first vertex the edge is adjacent to
     * @param w the second vertex the edge is adjacent to
     * @param max the maximum number of triangles to count to
     * @return the number of triangles the edge (v, w) is part of
     */
    public static int countEdgeTrianglesFromMap(RtcGraph graph, HashMap<Integer, HashSet<Triangle>> triangles, int v,
                                                int w, int max){
        if(!graph.contains(v) || triangles == null)
            return -1;

        if(!triangles.containsKey(v) || !triangles.containsKey(w)) // v or w was not part of any triangle at all
            return 0;

        int n = 0;
        for(Triangle t : triangles.get(v)){
            if(t.exists(graph) && t.contains(w)) // only count triangles that currently exist and have w
                n++;
            if(n >= max)
                break;
        }

        return n;
    }

    /**
     * Find vertices that are part of all triangles that also contain the specified vertex v.
     * Method for application of No Choice Rule 1.
     * @param v the specified vertex
     * @return a list of all vertices that are part of every triangle that v is part of or null if there are none
     */
    public List<Integer> getAllTriangleSharedVertices(int v){
        if(!algorithm.getTriangles().containsKey(v))
            return null;

        List<Integer> triangleVertices = null;
        for(Triangle t : algorithm.getTriangles().get(v)){
            // find fist triangle
            if(triangleVertices == null && t.exists(algorithm.getGraph())){
                triangleVertices = t.getVertices();
                continue;
            }
            // remove from list that are not in the first triangle
            if(triangleVertices != null && t.exists(algorithm.getGraph())){
                for(Integer w : t.getVertices()){
                    if(!triangleVertices.contains(w))
                        triangleVertices.remove(w);
                }
            }
        }

        return triangleVertices;
    }

    /**
     * Counts the common neighbors of two vertices up until a specified maximum.
     * @param v the first vertex
     * @param w the second vertex
     * @param max the maximum number to count to
     * @return the number of common neighbors between v and w
     */
    public int countCommonNeighbors(int v, int w, int max){
        RtcGraph g = algorithm.getGraph();
        if(!g.contains(v) || !g.contains(w) || max <= 0)
            return 0;
        if(v == w)
            return g.getNeighbors(v).size();

        int count = 0;
        for(int u : g.getNeighbors(v)){
            if(g.getNeighbors(w).contains(u)){
                if(++count >= max)
                    break;
            }
        }
        return count;
    }

    public boolean isLTriangleTwoClique(){
        RtcGraph graph = algorithm.getGraph();
        int l = graph.size();
        if(l < 3) // nothing to check for
            return true;

        if(!algorithm.containsOnlyTriangle(graph)){
            return false;
        }

        return(algorithm.isClique(graph));
    }

    /**
     * Checks if the current algorithm state is a valid solution to the vertex l-triangle-2-club problem.
     * @return true if the graph is a l-triangle-2-club
     */
    public boolean isLTriangleTwoClub(){
        RtcGraph graph = algorithm.getGraph();
        int l = algorithm.getL();

        if(l < 1) // nothing to check for
            return true;

        for(int v : graph.getVertices()){
            // is the distance between v and every other vertex at most two
            for(int w : graph.getVertices())
                if(v < w && !verticesAreCompatible(v, w))
                    return false;
        }

        for(int v : graph.getVertices()){
            // is v part of ell triangles?
            if(countVertexTrianglesFromMap(graph, algorithm.getTriangles(), v, l) < l)
                return false;
        }

        return true;
    }

    /**
     * Method to check if the given graph is a l - Triangle Two Club,
     * where each vertex is in at least l triangles
     * @return true when the graph is a Triangle Two Club, false else
     */
    public static boolean isLTriangleTwoClub(RtcGraph graph, HashMap<Integer, HashSet<Triangle>> triangles, int l){
        if(l < 1) // nothing to check for
            return true;

        for(int v : graph.getVertices()){
            // is v part of ell triangles?
            if(countVertexTrianglesFromMap(graph, triangles, v, l) < l)
                return false;

            // is the distance between v and every other vertex at most two
            for(int w : graph.getVertices())
                if(v < w && !verticesAreCompatible(graph, v, w))
                    return false;
        }
        return true;
    }

    /**
     * Method to check if the given graph is a l - Triangle Two Club,
     * where each vertex is in at least l triangles. If it is not a triangle 2 club, this method returns a pair of
     * incompatible vertices as an edge.
     * @param graph the graph
     * @return null if it is a valid triangle 2 club, an edge object containing a single vertex if this vertex has not
     *         enough triangles, or an edge containing two incompatible vertices
     */
    public static Edge getIncompatiblePair(RtcGraph graph){
        for(int v : graph.getVertices()){
            // is the distance between v and every other vertex at most two
            for(int w : graph.getVertices())
                if(v < w && !verticesAreCompatible(graph, v, w))
                    return new Edge(v, w);
        }
        return null;
    }

    /**
     * Checks if the current algorithm state is a valid solution to the edge l-triangle-2-club problem.
     * @return true if the graph is a edge l-triangle-2-club
     */
    public boolean isEdgeLTriangleTwoClub(){
        int l = algorithm.getL();
        RtcGraph graph = algorithm.getGraph();

        if(l < 1)
            return true;

        for(Integer v : graph.getVertices()) {
            // check if every edge is part of at least l triangles
            for (Integer w : graph.getNeighbors(v))
                if (v < w && countEdgeTrianglesFromMap(graph, algorithm.getTriangles(), v, w, l) < l)
                    return false;

            // check if every vertex is compatible to every other vertex
            for(Integer w : graph.getVertices())
                if(v < w && !verticesAreCompatible(v, w))
                    return false;
        }
        return true;
    }

    /**
     * Checks if the given graph with its triangles is an edge-l-triangle-2-club.
     * @param graph a graph
     * @param triangles the triangles of the graph
     * @param l the number of triangles each edge has to be part of
     * @return true if the graph is an edge-l-triangle-2-club
     */
    public static boolean isEdgeLTriangleTwoClub(RtcGraph graph, HashMap<Integer, HashSet<Triangle>> triangles, int l){
        if(l < 1)
            return true;

        for(Integer v : graph.getVertices()) {
            // check if every edge is part of at least l triangles
            for (Integer w : graph.getNeighbors(v))
                if (v < w && countEdgeTrianglesFromMap(graph, triangles, v, w, l) < l)
                    return false;

            // check if every vertex is compatible to every other vertex
            for(Integer w : graph.getVertices())
                if(v < w && !ConditionChecker.verticesAreCompatible(graph, v, w))
                    return false;
        }
        return true;
    }


    /**
     * Counts the number of vertices incompatible to a specified one.
     * @param v the vertex
     * @param currentWorst the maximum to count to
     * @return the number of incompatibilities v has
     */
    public int countIncompatibles(int v,int currentWorst){
        if(!algorithm.getGraph().contains(v))
            return 0;

        int n = algorithm.getGraph().size();
        if(currentWorst > n -1)
            currentWorst = n -1;

        int c = 0;
        for(int w : algorithm.getGraph().getVertices()){
            n--;
            if(!verticesAreCompatible(v, w))
                c++;
            if (c+n<= currentWorst){
                return 0;
            }
        }
        return c;
    }


    /**
     * Counts the number of vertices incompatible to a specified one.
     * @param v the vertex
     * @param currentWorst the maximum to count to
     * @return the number of incompatibilities v has
     */
    public int countIncompatiblesClique(int v,int currentWorst){
        if(!algorithm.getGraph().contains(v))
            return 0;
        int n = algorithm.getGraph().size();
        if(!vertexUsedInAnyActiveTriangle(v)){
            return(n);
        }
        if(currentWorst > n -1)
            currentWorst = n -1;

        int c = 0;
        for(int w : algorithm.getGraph().getVertices()){
            if(!algorithm.getOriginalGraph().getGraph().getTwoNeighbors(w).contains(v)){
                c++;
            }
            n--;
            if (c+n<= currentWorst){
                return 0;
            }
        }
        return c;
    }


    private boolean vertexUsedInAnyActiveTriangle(Integer v){
        Set<Set<Integer>> setsOfTriangle = algorithm.getSetsOfTriangle();
        Set<Integer> vertices = algorithm.getGraph().getVertices();
        for (Set<Integer> triangle : setsOfTriangle){
            if(triangle.contains(v) && vertices.containsAll(triangle)){
                return(true);
            }
        }
        return (false);
    }

    /**
     * Determines the number of incompatible vertices a given vertex has in the graph, using a BFS.
     * @param g a graph
     * @param v the vertex to count the incompatibles from
     * @return the number of vertices incompatible to v, -1 if the graph is null or does not contain v
     */
    public static int countIncompatibles(RtcGraph g, int v){
        if(g == null || !g.contains(v))
            return -1;

        HashSet<Integer> compatibles = new HashSet<>(); // add neighbors
        compatibles.add(v);
        g.getNeighbors(v).forEach(u -> compatibles.addAll(g.getNeighbors(u))); // add neighbors of neighbors
        compatibles.addAll(g.getNeighbors(v));

        return g.size() - compatibles.size();
    }



    /**
     * Counts the number of vertices incompatible to a specified one.
     * @param v the vertex
     * @return the number of incompatibilities v has
     */
    public int countIncompatibles(int v){
        if(!algorithm.getGraph().contains(v))
            return 0;

        int c = 0;
        for(int w : algorithm.getGraph().getVertices()){
            if(!verticesAreCompatible(v, w))
                c++;
        }
        return c;
    }

    /**
     * Counts the number of vertices incompatible to a specified one.
     * @param v the vertex
     * @param currentWorst the maximum to count to
     * @return the number of incompatibilities v has
     */
    public int countIncompatibles(int v,int currentWorst, RtcGraph inputGraph){
        if(!algorithm.getGraph().contains(v))
            return 0; //graph am arbeitsstand

        int n = algorithm.getGraph().size(); //graph jeweiliger arbeitsstand
        if(currentWorst > n -1)
            currentWorst = n -1;

        int c = 0;
        for(int w : algorithm.getGraph().getVertices()){ //anpassen

            n--;
            if(!verticesAreCompatible(inputGraph, v, w))
                c++;
            if (c+n<= currentWorst){
                return 0;
            }
        }
        return c;
    }

    /**
     * Finds the id of the vertex with the highest number of incompatibilities in the graph.
     * @return the vertex id with the most incompatibilities
     */
    public int getMostIncompatibleVertexID(){
        int worstID = -1;
        int worstIncompatibilityCount = 0;

        int c;
        for(int v : algorithm.getGraph().getVertices()){
            c = countIncompatibles(v,worstIncompatibilityCount);
            if(c > worstIncompatibilityCount){
                worstID = v;
                worstIncompatibilityCount = c;
            }
        }
        return worstID;
    }

    /**
     * Finds the id of the vertex with the highest number of incompatibilities in the graph.
     * @return the vertex id with the most incompatibilities
     */
    public int getMostIncompatibleVertexID(RtcGraph inputGraph){
        int worstID = -1;
        int worstIncompatibilityCount = 0;

        int c;
        for(int v : algorithm.getGraph().getVertices()){ // bearbeitet
            c = countIncompatibles(v,worstIncompatibilityCount, inputGraph);
            if(c > worstIncompatibilityCount){
                worstID = v;
                worstIncompatibilityCount = c;
            }
        }
        return worstID;
    }


    /**
     * Finds the id of the vertex with the highest number of incompatibilities in the graph.
     * @return the vertex id with the most incompatibilities
     */
    public int getMostIncompatibleVertexIDClique(){
        int worstID = -1;
        int worstIncompatibilityCount = 0;

        int c;
        for(int v : algorithm.getGraph().getVertices()){
            c = countIncompatiblesClique(v,worstIncompatibilityCount);
            if(c > worstIncompatibilityCount){
                worstID = v;
                worstIncompatibilityCount = c;
            }
        }
        return worstID;
    }

    /**
     * Finds all vertices that a given vertex shares triangles with.
     * @param v the vertex
     * @return a set of triangle vertices, null if there are none
     */
    public HashSet<Integer> getTriangleVertices(int v){
        if(!algorithm.getGraph().contains(v) || !algorithm.getTriangles().containsKey(v))
            return null;

        HashSet<Integer> vertices = new HashSet<>();
        for(Triangle t : algorithm.getTriangles().get(v))
            vertices.addAll(t.getVertices());
        return vertices;
    }


    /**
     * Determines if two vertices are compatible to each other with the new definition of compatible.
     *
     * Let T_v/T_u be the set of vertices that vertex v/u shares a triangle with.
     * While a vertex w in T_u has a distance of at most two to less than L vertices of T_v, remove w. (and vice versa)
     *
     * If |T_v| < L, |T_u| < L, or w was marked, then return no/false.
     *
     * @param v a vertex
     * @param u another vertex
     * @return true if v and u are compatible
     */
    public boolean verticesAreCompatibleV2(int v, int u){
        if(v == u)
            return true;

        HashSet<Integer> T_v = getTriangleVertices(v), T_u = getTriangleVertices(u);
        Integer w;

        // find w in T_u
        while(true){
            w = findW(T_u, T_v);
            if(w == null) break;
            T_u.remove(w);
            if(T_u.size() < algorithm.getL() || algorithm.getMarked().contains(w))
                return false;
        }

        // find w in T_v
        while(true){
            w = findW(T_v, T_u);
            if(w == null) break;
            T_v.remove(w);
            if(T_v.size() < algorithm.getL() || algorithm.getMarked().contains(w))
                return false;
        }

        return true;
    }

    /**
     * Finds a vertex w according to new compatibility definition (see method "verticesAreCompatibleV2(int v, int u)").
     * @param set1 the set w should be part of
     * @param set2 the other set
     * @return a vertex w, null if no vertex fits the definition
     */
    private Integer findW(HashSet<Integer> set1, HashSet<Integer> set2){
        if(set1 == null || set2 == null)
            return null;

        Integer vertexToRemove = null;
        int count;

        wLoop:
        for(Integer w : set1){
            count = 0;
            for(Integer u : set2){
                if(verticesAreCompatible(algorithm.getGraph(), w, u)) {
                    // if there are enough compatible vertices for w, then try another vertex
                    if (++count >= algorithm.getL())
                        continue wLoop;
                }
            }
            // if the inner loop did not find enough compatible vertices, then we found vertexToRemove
            vertexToRemove = w;
            break;
        }

        return vertexToRemove;
    }
}
