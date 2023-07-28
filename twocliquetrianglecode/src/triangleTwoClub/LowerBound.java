package triangleTwoClub;

import robustTwoClub.graph.Edge;
import robustTwoClub.graph.RollbackRtcGraph;
import robustTwoClub.graph.RtcGraph;
import robustTwoClub.graph.Triangle;
import triangleTwoClub.dataReduction.TwoNeighborhoodRule;
import triangleTwoClub.graph.ConflictGraph;
import triangleTwoClub.graph.T2CGraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

public class LowerBound {

    public enum LowerBoundVariant{BASIC, LB1, LB2, DISABLED}
    public enum LB2GreedyApproach{FIRST, WORST}
    static String lowerBoundVariantName = "LowerBoundVariant";
    static String lowerBoundSizeName = "LowerBoundSize";
    static String lowerBoundTimeName = "LowerBoundTime[s]";
    static String lowerBoundTimeLimitName = "LowerBoundTimeLimit[s]";
    static String lowerBoundTimeLimitAbortName = "LowerBoundAborted";
    static String LB1TimeForLB2 = "LB1Time[s]";
    static String LB1SizeForLB2 = "LB1SolutionSize";


    private ITriangleTwoClubAlgorithm algorithm;
    private Output output;
    private LowerBoundVariant variant;
    private LB2GreedyApproach lb2approach;
    private boolean aborted;
    private long timeLimit;
    private int vertexMinDegree;

    private RtcGraph conflict;

    public LowerBound(ITriangleTwoClubAlgorithm algorithm, Output output, LowerBoundVariant variant){
        this.algorithm = algorithm;
        this.output = output;
        this.variant = variant;
        this.timeLimit = 0;
        this.vertexMinDegree = 1;
        lb2approach = LB2GreedyApproach.FIRST;
    }

    public LowerBound(ITriangleTwoClubAlgorithm algorithm, LowerBoundVariant variant){
        this(algorithm, algorithm.getOutput(), variant);
    }

    /**
     * Computes the lower bound solution using the parameters of the constructor.
     * @return A lower bound solution as graph
     */
    public RtcGraph getLowerBoundSolution(){
        // add information to output
        if(output != null){
            if(variant == LowerBoundVariant.BASIC)
                output.addInformation(lowerBoundVariantName, "Basic");
            if(variant == LowerBoundVariant.LB1)
                output.addInformation(lowerBoundVariantName, "LB1");
            if(variant == LowerBoundVariant.LB2){
                output.addInformation(lowerBoundVariantName, "LB2");
                output.setFormat(LB1TimeForLB2, Output.DataFormat.MILLISECONDS);
                output.setFormat(LB1SizeForLB2, Output.DataFormat.INTEGER);
            }
            if(variant == LowerBoundVariant.DISABLED)
                output.addInformation(lowerBoundVariantName, "DISABLED");
            output.setFormat(lowerBoundTimeName, Output.DataFormat.MILLISECONDS);
            output.setFormat(lowerBoundSizeName, Output.DataFormat.INTEGER);
        }

        // compute minimal degree for a vertex
        if(algorithm.isEdgeConditionEnabled()) {
            this.vertexMinDegree = algorithm.getL() + 1;
        }else{
            double d = 0.5d + Math.sqrt(0.25d + 2 * algorithm.getL());
            this.vertexMinDegree = (int) Math.ceil(d);
        }

        // compute lowerbound
        RtcGraph lowerBoundSolution = null;
        aborted = false;
        long time = System.currentTimeMillis();
        if(variant == LowerBoundVariant.BASIC)
            lowerBoundSolution = getBasicLowerBoundSolution(algorithm.getGraph());
        if(variant == LowerBoundVariant.LB1)
            lowerBoundSolution = getLB1Solution(time);
        if(variant == LowerBoundVariant.LB2)
            lowerBoundSolution = getLB2Solution(time);
        time = System.currentTimeMillis() - time;

        // log
        if(output != null){
            output.addInformation(lowerBoundTimeLimitName, String.valueOf((int)(timeLimit / 1000)));
            output.addInformation(lowerBoundTimeLimitAbortName, String.valueOf(aborted));
            output.log(lowerBoundTimeName, time);
            output.log(lowerBoundSizeName, lowerBoundSolution == null ? 0 : lowerBoundSolution.size());
        }

        return lowerBoundSolution;
    }

    /**
     * Sets a time limit for the computation of the lower bound. Use a value < 1 to disable time limit.
     * Default is 0.
     * @param limit the time limit in ms
     */
    public void setTimeLimit(long limit){
        if(limit < 0)
            timeLimit = -1L;
        else
            timeLimit = limit;
    }

    /**
     * Sets the approach the LB2 uses to make a two club from a found solution.
     * @param a the approach to take
     */
    public void setLB2GreedyApproach(LB2GreedyApproach a){
        this.lb2approach = a;
    }

    /**
     * Calculates lower bound for a triangle 2-club.
     * Counts non-isolated vertices in the neighbourhood of every vertex to find a lower bound.
     * @param graph the graph
     * @return the lower bound solution
     */
    public static RtcGraph getBasicLowerBoundSolution(RtcGraph graph){
        // find best vertex
        int bestVertex = -1, best = 0;
        for(Integer v : graph.getVertices()){
            HashSet<Integer> neighbours = graph.getNeighbors(v);
            int isolatedCount = neighbours.size();
            wLoop:
            for(Integer w : neighbours){
                for(Integer x : neighbours){
                    if(!w.equals(x) && graph.adjacent(w, x)){
                        isolatedCount--;
                        continue wLoop;
                    }
                }
            }
            int size = neighbours.size() + 1 - isolatedCount; // neighbours + v itself - isolated neighbours
            if(size > 2 && size > best){
                best = size;
                bestVertex = v;
            }
        }

        if(bestVertex < 0) return null;

        // get lower bound solution from found vertex
        RtcGraph lowerBound = graph.getClone();
        HashSet<Integer> neighbours = graph.getNeighbors(bestVertex);
        HashSet<Integer> toDelete = new HashSet<>();
        for(int v : lowerBound.getVertices()){      // reduce to best vertex and it's neighbourhood
            if( v != bestVertex && !neighbours.contains(v)) toDelete.add(v);
        }
        for(Integer v : neighbours){
            boolean isIsolated = true;
            for(Integer w : neighbours) {
                if (!v.equals(w) && lowerBound.adjacent(v, w)) {
                    isIsolated = false;
                    break;
                }
            }
            if(isIsolated) toDelete.add(v);
        }
        toDelete.forEach(lowerBound::deleteVertex);
        return lowerBound;
    }

    private RtcGraph getLB1Solution(long startTime){
        RtcGraph bestSolution = null, g = algorithm.getGraph();
        HashSet<Integer> neighborhood;
        int bestSolutionSize;

        // find the best vertex
        for(Integer v : g.getVertices()){
            // check for abort due to time limit
            if(timeLimit > -1 && startTime + timeLimit <= System.currentTimeMillis()){
                aborted = true;
                break;
            }

            bestSolutionSize = bestSolution == null ? 0 : bestSolution.size();

            // create N[v]
            neighborhood = new HashSet<>();
            neighborhood.add(v);
            neighborhood.addAll(g.getNeighbors(v));
            if(neighborhood.size() <= bestSolutionSize)
                continue;

            RtcGraph solution = algorithm.getGraph().getSubgraph(neighborhood);
            if(!vertexDegreeReduction(solution, bestSolutionSize, null, v))
                continue;

            boolean res;
            // reduces the neighborhood, if it returns false, then abort v
            if(algorithm.isEdgeConditionEnabled()){
                res = doLBRemovalEdge(solution, v, bestSolutionSize);
            }else {
                // get subset of triangles in N[v]
                HashMap<Integer, HashSet<Triangle>> localTriangles = getNeighborhoodTriangles(algorithm.getTriangles(),
                        neighborhood);
                res = doLBRemovalVertex(solution, localTriangles, v, bestSolutionSize);
            }

            if(!res) // res is false if the removal did not generate a better or valid solution
                continue;
            bestSolution = solution;
        }

        // create solution for best neighborhood and return
        return bestSolution;
    }

    private boolean doLBRemovalVertex(RtcGraph graph, HashMap<Integer, HashSet<Triangle>> localTriangles,
                                      Integer centerVertex, int bestSize){
        // vertex LB2
        // find a vertex w that is part of at most l-1 triangles and remove it
        HashSet<Integer> affectedVertices;

        while(true){
            affectedVertices = new HashSet<>();
            // search for invalid vertices
            for(Integer u : graph.getVertices())
                if(!localTriangles.containsKey(u) || localTriangles.get(u).size() < algorithm.getL()){
                    affectedVertices.add(u);
                }

            // exit loop if there is no w
            if(affectedVertices.size() == 0)
                break;

            if(affectedVertices.contains(centerVertex) || graph.size() - affectedVertices.size() <= bestSize)
                return false;

            for(Integer u : affectedVertices){
                graph.deleteVertex(u);
                removeTrianglesOfVertex(localTriangles, u);
            }
            vertexDegreeReduction(graph, bestSize, localTriangles, centerVertex);

            if(graph.size() <= bestSize)
                return false;
        }
        return true;
    }

    private boolean doLBRemovalEdge(RtcGraph graph, Integer centerVertex, int bestSize){
        // edge LB2
        LinkedList<Edge> lowTriangleEdges;
        HashSet<Integer> affectedVertices;
        boolean removed;

        while(true){
            lowTriangleEdges = new LinkedList<>();
            affectedVertices = new HashSet<>();
            removed = false;

            // search for w
            for(Integer u : graph.getVertices())
                for(Integer w : graph.getNeighbors(u)){
                    if(u >= w)
                        continue;

                    if(graph.getCommonNeighbors(u, w).size() < algorithm.getL()){
                        lowTriangleEdges.add(new Edge(u, w));
                        affectedVertices.add(u);
                        affectedVertices.add(w);
                    }
                }


            // exit loop if there is all edges are part of enough triangles
            if(lowTriangleEdges.size() == 0)
                break;

            // remove bad edges from graph
            lowTriangleEdges.forEach( e -> {
                graph.getNeighbors(e.getV()).remove(e.getW());
                graph.getNeighbors(e.getW()).remove(e.getV());
            });

            // remove vertices with low degree
            for(Integer u : affectedVertices){
                if(graph.getNeighbors(u).size() > algorithm.getL())
                    continue;
                if(u.equals(centerVertex))
                    return false;
                graph.deleteVertex(u);
                removed = true;
            }

            if(removed)
                vertexDegreeReduction(graph, bestSize, null, centerVertex);

            if(graph.size() <= bestSize)
                return false;
        }
        return true;
    }


    private RtcGraph getLB2Solution(long startTime){
        RtcGraph g = algorithm.getGraph(), bestSolution;
        HashSet<Integer> neighborhood, toDelete = new HashSet<>();

        // get LB1 solution
        long lb1Time = System.currentTimeMillis();
        bestSolution = getLB1Solution(startTime);
        lb1Time = System.currentTimeMillis() - lb1Time;
        output.log(LB1TimeForLB2, lb1Time);
        int bestSolutionSize = bestSolution == null ? 0 : bestSolution.size();
        output.log(LB1SizeForLB2, bestSolutionSize);
        new TwoNeighborhoodRule(algorithm, bestSolutionSize).apply();

        // find the best vertex
        ArrayList<Integer> startPoints = getVerticesByAscendingDegree(g, true);
        for(int i = 0; i < startPoints.size(); i++){
            int v = startPoints.get(i);
            // check for abort due to time limit
            if(timeLimit > -1 && startTime + timeLimit <= System.currentTimeMillis()){
                aborted = true;
                break;
            }
            bestSolutionSize = bestSolution == null ? 0 : bestSolution.size();

            // create N_2[v], a copy is needed to not destroy the graphs adjacency
            neighborhood = g.getTwoNeighbors(v);
            neighborhood.add(v);
            neighborhood.removeAll(toDelete);
            if(neighborhood.size() <= bestSolutionSize){
                toDelete.add(v);
                continue;
            }

            conflict = new RtcGraph();

            boolean res;
            RtcGraph solution = g.getSubgraph(neighborhood);
            if(!vertexDegreeReduction(solution, bestSolutionSize, null, v)){
                toDelete.add(v);
                continue;
            }

            // build conflict graph
            for (int w : solution.getVertices()){
                conflict.addVertex(w);
            }
            for(int x : solution.getVertices()) {
                for (int w : solution.getVertices()) {
                    if (w <= x) continue;
                    if (!verticesAreCompatible(solution, x, w))
                        conflict.addEdge(x, w);
                }
            }

            // make graph satisfy the triangle condition
            HashMap<Integer, HashSet<Triangle>> localTriangles = null;
            if(algorithm.isEdgeConditionEnabled()){
                res = doLB2RemovalEdge(solution, v, bestSolutionSize);
            }else{
                // get subset of triangles in N_2[v]
                localTriangles = getNeighborhoodTriangles(algorithm.getTriangles(), neighborhood);
                res = doLB2RemovalVertex(solution, localTriangles, v, bestSolutionSize);
            }

            // continue if neighborhood is too small
            if(!res || solution.size() <= bestSolutionSize){
                toDelete.add(v);
                continue;
            }

            if(lb2approach == LB2GreedyApproach.FIRST) {
                if (!makeTwoClubGreedyFirst(solution, localTriangles, bestSolutionSize, algorithm.getL(), v))
                    continue;
            }else {
                makeTwoClubGreedyWorst(solution, localTriangles, bestSolutionSize, algorithm.getL(), v);
            }
            // check if remaining solution is better than current best
            if(bestSolutionSize < solution.size()){
                bestSolution = solution;
                //System.out.println(bestSolution.size());
            }
        }

        toDelete.forEach(v -> algorithm.deleteVertex(v));
        // create solution for best neighborhood
        return bestSolution;
    }

    private boolean doLB2RemovalVertex(RtcGraph graph, HashMap<Integer, HashSet<Triangle>> localTriangles,
                                       Integer centerVertex, int bestSize){
        // vertex LB2
        // find a vertex w that is part of at most l-1 triangles and remove it
        HashSet<Integer> affectedVertices;


        while(true){
            // low compatibility rule
            HashSet<Integer> toDelete2 = new HashSet<>();

            for(int v2 : graph.getVertices()){
                if(conflict.size() - conflict.getNeighbors(v2).size() <= bestSize){
                    toDelete2.add(v2);
                }
            }
            if (toDelete2.contains(centerVertex)){
                return false;
            }

            for (int w : toDelete2){
                HashSet<Integer> neigs = graph.getNeighbors(w);
                graph.deleteVertex(w);
                updateConflictsForVertexDeletionV1(w, graph, neigs);
                removeTrianglesOfVertex(localTriangles, w);
            }

            if(!vertexDegreeReduction(graph, bestSize, localTriangles, centerVertex)){
                return false;
            }

            if (graph.size() <= bestSize){
                return false;
            }

            if (toDelete2.size() >0){
                continue;
            }

            affectedVertices = new HashSet<>();
            // search for invalid vertices
            for(Integer u : graph.getVertices())
                if(!localTriangles.containsKey(u) || localTriangles.get(u).size() < algorithm.getL()){
                    affectedVertices.add(u);
                }

            // exit loop if there is no w
            if(affectedVertices.size() == 0)
                break;

            for(Integer u : affectedVertices){
                if(u.equals(centerVertex) || graph.size() - 1 <= bestSize )
                    return false;
                HashSet<Integer> neigs = graph.getNeighbors(u);
                graph.deleteVertex(u);
                updateConflictsForVertexDeletionV1(u, graph, neigs);
                removeTrianglesOfVertex(localTriangles, u);
            }
            vertexDegreeReduction(graph, bestSize, localTriangles, centerVertex);

            if(graph.size() <= bestSize)
                return false;
        }
        return true;
    }

    private boolean doLB2RemovalEdge(RtcGraph graph, Integer centerVertex, int bestSize){
        // edge LB2
        LinkedList<Edge> lowTriangleEdges;
        HashSet<Integer> affectedVertices;
        boolean removed;

        while(true){
            // low compatibility rule
            HashSet<Integer> toDelete2 = new HashSet<>();

            for(int v2 : graph.getVertices()){
                if(conflict.size() - conflict.getNeighbors(v2).size() <= bestSize){
                    toDelete2.add(v2);
                }
            }
            if (toDelete2.contains(centerVertex)){
                return false;
            }

            for (int w : toDelete2){
                HashSet<Integer> neigs = graph.getNeighbors(w);
                graph.deleteVertex(w);
                updateConflictsForVertexDeletionV1(w, graph, neigs);
            }

            if(!vertexDegreeReduction(graph, bestSize, null, centerVertex)){
                return false;
            }

            if (graph.size() <= bestSize){
                return false;
            }

            if (toDelete2.size() >0){
                continue;
            }

            lowTriangleEdges = new LinkedList<>();
            affectedVertices = new HashSet<>();
            removed = false;

            // search for w
            for(Integer u : graph.getVertices())
                for(Integer w : graph.getNeighbors(u)){
                    if(u >= w)
                        continue;

                    if(graph.getCommonNeighbors(u, w).size() < algorithm.getL()){
                        lowTriangleEdges.add(new Edge(u, w));
                        affectedVertices.add(u);
                        affectedVertices.add(w);
                    }
                }


            // exit loop if there is all edges are part of enough triangles
            if(lowTriangleEdges.size() == 0)
                break;

            for (Edge e : lowTriangleEdges){
                graph.getNeighbors(e.getV()).remove(e.getW());
                graph.getNeighbors(e.getW()).remove(e.getV());
                updateConflictsForEdgeDeletionV1(e.getV(), e.getW(), graph);
            }

            // remove vertices with low degree
            for(Integer u : affectedVertices){
                if(graph.getNeighbors(u).size() <= algorithm.getL()){
                    if(u.equals(centerVertex))
                        return false;
                    graph.deleteVertex(u);
                    removed = true;
                }
            }

            if(removed)
                vertexDegreeReduction(graph, bestSize, null, centerVertex);

            if(graph.size() <= bestSize)
                return false;
        }
        return true;
    }

    /**
     * Extracts triangles from a given map of triangles that are within a specified neighborhood.
     * @param triangles the complete set of triangles
     * @param neighborhood some vertices
     * @return only the triangles that are in the neighborhood
     */
    private static HashMap<Integer, HashSet<Triangle>> getNeighborhoodTriangles(
            HashMap<Integer, HashSet<Triangle>> triangles, HashSet<Integer> neighborhood){

        HashMap<Integer, HashSet<Triangle>> localTriangles = new HashMap<>();
        for(Integer v : neighborhood){
            // skip v if v has no triangles at all
            if(!triangles.containsKey(v))
                continue;

            // create a new set for v with only triangles from the neighborhood
            localTriangles.put(v, new HashSet<>());
            outerLoop:
            for(Triangle t : triangles.get(v)){
                // skip triangle if a vertex is not in the neighborhood
                for(Integer w : t.getVertices())
                    if(!neighborhood.contains(w))
                        continue outerLoop;

                // all vertices are in the neighborhood, add triangle
                localTriangles.get(v).add(t);
            }
        }

        return localTriangles;
    }

    /**
     * Removes the triangles of a specified vertex from a set of triangles.
     * @param triangles the set of triangles
     * @param v the vertex whose triangles needs to be removed
     */
    private static void removeTrianglesOfVertex(HashMap<Integer, HashSet<Triangle>> triangles, Integer v){
        HashSet<Triangle> vTriangles = triangles.remove(v);

        if(vTriangles == null)
            return;

        // remove all triangles v is part of without deleting the HashSets (for counting purpose)
        for(Triangle t : vTriangles)
            for(Integer w : t.getVertices())
                if(triangles.get(w) != null)
                    triangles.get(w).remove(t);
    }

    /**
     * Removes the triangles that edge e is part of
     * @param triangles the map of triangles
     * @param e an edge
     */
    private static void removeTrianglesOfEdge(HashMap<Integer, HashSet<Triangle>> triangles, Edge e){
        Integer v = e.getV(), w = e.getW();
        if(!triangles.containsKey(v))
            return;

        outer:
        while(true){
            for(Triangle t : triangles.get(v)){
                if(t.contains(w)){
                    for(Integer u : t.getVertices())
                        triangles.get(u).remove(t);
                    continue outer;
                }
            }
            break;
        }
    }

    public HashSet<Edge> matchingBoundVertexCover(RtcGraph working, int limit){
        // 2-approximation for vertex cover of conflict graph
        // Minimal vertex cover of conflict graph is lower bound on vertices to delete
        HashSet<Integer> marks = new HashSet<>();
        HashSet<Edge> edges = new HashSet<>();
        int size = 0;
        for (int v : conflict.getVertices())
            if (working.contains(v) && !marks.contains(v))
                for (int w : conflict.getNeighbors(v))
                    if (working.contains(w) && w > v)					// every remaining conflict is an edge
                        if (!marks.contains(w)) {	// if the edge is not covered
                            marks.add(v);
                            marks.add(w);
                            size++;
                            edges.add(new Edge(v, w)); // coverF edge with both end points
                            if (size > limit) return edges;				// if current size exceeds allowed limit
                            break;										// v us marked, so the condition of the if-statement is never met again.
                        }												// of deletions we do not need to carry on
        return edges;
    }

    // vertex v gets deleted from conflict graph
    private void updateConflictsForVertexDeletionV1(int v, RtcGraph solution, HashSet<Integer> toVisit){
        conflict.deleteVertex(v);

        for(Integer x : toVisit){
            for(Integer y : toVisit){
                if(y <= x) continue;
                if (!verticesAreCompatible(solution, x, y)){
                    conflict.addEdge(x, y);
                }
            }
        }
    }

    private void updateConflictsForEdgeDeletionV1(int v, int w, RtcGraph solution){

        // endpoint v with neighbors of w
        for(int u : solution.getNeighbors(w)){
            if(v != u && !ConditionChecker.verticesAreCompatible(solution, v, u))
                conflict.addEdge(v, u);
        }

        // endpoint w with neighbors of v
        for(int u : solution.getNeighbors(v)){
            if(w != u && !ConditionChecker.verticesAreCompatible(solution, w, u))
                conflict.addEdge(w, u);
        }

    }

    /**
     * Method used in LB2. The neighborhood of the LB2 algorithm may not be a 2-club. This method uses a greedy
     * approach to remove vertices from the solution until it is a valid 2-club.
     * This method deletes one of the incompatible vertices found first.
     * @param solution the LB2 solution
     * @param triangles the local triangles of the solution
     * @param bestSolutionSize the size of the current best solution (to abort)
     * @param l the algorithms value for L
     * @param v the center vertex
     */
    private boolean makeTwoClubGreedyFirst(RtcGraph solution, HashMap<Integer, HashSet<Triangle>> triangles,
                                           int bestSolutionSize, int l, int v){
        int uTriangles, wTriangles, toRemove;
        while (true){
            HashSet<Edge> matching = matchingBoundVertexCover(solution, conflict.size() - bestSolutionSize);
            if (matching.size() == 0){
                break;
            }
            if (solution.size() <= bestSolutionSize + matching.size()){
                return false;
            }
            for (Edge incompatible : matching){
                int u = incompatible.getV(), w = incompatible.getW();
                if(algorithm.isEdgeConditionEnabled()){
                    uTriangles = solution.degree(u);
                    wTriangles = solution.degree(w);
                }else{
                    uTriangles = triangles.containsKey(u) ? triangles.get(u).size() : 0;
                    wTriangles = triangles.containsKey(u) ? triangles.get(w).size() : 0;
                }
                toRemove = uTriangles < wTriangles ? u : w;
                if(toRemove == v) // switch if toRemove is center vertex
                    toRemove = uTriangles < wTriangles ? w : u;

                // remove vertex and restore triangle condition
                HashSet<Integer> neigs = solution.getNeighbors(toRemove);
                solution.deleteVertex(toRemove);
                updateConflictsForVertexDeletionV1(toRemove, solution, neigs);
                if(!algorithm.isEdgeConditionEnabled())
                    removeTrianglesOfVertex(triangles, toRemove);

            }
            if(!restoreTriangleCondition(solution, triangles, v, bestSolutionSize, conflict))
                return false;

        }
        //Edge incompatible;
        //while((incompatible = ConditionChecker.getLTriangleTwoClubIncompatible(solution, l)) != null
        //        && solution.size() > bestSolutionSize){
        //    if(solution.size() == bestSolutionSize + 1)
        //        return false;

        // select a vertex to remove
        //   int u = incompatible.getV(), w = incompatible.getW();
        //   if(algorithm.isEdgeConditionEnabled()){
        //       uTriangles = solution.degree(u);
        //       wTriangles = solution.degree(w);
        //   }else{
        //       uTriangles = triangles.containsKey(u) ? triangles.get(u).size() : 0;
        //       wTriangles = triangles.containsKey(u) ? triangles.get(w).size() : 0;
        //    }
        //    toRemove = uTriangles < wTriangles ? u : w;
        //    if(toRemove == v) // switch if toRemove is center vertex
        //        toRemove = uTriangles < wTriangles ? w : u;

        // remove vertex and restore triangle condition
        //    HashSet<Integer> neigs = solution.getNeighbors(toRemove);
        //    solution.deleteVertex(toRemove);
        //    updateConflictsForVertexDeletionV1(toRemove, solution, neigs);
        ///   if(!algorithm.isEdgeConditionEnabled())
        //        removeTrianglesOfVertex(triangles, toRemove);
        //    if(!restoreTriangleCondition(solution, triangles, v, bestSolutionSize, conflict))
        //        return false;
        // }
        return true;
    }

    /**
     * Method used in LB2. The neighborhood of the LB2 algorithm may not be a 2-club. This method uses a greedy
     * approach to remove vertices from the solution until it is a valid 2-club.
     * This method deletes the vertex with an incompatibility and the lowest number of triangles.
     * @param solution the LB2 solution
     * @param triangles the local triangles of the solution
     * @param bestSolutionSize the size of the current best solution (to abort)
     * @param l the algorithms value for L
     * @param v the center vertex
     */
    private void makeTwoClubGreedyWorst(RtcGraph solution, HashMap<Integer, HashSet<Triangle>> triangles,
                                        int bestSolutionSize, int l,  int v){
        while(((!algorithm.isEdgeConditionEnabled() && !ConditionChecker.isLTriangleTwoClub(solution, triangles, l))
                || (algorithm.isEdgeConditionEnabled()
                && !ConditionChecker.isEdgeLTriangleTwoClub(solution, triangles, l)))
                && solution.size() > bestSolutionSize){
            int worstID = -1, worstTriangleCount = -1;

            // find two incompatible vertices
            for(int i : solution.getVertices()){
                //check if i has still enough triangles
                if(ConditionChecker.countVertexTrianglesFromMap(solution, triangles, i, l) < l){
                    solution.deleteVertex(i);
                    break;
                }

                // find an incompatible vertex
                for(int j : solution.getVertices()){
                    if(j == i) continue;
                    if(!ConditionChecker.verticesAreCompatible(solution, i, j)){
                        int t_i = ConditionChecker.countVertexTrianglesFromMap(solution, triangles, i,
                                Integer.MAX_VALUE);
                        int t_j = ConditionChecker.countVertexTrianglesFromMap(solution, triangles, j,
                                Integer.MAX_VALUE);

                        if(t_i < worstTriangleCount){
                            worstID = i;
                            worstTriangleCount = t_i;
                        }
                        if(t_j < worstTriangleCount){
                            worstID = j;
                            worstTriangleCount = t_j;
                        }
                    }
                }
            }

            // delete worst
            solution.deleteVertex(worstID);
        }
    }

    private boolean restoreTriangleCondition(RtcGraph solution, HashMap<Integer, HashSet<Triangle>> triangles,
                                             int v, int currentBest, RtcGraph conflict){
        if(!vertexDegreeReduction(solution, currentBest, triangles, v))
            return false; // if false is returned, the graph is already too small

        if(algorithm.isEdgeConditionEnabled()){
            // edge
            return doLB2RemovalEdge(solution, v, currentBest);
        }else{
            // vertex
            return doLB2RemovalVertex(solution, triangles, v, currentBest);
        }
    }

    /**
     * Removes vertices with a to small degree. Aborts if the graph gets too small for a better solution.
     * This method also
     * @param graph a graph
     * @param currentBest the size of the current best solution
     * @param triangles the local triangle map, or null if no triangles were computed
     * @param centerVertex the center of the local neighborhood
     * @return true if the remaining graph is larger than the current solution, false else
     */
    private boolean vertexDegreeReduction(RtcGraph graph, int currentBest,
                                          HashMap<Integer, HashSet<Triangle>> triangles, int centerVertex){
        HashSet<Integer> neighbors, toVisit = new HashSet<>(graph.getVertices());
        boolean removeTriangles;
        while(!toVisit.isEmpty() && graph.size() > currentBest){
            removeTriangles = false;
            for(Integer v : toVisit){
                neighbors = graph.getNeighbors(v);
                if(neighbors.size() < this.vertexMinDegree
                        || !ConditionChecker.verticesAreCompatible(graph, v, centerVertex)
                        || (triangles != null
                        && (!triangles.containsKey(v) || triangles.get(v).size() < algorithm.getL()))){
                    // remove vertex v because of its degree
                    toVisit.addAll(neighbors);
                    HashSet<Integer> neigs = graph.getNeighbors(v);
                    graph.deleteVertex(v);
                    if (conflict != null && conflict.size() > 0)
                        //System.out.println("new test " + conflict.size());
                        updateConflictsForVertexDeletionV1(v, graph, neigs);
                    removeTriangles = triangles != null;
                }
                if(removeTriangles)
                    removeTrianglesOfVertex(triangles, v);

                toVisit.remove(v);
                break;
            }
        }
        return graph.size() > currentBest;
    }

    /**
     * Removes vertices with a to small degree. Aborts if the graph gets too small for a better solution.
     * This method also
     * @param neighborhood a graph neighborhood
     * @param adjacency the local adjacency in the neighborhood
     * @param currentBest the size of the current best solution
     * @param triangles the local triangle map, or null if no triangles were computed
     * @param centerVertex the center of the local neighborhood
     * @return true if the remaining graph is larger than the current solution, false else
     */
    private boolean vertexDegreeReduction(HashSet<Integer> neighborhood, HashMap<Integer, HashSet<Integer>> adjacency,
                                          int currentBest, HashMap<Integer, HashSet<Triangle>> triangles,
                                          int centerVertex){
        HashSet<Integer> neighbors, toVisit = new HashSet<>(neighborhood);
        boolean removeTriangles;
        while(!toVisit.isEmpty() && neighborhood.size() > currentBest){
            removeTriangles = false;
            for(Integer v : toVisit){
                neighbors = adjacency.get(v);
                if(neighbors.size() < this.vertexMinDegree
                        || !ConditionChecker.verticesAreCompatible(neighborhood, adjacency, centerVertex, v)
                        || (triangles != null
                        && (!triangles.containsKey(v) || triangles.get(v).size() < algorithm.getL()))){
                    // remove vertex v because of its degree
                    toVisit.addAll(neighbors);
                    adjacency.remove(v).forEach( u -> adjacency.get(u).remove(v));
                    neighborhood.remove(v);
                    removeTriangles = triangles != null;
                }
                if(removeTriangles)
                    removeTrianglesOfVertex(triangles, v);

                toVisit.remove(v);
                break;
            }
        }
        return neighborhood.size() > currentBest;
    }

    /**
     * Method used in LB2 for removing vertices successively. A single vertex removal also removes some of the
     * triangles and therefore may make other vertices invalid.
     * @param graph the neighborhood of vertices
     * @param triangles the map of triangles
     * @param initialDelete the vertex to start the removing with
     * @param l the algorithms setting for L
     * @param min the minimum size the neighborhood should have (aborts of reached)
     */
    private static void removeVerticesSuccessively(RtcGraph graph, int initialDelete, int l, int min,
                                                   HashMap<Integer, HashSet<Triangle>> triangles){
        if(!graph.contains(initialDelete))
            return; // initial delete vertex has already been deleted

        graph.deleteVertex(initialDelete);
        HashSet<Integer> toVisit = new HashSet<>(); // all vertices that lose a triangle in a deletion
        if(triangles.containsKey(initialDelete))
            triangles.get(initialDelete).forEach( t -> toVisit.addAll(t.getVertices()));
        toVisit.remove(initialDelete);
        removeTrianglesOfVertex(triangles, initialDelete);

        while(!toVisit.isEmpty() && graph.size() > min){
            for(Integer v : toVisit){
                int triangleCount = triangles.containsKey(v) ? triangles.get(v).size() : 0;
                if(triangleCount < l){
                    // v needs to get removed
                    graph.deleteVertex(v);
                    if(triangleCount > 0)
                        triangles.get(v).forEach( t -> toVisit.addAll(t.getVertices()));
                    removeTrianglesOfVertex(triangles, v);
                }
                toVisit.remove(v);
                break;
            }
        }
    }

    private static boolean hasLocalEdgeTriangle(HashMap<Integer, HashSet<Triangle>> triangles, int v, int w){
        if(!triangles.containsKey(v))
            return false;

        for(Triangle t : triangles.get(v))
            if(t.contains(w))
                return true;

        return false;
    }

    private static void removeLowTriangleEdges(RtcGraph graph, int l){
        for(int v : graph.getVertices()){
            for(int w : graph.getNeighbors(v)){
                if(v < w && graph.adjacent(v, w)
                        && graph.getCommonNeighbors(v, w).size() < l){
                    graph.getNeighbors(v).remove(w);
                    graph.getNeighbors(w).remove(v);
                }
            }
        }
    }

    private static boolean haveCommonNeighbor(HashMap<Integer, HashSet<Integer>> adjacency, int v, int w){
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

    public static boolean verticesAreCompatible(RtcGraph graph, int v, int w){
        if(v == w) return true;
        if(graph == null) return false;

        if(graph.getNeighbors(v).contains(w))
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

    public static int countCompatibles2(RtcGraph graph, int v, int max){
        if(graph == null || !graph.contains(v))
            return 0;

        int compatibles = 0;
        int n = graph.size();
        for(int w : graph.getVertices()){
            n--;
            if(w!= v && verticesAreCompatible(graph, v, w)){
                if(++compatibles > max)
                    break;
                if (compatibles + n <= max)
                    break;
            }
        }

        return compatibles;
    }

    private static ArrayList<Integer> getVerticesByAscendingDegree(RtcGraph g, boolean twoNeighborhood){
        ArrayList<Integer> startPoints = new ArrayList<>(g.size());
        HashMap<Integer, Integer> degrees = new HashMap<>();
        int vDeg = 0;

        for(int v : g.getVertices()){
            startPoints.add(v);
            vDeg = twoNeighborhood ? g.getTwoNeighbors(v).size() : g.getNeighbors(v).size();
            degrees.put(v, vDeg);
        }

        return degreeQuickSort(startPoints, degrees);
    }

    private static ArrayList<Integer> degreeQuickSort(ArrayList<Integer> toSort, HashMap<Integer, Integer> degrees){
        if(toSort.isEmpty())
            return toSort;
        int pivotElem = toSort.get((int) (Math.random() * (toSort.size())));
        int pivotVal = degrees.get(pivotElem);
        ArrayList<Integer> equalsList = new ArrayList<Integer>();
        ArrayList<Integer> smallerList = new ArrayList<Integer>();
        ArrayList<Integer> greaterList = new ArrayList<Integer>();
        for (int elem : toSort) {
            if (degrees.get(elem) < pivotVal) smallerList.add(elem);
            else if (degrees.get(elem) > pivotVal) greaterList.add(elem);
            else equalsList.add(elem);
        }
        ArrayList<Integer> solution = new ArrayList<Integer>();
        if (!smallerList.isEmpty()) solution.addAll(degreeQuickSort(smallerList, degrees));
        solution.addAll(equalsList);
        if (!greaterList.isEmpty()) solution.addAll(degreeQuickSort(greaterList, degrees));
        return solution;
    }

    /*
    // This one may be used instead of a global conflict graph
    private ITriangleTwoClubAlgorithm getAlgorithmIntefrace(RtcGraph graph, int bestSize){
        final T2CGraph g;
        Output o = new Output();
        ITriangleTwoClubAlgorithm a = new ITriangleTwoClubAlgorithm() {
            @Override
            public RtcGraph getGraph() {
                return g.getGraph();
            }

            @Override
            public RtcGraph getConflictGraph() {
                return g.getConflictGraph();
            }

            @Override
            public HashSet<Integer> getMarked() {
                return g.getMarkedVertices();
            }

            @Override
            public HashMap<Integer, HashSet<Triangle>> getTriangles() {
                return null;
            }

            @Override
            public ConditionChecker getConditionChecker() {
                return new ConditionChecker(this);
            }

            @Override
            public Output getOutput() {
                return o;
            }

            @Override
            public int getL() {
                return algorithm.getL();
            }

            @Override
            public int getBestSolutionSize() {
                return bestSize;
            }

            @Override
            public boolean isEdgeConditionEnabled() {
                return algorithm.isEdgeConditionEnabled();
            }

            @Override
            public void abortBranch() {

            }

            @Override
            public void deleteVertex(int v) {
                g.deleteVertex(v);
            }

            @Override
            public void deleteEdge(int v, int w) {
                g.deleteEdge(v, w);
            }

            @Override
            public void markVertex(int v) {
                g.markVertex(v);
            }
        };
        g = new T2CGraph(graph, a);
        return a;
    }
    */
}
