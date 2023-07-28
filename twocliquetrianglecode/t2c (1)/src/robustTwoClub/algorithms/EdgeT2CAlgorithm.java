package robustTwoClub.algorithms;

import robustTwoClub.graph.Edge;
import robustTwoClub.graph.RollbackRtcGraph;
import robustTwoClub.graph.RtcGraph;
import robustTwoClub.graph.Triangle;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.*;

public class EdgeT2CAlgorithm{

    public enum Branching {SINGLE_VERTEX, INCOMPATIBLE_VERTICES}

    private String graphName;        // name/filename of the graph
    private boolean printDetails, currentlyReducing, abortBranching, useVertexCoverRule = true,
            useNoChoiceRule = true, currentlyBranching;

    private final DecimalFormat tf = new DecimalFormat("0.000");		// Format for time output
    private static final DecimalFormat df = new DecimalFormat("00.00");		// Format for progress output

    private ArrayList<Integer> startPointList;
    private int[] startPointSizes;
    private int currentStartPoint;
    private RollbackRtcGraph inputGraph, workingGraph, bestSolution;
    private RtcGraph conflictGraph;
    private Stack<Integer> markRollbackStack;    // stack to undo delete operations in right order
    private HashSet<Integer> marked;
    private HashMap<Object, List<Edge>> deletionConflicts;
    private HashMap<Integer, HashSet<Triangle>> triangles;

    private long	rollbackTime, startPointBuildTime, triangleCheckTime,
            branchReductionTime, compatibilityChecks, initialReductionTime, lastStatusPrintTime,
            branchReductionDeletions, lowerBoundCalcTime, runStartTime, timeLimitMS = 0,
            vertexCoverRuleTime, noChoiceRuleTime;
    private int initialGraphSize, initialEdgeCount, r3InitialRemoves, totalBranches,
            maxBranchDepth, startPointSkips, startPointCounter, branchReductionSteps = 1, branchDeptLimit,
            bestBranchSolution, lowerBound, upperBoundSkips, noChoiceRuleMarked, vertexCoverRuleSteps = 1;
    private PrintWriter fileWriter = null;
    private TriangleTwoClubAlgorithm.Branching branchingMethod = TriangleTwoClubAlgorithm.Branching.SINGLE_VERTEX;

    public EdgeT2CAlgorithm(){
        resetCounters();
    }

    public RollbackRtcGraph run(RollbackRtcGraph graph, String name){
        inputGraph = graph;
        workingGraph = graph.getClone();
        graphName = name;
        initialGraphSize = TriangleTwoClubAlgorithm.getHighestVertexID(graph);
        initialEdgeCount = graph.getEdgeCount();
        branchDeptLimit = inputGraph.size() + 10;

        runStartTime = System.currentTimeMillis();
        lastStatusPrintTime = runStartTime;


        // initial R3 reduction (remove vertices with low degree or edges with no triangle)
        long time = System.nanoTime();
        r3InitialRemoves += removeNonTriangleEdges(true);
        r3InitialRemoves += TriangleTwoClubAlgorithm.removeLowDegreeVertices(workingGraph);
        initialReductionTime += System.nanoTime() - time;

        // count triangles
        triangles = Triangle.getTriangleHashMap(workingGraph);

        // prepare start points
        time = System.nanoTime();
        startPointSizes = new int[initialGraphSize];
        ArrayList<Integer> unsortedIDs = new ArrayList<>();
        for(int id: workingGraph.getVertices()){
            startPointSizes[id] = workingGraph.sizeOfTwoNeighborhood(id, true);
            unsortedIDs.add(id);
        }
        startPointList = kernelQuickSort(unsortedIDs);
        startPointBuildTime += System.nanoTime() - time;

        // determine lower bound
        time = System.nanoTime();
        bestSolution = (RollbackRtcGraph) TriangleTwoClubAlgorithm.getLowerBound(workingGraph);
        lowerBound = bestSolution == null ? 0 : bestSolution.size();
        lowerBoundCalcTime += System.nanoTime() - time;
        printOutputLine("initial lower bound: " + lowerBound,1);

        // search for triangle two clubs
        // start points are sorted ascending based on size of their 2-neighbourhood
        printOutputLine("branching started", 1);
        for(int i = startPointList.size() - 1; i >= 0 ; i--){
            if(abortBranching) break;
            int upperBound = TriangleTwoClubAlgorithm.calculateLocalUpperBound(workingGraph, startPointList.get(i));
            if(bestSolution != null){
                if(bestSolution.size() >= startPointSizes[i]){
                    startPointSkips++;
                    continue;
                }
                if(upperBound < lowerBound){
                    startPointSkips++;
                    upperBoundSkips++;
                    continue;
                }
            }
            searchAtStartPoint(i);
            startPointCounter++;
        }

        // end
        time = System.currentTimeMillis() - runStartTime;
        printInfo();
        double seconds = (double)time / 1000.0;
        int minutes = (int)(seconds / 60);
        printOutputLine("Finished after " + minutes + "m " + tf.format(seconds%60)
                + "s" + " (" + tf.format(seconds) + "s)");
        printOutputLine("Result size: " + (bestSolution == null ? "no solution" : bestSolution.size() ));

        if(fileWriter != null){
            fileWriter.flush();
            fileWriter.close();
            fileWriter = null; // no multiple runs in a single file
        }
        return bestSolution;
    }

    private void searchAtStartPoint(int startPointID){
        currentStartPoint = startPointList.get(startPointID);
        marked = new HashSet<>();
        markRollbackStack = new Stack<>();
        markVertex(currentStartPoint);

        Edge rollbackPoint = workingGraph.getRollbackPoint();

        // data reduction
        HashSet<Integer> deleteSet = new HashSet<>();
        long time = System.nanoTime();
        // remove the vertices incompatible to marked ones
        for(int v : marked){
            for(int w : workingGraph.getVertices()){
                if(!TriangleTwoClubAlgorithm.verticesAreCompatible(workingGraph, v, w)) deleteSet.add(w);
            }
        }
        branchReductionDeletions += deleteSet.size();
        for(int v : deleteSet) deleteVertex(v);
        branchReductionDeletions += removeNonTriangleEdges();
        branchReductionDeletions += TriangleTwoClubAlgorithm.removeLowDegreeVertices(workingGraph);
        branchReductionTime += System.nanoTime() - time;

        // build conflict graph
        if(useVertexCoverRule){
            RtcGraph con = TriangleTwoClubAlgorithm.buildConflictGraph(workingGraph);
            conflictGraph = new RollbackRtcGraph(con);
            deletionConflicts = new HashMap<>();
        }

        // branching
        bestBranchSolution = 0;
        currentlyBranching = true;
        branch(1);
        currentlyBranching = false;

        conflictGraph = null;
        rollback(rollbackPoint);
        marked.remove(currentStartPoint);
        deleteVertex(currentStartPoint); // old center points can not be part of another
    }

    public void branch(int depth){
        totalBranches++;
        if(abortBranching) return;

        if(depth > branchDeptLimit){
            printOutputLine("Branch depth limit of " + branchDeptLimit + " reached", 1);
            return;
        }

        if(maxBranchDepth < depth) maxBranchDepth = depth;
        if(timeLimitMS > 0 && System.currentTimeMillis() > runStartTime + timeLimitMS){
            abortBranching = true;
            printOutputLine("Time limit reached, abort branching ...", 1);
            return;
        }
        if(printDetails && System.currentTimeMillis() > lastStatusPrintTime + 10000){
            printOutputLine("startpoint " + startPointCounter + "/" + startPointList.size() + ", total branches: "
                    + totalBranches + ", max depth: " + maxBranchDepth + ", current best: "
                    + (bestSolution == null ? 0 : bestSolution.size()), 1);
            lastStatusPrintTime = System.currentTimeMillis();
        }

        // marking rule (no choice rule)
        if(useNoChoiceRule)
            applyNoChoiceRule();
        int markRollbackPoint = -1;
        if (!markRollbackStack.isEmpty()) markRollbackPoint = markRollbackStack.peek();

        Edge branchRollbackPoint = workingGraph.getRollbackPoint();

        long time = 0L;
        currentlyReducing = depth % branchReductionSteps == 0;
        if(currentlyReducing){
            // data reduction
            HashSet<Integer> deleteSet = new HashSet<>();
            time = System.nanoTime();
            // remove the vertices incompatible to marked ones
            for(int v : marked){
                for(int w : workingGraph.getVertices()){
                    if(v == w || marked.contains(w)) continue;
                    if(!TriangleTwoClubAlgorithm.verticesAreCompatible(workingGraph, v, w)) deleteSet.add(w);
                }
            }
            branchReductionDeletions += deleteSet.size();
            for(int v : deleteSet) deleteVertex(v);
            branchReductionDeletions += TriangleTwoClubAlgorithm.removeLowDegreeVertices(workingGraph); // TODO: change
            branchReductionDeletions += removeNonTriangleEdges();
            branchReductionTime += System.nanoTime() - time;
        }

        // check if all marked vertices are still valid
        for(int i : marked) {
            if (!TriangleTwoClubAlgorithm.isPartOfTriangle(workingGraph, i))
                return;
            for (int j : marked) {
                if (i < j) continue;
                if((useVertexCoverRule && conflictGraph.adjacent(i, j))
                        || !TriangleTwoClubAlgorithm.verticesAreCompatible(workingGraph, i, j))
                    return; // two marked vertices are incompatible
            }
        }


        // if current working branch is smaller than the best current solution we can abort this branch
        if(workingGraph.size() < 2
                || (bestSolution != null
                && (workingGraph.size() <= bestSolution.size() || workingGraph.size() < lowerBound))) return;

        // if current working tree is a triangle two club this branch is finished
        time = System.nanoTime();
        if(isEdgeTriangleTwoClub(workingGraph)){
            bestBranchSolution = bestBranchSolution < workingGraph.size() ? workingGraph.size() : bestBranchSolution;
            if(bestSolution == null || workingGraph.size() > bestSolution.size()){
                // new best solution
                bestSolution = workingGraph.getClone();
            }
            return;
        }
        triangleCheckTime += System.nanoTime() - time;

        // vertex cover rule
        if(useVertexCoverRule && depth % vertexCoverRuleSteps == 0){
            time = System.nanoTime();
            int limit = workingGraph.size() - lowerBound - 1;
            int vc = TriangleTwoClubAlgorithm.matchingBoundVertexCover(workingGraph, conflictGraph, limit);
            vertexCoverRuleTime += System.nanoTime() - time;
            if(workingGraph.size() - vc <= lowerBound) return;
        }

        // find any two incompatible nodes
        int v = getMostIncompatibleVertex(workingGraph);
        int w = -1;
        if(branchingMethod == TriangleTwoClubAlgorithm.Branching.INCOMPATIBLE_VERTICES){
            // find a vertex incompatible with v
            if(useVertexCoverRule)
                // use conflict graph if possible
                w = conflictGraph.getNeighbors(v).iterator().next();
            else
                for(int i : workingGraph.getVertices()){
                    if(!TriangleTwoClubAlgorithm.verticesAreCompatible(workingGraph, v, i)){
                        w = i;
                        break;
                    }
                }
        }

        // branching
        Edge rbp = workingGraph.getRollbackPoint();

        if(branchingMethod == TriangleTwoClubAlgorithm.Branching.INCOMPATIBLE_VERTICES){
            // branch without w
            deleteVertex(v);
            branch(depth + 1);
            rollback(rbp);

            // branch with marked w and deleted v
            markVertex(v);
            deleteVertex(w);
            branch(depth + 1);
            rollbackMarkedVertices(markRollbackPoint);
            rollback(rbp);
        }else{
            // single vertex branching
            // branch without v
            deleteVertex(v);
            branch(depth + 1);
            rollback(rbp);
            rollbackMarkedVertices(markRollbackPoint);

            // branch with v
            markVertex(v);
            branch(depth + 1);
            rollback(rbp);
            rollbackMarkedVertices(markRollbackPoint);
        }

        if(currentlyReducing){
            rollback(branchRollbackPoint);
        }
    }

    private void applyNoChoiceRule(){
        long time = System.nanoTime();
        LinkedList<Integer> toMark = new LinkedList<>();

        // single vertex between two marked or multiple connected ones
        for(int v : marked)
            for(int w : marked)
                if(v < w){
                    HashSet<Integer> commonNeighbours = workingGraph.getCommonNeighbors(v, w);
                    if(commonNeighbours.size() == 1) {
                        toMark.addAll(commonNeighbours); // single vertex between two marked
                    }
                }

        toMark.forEach(this::markVertex);
        noChoiceRuleTime += System.nanoTime() - time;
        noChoiceRuleMarked += toMark.size();
    }

    public int removeNonTriangleEdges(){
        HashSet<Integer> seen = new HashSet<>();
        List<Edge> toDelete = new LinkedList<>();
        for(Integer i : workingGraph.getVertices()){
            if(seen.contains(i)) continue;
            for(Integer j : workingGraph.getNeighbors(i)){
                if(seen.contains(j)) continue;
                if(!isPartOfEdgeTriangle(workingGraph, i, j)){
                    toDelete.add(new Edge(i, j));
                }
            }
            seen.add(i);
        }
        for(Edge e : toDelete){
            deleteEdge(e);
        }
        return toDelete.size();
    }

    public int removeNonTriangleEdges(boolean toExhaustion){
        if(!toExhaustion)
            return removeNonTriangleEdges();
        int old_i, i = 0;
        do{
            old_i = i;
            i += removeNonTriangleEdges();
        }while(old_i < i);
        return i;
    }

    private void markVertex(int v){
        if(workingGraph != null && workingGraph.contains(v)){
            marked.add(v);
            markRollbackStack.push(v);
        }else{
            if(workingGraph == null)
                printOutputLine("bad mark: no working graph", 1);
            else
                printOutputLine("bad mark, no vertex " + v + " in working graph", 1);
        }
    }

    private void rollbackMarkedVertices(int rollbackPoint){
        if(rollbackPoint < -1) return;
        long time = System.nanoTime();
        while(!markRollbackStack.empty()
                && markRollbackStack.peek() != rollbackPoint
                && markRollbackStack.peek() != currentStartPoint){
            marked.remove(markRollbackStack.pop());
        }
        rollbackTime += System.nanoTime() - time;
    }

    private void rollback(Edge rollbackPoint){
        long time = System.nanoTime();
        workingGraph.rollback(rollbackPoint);
        if(useVertexCoverRule && currentlyBranching)
            rollbackConflictGraph(rollbackPoint);
        rollbackTime += System.nanoTime() - time;
    }

    private void deleteVertex(int v){
        if(workingGraph != null && workingGraph.contains(v)){
            if(useVertexCoverRule && currentlyBranching)
                deleteInConflictGraph(v);
            workingGraph.deleteVertex(v);
            if(marked != null && marked.contains(v)){
                printOutputLine("Deleted marked vertex " + v + "!", 1);
                marked.remove(v);
            }
        }else{
            if(workingGraph == null)
                printOutputLine("bad delete: no working graph", 1);
            else
                printOutputLine("bad delete, no vertex " + v + " in working graph", 1);
        }
    }

    private void deleteEdge(Edge e){
        if(workingGraph != null && workingGraph.contains(e.getV()) && workingGraph.contains(e.getW())){
            if(useVertexCoverRule && currentlyBranching)
                deleteEdgeInConflictGraph(e);
            workingGraph.deleteEdge(e);
            if(marked != null && marked.contains(e.getV())){
                printOutputLine("Deleted marked vertex " + e.getV() + "!", 1);
                marked.remove(e.getV());
            }
            if(marked != null && marked.contains(e.getW())){
                printOutputLine("Deleted marked vertex " + e.getW() + "!", 1);
                marked.remove(e.getW());
            }
        }else{
            if(workingGraph == null)
                printOutputLine("bad delete: no working graph", 1);
            else
                printOutputLine("bad delete, a vertex of the edge is not in working graph", 1);
        }
    }

    private void deleteInConflictGraph(int v){
        if(conflictGraph == null){
            printOutputLine("algorithm should use VCR but the conflict graph is null", 1);
            return;
        }
        HashSet<Integer> toVisit = workingGraph.getTwoNeighbors(v);
        Edge rbp = workingGraph.getRollbackPoint();
        workingGraph.deleteVertex(v);
        for(Integer x : toVisit){
            for(Integer y : toVisit){
                if(y < x) continue;
                if(!conflictGraph.adjacent(x, y) && !TriangleTwoClubAlgorithm.verticesAreCloseEnough(workingGraph, x, y)){
                    // x and y are in conflict after deleting v
                    deletionConflicts.computeIfAbsent(v, k -> new LinkedList<>());
                    deletionConflicts.get(v).add(new Edge(x, y));
                    conflictGraph.addEdge(x, y);
                }
            }
        }
        conflictGraph.deleteVertex(v);
        workingGraph.rollback(rbp); // v gets deleted by deleteVertex(int v)
    }

    private void deleteEdgeInConflictGraph(Edge e){
        if(conflictGraph == null){
            printOutputLine("algorithm should use VCR but the conflict graph is null", 1);
            return;
        }
        HashSet<Integer> toVisit = workingGraph.getTwoNeighbors(e.getV());
        toVisit.addAll(workingGraph.getTwoNeighbors(e.getW()));
        Edge rbp = workingGraph.getRollbackPoint();
        workingGraph.deleteEdge(e.getV(), e.getW());
        for(Integer x : toVisit){
            for(Integer y : toVisit){
                if(y < x) continue;
                if(!conflictGraph.adjacent(x, y) && !TriangleTwoClubAlgorithm.verticesAreCloseEnough(workingGraph, x, y)){
                    // x and y are in conflict after deleting v
                    deletionConflicts.computeIfAbsent(e, k -> new LinkedList<>());
                    deletionConflicts.get(e).add(new Edge(x, y));
                    conflictGraph.addEdge(x, y);
                }
            }
        }
        workingGraph.rollback(rbp); // add removed edge again
    }

    private void rollbackConflictGraph(Edge rollbackPoint){
        if(conflictGraph == null){
            printOutputLine("can not rollback conflict graph: it is null",1);
            return;
        }
        if(deletionConflicts == null){
            printOutputLine("can not rollback conflict graph: conflict list is null",1);
            return;
        }
        Stack<Edge> rb = workingGraph.getRollbacksTo(rollbackPoint);
        Edge edge;
        while(!rb.empty()){
            edge = rb.pop();
            List<Edge> conflictsToRemove;

            if(edge.hasNull())
                conflictsToRemove = deletionConflicts.get(edge.getNonNullVertex());
            else
                conflictsToRemove = deletionConflicts.get(edge);

            for(Edge e : conflictsToRemove){
                if(e.hasNull()){
                    // it's a vertex
                    conflictGraph.undeleteVertex(e.getNonNullVertex());
                }else{
                    // it's an edge
                    conflictGraph.getNeighbors(e.getV()).remove(e.getW());
                    conflictGraph.getNeighbors(e.getW()).remove(e.getV());
                }
            }
        }
        deletionConflicts.remove(rollbackPoint);
    }

    /**
     * naive approach to determine if an edge is part of a edge triangle
     * @param graph a graph
     * @param v id of the vertex at one end of an edge
     * @param w id of the vertex on the oder side of the edge
     * @return true when edge between v and w is part of an edge triangle in the graph, false else
     */
    public static boolean isPartOfEdgeTriangle(RtcGraph graph, int v, int w){
        if(!graph.adjacent(v, w)) return false;
        return graph.getCommonNeighbors(v, w).size() > 0; // a common neighbor
    }

    /**
     * naive method to test if a graph is an edge triangle two club
     * @param graph the input graph
     * @return true when the graph is an edge triangle two club, false else
     */
    public static boolean isEdgeTriangleTwoClub(RtcGraph graph){
        for(Integer v : graph.getVertices()) {
            for (Integer w : graph.getNeighbors(v))
                if (v < w && !isPartOfEdgeTriangle(graph, v, w))
                    return false;

            for(Integer w : graph.getVertices())
                if(v < w && !TriangleTwoClubAlgorithm.verticesAreCompatible(graph, v ,w))
                    return false;

        }
        return true;
    }

    /**
     * Finds the vertex that has the most incompatibilities to other vertices in the graph
     * @param graph the graph
     * @return the vertex id, -1 when no vertices are in the graph or graph is null
     */
    public int getMostIncompatibleVertex(RtcGraph graph){
        if(conflictGraph != null){
            int worstVertex = -1;
            int incompatibilities = 0;
            for(Integer v : conflictGraph.getVertices()){
                int i = conflictGraph.getNeighbors(v).size();
                if(i > incompatibilities && !marked.contains(v) && workingGraph.contains(v)){
                    incompatibilities = i;
                    worstVertex = v;
                }
            }
            if(worstVertex != -1) return worstVertex;
        }
        if(graph == null)
            return -1;
        int worstVertex = -1;
        int incompatibilities = 0;
        for(int i : graph.getVertices()){
            if(marked.contains(i)) continue;
            int currentCount = 0;
            for(int j : graph.getVertices()){
                if(!TriangleTwoClubAlgorithm.verticesAreCloseEnough(graph, i, j))
                    currentCount++;
            }
            if(currentCount > incompatibilities){
                worstVertex = i;
                incompatibilities = currentCount;
            }
        }
        return worstVertex;
    }

    private ArrayList<Integer> kernelQuickSort(ArrayList<Integer> toSort){
        if(toSort.isEmpty())
            return toSort;
        int pivotElem = toSort.get((int) (Math.random() * (toSort.size())));
        int pivotVal = startPointSizes[pivotElem];
        ArrayList<Integer> equalsList = new ArrayList<Integer>();
        ArrayList<Integer> smallerList = new ArrayList<Integer>();
        ArrayList<Integer> greaterList = new ArrayList<Integer>();
        for (int elem : toSort) {
            if (startPointSizes[elem] < pivotVal) smallerList.add(elem);
            else if (startPointSizes[elem] > pivotVal) greaterList.add(elem);
            else equalsList.add(elem);
        }
        ArrayList<Integer> solution = new ArrayList<Integer>();
        if (!smallerList.isEmpty()) solution.addAll(kernelQuickSort(smallerList));
        solution.addAll(equalsList);
        if (!greaterList.isEmpty()) solution.addAll(kernelQuickSort(greaterList));
        return solution;
    }

    /**
     * Gets the current number of triangles an edge is part of based on the working graph
     * @param v the id of the first vertex of the edge
     * @param w the id of the second vertex of the edge
     * @return the number of triangles in v is part of in the working graph,
     * -1 if v does not exist or the triangle map is not initialized
     */
    private int countTrianglesFromMap(int v, int w){
        return countTrianglesFromMap(v, w, Integer.MAX_VALUE);
    }

    /**
     * Gets the current number of triangles an edge is part of based on the working graph
     * @param v the id of the first vertex of the edge
     * @param w the id of the second vertex of the edge
     * @param max the maximum number od triangles to count to
     * @return the number of triangles in v is part of in the working graph,
     * -1 if v does not exist or the triangle map is not initialized
     */
    private int countTrianglesFromMap(int v, int w, int max){
        if(!workingGraph.contains(v) || triangles == null)
            return -1;

        if(!triangles.containsKey(v) || !triangles.containsKey(w)) // v or w was not part of any triangle at all
            return 0;

        int n = 0;
        for(Triangle t : triangles.get(v)){
            if(t.exists(workingGraph) && t.contains(w)) // only count triangles that currently exist and have w
                n++;
            if(n >= max)
                break;
        }

        return n;
    }


    /**
     * Enables or disables console output
     * @param b if true the results will be printed on console
     */
    public void printDetails(boolean b){
        printDetails = b;
    }

    public void setOutputFile(String path){
        File f = new File(path);
        try {
            fileWriter = new PrintWriter(new FileWriter(f), true);
        } catch (IOException e) {
            printOutputLine("Failed to initialize file writing", 1);
        }
    }

    /**
     * Method for testing purposes.
     * @param graph the new working graph
     */
    public void setWorkingGraph(RollbackRtcGraph graph){
        workingGraph = graph;
    }

    /**
     * Method for testing purposes.
     */
    public RollbackRtcGraph getWorkingGraph(){
        return workingGraph;
    }

    /**
     * Method for testing purposes.
     * @param m the new set of marked vertices
     */
    public void setMarked(HashSet<Integer> m){
        marked = m;
    }

    /**
     * Writes output to console and/or file
     * @param line the line to output
     * @param mode 0 = file + console, 1 = console only, 2 = file only
     */
    private void printOutputLine(String line, int mode){
        if(printDetails && mode < 2) System.out.println(line);
        if(fileWriter != null && (mode == 0 || mode == 2)) fileWriter.println(line);
    }

    private void printOutputLine(String line){
        printOutputLine(line, 0);
    }

    private void printInfo(){
        // add used settings to file
        if(fileWriter != null){
            printOutputLine("Settings: ", 2);
            printOutputLine("Branching Method: "
                    + (branchingMethod == TriangleTwoClubAlgorithm.Branching.SINGLE_VERTEX ? "Single Vertex" : "Incompatible Vertices"),2);
            printOutputLine("Branch reduction steps: " + branchReductionSteps,2);
            printOutputLine("Vertex Cover Rule used: " + useVertexCoverRule, 2);
            if(useVertexCoverRule)
                printOutputLine("Vertex Cover Rule steps: " + vertexCoverRuleSteps,2);
            printOutputLine("No Choice Rule used: " + useNoChoiceRule, 2);
        }
        printOutputLine("");

        printOutputLine("Graph name: " + graphName);
        printOutputLine("Vertices: " + initialGraphSize + "\t Edges: " + initialEdgeCount);
        printOutputLine("Initial Reduction: " + r3InitialRemoves + " deleted ("
                + df.format((double) r3InitialRemoves / initialGraphSize * 100.0) + "%) in "
                + tf.format(initialReductionTime / 1000000000.0) + "s");

        int kernelSum = 0, kernelMin = initialGraphSize, kernelMax = -1;
        for(int i : startPointSizes){
            kernelSum += i;
            if(i < kernelMin && i > 0) kernelMin = i; // do not consider empty kernels
            if(i > kernelMax) kernelMax = i;
        }
        printOutputLine("Startpoints: " + startPointList.size() + "\t minSize: " + kernelMin + "\t maxSize: "
                + kernelMax + "\t avgSize: " + df.format((double) kernelSum / startPointList.size()));
        printOutputLine("Startpoints visited: " + startPointCounter);
        printOutputLine("Initial lower bound: " + lowerBound);
        printOutputLine("Time spent sorting start points: "
                + tf.format(startPointBuildTime / 1000000000.0) + "s");
        printOutputLine("Time spent calculating lower bound: "
                + tf.format(lowerBoundCalcTime / 1000000000.0) + "s");
        printOutputLine("Skipped start points: " + startPointSkips + " ("
                + df.format((double) startPointSkips / startPointList.size() * 100.0) + "%)");
        printOutputLine("Upper bound skips: " + upperBoundSkips + " ("
                + df.format((double) upperBoundSkips / startPointSkips * 100.0) + "% of skipped points)");
        printOutputLine("Total number of branches: " + totalBranches);
        printOutputLine("Maximum branch depth: " + maxBranchDepth);
        printOutputLine("Branch reduction time: " + tf.format(branchReductionTime / 1000000000.0) + "s");
        printOutputLine("Branch reduction deletions: " + branchReductionDeletions + " ("
                + String.format("%.2f", (double)branchReductionDeletions / totalBranches) + " per branch avg)");
        printOutputLine("ttclub check time: " + tf.format(triangleCheckTime / 1000000000.0) + "s");
        if(useVertexCoverRule)
            printOutputLine("Vertex cover rule time: " + tf.format(vertexCoverRuleTime / 1000000000.0) + "s");
        if(useNoChoiceRule)
            printOutputLine("No-Choice rule: marked " + noChoiceRuleMarked + " in "
                    + tf.format(noChoiceRuleTime / 1000000000.0) + "s");
        printOutputLine("Rollback time: " + tf.format(rollbackTime / 1000000000.0) + "s");
        printOutputLine("Number of compatibility checks: " + compatibilityChecks);
        if(timeLimitMS > 0){
            printOutputLine("Time limit: " + df.format(timeLimitMS/1000.0) + "s" );
            printOutputLine("Time limit reached: " + abortBranching);
        }
    }

    public void resetCounters(){
        // initialize
        inputGraph = null;
        workingGraph = null; // g will be our copy to work with
        bestSolution = null;
        conflictGraph = null;
        graphName = "";
        initialGraphSize = -1;
        initialEdgeCount = -1;
        printDetails = false;
        currentlyReducing = false;
        abortBranching = false;
        currentlyBranching = false;
        branchDeptLimit = -1;

        // reset counts
        totalBranches = 0;
        maxBranchDepth = 0;
        lowerBound = 0;
        compatibilityChecks = 0;
        branchReductionTime = 0L;
        initialReductionTime = 0L;
        lowerBoundCalcTime = 0L;
        vertexCoverRuleTime = 0L;
        noChoiceRuleTime = 0L;
        rollbackTime = 0L;
        triangleCheckTime = 0L;
        r3InitialRemoves = 0;
        startPointSkips = 0;
        upperBoundSkips = 0;
        noChoiceRuleMarked = 0;
        startPointCounter = 1;
        branchReductionDeletions = 0L;
    }
}
