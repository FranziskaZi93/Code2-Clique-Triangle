package robustTwoClub.algorithms;

import robustTwoClub.graph.Edge;
import robustTwoClub.graph.RtcGraph;
import robustTwoClub.graph.Triangle;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class TriangleTwoClubAlgorithm {

    public enum Branching {SINGLE_VERTEX, INCOMPATIBLE_VERTICES}

    private static String graphName;        // name/filename of the graph
    private static boolean printDetails, currentlyReducing, abortBranching, useVertexCoverRule = true,
            useNoChoiceRule = true;

    private static final DecimalFormat tf = new DecimalFormat("0.000");		// Format for time output
    private static final DecimalFormat df = new DecimalFormat("00.00");		// Format for progress output

    private static ArrayList<Integer> startPointList;
    private static int[] startPointSizes;
    private static int currentStartPoint, ell = 1;
    private static RtcGraph inputGraph, workingGraph, bestSolution, conflictGraph;
    private static Stack<Integer> rollbackStack, markRollbackStack;    // stack to undo delete operations in right order
    private static HashSet<Integer> marked;
    private static HashMap<Integer, List<Edge>> deletionConflicts;
    private static HashMap<Integer, HashSet<Triangle>> triangles;

    private static long	rollbackTime, startPointBuildTime, triangleCheckTime,
            branchReductionTime, compatibilityChecks, initialReductionTime, lastStatusPrintTime,
            branchReductionDeletions, lowerBoundCalcTime, runStartTime, timeLimitMS = 0,
            vertexCoverRuleTime, noChoiceRuleTime, triangleComputationTime;
    private static int initialGraphSize, initialEdgeCount, r3InitialRemoves, totalBranches,
            maxBranchDepth, startPointSkips, startPointCounter, branchReductionSteps = 1, branchDeptLimit,
            bestBranchSolution, lowerBound, upperBoundSkips, noChoiceRuleMarked, vertexCoverRuleSteps = 1;
    private static PrintWriter fileWriter = null;
    private static Branching branchingMethod = Branching.SINGLE_VERTEX;

    /**
     * Finds the largest Triangle Two Club in a graph
     * @param graph input graph
     * @param name the name of the graph (shown in results)
     * @param printDetails if true the results will be printed
     * @return the subgraph that is the largest TTC or null if there is no TTC
     */
    public static RtcGraph run(RtcGraph graph, String name, boolean printDetails){
        // initialize
        inputGraph = graph;
        workingGraph = graph.getClone(); // g will be our copy to work with
        bestSolution = null;
        conflictGraph = null;
        graphName = name;
        initialGraphSize = getHighestVertexID(graph);
        initialEdgeCount = graph.getEdgeCount();
        TriangleTwoClubAlgorithm.printDetails = printDetails;
        rollbackStack = new Stack<>();
        currentlyReducing = false;
        abortBranching = false;
        branchDeptLimit = inputGraph.size() + 10;

        // reset counts
        resetCounts();

        runStartTime = System.currentTimeMillis();
        lastStatusPrintTime = runStartTime;

        // compute all triangles in the graph
        long time = System.nanoTime();
        triangles = Triangle.getTriangleHashMap(workingGraph);
        triangleComputationTime += System.nanoTime() - time;

        // initial R3 reduction (remove vertices with low degree or no triangle)
        time = System.nanoTime();
        r3InitialRemoves += removeLowDegreeVertices();
        r3InitialRemoves += removeLNonTriangleVertices();
        initialReductionTime += System.nanoTime() - time;

        removeNonExistingTriangles();

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
        bestSolution = getLowerBound(workingGraph);
        lowerBound = bestSolution == null ? 0 : bestSolution.size();
        lowerBoundCalcTime += System.nanoTime() - time;
        printOutputLine("initial lower bound: " + lowerBound,1);

        // search for triangle two clubs
        // start points are sorted ascending based on size of their 2-neighbourhood
        printOutputLine("branching started", 1);
        for(int i = startPointList.size() - 1; i >= 0 ; i--){
            if(abortBranching) break;
            int upperBound = calculateLocalUpperBound(workingGraph, startPointList.get(i));
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

    private static void searchAtStartPoint(int startPointID){
        currentStartPoint = startPointList.get(startPointID);
        marked = new HashSet<>();
        markRollbackStack = new Stack<>();
        markVertex(currentStartPoint);

        int rollbackPoint = -1;
        if (!rollbackStack.isEmpty()) rollbackPoint = rollbackStack.peek();

        // data reduction
        HashSet<Integer> deleteSet = new HashSet<>();
        long time = System.nanoTime();
        // remove the vertices incompatible to marked ones
        for(int v : marked){
            for(int w : workingGraph.getVertices()){
                if(!verticesAreCompatible(workingGraph, v, w)) deleteSet.add(w);
            }
        }
        branchReductionDeletions += deleteSet.size();
        for(int v : deleteSet) deleteVertex(v);
        branchReductionDeletions += removeLowDegreeVertices();
        branchReductionDeletions += removeNonTriangleVertices();
        branchReductionTime += System.nanoTime() - time;

        // build conflict graph
        if(useVertexCoverRule){
            conflictGraph = buildConflictGraph(workingGraph);
            deletionConflicts = new HashMap<>();
        }

        // branching
        bestBranchSolution = 0;
        branch(1);

        conflictGraph = null;
        rollback(rollbackPoint);
        marked.remove(currentStartPoint);
        deleteVertex(currentStartPoint); // old center points can not be part of another
    }

    private static void branch(int depth){
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

        int branchRollbackPoint = -2;
        long time = 0L;
        currentlyReducing = depth % branchReductionSteps == 0;
        if(currentlyReducing){
            // init rollback point
            if (!rollbackStack.isEmpty()) branchRollbackPoint = rollbackStack.peek();

            // data reduction
            HashSet<Integer> deleteSet = new HashSet<>();
            time = System.nanoTime();
            // remove the vertices incompatible to marked ones
            for(int v : marked){
                for(int w : workingGraph.getVertices()){
                    if(v == w || marked.contains(w)) continue;
                    if(!verticesAreCompatible(workingGraph, v, w)) deleteSet.add(w);
                }
            }
            branchReductionDeletions += deleteSet.size();
            for(int v : deleteSet) deleteVertex(v);
            branchReductionDeletions += removeLowDegreeVertices();
            branchReductionDeletions += removeLNonTriangleVertices();
            branchReductionTime += System.nanoTime() - time;
        }

        // check if all marked vertices are still valid
        for(int i : marked) {
            if (countTrianglesFromMap(i, ell) < ell)
                return;
            for (int j : marked) {
                if (i < j) continue;
                if((useVertexCoverRule && conflictGraph.adjacent(i, j)) || !verticesAreCompatible(workingGraph, i, j))
                    return; // two marked vertices are incompatible
            }
        }
        // if current working branch is smaller than the best current solution we can abort this branch
        if(workingGraph.size() < 2
                || (bestSolution != null
                    && (workingGraph.size() <= bestSolution.size() || workingGraph.size() < lowerBound))) return;

        // if current working tree is a triangle two club this branch is finished
        time = System.nanoTime();
        if(isLTriangleTwoClub()){
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
            int vc = matchingBoundVertexCover(workingGraph, conflictGraph, limit);
            vertexCoverRuleTime += System.nanoTime() - time;
            if(workingGraph.size() - vc <= lowerBound) return;
        }

        // find any two incompatible nodes
        int v = getMostIncompatibleVertex(workingGraph);
        int w = -1;
        if(branchingMethod == Branching.INCOMPATIBLE_VERTICES){
            // find a vertex incompatible with v
            if(useVertexCoverRule)
                // use conflict graph if possible
                w = conflictGraph.getNeighbors(v).iterator().next();
            else
                for(int i : workingGraph.getVertices()){
                    if(!verticesAreCompatible(workingGraph, v, i)){
                        w = i;
                        break;
                    }
                }
        }

        // branching
        int rollbackPoint = getRollbackPoint();

        if(branchingMethod == Branching.INCOMPATIBLE_VERTICES){
            // branch without w
            deleteVertex(v);
            branch(depth + 1);
            rollback(rollbackPoint);

            // branch with marked w and deleted v
            markVertex(v);
            deleteVertex(w);
            branch(depth + 1);
            rollbackMarkedVertices(markRollbackPoint);
            rollback(rollbackPoint);
        }else{
            // single vertex branching
            // branch without v
            deleteVertex(v);
            branch(depth + 1);
            rollback(rollbackPoint);
            rollbackMarkedVertices(markRollbackPoint);

            // branch with v
            markVertex(v);
            branch(depth + 1);
            rollback(rollbackPoint);
            rollbackMarkedVertices(markRollbackPoint);
        }

        if(currentlyReducing){
            rollback(branchRollbackPoint);
        }
    }

    /**
     *  removes every vertex from the input graph that is not part of a triangle
     * @param graph a graph
     * @return the number of vertices deleted
     */
    public static int removeNonTriangleVertices(RtcGraph graph){
        RtcGraph backup = workingGraph;
        RtcGraph conflictBackup = conflictGraph;
        Stack<Integer> backupStack = rollbackStack;
        workingGraph = graph;
        conflictGraph = null;
        rollbackStack = new Stack<>();
        int i = removeNonTriangleVertices();
        workingGraph = backup;
        conflictGraph = conflictBackup;
        rollbackStack = backupStack;
        return i;
    }

    private static int removeNonTriangleVertices(){
        HashSet<Integer> seen = new HashSet<>();
        List<Integer> toDelete = new LinkedList<>();
        outerLoop:
        for(int id : workingGraph.getVertices()){
            if(seen.contains(id)) continue;
            seen.add(id);
            HashSet<Integer> neighbors = workingGraph.getNeighbors(id);
            for(Integer i : neighbors){
                for (Integer j : neighbors){
                    if (!i.equals(j) && workingGraph.adjacent(i, j)) {
                        seen.add(i);
                        seen.add(j);
                        continue outerLoop;
                    }
                }
            }
            if(marked == null || !marked.contains(id))
                toDelete.add(id);
        }
        if(marked != null) toDelete.removeAll(marked);
        for(int id : toDelete) deleteVertex(id);
        return toDelete.size();
    }

    private static int removeLNonTriangleVertices(){
        List<Integer> toDelete = new LinkedList<>();
        for(int v : workingGraph.getVertices()){
            if(countTrianglesFromMap(v, ell) < ell)
                toDelete.add(v);
        }

        for(int v : toDelete)
            deleteVertex(v);
        return toDelete.size();
    }

    /**
     * deletes every vertex with degree lower than two from the graph
     * @param graph the input graph
     * @return the number of vertices deleted
     */
    public static int removeLowDegreeVertices(RtcGraph graph){
        RtcGraph backup = workingGraph;
        RtcGraph conflictBackup = conflictGraph;
        Stack<Integer> backupStack = rollbackStack;
        workingGraph = graph;
        conflictGraph = null;
        rollbackStack = new Stack<>();
        int i = removeLowDegreeVertices();
        workingGraph = backup;
        conflictGraph = conflictBackup;
        rollbackStack = backupStack;
        return i;
    }

    private static int removeLowDegreeVertices(){
        HashSet<Integer> seen = new HashSet<>();
        List<Integer> toDelete = new LinkedList<>();
        for(Integer v : workingGraph.getVertices())                // find vertices to delete
            removeRecursive(seen, toDelete, v);
        if(marked != null) toDelete.removeAll(marked);
        for(Integer v : toDelete) deleteVertex(v);// remove vertices from graph
        return toDelete.size();
    }

    private static void removeRecursive(HashSet<Integer> seen, List<Integer> toDelete, int v){
        if(seen.contains(v)) return;
        seen.add(v);
        if(workingGraph.getNeighbors(v).size() < 2){
            if(marked != null && marked.contains(v))
                return;
            toDelete.add(v);
            for(int i: workingGraph.getNeighbors(v))
                removeRecursive(seen, toDelete, i);
        }
    }

    /**
     * naive approach to determine if a single vertex is part of a triangle
     * @param graph a graph
     * @param v id of a vertex in graph
     * @return true when v is part of a triangle in the graph, false else
     */
    public static boolean isPartOfTriangle(RtcGraph graph, int v){
        HashSet<Integer> neighbors = graph.getNeighbors(v);
        for(Integer i : neighbors)
            for(Integer j : neighbors)
                if(!i.equals(j) && graph.adjacent(i, j))
                    return true;
        return false;
    }

    /**
     * @param graph the graph
     * @param v first vertex
     * @param w second vertex
     * @return true when v and w can both be part of a triangle-two-club, false else
     */
    public static boolean verticesAreCompatible(RtcGraph graph, int v, int w){
        compatibilityChecks++;
        if(conflictGraph != null && graph == workingGraph && conflictGraph.contains(v) && conflictGraph.contains(w))
            return !conflictGraph.adjacent(v, w);   // not adjacent in conflict graph -> compatible
        return verticesAreCloseEnough(graph, v, w);
    }

    /**
     * @param graph a graph
     * @param v first vertex
     * @param w second vertex
     * @return true when both vertices are part of the same triangle
     */
    public static boolean verticesShareTriangle(RtcGraph graph, int v, int w){
        if(v == w)
            return isPartOfTriangle(graph, v);
        if(!graph.adjacent(v, w))
            return false;
        HashSet<Integer> neighborsV = graph.getNeighbors(v);
        HashSet<Integer> neighborsW = graph.getNeighbors(w);
        for(Integer i : neighborsV)
            if(neighborsW.contains(i))
                return true;
        return false;
    }

    /**
     * naive method to test if a graph is a Triangle Two Club
     * @param graph the input graph
     * @return true when the graph is a Triangle Two Club, false else
     */
    public static boolean isTriangleTwoClub(RtcGraph graph){
        outerLoop:
        for(int id: graph.getVertices()){
            // distance to every other vertex max. 2
            for(int j: graph.getVertices())
                if(id != j && !verticesAreCloseEnough(graph, id, j))
                    return false;
            // shares a triangle with at least one neighbor (is not connected over "green" edge)
            for(int j: graph.getNeighbors(id))
                if(verticesShareTriangle(graph, id, j))
                    continue outerLoop;
            return false;
        }
        return true;
    }

    /**
     * method to test if the current working graph is a l - Triangle Two Club,
     * where each vertex is in at least l triangles
     * @return true when the graph is a Triangle Two Club, false else
     */
    public static boolean isLTriangleTwoClub(){
        if(ell < 1) // nothing to check for
            return true;

        for(int v : workingGraph.getVertices()){
            // is v part of ell triangles?
            if(countTrianglesFromMap(v, ell) < ell)
                return false;

            // is the distance between v and every other vertex at most two
            for(int w : workingGraph.getVertices())
                if(!verticesAreCloseEnough(workingGraph, v, w))
                    return false;
        }
        return true;
    }

    /**
     * Sets the steps in which data reduction is applied while branching.
     * For given i, only on branches with depth % i == 0, reduction will be allied
     * Default is 1.
     * @param i
     */
    public static void setBranchReductionSteps(int i){
        if(i > 0 && i < 10)
            TriangleTwoClubAlgorithm.branchReductionSteps = i;
        else
            printOutputLine("Invalid branching step size: " + i, 1);
    }

    private static void markVertex(int v){
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

    private static void rollbackMarkedVertices(int rollbackPoint){
        if(rollbackPoint < -1) return;
        long time = System.nanoTime();
        while(!markRollbackStack.empty()
                && markRollbackStack.peek() != rollbackPoint
                && markRollbackStack.peek() != currentStartPoint){
            marked.remove(markRollbackStack.pop());
        }
        rollbackTime += System.nanoTime() - time;
    }

    private static void rollback(int rollbackPoint){
        if(rollbackPoint < -1) return;
        long time = System.nanoTime();
        while(!rollbackStack.empty() && rollbackStack.peek() != rollbackPoint){
            int v = rollbackStack.pop();
            workingGraph.undeleteVertex(v);
            rollbackConflictGraph(v);
        }
        rollbackTime += System.nanoTime() - time;
    }

    private static void deleteVertex(int v){
        if(workingGraph != null && workingGraph.contains(v)){
            rollbackStack.push(v);
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

    private static void deleteInConflictGraph(int v){
        if(conflictGraph == null)
            return;
        HashSet<Integer> toVisit = workingGraph.getTwoNeighbors(v);
        workingGraph.deleteVertex(v);
        for(Integer x : toVisit){
            for(Integer y : toVisit){
                if(y < x) continue;
                if(!conflictGraph.adjacent(x, y) && !verticesAreCloseEnough(workingGraph, x, y)){
                    // x and y are in conflict after deleting v
                    deletionConflicts.computeIfAbsent(v, k -> new LinkedList<>());
                    deletionConflicts.get(v).add(new Edge(x, y));
                    conflictGraph.addEdge(x, y);
                }
            }
        }
        conflictGraph.deleteVertex(v);
        workingGraph.undeleteVertex(v); // v gets deleted by deleteVertex(int v)
    }

    private static void rollbackConflictGraph(int v){
        if(conflictGraph == null)
            return;
        conflictGraph.undeleteVertex(v);
        List<Edge> conflictsToRemove = deletionConflicts.get(v);
        if(conflictsToRemove == null) return;
        for(Edge e : conflictsToRemove){
            // remove conflict edge of pair
            Integer x = e.getV();
            Integer y = e.getW();
            conflictGraph.getNeighbors(x).remove(y);
            conflictGraph.getNeighbors(y).remove(x);
        }
        deletionConflicts.remove(v);
    }

    public static boolean verticesAreCloseEnough(RtcGraph graph, int v, int w){
        // vertices are close enough when they are adjacent or have a common neighbour
        if(v == w)
            return true;
        if(graph == null || !graph.contains(v) || !graph.contains(w))
            return false;

        // TODO: use conflict graph if possible
        /*
        if(useVertexCoverRule && conflictGraph != null && graph == workingGraph)
            return !conflictGraph.adjacent(v, w);
        */

        // the shortest path can be over two edges/three vertices max
        // therefore when v and w are adjacent or share a neighbour
        if(graph.adjacent(v, w)){
            return true;
        }else{
            for(int i : graph.getNeighbors(v))
                if(graph.getNeighbors(i).contains(w))
                    return true;
        }
        return false;
    }

    private static void printInfo(){
        // add used settings to file
        if(fileWriter != null){
            printOutputLine("Settings: ", 2);
            printOutputLine("Branching Method: "
                    + (branchingMethod == Branching.SINGLE_VERTEX ? "Single Vertex" : "Incompatible Vertices"),2);
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

    private static ArrayList<Integer> kernelQuickSort(ArrayList<Integer> toSort){
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

    public static int getHighestVertexID(RtcGraph g){
        int max = 0;
        try{
            for(int i : g.getVertices())
                if(Integer.parseInt(g.getVertexName(i)) > max)
                    max = i;
            max++;
        }catch(Exception e){
            max = g.getVertexNames().size();
        }

        return max > g.size() ? max : g.size();
    }

    private static void applyNoChoiceRule(){
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
        /*
                // marked vertex in single triangle
        for(int v : marked){
            if(isPartOfSingleTriangle(workingGraph, v)){
                for(int i : getTriangleVertices(workingGraph, v)) {
                    if(!marked.contains(i))
                        toMark.add(i);
                }
            }
        }
        */
        for(int v : marked){
            Set<Integer> vertices = getAllTriangleSharedVertices(workingGraph, v);
            if(vertices != null) toMark.addAll(vertices);
        }

        toMark.forEach(TriangleTwoClubAlgorithm::markVertex);
        noChoiceRuleTime += System.nanoTime() - time;
        noChoiceRuleMarked += toMark.size();
    }

    /**
     * Finds the vertex that has the most incompatibilities to other vertices in the graph
     * @param graph the graph
     * @return the vertex id, -1 when no vertices are in the graph or graph is null
     */
    public static int getMostIncompatibleVertex(RtcGraph graph){
        if(graph == null)
            return -1;
        if(conflictGraph != null && graph == workingGraph){
            int worstVertex = -1;
            int incompatibilities = 0;
            for(Integer v : conflictGraph.getVertices()){
                int i = conflictGraph.getNeighbors(v).size();
                if(i > incompatibilities && !marked.contains(v)){
                    incompatibilities = i;
                    worstVertex = v;
                }
            }
            if(worstVertex != -1) return worstVertex;
        }
        int worstVertex = -1;
        int incompatibilities = 0;
        for(int i : graph.getVertices()){
            if(marked.contains(i)) continue;
            int currentCount = 0;
            for(int j : graph.getVertices()){
                if(!verticesAreCloseEnough(graph, i, j))
                    currentCount++;
            }
            if(currentCount > incompatibilities){
                worstVertex = i;
                incompatibilities = currentCount;
            }
        }
        return worstVertex;
    }

    /** Determines if given v is part of exactly one triangle
     * @param graph the graph
     * @param v the vertex id in the graph
     * @return true when v is part of exactly one triangle
     */
    public static boolean isPartOfSingleTriangle(RtcGraph graph, int v){
        return countTriangles(graph, v, 2) == 1;
    }

    /**
     * Counts the number of vertex triangles vertex v is part of. Choose Integer.MAX_VALUE to count all.
     * @param graph a graph
     * @param v the id of a vertex
     * @param max the maximum number to count to
     * @return the number of triangles, at most max
     */
    public static int countTriangles(RtcGraph graph, int v, int max){
        if(graph == null || max < 1 || !graph.contains(v))
            return -1;

        int triangles = 0;
        HashSet<Integer> seen = new HashSet<>();
        for(int w : graph.getNeighbors(v)){
            for(int x : graph.getNeighbors(w)){
                if(!graph.getNeighbors(x).contains(v) || seen.contains(x)) continue; // no new triangle
                triangles++;
            }
            if(triangles >= max) return triangles;
            seen.add(w);
        }
        return triangles;
    }

    /**
     *  Finds vertex id's of a triangle containing v
     * @param graph the graph
     * @param v the vertex that needs to be part of the triangle
     * @return int[] of length 3, containing v and it's triangle neighbours (if there are any)
     */
    private static int[] getTriangleVertices(RtcGraph graph, int v){
        int[] triangle = new int[3];
        triangle[0] = v;
        outer:
        for(int i : graph.getNeighbors(v)){
            for(int j : graph.getNeighbors(i)){
                if(graph.adjacent(j, v)){
                    triangle[1] = i;
                    triangle[2] = j;
                    break outer;
                }
            }
        }
        return triangle;
    }

    /**
     * Approximates the size of a vertex cover. Used with conflict graph.
     * @param g the conflict graph
     * @return the minimum number of vertices you need to delete for a triangle 2 club
     */
    public static int approximateVertexCoverSize(RtcGraph g, int limit){
        Stack<Integer> vcRollbackStack = new Stack<>();
        int count = 0;
        while(g.getEdgeCount() > 0){
            if(count > limit) break; // break if limit reached
            int toDelete = -1, best = 0;
            for(Integer v : g.getVertices()){
                if(g.getNeighbors(v).size() > best){
                    toDelete = v;
                    best = g.getNeighbors(v).size();
                }
            }
            vcRollbackStack.push(toDelete);
            g.deleteVertex(toDelete);
            count++;
        }
        while(!vcRollbackStack.empty()) g.undeleteVertex(vcRollbackStack.pop());
        return count;
    }

    /**
     * Uses the matching bound approximation for vertex cover.
     * @param working the working graph
     * @param conflict the conflict graph
     * @return the minimum number of vertices you need to delete for a triangle 2 club
     */
    public static int matchingBoundVertexCover(RtcGraph working, RtcGraph conflict, int limit){
        // 2-approximation for vertex cover of conflict graph
        // Minimal vertex cover of conflict graph is lower bound on vertices to delete
        HashSet<Integer> marks = new HashSet<Integer>();
        int size = 0;
        for (int v : conflict.getVertices())
            if (working.contains(v) && !marks.contains(v))
                for (int w : conflict.getNeighbors(v))
                    if (working.contains(w) && w > v)					// every remaining conflict is an edge
                        if (!marks.contains(v) && !marks.contains(w)) {	// if the edge is not covered
                            marks.add(v); marks.add(w); size++;			// cover edge with both end points
                            if (size > limit) return size;				// if current size exceeds allowed limit
                            break;										// v us marked, so the condition of the if-statement is never met again.
                        }												// of deletions we do not need to carry on
        return size;
    }

    /**
     * Calculates lower bound for a triangle 2-club.
     * Counts non-isolated vertices in the neighbourhood of every vertex to find a lower bound.
     * @param graph the graph
     * @return the lower bound solution
     */
    public static RtcGraph getLowerBound(RtcGraph graph){
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


    /**
     * Gets the current number of triangles a vertex is part of based on the working graph
     * @param v the id of a vertex
     * @return the number of triangles in v is part of in the working graph,
     * -1 if v does not exist or the triangle map is not initialized
     */
    private static int countTrianglesFromMap(int v){
        return countTrianglesFromMap(v, Integer.MAX_VALUE);
    }

    /**
     * Gets the current number of triangles a vertex is part of based on the working graph
     * @param v the id of a vertex
     * @return the number of triangles in v is part of in the working graph,
     * -1 if v does not exist or the triangle map is not initialized
     */
    private static int countTrianglesFromMap(int v, int max){
        if(!workingGraph.contains(v) || triangles == null)
            return -1;

        if(!triangles.containsKey(v)) // v was not part of any triangle at all
            return 0;

        int n = 0;
        for(Triangle t : triangles.get(v)){
            if(t.exists(workingGraph)) // only count triangles that currently exist
                n++;
            if(n >= max)
                break;
        }

        return n;
    }



    /**
     * Calculates upper bound for a triangle 2-club. Uses method 'int calculateLocalUpperBound(RtcGraph, int)' on every
     * vertex in the graph for global maximum.
     * @param graph the graph
     * @return upper bound for a triangle 2-club
     */
    public static int calculateUpperBound(RtcGraph graph){
        int best = 0, current;
        for(Integer v : graph.getVertices()){
            current = calculateLocalUpperBound(graph, v);
            if(current > best) best = current;
        }
        return best;
    }

    /**
     * Calculates upper bound for a triangle 2-club containing v.
     * Returns the sum of the degrees of every neighbour of v.
     * @param graph the graph
     * @param v id of a vertex in the graph
     * @return upper bound for a triangle 2-club
     */
    public static int calculateLocalUpperBound(RtcGraph graph, int v){
        int sum = 0;
        for(Integer w : graph.getNeighbors(v)) sum += graph.getNeighbors(w).size();
        return sum + 1;
    }

    /**
     * Calculates the intersection of every vertices that share a triangle with v
     * @param g the graph
     * @param v center vertex
     * @return a list of intersecting triangle vertices not including v,
     *          list may be empty or null if v is not part of a triangle
     */
    public static Set<Integer> getAllTriangleSharedVertices(RtcGraph g, Integer v){
        if(g == null || v < 0 || !g.contains(v) || !isPartOfTriangle(g, v))
            return null;
        Set<Integer> intersection = new HashSet<>(g.getNeighbors(v)); // starting will all
        // remove vertices that don't apply to condition
        for(int u : g.getNeighbors(v)){
            // checking if u shares every triangle with v
            boolean removeU = false;
            if(!isPartOfTriangle(g, u)){
                removeU = true;
            }else{
                for(int w : g.getNeighbors(v)){
                    if(u == w) continue;
                    if(verticesShareTriangle(g, v, w) && !verticesShareTriangle(g, u, w)){
                        // if v is in triangle with w but not with u, u is not part of every triangle of v
                        removeU = true;
                        break;
                    }
                }
            }
            if(removeU) intersection.remove(u);
        }
        return intersection;
    }

    /**
     * Removes the triangle of removed vertices, a removed vertex is part of the inputGraph but not of the workingGraph.
     * This aims to reduce the number of triangles in sparse graphs.
     */
    private static void removeNonExistingTriangles(){
        for(int v : inputGraph.getVertices()){
            if(workingGraph.contains(v) || !triangles.containsKey(v))
                continue;

            HashSet<Triangle> tri = triangles.remove(v);
            for(Triangle t : tri){
                for(int w : t.getVertices()){
                    if(triangles.containsKey(w))
                        triangles.get(w).remove(t);
                }
            }
        }
    }

    /**
     * Builds the conflict graph for vertex cover rule.
     * The conflict graph contains the same vertices as the input graph. The conflict graph contains an edge between
     * any vertices v and w if and only if v and w are incompatible.
     * @param graph the input graph
     * @return the corresponding conflict graph
     */
    public static RtcGraph buildConflictGraph(RtcGraph graph){
        RtcGraph conGraph = new RtcGraph();
        HashSet<Integer> vertices = graph.getVertices();
        vertices.forEach( v -> conGraph.addVertex(v, graph.getVertexName(v)));
        for(Integer v : vertices)
            for(Integer w : vertices)
                if( !v.equals(w) && !verticesAreCompatible(graph, v, w))
                    conGraph.addEdge(v, w);
        return conGraph;
    }

    /**
     * Sets a time limit for the algorithm. Branching will stop and the best current result will be returned by the
     * run() method. Time set to =< 0 will disable the the time limit.
     * Default is 0.
     * @param seconds time in seconds after which the run method will abort
     */
    public static void setTimeLimit(float seconds){
        timeLimitMS = (long) (seconds * 1000.0);
    }

    public static void setOutputFile(String path){
        File f = new File(path);
        try {
            fileWriter = new PrintWriter(new FileWriter(f), true);
        } catch (IOException e) {
            printOutputLine("Failed to initialize file writing", 1);
        }
    }

    /**
     * Sets a mode for branching
     * @param b the wanted method
     */
    public static void setBranchingMethod(Branching b){
        if(b == null)
            printOutputLine("null is not a valid branching method", 1);
        else
            branchingMethod = b;
    }

    /**
     * @param b If true, vertex cover rule will be used while branching
     */
    public static void useVertexCoverRule(boolean b){
        useVertexCoverRule = b;
    }

    /**
     * Sets frequency in which the vertex cover rule is used while branching
     * @param steps use vc rule every [steps] branching depth
     */
    public static void setVertexCoverRuleSteps(int steps){
        if(steps > 0 && steps < 10)
            vertexCoverRuleSteps = steps;
    }

    /**
     * @param b If true, vertex cover rule will be used while branching
     */
    public static void useNoChoiceRule(boolean b){
        useNoChoiceRule = b;
    }

    /**
     * l is the number of triangles each vertex has to be part of
     * @param l
     */
    public static void setL(int l){
        if(l > 0){
            ell = l;
        }else{
            printOutputLine("Invalid value for l: " + l + ", using " + ell + "instead", 1);
        }
    }

    private static void resetCounts(){
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
        triangleComputationTime = 0L;
    }

    private static int getRollbackPoint(){
        if(rollbackStack == null || rollbackStack.empty())
            return -1;
        else
            return rollbackStack.peek();
    }

    /**
     * Writes output to console and/or file
     * @param line the line to output
     * @param mode 0 = file + console, 1 = console only, 2 = file only
     */
    private static void printOutputLine(String line, int mode){
        if(printDetails && mode < 2) System.out.println(line);
        if(fileWriter != null && (mode == 0 || mode == 2)) fileWriter.println(line);
    }

    private static void printOutputLine(String line){
        printOutputLine(line, 0);
    }
}
