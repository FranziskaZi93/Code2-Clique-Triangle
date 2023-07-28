package triangleTwoClub.graph;

//import jdk.jshell.spi.ExecutionControl;
import robustTwoClub.graph.Edge;
import robustTwoClub.graph.RollbackRtcGraph;
import robustTwoClub.graph.RtcGraph;
import robustTwoClub.graph.Triangle;
import triangleTwoClub.ConditionChecker;
import triangleTwoClub.ITriangleTwoClubAlgorithm;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Stack;

/**
 * A graph specifically for the conflict graph. Removing vertices and edges in a graph can only lead to adding edges in
 * the conflict graph. This graph class logs edges that have been added and is able to rollback the edges to a given
 * point.
 */
public class ConflictGraph extends RtcGraph {

    private RollbackRtcGraph original;
    private Stack<Edge> rollbackStack;
    private ITriangleTwoClubAlgorithm algorithm; // required for new conflict definition
    private boolean useCompatibilityDefinitionV2;

    public ConflictGraph(RollbackRtcGraph graph){
        super();
        this.original = graph;

        // clone vertices
        for(int v : graph.getVertices()){
            super.addVertex(v, graph.getVertexName(v));
        }

        useCompatibilityDefinitionV2 = false;
        initRollbackStack();
    }

    @Override
    public void addVertex(Integer v, String name) {
        System.out.println("Adding vertices is not allowed in the conflict graph");
    }

    @Override
    public void addEdge(Integer v, Integer w) {
        super.addEdge(v, w);
        Edge e = new Edge(v, w);
        rollbackStack.push(e);
    }

    /**
     * Get the current state of the graph as rollback point. Revert changes to the graph later by using
     * rollback()-Method with this point
     * @return the rollback point or null if there is nothing to rollback
     */
    public Edge getRollbackPoint(){
        return rollbackStack.empty() ? null : rollbackStack.peek();
    }

    /**
     * Rollback to a specific point, restores all deleted vertices and edges to this point
     * @param rollbackPoint the point to rollback too
     */
    public void rollback(Edge rollbackPoint){
        Edge edge;
        while((edge = rollbackStack.peek()) != null){
            if(edge.equals(rollbackPoint)) break; // we reached the requested rollback point
            removeEdge(rollbackStack.pop());
        }
    }

    public void initRollbackStack(){
        rollbackStack = new Stack<>();
        rollbackStack.push(null);
    }

    /**
     * Creates a conflict graph based on the current graphs state.
     */
    public void buildConflictGraph(){
        for(int v : original.getVertices())
            for(int w : original.getVertices()){
                if (w <= v) continue;
                if ((useCompatibilityDefinitionV2 && !algorithm.getConditionChecker().verticesAreCompatibleV2(v, w))
                        || !ConditionChecker.verticesAreCompatible(original, v, w))
                    addEdge(v, w);
            }

        initRollbackStack();
    }

    /**
     * Updates the conflict graph for a deletion.
     * @param v the first vertex id
     * @param w the second vertex id if it's an edge, or null if vertex v shall get deleted
     */
    public void updateConflictsForDeletion(Integer v, Integer w){
        if(useCompatibilityDefinitionV2){
            // use new methods
            if(w == null)
                updateConflictsForVertexDeletionV2(v);
            else
                updateConflictsForEdgeDeletionV2(v, w);
        }else{
            // use old methods
            if(w == null)
                updateConflictsForVertexDeletionV1(v);
            else
                updateConflictsForEdgeDeletionV1(v, w);
        }
    }

    /**
     * Sets the algorithm for computation of V2 conflict graph. Triangles and L is required from the algorithm.
     * @param algorithm the triangles of the original graph
     */
    public void setAlgorithm(ITriangleTwoClubAlgorithm algorithm){
        this.algorithm = algorithm;
    }

    /**
     * Enables/disables the use of the better compatibility definition. Requires an algorithm to be set.
     * @param flag true if the better definition should be used
     */
    public void useBetterCompatibility(boolean flag){
        if(flag && (algorithm == null || algorithm.getConditionChecker() == null)){
            throw new RuntimeException("Can not use better compatibility without an algorithm object");
        }
        useCompatibilityDefinitionV2 = flag;
    }

    private void removeEdge(Edge e){
        getNeighbors(e.getV()).remove(e.getW());
        getNeighbors(e.getW()).remove(e.getV());
    }


    private void updateConflictsForVertexDeletionV1(int v){
        HashSet<Integer> toVisit = new HashSet<>(original.getNeighbors(v));

        Edge rollbackPoint = original.getRollbackPoint();
        original.deleteVertex(v);


        for(Integer x : toVisit){
            for(Integer y : toVisit){
                if(y <= x) continue;
                if(!adjacent(x, y) && !ConditionChecker.verticesAreCompatible(original, x, y)){
                    // x and y are in conflict after deleting v
                    addEdge(x, y);
                }
            }
        }
        original.rollback(rollbackPoint);
    }

    private void updateConflictsForEdgeDeletionV1(int v, int w){
        Edge rollbackPoint = original.getRollbackPoint();
        original.deleteEdge(v, w);

        // endpoint v with neighbors of w
        for(int u : original.getNeighbors(w)){
            if(v != u && !ConditionChecker.verticesAreCompatible(original, v, u))
                addEdge(v, u);
        }

        // endpoint w with neighbors of v
        for(int u : original.getNeighbors(v)){
            if(w != u && !ConditionChecker.verticesAreCompatible(original, w, u))
                addEdge(w, u);
        }

        original.rollback(rollbackPoint);
    }

    private void updateConflictsForVertexDeletionV2(int v){
        HashSet<Integer> toVisit = original.getNeighbors(v);

        Edge rollbackPoint = original.getRollbackPoint();
        original.deleteVertex(v);

        for(Integer x : toVisit){
            // check pairs from N(v)
            for(Integer y : toVisit){
                if(y <= x) continue;
                if(!adjacent(x, y) && !algorithm.getConditionChecker().verticesAreCompatibleV2(x, y)){
                    // x and y are in conflict after deleting v
                    addEdge(x, y);
                }
            }

            // check neighbors x \in N(v) with (N_2(v) \cap N_2(x))
            for(Integer y : original.getCommonTwoNeighborhood(v, x)){
                if(y.equals(x)) continue;
                if(!adjacent(x, y) && !algorithm.getConditionChecker().verticesAreCompatibleV2(x, y)){
                    // x and y are in conflict after deleting v
                    addEdge(x, y);
                }
            }
        }

        original.rollback(rollbackPoint);
        throw new RuntimeException("Not implemented yet");
    }

    private void updateConflictsForEdgeDeletionV2(int v, int w){
        Edge rollbackPoint = original.getRollbackPoint();
        original.deleteEdge(v, w);

        HashSet<Integer> commonNeighborhood = original.getCommonTwoNeighborhood(v, w);

        // edge (v, w): check v with N(w) and u with N(v), and


        // endpoint v with neighbors of w
        for(int u : original.getNeighbors(w)){
            if(v != u && !algorithm.getConditionChecker().verticesAreCompatibleV2(v, u))
                addEdge(v, u);
        }

        // endpoint w with neighbors of v
        for(int u : original.getNeighbors(v)){
            if(w != u && !algorithm.getConditionChecker().verticesAreCompatibleV2(w, u))
                addEdge(w, u);
        }

        // check v and w with (N_2(v) \cap N_2(w))
        for(int u : original.getCommonTwoNeighborhood(v, w)){
            if(v != u && !algorithm.getConditionChecker().verticesAreCompatibleV2(v, u))
                addEdge(v, u);
            if(w != u && !algorithm.getConditionChecker().verticesAreCompatibleV2(w, u))
                addEdge(w, u);
        }


        original.rollback(rollbackPoint);
        throw new RuntimeException("Not implemented yet");
    }
}
