package robustTwoClub.graph;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.*;

public class RollbackRtcGraph extends RtcGraph{

    private Stack<Edge> rollbackStack;

    public RollbackRtcGraph(){
        super();
        initRollbackStack();
    }

    public RollbackRtcGraph(String filename, int type){
        super(filename, type);
        initRollbackStack();
    }

    /**
     * Private constructor used to clone an instance of this class
     * @param g the graph to clone
     */
    public RollbackRtcGraph(RtcGraph g){
        this();

        // clone vertices + vertex names
        for (int v : g.getVertices())
            addVertex(v, g.getVertexName(v));

        // clone edges
        for (int v : g.getVertices())
            for (int w : g.getNeighbors(v))
                if (w > v)
                    addEdge(v, w);
    }

    /**
     * Removes a vertex and its edges from the graph. A rollbackpoint made before deleting a vertex also restores
     * the removed edges.
     * @param vertex the id of the vertex that shall be deleted
     */
    @Override
    public void deleteVertex(int vertex) {
        super.deleteVertex(vertex);
        rollbackStack.push(new Edge(vertex, null));
        HashSet<Integer> neighbors = (HashSet<Integer>) getNeighbors(vertex).clone();
        for(Integer i : neighbors)
            deleteEdge(vertex, i);
    }

    /**
     * Clones the current graph.
     * The rollback stack is not cloned. Undelete operations cannot undelete vertices/edges that were removed
     * before cloning.
     * @return a new graph with the same vertices and edges as the current one
     */
    @Override
    public RollbackRtcGraph getClone() {
        return new RollbackRtcGraph(this);
    }

    /**
     * Rollback to a specific point, restores all deleted vertices and edges to this point
     * @param rollbackPoint the point to rollback too
     */
    public void rollback(Edge rollbackPoint){
        Edge edge;
        while((edge = rollbackStack.peek()) != null){
            if(edge.equals(rollbackPoint)) break; // we reached the requested rollback point
            rollbackStack.pop();
            if(edge.hasNull()){
                // restore a vertex
                Integer v = edge.getV();
                if(v == null) v = edge.getW(); // get the not null field
                undeleteVertex(v);
            }else{
                // restore an edge
                addEdge(edge.getV(), edge.getW());
            }
        }
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
     * Returns the rollback stack up to a certain point
     * @param rollbackPoint the last point you want to rollback to
     * @return the top of the actual rollback stack
     */
    public Stack<Edge> getRollbacksTo(Edge rollbackPoint){
        Stack<Edge> copy = new Stack<>(), movingStack = new Stack<>();
        Edge e;
        // move relevant part to temporary stack
        while((e = rollbackStack.peek()) != null){
            if(e.equals(rollbackPoint)) break; // we reached the requested rollback point
            rollbackStack.pop();
            movingStack.push(e);
        }
        // copy stack items
        while(!movingStack.empty()){
            e = movingStack.pop();
            rollbackStack.push(e);
            copy.push(e);
        }

        return copy;
    }

    /**
     * Removes an edge between two vertices.
     * @param v first vertex
     * @param w second vertex
     */
    public void deleteEdge(Integer v, Integer w){
        getNeighbors(v).remove(w);
        getNeighbors(w).remove(v);
        rollbackStack.push(new Edge(v, w));
    }

    /**
     * Removes an edge of the graph.
     * @param e the edge
     */
    public void deleteEdge(Edge e){
        deleteEdge(e.getV(), e.getW());
    }

    private void initRollbackStack(){
        rollbackStack = new Stack<>();
        rollbackStack.push(null);
    }
}
