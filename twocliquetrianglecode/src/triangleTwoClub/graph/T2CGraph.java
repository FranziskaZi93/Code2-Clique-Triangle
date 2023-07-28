package triangleTwoClub.graph;

import robustTwoClub.graph.Edge;
import robustTwoClub.graph.RollbackRtcGraph;
import robustTwoClub.graph.RtcGraph;
import triangleTwoClub.ConditionChecker;
import triangleTwoClub.ITriangleTwoClubAlgorithm;
import triangleTwoClub.Output;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Stack;

/**
 * Implements a RollbackRtcGraph with a corresponding conflict graph
 */
public class T2CGraph {

    public enum ConflictGraphType{DISABLED, SIMPLE, ADVANCED}

    private RollbackRtcGraph graph;
    private ConflictGraph conflictGraph;
    private HashSet<Integer> markedVertices;
    private Stack<Integer> markRollbackStack;
    private HashMap<Edge, RollbackPoint> rollbackMap;
    private Output output;
    private ITriangleTwoClubAlgorithm algorithm;
    private ConflictGraphType conflictType;

    private static String conflictGraphTimeName = "ConflictGraphTime[s]";
    private static String conflictGraphTypeName = "ConflictGraphType";

    public T2CGraph(RtcGraph g, ITriangleTwoClubAlgorithm algorithm){
        if(g.getClass().equals(RollbackRtcGraph.class))
            this.graph = (RollbackRtcGraph) g;
        else
            this.graph = new RollbackRtcGraph(g);

        this.algorithm = algorithm;
        this.output = algorithm.getOutput();
        this.conflictType = ConflictGraphType.SIMPLE;
        output.setFormat(conflictGraphTimeName, Output.DataFormat.MILLISECONDS);
        resetConflictGraphAndMarks();
        rollbackMap = new HashMap<>();
    }

    public void deleteVertex(int v){
        if(!graph.contains(v))
            return;

     //   updateConflictsForDeletion(v, null);
        graph.deleteVertex(v);
        markedVertices.remove(v);
    }

    public void deleteEdge(int v, int w){
        if(!graph.contains(v) || !graph.contains(w) || !graph.adjacent(v, w))
            return;

        updateConflictsForDeletion(v, w);
        graph.deleteEdge(v, w);
    }

    public void markVertex(int v){
        markedVertices.add(v);
        markRollbackStack.push(v);
    }

    public Edge getRollbackPoint(){
        Edge rbp = graph.getRollbackPoint();
        Edge conflictRBP = conflictGraph == null ? null : conflictGraph.getRollbackPoint();
        Integer marksRBP = markRollbackStack.peek();
        rollbackMap.put(rbp, new RollbackPoint(conflictRBP, marksRBP));
        return rbp;
    }

    public void rollback(Edge rollbackPoint){
        RollbackPoint p = rollbackMap.get(rollbackPoint);
        if(p != null){
            if(conflictGraph != null)
                conflictGraph.rollback(p.getConflictRollbackPoint());
            rollbackMarked(p.getMarkRollbackPoint());
        }
        graph.rollback(rollbackPoint);
    }

    public void resetMarks(){
        this.markedVertices = new HashSet<>();
        markRollbackStack = new Stack<>();
        markRollbackStack.push(null); // use null ase base to avoid try/catch for EmptyStackException
    }

    public RollbackRtcGraph getGraph(){
        return this.graph;
    }

    public RtcGraph getConflictGraph(){
        return conflictGraph;
    }

    public HashSet<Integer> getMarkedVertices(){return  this.markedVertices;}

    private void updateConflictsForDeletion(Integer v, Integer w){
        if(conflictGraph == null || v == null || conflictType == ConflictGraphType.DISABLED)
            return;

        long time = System.currentTimeMillis();
       // conflictGraph.updateConflictsForDeletion(v, w);
        output.log(conflictGraphTimeName, System.currentTimeMillis() - time);
    }

    private void rollbackMarked(Integer rollbackPoint){
        Integer p;
        while(!markRollbackStack.isEmpty()){
            p = markRollbackStack.peek();
            if(p == null || p.equals(rollbackPoint)){
                break;
            }else
                markedVertices.remove(markRollbackStack.pop());
            }
    }

    /**
     * Generated the conflict graph for the current graph.
     */
    public void buildConflictGraph(){
        if(conflictType == ConflictGraphType.DISABLED)
            return;

        long time = System.currentTimeMillis();
        conflictGraph = new ConflictGraph(graph);
        conflictGraph.setAlgorithm(algorithm);
        conflictGraph.useBetterCompatibility(conflictType == ConflictGraphType.ADVANCED);
        conflictGraph.buildConflictGraph();
        output.log(conflictGraphTimeName, System.currentTimeMillis() - time);
    }

    /**
     * Resets the conflict graph and the marked vertices.
     */
    public void resetConflictGraphAndMarks(){
        conflictGraph = null;
        resetMarks();
        rollbackMap = new HashMap<>();
    }

    /**
     * Sets the conflict type used.
     * @param t the type
     */
    public void setConflictType(ConflictGraphType t){
        this.conflictType = t;
        if(conflictGraph != null)
            conflictGraph.useBetterCompatibility(t == ConflictGraphType.ADVANCED);
        output.addInformation(conflictGraphTypeName, t.toString());
    }

    /**
     * Represents a point to rollback the conflict graph and the marked vertices to.
     */
    private class RollbackPoint{
        private Edge conflictRollbackPoint;
        private Integer markRollbackPoint;

        private RollbackPoint(Edge conflictRollbackPoint, Integer markRollbackPoint){
            this.conflictRollbackPoint = conflictRollbackPoint;
            this.markRollbackPoint = markRollbackPoint;
        }

        private Edge getConflictRollbackPoint(){
            return conflictRollbackPoint;
        }

        private Integer getMarkRollbackPoint(){
            return markRollbackPoint;
        }
    }
}
