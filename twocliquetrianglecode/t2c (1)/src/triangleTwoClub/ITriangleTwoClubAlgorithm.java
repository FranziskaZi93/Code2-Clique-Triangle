package triangleTwoClub;

import robustTwoClub.graph.RtcGraph;
import robustTwoClub.graph.Triangle;

import java.util.HashMap;
import java.util.HashSet;

public interface ITriangleTwoClubAlgorithm {

    public RtcGraph getGraph();
    public RtcGraph getConflictGraph();
    public HashSet<Integer> getMarked();
    public HashMap<Integer, HashSet<Triangle>> getTriangles();
    public ConditionChecker getConditionChecker();
    public Output getOutput();
    public int getL();
    public int getBestSolutionSize();
    public boolean isEdgeConditionEnabled();

    public void abortBranch();
    public void deleteVertex(int v);
    public void deleteEdge(int v, int w);
    public void markVertex(int v);

}
