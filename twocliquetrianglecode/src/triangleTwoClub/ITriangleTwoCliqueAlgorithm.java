package triangleTwoClub;

import robustTwoClub.graph.RtcGraph;
import triangleTwoClub.graph.T2CGraph;

import java.util.Map;
import java.util.Set;

public interface ITriangleTwoCliqueAlgorithm extends ITriangleTwoClubAlgorithm{

    Map<Integer, Set<Integer>> identifyVerticesWithDistanceMaxTwo(T2CGraph graph);

    public Set<Set<Integer>> getSetsOfTriangle();

    public Map<Integer,Set<Set<Integer>>> mapRemainingTrianglesInGraph(T2CGraph globalGraph);

    public boolean containsOnlyTriangle(RtcGraph graph);

    public boolean isClique(RtcGraph graph);

    Set getElementsOfSetOneNotInSetTwo(Set setOne, Set setTwo);

    Set overlapOfTwoSets(Set setA, Set setB);

    T2CGraph getOriginalGraph();
}
