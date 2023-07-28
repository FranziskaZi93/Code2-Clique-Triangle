package triangleTwoClub.dataReduction;

import robustTwoClub.graph.RtcGraph;
import triangleTwoClub.ConditionChecker;
import triangleTwoClub.ITriangleTwoClubAlgorithm;

import java.util.HashSet;

/**
 * Reduction Rule 3 - Triangle Rule (Vertex)
 * Removes every vertex that is not part of at least L triangles.
 */
public class VertexTriangleRule extends DeleteRule {

    public VertexTriangleRule(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    @Override
    int reduce() {
        RtcGraph graph = algorithm.getGraph();
        int l = algorithm.getL();
        HashSet<Integer> toDelete = new HashSet<>();

        for(int v : graph.getVertices()){
            if(ConditionChecker.countVertexTrianglesFromMap(graph, algorithm.getTriangles(), v, l) < l)
                toDelete.add(v);
        }

        toDelete.forEach(v -> algorithm.deleteVertex(v));
        return toDelete.size();
    }

    @Override
    public String getRuleName() {
        return "VertexTriangleRule";
    }

    @Override
    boolean useUntilExhaustion() {
        return true;
    }
}
