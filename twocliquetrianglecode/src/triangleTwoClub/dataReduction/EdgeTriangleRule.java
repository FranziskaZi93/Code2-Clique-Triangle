package triangleTwoClub.dataReduction;

import robustTwoClub.graph.Edge;
import triangleTwoClub.ITriangleTwoClubAlgorithm;

import java.util.HashSet;
import java.util.LinkedList;

/**
 * Reduction Rule 3 - Triangle Rule (Edge)
 * Removes every edge that is not part of at least L triangles.
 */
public  class EdgeTriangleRule extends DeleteRule {

    public EdgeTriangleRule(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    @Override
       int reduce() {
      HashSet<Integer> seen = new HashSet<>();
        LinkedList<Edge> toDelete = new LinkedList<>();

        for(Integer i : algorithm.getGraph().getVertices()){
            if(seen.contains(i)) continue;
            for(Integer j : algorithm.getGraph().getNeighbors(i)){
                if(seen.contains(j)) continue;
                if(algorithm.getConditionChecker().countEdgeTriangles(i, j, algorithm.getL()) < algorithm.getL()){
                    toDelete.add(new Edge(i, j));
                }
            }
            seen.add(i);
        }

        toDelete.forEach(e -> algorithm.deleteEdge(e.getV(), e.getW()));
        return toDelete.size();
    }

    @Override
    public String getRuleName() {
        return "EdgeTriangleRule";
    }

    @Override
    boolean useUntilExhaustion() {
        return true;
    }
}
